MODULE tradmp
   !!======================================================================
   !!                       ***  MODULE  tradmp  ***
   !! Ocean physics: internal restoring trend on active tracers (T and S)
   !!======================================================================
   !! History :  OPA  ! 1991-03  (O. Marti, G. Madec)  Original code
   !!                 ! 1992-06  (M. Imbard)  doctor norme
   !!                 ! 1998-07  (M. Imbard, G. Madec) ORCA version
   !!            7.0  ! 2001-02  (M. Imbard)  add distance to coast, Original code
   !!            8.1  ! 2001-02  (G. Madec, E. Durand)  cleaning
   !!  NEMO      1.0  ! 2002-08  (G. Madec, E. Durand)  free form + modules
   !!            3.2  ! 2009-08  (G. Madec, C. Talandier)  DOCTOR norm for namelist parameter
   !!            3.3  ! 2010-06  (C. Ethe, G. Madec) merge TRA-TRC 
   !!            3.4  ! 2011-04  (G. Madec, C. Ethe) Merge of dtatem and dtasal + suppression of CPP keys
   !!            3.6  ! 2015-06  (T. Graham)  read restoring coefficient in a file
   !!            3.7  ! 2015-10  (G. Madec)  remove useless trends arrays
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_dmp_alloc : allocate tradmp arrays
   !!   tra_dmp       : update the tracer trend with the internal damping
   !!   tra_dmp_init  : initialization, namlist read, parameters control
   !!----------------------------------------------------------------------
   USE oce            ! ocean: variables
   USE dom_oce        ! ocean: domain variables
   USE c1d            ! 1D vertical configuration
   USE trd_oce        ! trends: ocean variables
   USE trdtra         ! trends manager: tracers 
   USE zdf_oce        ! ocean: vertical physics
   USE phycst         ! physical constants
   USE dtatsd         ! data: temperature & salinity
   USE zdfmxl         ! vertical physics: mixed layer depth
   !
   USE in_out_manager ! I/O manager
   USE iom            ! XIOS
   USE lib_mpp        ! MPP library
   USE prtctl         ! Print control
   USE timing         ! Timing
   USE stopack        ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_dmp        ! called by step.F90
   PUBLIC   tra_dmp_init   ! called by nemogcm.F90

   !                                           !!* Namelist namtra_dmp : T & S newtonian damping *
   LOGICAL            , PUBLIC ::   ln_tradmp   !: internal damping flag
   INTEGER            , PUBLIC ::   nn_zdmp     !: = 0/1/2 flag for damping in the mixed layer
   INTEGER            , PUBLIC ::   nn_shiter   !: Shapiro filter iterations (<=0 for no filter)
   INTEGER            , PUBLIC ::   nn_dfreq    !: Frequency of diff calculation
   CHARACTER(LEN=200) , PUBLIC ::   cn_resto    !: name of netcdf file containing restoration coefficient field
   !
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   resto    !: restoring coeff. on T and S (s-1)
   REAL(wp),         ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   z2da     !: restoring coeff. on T and S (s-1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   ztsdiff  !: restoring coeff. on T and S (s-1)

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: tradmp.F90 11536 2019-09-11 13:54:18Z smasson $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION tra_dmp_alloc()
      !!----------------------------------------------------------------------
      !!                ***  FUNCTION tra_dmp_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( resto(jpi,jpj,jpk), STAT= tra_dmp_alloc )
      ALLOCATE( ztsdiff(jpi,jpj,jpk,jpts), z2da(jpi,jpj), STAT= tra_dmp_alloc )
      CALL mpp_sum ( 'tradmp', tra_dmp_alloc )
      IF( tra_dmp_alloc > 0 )   CALL ctl_warn('tra_dmp_alloc: allocation of arrays failed')
      z2da(:,:) = 1._wp
      !
   END FUNCTION tra_dmp_alloc


   SUBROUTINE tra_dmp( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE tra_dmp  ***
      !!                  
      !! ** Purpose :   Compute the tracer trend due to a newtonian damping
      !!      of the tracer field towards given data field and add it to the
      !!      general tracer trends.
      !!
      !! ** Method  :   Newtonian damping towards t_dta and s_dta computed 
      !!      and add to the general tracer trends:
      !!                     ta = ta + resto * (t_dta - tb)
      !!                     sa = sa + resto * (s_dta - sb)
      !!         The trend is computed either throughout the water column
      !!      (nlmdmp=0) or in area of weak vertical mixing (nlmdmp=1) or
      !!      below the well mixed layer (nlmdmp=2)
      !!
      !! ** Action  : - tsa: tracer trends updated with the damping trend
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !
      INTEGER ::   ji, jj, jk, jn   ! dummy loop indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts)     ::  zts_dta
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE     ::  zzt
      REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  ztrdts
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('tra_dmp')
      !
      IF( l_trdtra .OR. iom_use('hflx_dmp_cea') .OR. iom_use('sflx_dmp_cea') ) THEN   !* Save ta and sa trends
         ALLOCATE( ztrdts(jpi,jpj,jpk,jpts) )
         ztrdts(:,:,:,:) = tsa(:,:,:,:)
      ENDIF

      !                           !==  input T-S data at kt  ==!
      !
      IF( MOD(kt,nn_dfreq) .EQ. 0 .OR. kt == nit000 ) THEN 
         CALL dta_tsd( kt, zts_dta )            ! read and interpolates T-S data at kt
         IF( nn_shiter > 0 ) ALLOCATE(zzt(jpi,jpj,jpk))
         DO jn=1,jpts
            ztsdiff(:,:,:,jn) = ( zts_dta(:,:,:,jn) - tsb(:,:,:,jn) ) * tmask(:,:,:)
            IF( nn_shiter > 0 ) THEN
!               DO jk=1,jpkm1
!                  zzt = ztsdiff(:,:,jk,jn)
!                  CALL Shapiro_1D( zzt, nn_shiter, ztsdiff(:,:,jk,jn), jk)
!               ENDDO
                zzt (:,:,:)       = ztsdiff(:,:,:,jn)
                CALL Shapiro_3D( zzt, nn_shiter, ztsdiff(:,:,:,jn) )
            ENDIF
         ENDDO
         IF( nn_shiter > 0 ) DEALLOCATE(zzt)
      ENDIF
      !
      IF( ln_stopack .AND. nn_spp_tdmp > 0 ) THEN
           z2da(:,:) = 1._wp
           CALL spp_gen(kt, z2da, nn_spp_tdmp, rn_tdmp_sd, jk_spp_tdmp )
      ENDIF
      !
      SELECT CASE ( nn_zdmp )     !==  type of damping  ==!
      !
      CASE( 0 )                        !*  newtonian damping throughout the water column  *!
         DO jn = 1, jpts
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1   ! vector opt.
                     tsa(ji,jj,jk,jn) = tsa(ji,jj,jk,jn) + z2da(ji,jj)*resto(ji,jj,jk) * ztsdiff(ji,jj,jk,jn) 
                  END DO
               END DO
            END DO
         END DO
         !
      CASE ( 1 )                       !*  no damping in the turbocline (avt > 5 cm2/s)  *!
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  IF( avt(ji,jj,jk) <= avt_c ) THEN
                    tsa(ji,jj,jk,jp_tem)=tsa(ji,jj,jk,jp_tem)+z2da(ji,jj)*resto(ji,jj,jk)*ztsdiff(ji,jj,jk,jp_tem)
                    tsa(ji,jj,jk,jp_sal)=tsa(ji,jj,jk,jp_sal)+z2da(ji,jj)*resto(ji,jj,jk)*ztsdiff(ji,jj,jk,jp_sal)
                  ENDIF
               END DO
            END DO
         END DO
         !
      CASE ( 2 )                       !*  no damping in the mixed layer   *!
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  IF( gdept_n(ji,jj,jk) >= hmlp (ji,jj) ) THEN
                     tsa(ji,jj,jk,jp_tem)=tsa(ji,jj,jk,jp_tem)+z2da(ji,jj)*resto(ji,jj,jk)*ztsdiff(ji,jj,jk,jp_tem)
                     tsa(ji,jj,jk,jp_sal)=tsa(ji,jj,jk,jp_sal)+z2da(ji,jj)*resto(ji,jj,jk)*ztsdiff(ji,jj,jk,jp_sal)
                  ENDIF
               END DO
            END DO
         END DO
         !
      END SELECT
      !
      IF( iom_use('hflx_dmp_cea') ) &
         & CALL iom_put('hflx_dmp_cea', SUM( ( tsa(:,:,:,jp_tem) - ztrdts(:,:,:,jp_tem) ) * e3t_n(:,:,:), dim=3 ) * rcp * rau0 ) ! W/m2
      IF( iom_use('sflx_dmp_cea') ) &
         & CALL iom_put('sflx_dmp_cea', SUM( ( tsa(:,:,:,jp_sal) - ztrdts(:,:,:,jp_sal) ) * e3t_n(:,:,:), dim=3 ) * rau0 )       ! g/m2/s

      IF( l_trdtra )   THEN       ! trend diagnostic
         ztrdts(:,:,:,:) = tsa(:,:,:,:) - ztrdts(:,:,:,:)
         CALL trd_tra( kt, 'TRA', jp_tem, jptra_dmp, ztrdts(:,:,:,jp_tem) )
         CALL trd_tra( kt, 'TRA', jp_sal, jptra_dmp, ztrdts(:,:,:,jp_sal) )
         DEALLOCATE( ztrdts ) 
      ENDIF
      !                           ! Control print
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=tsa(:,:,:,jp_tem), clinfo1=' dmp  - Ta: ', mask1=tmask,   &
         &                       tab3d_2=tsa(:,:,:,jp_sal), clinfo2=       ' Sa: ', mask2=tmask, clinfo3='tra' )
      !
      IF( ln_timing )   CALL timing_stop('tra_dmp')
      !
   END SUBROUTINE tra_dmp


   SUBROUTINE tra_dmp_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_dmp_init  ***
      !! 
      !! ** Purpose :   Initialization for the newtonian damping 
      !!
      !! ** Method  :   read the namtra_dmp namelist and check the parameters
      !!----------------------------------------------------------------------
      INTEGER ::   ios, imask   ! local integers 
      !
      NAMELIST/namtra_dmp/ ln_tradmp, nn_zdmp, cn_resto, nn_shiter, nn_dfreq
      !!----------------------------------------------------------------------
      !
      REWIND( numnam_ref )   ! Namelist namtra_dmp in reference namelist : T & S relaxation
      READ  ( numnam_ref, namtra_dmp, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namtra_dmp in reference namelist' )
      !
      REWIND( numnam_cfg )   ! Namelist namtra_dmp in configuration namelist : T & S relaxation
      READ  ( numnam_cfg, namtra_dmp, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namtra_dmp in configuration namelist' )
      IF(lwm) WRITE ( numond, namtra_dmp )
      !
      IF(lwp) THEN                  ! Namelist print
         WRITE(numout,*)
         WRITE(numout,*) 'tra_dmp_init : T and S newtonian relaxation'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namtra_dmp : set relaxation parameters'
         WRITE(numout,*) '      Apply relaxation   or not       ln_tradmp   = ', ln_tradmp
         WRITE(numout,*) '         mixed layer damping option      nn_zdmp  = ', nn_zdmp
         WRITE(numout,*) '         Damping file name               cn_resto = ', cn_resto
         WRITE(numout,*) '         Shapiro filter iterations       nn_shiter= ', nn_shiter
         WRITE(numout,*) '         Frequency of calculations       nn_dfreq = ', nn_dfreq
         WRITE(numout,*)
      ENDIF
      !
      IF( ln_tradmp ) THEN
         !                          ! Allocate arrays
         IF( tra_dmp_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'tra_dmp_init: unable to allocate arrays' )
         !
         SELECT CASE (nn_zdmp)      ! Check values of nn_zdmp
         CASE ( 0 )   ;   IF(lwp) WRITE(numout,*) '   tracer damping as specified by mask'
         CASE ( 1 )   ;   IF(lwp) WRITE(numout,*) '   no tracer damping in the mixing layer (kz > 5 cm2/s)'
         CASE ( 2 )   ;   IF(lwp) WRITE(numout,*) '   no tracer damping in the mixed  layer'
         CASE DEFAULT
            CALL ctl_stop('tra_dmp_init : wrong value of nn_zdmp')
         END SELECT
         !
         IF( .NOT.ln_tsd_dmp ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout, *)  '   read T-S data not initialized, we force ln_tsd_dmp=T'
            CALL dta_tsd_init( ld_tradmp=ln_tradmp )        ! forces the initialisation of T-S data
         ENDIF
         !                          ! Read in mask from file
         CALL iom_open ( cn_resto, imask)
         CALL iom_get  ( imask, jpdom_autoglo, 'resto', resto )
         CALL iom_close( imask )
      ENDIF
      !
   END SUBROUTINE tra_dmp_init

   SUBROUTINE Shapiro_1D(rla_varin,id_np, rlpa_varout, nk) !GIG
      IMPLICIT NONE
      INTEGER, INTENT(IN)                       :: id_np, nk
      REAL(wp), DIMENSION(jpi,jpj), INTENT(IN)  :: rla_varin !GIG
      REAL(wp), DIMENSION(jpi,jpj), INTENT(OUT) :: rlpa_varout !GIG

      REAL(wp), DIMENSION(jpi,jpj)              :: rlpa_varout_tmp
      REAL, PARAMETER                           :: rl_alpha = 1./2.    ! fixed stability coefficient (isotrope case)
      REAL, parameter                           :: rap_aniso_diff_XY=2.25 !  anisotrope case
      REAL                                      :: alphax,alphay, znum, zden,test
      INTEGER                                   :: ji, jj, jn, nn
!
!------------------------------------------------------------------------------
!
! Loop on several filter iterations
       rlpa_varout(:,:) = rla_varin(:,:)
       rlpa_varout_tmp(:,:) = rlpa_varout(:,:)
!
       alphax=1./2.
       alphay=1./2.
!  Dx/Dy=rap_aniso_diff_XY  , D_ = vitesse de diffusion
!  140 passes du fitre, Lx/Ly=1.5, le rap_aniso_diff_XY correspondant est:
       IF ( rap_aniso_diff_XY .GE. 1. ) alphay=alphay/rap_aniso_diff_XY
       IF ( rap_aniso_diff_XY .LT. 1. ) alphax=alphax*rap_aniso_diff_XY
        DO jn = 1,id_np   ! number of passes of the filter
            DO jj = 2,jpjm1
               DO ji = 2,jpim1
                  ! We crop on the coast
                   znum = rlpa_varout_tmp(ji,jj)   &
                          + 0.25*alphax*(rlpa_varout_tmp(ji-1,jj)-rlpa_varout_tmp(ji,jj))*tmask(ji-1,jj ,nk)  &
                          + 0.25*alphax*(rlpa_varout_tmp(ji+1,jj)-rlpa_varout_tmp(ji,jj))*tmask(ji+1,jj ,nk)  &
                          + 0.25*alphay*(rlpa_varout_tmp(ji ,jj-1)-rlpa_varout_tmp(ji,jj))*tmask(ji  ,jj-1,nk)  &
                          + 0.25*alphay*(rlpa_varout_tmp(ji ,jj+1)-rlpa_varout_tmp(ji,jj))*tmask(ji  ,jj+1,nk)
                   rlpa_varout(ji,jj)=znum*tmask(ji,jj,nk)+rla_varin(ji,jj)*(1.-tmask(ji,jj,nk))
                ENDDO  ! end loop ji
            ENDDO  ! end loop jj
!
!
!           Periodical condition in case of cd_overlap (global ocean)
!           - on a mercator projection grid we consider that singular point at
!           poles
!             are a mean of the values at points of the previous latitude
!           - on ORCA and regular grid we copy the values at points of the
!           previous latitude
            call lbc_lnk('tradmp_shp', rlpa_varout, 'T', 1.) ! Boundary condition
            rlpa_varout_tmp(:,:) = rlpa_varout(:,:)
         ENDDO  ! end loop jn
!
END SUBROUTINE Shapiro_1D
!---
SUBROUTINE Shapiro_3D(rla_varin,id_np, rlpa_varout) !GIG
      IMPLICIT NONE
      INTEGER, INTENT(IN)                       :: id_np
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(IN)  :: rla_varin !GIG
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(OUT) :: rlpa_varout !GIG

      REAL(wp), DIMENSION(jpi,jpj,jpk)          :: rlpa_varout_tmp
      REAL, PARAMETER                           :: rl_alpha = 1./2.    ! fixed stability coefficient (isotrope case)
      REAL, parameter                           :: rap_aniso_diff_XY=2.25 !  anisotrope case
      REAL                                      :: alphax,alphay, znum(jpk), zden,test
      INTEGER                                   :: ji, jj, jn, nn
!
!------------------------------------------------------------------------------
!
! Loop on several filter iterations
       rlpa_varout(:,:,:) = rla_varin(:,:,:)
       rlpa_varout_tmp(:,:,:) = rlpa_varout(:,:,:)
!
       alphax=1./2.
       alphay=1./2.
!  Dx/Dy=rap_aniso_diff_XY  , D_ = vitesse de diffusion
!  140 passes du fitre, Lx/Ly=1.5, le rap_aniso_diff_XY correspondant est:
       IF ( rap_aniso_diff_XY .GE. 1. ) alphay=alphay/rap_aniso_diff_XY
       IF ( rap_aniso_diff_XY .LT. 1. ) alphax=alphax*rap_aniso_diff_XY
        DO jn = 1,id_np   ! number of passes of the filter
            DO jj = 2,jpjm1
               DO ji = 2,jpim1
                  ! We crop on the coast
                   znum(:) = rlpa_varout_tmp(ji,jj,:)   &
                          + 0.25*alphax*(rlpa_varout_tmp(ji-1,jj,:)-rlpa_varout_tmp(ji,jj,:))*tmask(ji-1,jj ,:)  &
                          + 0.25*alphax*(rlpa_varout_tmp(ji+1,jj,:)-rlpa_varout_tmp(ji,jj,:))*tmask(ji+1,jj ,:)  &
                          + 0.25*alphay*(rlpa_varout_tmp(ji ,jj-1,:)-rlpa_varout_tmp(ji,jj,:))*tmask(ji  ,jj-1,:)  &
                          + 0.25*alphay*(rlpa_varout_tmp(ji ,jj+1,:)-rlpa_varout_tmp(ji,jj,:))*tmask(ji  ,jj+1,:)
                   rlpa_varout(ji,jj,:)=znum*tmask(ji,jj,:)+rla_varin(ji,jj,:)*(1.-tmask(ji,jj,:))
                ENDDO  ! end loop ji
            ENDDO  ! end loop jj
!
!
!           Periodical condition in case of cd_overlap (global ocean)
!           - on a mercator projection grid we consider that singular point at
!           poles
!             are a mean of the values at points of the previous latitude
!           - on ORCA and regular grid we copy the values at points of the
!           previous latitude
            call lbc_lnk('tradmp_shp', rlpa_varout, 'T', 1.) ! Boundary condition
            rlpa_varout_tmp(:,:,:) = rlpa_varout(:,:,:)
         ENDDO  ! end loop jn
!
END SUBROUTINE Shapiro_3D
    !!======================================================================
END MODULE tradmp
