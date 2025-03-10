MODULE TRAMEFT
   !!======================================================================
   !!                       ***  MODULE trameft ***
   !!        Model error forcing term (estimated offline or online)
   !! 
   !!======================================================================
   !! History :       ! 2023-07  (A. Storto) Initial version (NEMO 4.0.7)
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_meft_init  : Initialize the model errors
   !!   tra_meft       : Apply the tracer (T and S) increments
   !!----------------------------------------------------------------------

   USE oce             ! Dynamics and active tracers defined in memory
   USE par_oce         ! Ocean space and time domain variables
   USE dom_oce         ! Ocean space and time domain
   !
   USE in_out_manager  ! I/O manager
   USE iom             ! Library to read input files
   USE fldread        ! read input fields
   USE lib_mpp         ! MPP library
   !
#if defined key_si3
   USE ice
#endif
   !
   IMPLICIT NONE
   PRIVATE

   PUBLIC TRA_MEFT

   LOGICAL, PUBLIC, SAVE :: ln_meft   = .TRUE.
   LOGICAL, SAVE :: ln_icemask
   REAL(wp), SAVE :: zincwgt
   REAL(wp), ALLOCATABLE, SAVE :: zicems(:,:)

   LOGICAL, SAVE :: ln_online = .FALSE.
   INTEGER :: nn_online, nn_features, nn_outputs, nn_shiter, nn_meft
   REAL(wp), ALLOCATABLE, SAVE :: wrk_array(:,:,:)

   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_bct
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_bcs

CONTAINS

   SUBROUTINE tra_meft_init
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE tra_meft_init  ***
      !!
      !! ** Purpose : Initialize the model error forcing term
      !!
      !! ** Action  :
      !!----------------------------------------------------------------------
      INTEGER :: inum, ios, ierror
      INTEGER :: nn_days
      LOGICAL :: ln_trameft
      !
      TYPE(FLD_N) ::   sn_bct, sn_bcs
      CHARACTER(len=100) ::  cn_dir
      !
      NAMELIST/nam_trameft/ ln_trameft, nn_days, ln_icemask, &
                          & ln_online, nn_online, nn_shiter, &
                          & sn_bct, sn_bcs, nn_meft, cn_dir
      !-----------------------------------------------------------------------
      ln_trameft     = .FALSE.
      ln_online      = .FALSE.
      ln_icemask     = .FALSE.
      nn_online      = 0
      nn_shiter      = 0
      nn_meft        = 1

      REWIND( numnam_ref )              ! Namelist nam_trameft in reference namelist : Assimilation increment
      READ  ( numnam_ref, nam_trameft, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'nam_trameft in reference namelist' )
      REWIND( numnam_cfg )              ! Namelist nam_trameft in configuration namelist : Assimilation increment
      READ  ( numnam_cfg, nam_trameft, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'nam_trameft in configuration namelist' )
      IF(lwm) WRITE ( numond, nam_trameft )

      ! Control print
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'tra_meft_init : Model error forcing term'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist nam_trameft'
         WRITE(numout,*) '      Logical switch for applying model error forcing term    ln_trameft = ', ln_trameft
         WRITE(numout,*) '      Timescale of the application (in days)                  nn_days    = ', nn_days
         WRITE(numout,*) '      Frequency of reading files                              nn_meft    = ', nn_meft
         WRITE(numout,*) '      Online calculation                                      ln_online  = ', ln_online
         WRITE(numout,*) '      Type of Online calculation                              nn_online  = ', nn_online
         WRITE(numout,*) '      Number of Shapiro filter iterations                     nn_shiter  = ', nn_shiter
         WRITE(numout,*)
       ENDIF

       ln_meft = ln_trameft

       IF(.NOT.ln_trameft) RETURN

       zincwgt = 1._wp / (86400.* nn_days)

       ALLOCATE( zicems(jpi,jpj) )      ; zicems(:,:)   = 1._wp

       IF( .NOT. ln_online ) THEN

            ALLOCATE( sf_bct(1), STAT=ierror )
            IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'trameft: unable to allocate sf_bct structure' )
            ALLOCATE( sf_bct(1)%fnow(jpi,jpj,jpk), STAT=ierror )
            IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'trameft: unable to allocate sf_bct now array' )
            CALL fld_fill( sf_bct, (/ sn_bct /), cn_dir, 'trameft', 'TRAMEFT temperature term', 'nam_trameft', no_print )
            IF( sf_bct(1)%ln_tint )   ALLOCATE( sf_bct(1)%fdta(jpi,jpj,jpk,2), STAT=ierror )
            IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'trameft: unable to allocate sf_bct data array' )
     
            ALLOCATE( sf_bcs(1), STAT=ierror )
            IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'trameft: unable to allocate sf_bcs structure' )
            ALLOCATE( sf_bcs(1)%fnow(jpi,jpj,jpk), STAT=ierror )
            IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'trameft: unable to allocate sf_bcs now array' )
            CALL fld_fill( sf_bcs, (/ sn_bcs /), cn_dir, 'trameft', 'TRAMEFT temperature term', 'nam_trameft', no_print )
            IF( sf_bcs(1)%ln_tint )   ALLOCATE( sf_bcs(1)%fdta(jpi,jpj,jpk,2), STAT=ierror )
            IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'trameft: unable to allocate sf_bcs data array' )

       ELSE
            CALL ctl_stop('trameft : ln_online not supported yet')
            CALL online_init
       ENDIF

   END SUBROUTINE tra_meft_init

   SUBROUTINE tra_meft(kt)
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE tra_meft  ***
      !!
      !! ** Purpose : Apply the model error forcing term
      !!
      !! ** Method  : Incremental Analysis Updating on tendencies
      !!
      !! ** Action  :
      !!----------------------------------------------------------------------
      INTEGER, INTENT(IN) ::   kt   ! Current time step
      !
      INTEGER  :: ji, jj, jk
      INTEGER  :: it
      REAL(wp), ALLOCATABLE :: zzt(:,:,:)
      !
      IF( kt .eq. nit000 ) THEN
            CALL tra_meft_init
            IF(.NOT.ln_meft) RETURN
      ENDIF
      !
      IF( MOD( kt-1, nn_meft ) == 0 ) THEN
          CALL fld_read( kt, nn_meft, sf_bct )
          CALL fld_read( kt, nn_meft, sf_bcs )
          IF( nn_shiter > 0 ) THEN
                IF(lwp) WRITE(numout,*)  'Applying filtering, iterations:',nn_shiter
                ALLOCATE ( zzt(jpi,jpj,jpk) )
                zzt (:,:,:)       = sf_bct(1)%fnow(:,:,:)
                CALL Shapiro_3D( zzt, nn_shiter, sf_bct(1)%fnow(:,:,:) )
                zzt (:,:,:)       = sf_bcs(1)%fnow(:,:,:)
                CALL Shapiro_3D( zzt, nn_shiter, sf_bcs(1)%fnow(:,:,:) )
                DEALLOCATE ( zzt )
          ENDIF
      ENDIF
      !
      zicems(:,:) = 1._wp
#if defined key_si3
      IF( ln_icemask ) THEN
              zicems(:,:) = 1._wp - SUM( a_i , dim=3 )
      ENDIF
#endif
      !
      ! Update the tracer tendencies
      !
      DO jk = 1, jpkm1
               tsa(:,:,jk,jp_tem) = tsa(:,:,jk,jp_tem) + sf_bct(1)%fnow(:,:,jk) * zicems(:,:) * zincwgt
               tsa(:,:,jk,jp_sal) = tsa(:,:,jk,jp_sal) + sf_bcs(1)%fnow(:,:,jk) * zicems(:,:) * zincwgt
      ENDDO
      !
   END SUBROUTINE tra_meft

   SUBROUTINE online_init
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE online_init  ***
      !!
      !! ** Purpose : Initialize in online mode
      !!
      !! ** Action  :
      !!----------------------------------------------------------------------

      SELECT CASE(nn_online)
                CASE(1)
                        nn_features = 177
                        nn_outputs  = jpkm1*2
                CASE DEFAULT
                        CALL ctl_stop('trameft : nn_online not supported')
      END SELECT

      ALLOCATE ( wrk_array(jpi,jpj,nn_features) )

   END SUBROUTINE online_init

#ifdef ONLINE_INFERENCE
   SUBROUTINE collect(kt)
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE collect  ***
      !!
      !! ** Purpose : Collect predictors
      !!
      !! ** Action  :
      !!----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: kt
      INTEGER :: kp
      REAL(wp) :: jul

      ! 3D fields
      wrk_array(:,:,1:jpkm1)             = tsn(:,:,1:jpkm1,jp_tem)
      wrk_array(:,:,(jpkm1+1):(2*jpkm1)) = tsn(:,:,1:jpkm1,jp_sal)
      kp=2*jpkm1

      ! 2D fields
      wrk_array(:,:,kp+1)  =  sshb(:,)
      wrk_array(:,:,kp+2)  =  hmlp(:,)
      wrk_array(:,:,kp+3)  =  emp - rnf
      wrk_array(:,:,kp+4)  =  qsr
      wrk_array(:,:,kp+5)  =  qns
      wrk_array(:,:,kp+6)  =  qns+qsr
      wrk_array(:,:,kp+7)  =  qns-zqlw+zqsb+zqla
      wrk_array(:,:,kp+8)  =  sfx
      wrk_array(:,:,kp+9)  =  taum
      wrk_array(:,:,kp+10) =  ub(:,:,1)
      wrk_array(:,:,kp+11) =  utau
      wrk_array(:,:,kp+12) =  0.5*rcp * z2d 
      wrk_array(:,:,kp+13) =  0.5 * z2d 
      wrk_array(:,:,kp+14) =  vtau
      wrk_array(:,:,kp+15) =  0.5*rcp * z2d 
      wrk_array(:,:,kp+16) =  0.5 * z2d 
      wrk_array(:,:,kp+17) =  sshb/e1t
      wrk_array(:,:,kp+18) =  sshb/e2t
      wrk_array(:,:,kp+19) =  vt_s * zmsksn
      wrk_array(:,:,kp+20) =  hm_i * zmsk00
      wrk_array(:,:,kp+21) =  at_i * zmsk00
      kp=kp+21

      ! Static fields
      deg2rad=2._wp*rpi/360._wp
      jul2rad=2._wp*rpi/365.25_wp
      jul = REAL(nday_year, wp)
      wrk_array(:,:,kp+1)  =  k_bot
      wrk_array(:,:,kp+2)  =  distc
      wrk_array(:,:,kp+3)  =  gphit
      wrk_array(:,:,kp+4)  =  sin(glamt*deg2rad)
      wrk_array(:,:,kp+5)  =  cos(glamt*deg2rad)
      wrk_array(:,:,kp+6)  =  sin(jul*jul2rad)
      wrk_array(:,:,kp+7)  =  cos(jul*jul2rad)
      kp=kp+7

      IF( kp .ne. nn_features ) CALL ctl_stop('trameft : inconsistent numb of features')

   END SUBROUTINE collect
#endif

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
            call lbc_lnk('trameft_shp', rlpa_varout, 'T', 1.) ! Boundary condition
            rlpa_varout_tmp(:,:,:) = rlpa_varout(:,:,:)
         ENDDO  ! end loop jn
!
END SUBROUTINE Shapiro_3D

END MODULE TRAMEFT
