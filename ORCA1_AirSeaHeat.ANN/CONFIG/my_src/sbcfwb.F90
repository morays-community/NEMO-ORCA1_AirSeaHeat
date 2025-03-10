MODULE sbcfwb
   !!======================================================================
   !!                       ***  MODULE  sbcfwb  ***
   !! Ocean fluxes   : domain averaged freshwater budget
   !!======================================================================
   !! History :  OPA  ! 2001-02  (E. Durand)  Original code
   !!   NEMO     1.0  ! 2002-06  (G. Madec)  F90: Free form and module
   !!            3.0  ! 2006-08  (G. Madec)  Surface module
   !!            3.2  ! 2009-07  (C. Talandier) emp mean s spread over erp area 
   !!            3.6  ! 2014-11  (P. Mathiot  ) add ice shelf melting
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_fwb       : freshwater budget for global ocean configurations (free surface & forced mode)
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE sbc_oce        ! surface ocean boundary condition
   USE sbc_ice , ONLY : snwice_mass, snwice_mass_b, snwice_fmass
   USE phycst         ! physical constants
   USE sbcrnf         ! ocean runoffs
   USE sbcisf         ! ice shelf melting contribution
   USE sbcssr         ! Sea-Surface damping terms
   !
   USE in_out_manager ! I/O manager
   USE iom            ! IOM
   USE lib_mpp        ! distribued memory computing library
   USE timing         ! Timing
   USE lbclnk         ! ocean lateral boundary conditions
   USE lib_fortran    ! 
   USE fldread        ! read input fields

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_fwb    ! routine called by step

   REAL(wp) ::   a_fwb_b   ! annual domain averaged freshwater budget
   REAL(wp) ::   a_fwb     ! for 2 year before (_b) and before year.
   REAL(wp) ::   fwfold    ! fwfold to be suppressed
   REAL(wp) ::   area      ! global mean ocean surface (interior domain)

   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_bsl   ! structure of input SST (file informations, fields read)
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_pmc

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: sbcfwb.F90 13581 2020-10-09 11:49:08Z mathiot $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sbc_fwb( kt, kn_fwb, kn_fsbc )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_fwb  ***
      !!
      !! ** Purpose :   Control the mean sea surface drift
      !!
      !! ** Method  :   several ways  depending on kn_fwb
      !!                =0 no control 
      !!                =1 global mean of emp set to zero at each nn_fsbc time step
      !!                =2 annual global mean corrected from previous year
      !!                =3 global mean of emp set to zero at each nn_fsbc time step
      !!                   & spread out over erp area depending its sign
      !! Note: if sea ice is embedded it is taken into account when computing the budget 
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt       ! ocean time-step index
      INTEGER, INTENT( in ) ::   kn_fsbc  ! 
      INTEGER, INTENT( in ) ::   kn_fwb   ! ocean time-step index
      !
      INTEGER  ::   inum, ikty, iyear     ! local integers
      INTEGER  ::   ierror
      REAL(wp) ::   z_fwf, z_fwf_nsrf, zsum_fwf, zsum_erp                ! local scalars
      REAL(wp) ::   zsurf_neg, zsurf_pos, zsurf_tospread, zcoef          !   -      -
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   ztmsk_neg, ztmsk_pos, z_wgt ! 2D workspaces
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   ztmsk_tospread, zerp_cor    !   -      -
      REAL(wp), DIMENSION(jpi,jpj) :: zdiff, zdiffo
      REAL(wp)   ,DIMENSION(1) ::   z_fwfprv  , ztcorr
      COMPLEX(wp),DIMENSION(1) ::   y_fwfnow  
      TYPE(FLD_N) ::   sn_bsl
      TYPE(FLD_N) ::   sn_pmc
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 ) THEN
         IF(lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'sbc_fwb : FreshWater Budget correction'
            WRITE(numout,*) '~~~~~~~'
            IF( kn_fwb == -3)   WRITE(numout,*) '          instantaneously set to external barystatic sea level with EMP mask (2 steps)'
            IF( kn_fwb == -2)   WRITE(numout,*) '          instantaneously set to external barystatic sea level with EMP mask'
            IF( kn_fwb == -1)   WRITE(numout,*) '          instantaneously set to external barystatic sea level'
            IF( kn_fwb == 1 )   WRITE(numout,*) '          instantaneously set to zero'
            IF( kn_fwb == 2 )   WRITE(numout,*) '          adjusted from previous year budget'
            IF( kn_fwb == 3 )   WRITE(numout,*) '          fwf set to zero and spread out over erp area'
         ENDIF
         !
         IF( kn_fwb == 3 .AND. nn_sssr /= 2 )   CALL ctl_stop( 'sbc_fwb: nn_fwb = 3 requires nn_sssr = 2, we stop ' )
         IF( kn_fwb == 3 .AND. ln_isfcav    )   CALL ctl_stop( 'sbc_fwb: nn_fwb = 3 with ln_isfcav = .TRUE. not working, we stop ' )
         !
         area = glob_sum( 'sbcfwb', e1e2t(:,:) * tmask(:,:,1))           ! interior global domain surface
         ! isf cavities are excluded because it can feedback to the melting with generation of inhibition of plumes
         ! and in case of no melt, it can generate HSSW.
         !
#if ! defined key_si3 && ! defined key_cice
         snwice_mass_b(:,:) = 0.e0               ! no sea-ice model is being used : no snow+ice mass
         snwice_mass  (:,:) = 0.e0
         snwice_fmass (:,:) = 0.e0
#endif
         IF( kn_fwb < 0 ) THEN
             sn_bsl = FLD_N( 'barystatic_sealevel', -1. ,  'emp' , .true. ,.false., 'yearly' , '', '', '')
             ALLOCATE( sf_bsl(1), STAT=ierror )
             IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_fwb: unable to allocate sf_bsl structure' )
             ALLOCATE( sf_bsl(1)%fnow(jpi,jpj,1), STAT=ierror )
             IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_fwb: unable to allocate sf_bsl now array' )
             CALL fld_fill( sf_bsl, (/ sn_bsl /), './', 'sbc_fwb', 'Barystatic sea level', 'namsbc', no_print )
             IF( sf_bsl(1)%ln_tint )   ALLOCATE( sf_bsl(1)%fdta(jpi,jpj,1,2), STAT=ierror )
             IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_fwb: unable to allocate sf_bsl data' )
         ENDIF
         IF( kn_fwb < -1 ) THEN
             sn_pmc = FLD_N( 'pme_clim', -1. ,  'PME0' , .true. ,.true., 'yearly' , '', '', '')
             ALLOCATE( sf_pmc(1), STAT=ierror )
             IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_fwb: unable to allocate sf_pmc structure' )
             ALLOCATE( sf_pmc(1)%fnow(jpi,jpj,1), STAT=ierror )
             IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_fwb: unable to allocate sf_pmc now array' )
             CALL fld_fill( sf_pmc, (/ sn_pmc/), './', 'sbc_fwb', 'Climatological EMP', 'namsbc', no_print )
             IF( sf_pmc(1)%ln_tint )   ALLOCATE( sf_pmc(1)%fdta(jpi,jpj,1,2), STAT=ierror )
             IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_fwb: unable to allocate sf_pmc data' )
         ENDIF
          !
      ENDIF

      SELECT CASE ( kn_fwb )
      !
      CASE ( -3 )                             !==  global mean fwf set to barystatic data  ==!
         !
         IF( MOD( kt-1, kn_fsbc ) == 0 ) THEN
            CALL fld_read( kt, nn_fsbc, sf_bsl )   ! Read data
            CALL fld_read( kt, nn_fsbc, sf_pmc )   ! Read data

            ! step 1
            z_fwfprv(1) = glob_sum( 'sbcfwb', e1e2t(:,:) * ( emp(:,:) - rnf(:,:) + fwfisf(:,:) - snwice_fmass(:,:) ) ) / area
            zdiff(:,:) = e1e2t(:,:) * ( emp(:,:) - rnf(:,:) + fwfisf(:,:) - snwice_fmass(:,:) + sf_pmc(1)%fnow(:,:,1) ) 
            CALL shapiro_1D( zdiff, 360, zdiffo )
            zdiff(:,:) = zdiffo(:,:) / ( glob_sum( 'sbcfwb', zdiffo(:,:)*e1e2t(:,:) * tmask(:,:,1) ) / area )

            ztcorr(1) = sf_bsl(1)%fnow(jpi/2,jpj/2,1)

            emp(:,:) = emp(:,:) -(z_fwfprv(1)*zdiff(:,:))*tmask(:,:,1)
            qns(:,:) = qns(:,:) +(z_fwfprv(1)*zdiff(:,:)) * rcp * sst_m(:,:) * tmask(:,:,1) ! account for change to the heat

            ! step 2
            z_fwfprv(1) = glob_sum( 'sbcfwb', e1e2t(:,:) * ( emp(:,:) - rnf(:,:) + fwfisf(:,:) - snwice_fmass(:,:) ) ) / area
            ztcorr(1) = sf_bsl(1)%fnow(jpi/2,jpj/2,1)
            zcoef = (z_fwfprv(1)-ztcorr(1)) * rcp
            emp(:,:) = emp(:,:) - z_fwfprv(1)        * tmask(:,:,1) + ztcorr(1)*tmask(:,:,1)
            qns(:,:) = qns(:,:) + zcoef * sst_m(:,:) * tmask(:,:,1) ! account for change to the heat budget due to fw correction


            IF(lwp) WRITE(numout,'(X,A,I9,A,2E12.3)') ' Corr at kt = ',kt,'  >>>  ',z_fwfprv(1),ztcorr
            IF( iom_use('fwbcorr') ) CALL iom_put( 'fwbcorr', (ztcorr(1)-z_fwfprv(1))*tmask(:,:,1)*zdiff(:,:) )

         ENDIF
      !
      CASE ( -2 )                             !==  global mean fwf set to barystatic data  ==!
         !
         IF( MOD( kt-1, kn_fsbc ) == 0 ) THEN
            CALL fld_read( kt, nn_fsbc, sf_bsl )   ! Read data
            CALL fld_read( kt, nn_fsbc, sf_pmc )   ! Read data

            z_fwfprv(1) = glob_sum( 'sbcfwb', e1e2t(:,:) * ( emp(:,:) - rnf(:,:) + fwfisf(:,:) - snwice_fmass(:,:) ) ) / area
            zdiff(:,:) = e1e2t(:,:) * ( emp(:,:) - rnf(:,:) + fwfisf(:,:) - snwice_fmass(:,:) + sf_pmc(1)%fnow(:,:,1) ) 
            CALL shapiro_1D( zdiff, 360, zdiffo )
            zdiff(:,:) = zdiffo(:,:) / ( glob_sum( 'sbcfwb', zdiffo(:,:)*e1e2t(:,:) * tmask(:,:,1) ) / area )

            ztcorr(1) = sf_bsl(1)%fnow(jpi/2,jpj/2,1)

            emp(:,:) = emp(:,:) -(z_fwfprv(1)*zdiff(:,:)+ztcorr(1))*tmask(:,:,1)
            qns(:,:) = qns(:,:) +(z_fwfprv(1)*zdiff(:,:)+ztcorr(1)) * rcp * sst_m(:,:) * tmask(:,:,1) ! account for change to the heat
            IF(lwp) WRITE(numout,'(X,A,I9,A,2E12.3)') ' Corr at kt = ',kt,'  >>>  ',z_fwfprv(1),ztcorr
            IF( iom_use('fwbcorr') ) CALL iom_put( 'fwbcorr', (ztcorr(1)-z_fwfprv(1))*tmask(:,:,1)*zdiff(:,:) )

         ENDIF
         !
       CASE ( -1 )                             !==  global mean fwf set to barystatic data  ==!
         !
         IF( MOD( kt-1, kn_fsbc ) == 0 ) THEN
            CALL fld_read( kt, nn_fsbc, sf_bsl )   ! Read data
            z_fwfprv(1) = glob_sum( 'sbcfwb', e1e2t(:,:) * ( emp(:,:) - rnf(:,:) + fwfisf(:,:) - snwice_fmass(:,:) ) ) / area
            ztcorr(1) = sf_bsl(1)%fnow(jpi/2,jpj/2,1)
            zcoef = (z_fwfprv(1)-ztcorr(1)) * rcp
            emp(:,:) = emp(:,:) - z_fwfprv(1)        * tmask(:,:,1) + ztcorr(1)*tmask(:,:,1)
            qns(:,:) = qns(:,:) + zcoef * sst_m(:,:) * tmask(:,:,1) ! account for change to the heat budget due to fw correction
            IF(lwp) WRITE(numout,'(X,A,I9,A,2E12.3)') ' Corr at kt = ',kt,'  >>>  ',z_fwfprv(1),ztcorr
            IF( iom_use('fwbcorr') ) CALL iom_put( 'fwbcorr', (ztcorr(1)-z_fwfprv(1))*tmask(:,:,1) )
         ENDIF
         !
       CASE ( 1 )                             !==  global mean fwf set to zero  ==!
         !
         IF( MOD( kt-1, kn_fsbc ) == 0 ) THEN
            y_fwfnow(1) = local_sum( e1e2t(:,:) * ( emp(:,:) - rnf(:,:) + fwfisf(:,:) - snwice_fmass(:,:) ) )
            CALL mpp_delay_sum( 'sbcfwb', 'fwb', y_fwfnow(:), z_fwfprv(:), kt == nitend - nn_fsbc + 1 )
            z_fwfprv(1) = z_fwfprv(1) / area
            zcoef = z_fwfprv(1) * rcp
            emp(:,:) = emp(:,:) - z_fwfprv(1)        * tmask(:,:,1)
            qns(:,:) = qns(:,:) + zcoef * sst_m(:,:) * tmask(:,:,1) ! account for change to the heat budget due to fw correction
         ENDIF
         !
      CASE ( 2 )                             !==  fwf budget adjusted from the previous year  ==!
         !
         IF( kt == nit000 ) THEN                      ! initialisation
            !                                         ! Read the corrective factor on precipitations (fwfold)
            IF ( ln_rstart .AND. iom_varid( numror, 'a_fwb_b', ldstop = .FALSE. ) > 0     &
               &           .AND. iom_varid( numror, 'a_fwb',   ldstop = .FALSE. ) > 0 ) THEN
               IF(lwp) WRITE(numout,*) 'sbc_fwb : reading FW-budget adjustment from restart file'
               CALL iom_get( numror, 'a_fwb_b', a_fwb_b, ldxios = lrxios )
               CALL iom_get( numror, 'a_fwb',   a_fwb,   ldxios = lrxios )
            ELSE
               CALL ctl_opn( inum, 'EMPave_old.dat'//TRIM(c_enss), 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
               READ ( inum, "(24X,I8,2ES24.16)" ) iyear, a_fwb_b, a_fwb
               CLOSE( inum )
            END IF
            fwfold = a_fwb                            ! current year freshwater budget correction
            !                                         ! estimate from the previous year budget
            IF(lwp)WRITE(numout,*)
            IF(lwp)WRITE(numout,*)'sbc_fwb : year = ',iyear  , ' freshwater budget correction = ', fwfold
            IF(lwp)WRITE(numout,*)'          year = ',iyear-1, ' freshwater budget read       = ', a_fwb
            IF(lwp)WRITE(numout,*)'          year = ',iyear-2, ' freshwater budget read       = ', a_fwb_b
            !
            IF( lwxios ) THEN                         ! Activate output of restart variables
               CALL iom_set_rstw_var_active( 'a_fwb_b' )
               CALL iom_set_rstw_var_active( 'a_fwb'   )
            END IF
         ENDIF   
         !                                         ! Update fwfold if new year start
         ikty = 365 * 86400 / rdt                  !!bug  use of 365 days leap year or 360d year !!!!!!!
         IF( MOD( kt, ikty ) == 0 ) THEN
            a_fwb_b = a_fwb                           ! mean sea level taking into account the ice+snow
                                                      ! sum over the global domain
            a_fwb   = glob_sum( 'sbcfwb', e1e2t(:,:) * ( sshn(:,:) + snwice_mass(:,:) * r1_rau0 ) )
            a_fwb   = a_fwb * 1.e+3 / ( area * rday * 365. )     ! convert in Kg/m3/s = mm/s
!!gm        !                                                      !!bug 365d year 
            fwfold =  a_fwb                           ! current year freshwater budget correction
            !                                         ! estimate from the previous year budget
         ENDIF
         ! 
         IF( MOD( kt-1, kn_fsbc ) == 0 ) THEN         ! correct the freshwater fluxes
            zcoef = fwfold * rcp
            emp(:,:) = emp(:,:) + fwfold             * tmask(:,:,1)
            qns(:,:) = qns(:,:) - zcoef * sst_m(:,:) * tmask(:,:,1) ! account for change to the heat budget due to fw correction
         ENDIF
         ! Output restart information
         IF( lrst_oce ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'sbc_fwb : writing FW-budget adjustment to ocean restart file at it = ', kt
            IF(lwp) WRITE(numout,*) '~~~~'
            IF( lwxios ) CALL iom_swap( cwxios_context )
            CALL iom_rstput( kt, nitrst, numrow, 'a_fwb_b', a_fwb_b, ldxios = lwxios )
            CALL iom_rstput( kt, nitrst, numrow, 'a_fwb',   a_fwb,   ldxios = lwxios )
            IF( lwxios ) CALL iom_swap( cxios_context  )
         END IF
         !
         IF( kt == nitend .AND. lwm ) THEN            ! save fwfold value in a file (only one required)
            CALL ctl_opn( inum, 'EMPave.dat'//TRIM(c_enss), 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE., narea )
            WRITE( inum, "(24X,I8,2ES24.16)" ) nyear, a_fwb_b, a_fwb
            CLOSE( inum )
         ENDIF
         !
      CASE ( 3 )                             !==  global fwf set to zero and spread out over erp area  ==!
         !
         ALLOCATE( ztmsk_neg(jpi,jpj) , ztmsk_pos(jpi,jpj) , ztmsk_tospread(jpi,jpj) , z_wgt(jpi,jpj) , zerp_cor(jpi,jpj) )
         !
         IF( MOD( kt-1, kn_fsbc ) == 0 ) THEN
            ztmsk_pos(:,:) = tmask_i(:,:)                      ! Select <0 and >0 area of erp
            WHERE( erp < 0._wp )   ztmsk_pos = 0._wp
            ztmsk_neg(:,:) = tmask_i(:,:) - ztmsk_pos(:,:)
            !                                                  ! fwf global mean (excluding ocean to ice/snow exchanges) 
            z_fwf     = glob_sum( 'sbcfwb', e1e2t(:,:) * ( emp(:,:) - rnf(:,:) + fwfisf(:,:) - snwice_fmass(:,:) ) ) / area
            !            
            IF( z_fwf < 0._wp ) THEN         ! spread out over >0 erp area to increase evaporation
               zsurf_pos = glob_sum( 'sbcfwb', e1e2t(:,:)*ztmsk_pos(:,:) )
               zsurf_tospread      = zsurf_pos
               ztmsk_tospread(:,:) = ztmsk_pos(:,:)
            ELSE                             ! spread out over <0 erp area to increase precipitation
               zsurf_neg = glob_sum( 'sbcfwb', e1e2t(:,:)*ztmsk_neg(:,:) )  ! Area filled by <0 and >0 erp 
               zsurf_tospread      = zsurf_neg
               ztmsk_tospread(:,:) = ztmsk_neg(:,:)
            ENDIF
            !
            zsum_fwf   = glob_sum( 'sbcfwb', e1e2t(:,:) * z_fwf )         ! fwf global mean over <0 or >0 erp area
!!gm :  zsum_fwf   = z_fwf * area   ???  it is right?  I think so....
            z_fwf_nsrf =  zsum_fwf / ( zsurf_tospread + rsmall )
            !                                                  ! weight to respect erp field 2D structure 
            zsum_erp   = glob_sum( 'sbcfwb', ztmsk_tospread(:,:) * erp(:,:) * e1e2t(:,:) )
            z_wgt(:,:) = ztmsk_tospread(:,:) * erp(:,:) / ( zsum_erp + rsmall )
            !                                                  ! final correction term to apply
            zerp_cor(:,:) = -1. * z_fwf_nsrf * zsurf_tospread * z_wgt(:,:)
            !
!!gm   ===>>>>  lbc_lnk should be useless as all the computation is done over the whole domain !
            CALL lbc_lnk( 'sbcfwb', zerp_cor, 'T', 1. )
            !
            emp(:,:) = emp(:,:) + zerp_cor(:,:)
            qns(:,:) = qns(:,:) - zerp_cor(:,:) * rcp * sst_m(:,:)  ! account for change to the heat budget due to fw correction
            erp(:,:) = erp(:,:) + zerp_cor(:,:)
            !
            IF( nprint == 1 .AND. lwp ) THEN                   ! control print
               IF( z_fwf < 0._wp ) THEN
                  WRITE(numout,*)'   z_fwf < 0'
                  WRITE(numout,*)'   SUM(erp+)     = ', SUM( ztmsk_tospread(:,:)*erp(:,:)*e1e2t(:,:) )*1.e-9,' Sv'
               ELSE
                  WRITE(numout,*)'   z_fwf >= 0'
                  WRITE(numout,*)'   SUM(erp-)     = ', SUM( ztmsk_tospread(:,:)*erp(:,:)*e1e2t(:,:) )*1.e-9,' Sv'
               ENDIF
               WRITE(numout,*)'   SUM(empG)     = ', SUM( z_fwf*e1e2t(:,:) )*1.e-9,' Sv'
               WRITE(numout,*)'   z_fwf         = ', z_fwf      ,' Kg/m2/s'
               WRITE(numout,*)'   z_fwf_nsrf    = ', z_fwf_nsrf ,' Kg/m2/s'
               WRITE(numout,*)'   MIN(zerp_cor) = ', MINVAL(zerp_cor) 
               WRITE(numout,*)'   MAX(zerp_cor) = ', MAXVAL(zerp_cor) 
            ENDIF
         ENDIF
         DEALLOCATE( ztmsk_neg , ztmsk_pos , ztmsk_tospread , z_wgt , zerp_cor )
         !
      CASE DEFAULT                           !==  you should never be there  ==!
         CALL ctl_stop( 'sbc_fwb : wrong nn_fwb value for the FreshWater Budget correction, choose either 1, 2 or 3' )
         !
      END SELECT
      !
   END SUBROUTINE sbc_fwb

   SUBROUTINE Shapiro_1D(rla_varin,id_np, rlpa_varout) !GIG
      IMPLICIT NONE
      INTEGER, INTENT(IN)                       :: id_np
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
                          + 0.25*alphax*(rlpa_varout_tmp(ji-1,jj)-rlpa_varout_tmp(ji,jj))*tmask(ji-1,jj ,1)  &
                          + 0.25*alphax*(rlpa_varout_tmp(ji+1,jj)-rlpa_varout_tmp(ji,jj))*tmask(ji+1,jj ,1)  &
                          + 0.25*alphay*(rlpa_varout_tmp(ji ,jj-1)-rlpa_varout_tmp(ji,jj))*tmask(ji  ,jj-1,1)  &
                          + 0.25*alphay*(rlpa_varout_tmp(ji ,jj+1)-rlpa_varout_tmp(ji,jj))*tmask(ji  ,jj+1,1)
                   rlpa_varout(ji,jj)=znum*tmask(ji,jj,1)+rla_varin(ji,jj)*(1.-tmask(ji,jj,1))
                ENDDO  ! end loop ji
            ENDDO  ! end loop jj
            call lbc_lnk('sbcfwb_shp', rlpa_varout, 'T', 1.) ! Boundary condition
            rlpa_varout_tmp(:,:) = rlpa_varout(:,:)
         ENDDO  ! end loop jn
!
END SUBROUTINE Shapiro_1D
   !!======================================================================
END MODULE sbcfwb
