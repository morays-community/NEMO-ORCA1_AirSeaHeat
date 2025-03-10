MODULE dnnqcorr
#if defined key_annif
   !!======================================================================
   !!                       ***  MODULE dnnqcorr ***
   !!        Surface model error forcing term (estimated online with ANN)
   !! 
   !!======================================================================
   !! History :       ! 2024-04  (A. Storto) Initial version (NEMO 4.0.7)
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dnn_qcorr_init : Initialize the model errors
   !!   dnn_qcorr      : Apply the qrp increments
   !!----------------------------------------------------------------------

   USE oce             ! Dynamics and active tracers defined in memory
   USE par_oce         ! Ocean space and time domain variables
   USE dom_oce         ! Ocean space and time domain
   !
   USE in_out_manager  ! I/O manager
   USE iom             ! Library to read input files
   USE fldread         ! read input fields
   USE lib_mpp         ! MPP library
   USE zdfmxl
   USE sbc_oce
   USE sbc_ice
   USE sbcssr
   USE sbcblk
   USE annif
   USE eosbn2          ! equation of state                (eos_bn2 routine)
   USE nemo2infero
   USE infmod
   !
#if defined key_si3
   USE ice
#endif

   IMPLICIT NONE
   PRIVATE

   INTEGER,PARAMETER :: jpfea = 24
   LOGICAL,  SAVE :: ln_spinup = .FALSE.
   LOGICAL,  SAVE :: ln_infero = .FALSE.
   INTEGER,  SAVE :: np, ndays, kdays, kdd, nn_fdnq, nn_secave, ktt, kd1
   INTEGER,  SAVE :: knf
   REAL(wp), SAVE :: zdvd, zdv1, rn_qampl
   LOGICAL        :: ln_dnnqcorr = .true.
   LOGICAL        :: ln_subdaily = .false.
   LOGICAL        :: ln_climcorr = .false.
   LOGICAL        :: ln_eophis   = .false.

   REAL(wp), ALLOCATABLE :: zifld(:,:,:), zimea(:,:,:), ziuse(:,:,:,:), zofld(:,:)
   REAL(wp), ALLOCATABLE :: zi1d (:,:)  , zo1d(:,:)
   REAL(wp), ALLOCATABLE :: mxld (:,:)

   TYPE(ANNIF_TYPE), ALLOCATABLE :: NMOD(:)

   PUBLIC :: dnn_qcorr
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_qcl
#if defined key_infero
   TYPE(infero_instance) :: dnninfmod
#endif

   CONTAINS

   SUBROUTINE dnn_qcorr_init

      IMPLICIT NONE
      INTEGER  :: ios, kvrb, ierror, ic
      CHARACTER(LEN=299) cfnmod
      TYPE(FLD_N) ::   sn_qcl
      CHARACTER(LEN=128) :: input_layer_name, output_layer_name

      NAMELIST/namdnn_qrp/ln_dnnqcorr, ndays, cfnmod, nn_fdnq, &
      & nn_secave, rn_qampl, ln_subdaily, ln_climcorr, ln_infero, &
      & input_layer_name, output_layer_name, ln_eophis

      input_layer_name = 'serving_default_dense_input'
      output_layer_name = 'StatefulPartitionedCall'

      ndays = 30
      nn_fdnq = nn_fsbc
      nn_secave = 86400
      cfnmod    = 'model.nc'
      rn_qampl  = 1._wp

      REWIND( numnam_ref )
      READ  ( numnam_ref, namdnn_qrp, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namdnn_qrp in reference namelist')

      REWIND( numnam_cfg )
      READ  ( numnam_cfg, namdnn_qrp, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namdnn_qrp in configuration namelist')
      IF(lwm) WRITE ( numond, namdnn_qrp )

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'init_stopack : Stochastic physics package'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namdnn_qrp:'
         WRITE(numout,*)
         WRITE(numout,*) '       Switch on dnn_qrp                                   ln_dnnqcorr : ',ln_dnnqcorr
         WRITE(numout,*) '       Number of days for calculating mean                       ndays :', ndays
         WRITE(numout,*) '       Frequency of qns update                                 nn_fdnq :', nn_fdnq
         WRITE(numout,*) '       Period for averaging (in seconds)                     nn_secave :', nn_secave
         WRITE(numout,*) '       Amplitude of correction                                rn_qampl :', rn_qampl
         WRITE(numout,*) '       Neural network file                                      cfnmod :', cfnmod
         WRITE(numout,*) '       Subdaily update                                     ln_subdaily :', ln_subdaily
         WRITE(numout,*) '       Climatological correction                           ln_climcorr :', ln_climcorr
         WRITE(numout,*) '       Use INFERO for inference instead of ANNIF             ln_infero :', ln_infero
         WRITE(numout,*) '       Eophis for OASIS-coupled inference instead of ANNIF   ln_eophis :', ln_eophis
         WRITE(numout,*) '       TF Layer names for Infero                           layer_names :', &
         & TRIM(input_layer_name), TRIM(output_layer_name)
      ENDIF

      IF(.NOT. ln_dnnqcorr ) RETURN

      IF( nn_sstr .GT. 0 .AND. ln_ssr ) CALL ctl_stop( 'STOP', 'dnnqcorr : cannot work with sstr' )

      IF( ln_climcorr ) THEN
             sn_qcl = FLD_N( 'qrp_climcorr', -1. ,  'qrp' , .true. ,.true., 'yearly' , '', '', '')
             ALLOCATE( sf_qcl(1), STAT=ierror )
             IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'dnn_qrp: unable to allocate sf_qcl structure' )
             ALLOCATE( sf_qcl(1)%fnow(jpi,jpj,1), STAT=ierror )
             IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'dnn_qrp: unable to allocate sf_qcl now array' )
             CALL fld_fill( sf_qcl, (/ sn_qcl /), './', 'dnn_qrp', 'qrp correction', 'namdnn', no_print )
             IF( sf_qcl(1)%ln_tint )   ALLOCATE( sf_qcl(1)%fdta(jpi,jpj,1,2), STAT=ierror )
             IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'dnn_qrp: unable to allocate sf_qcl data' )

             RETURN

      ENDIF

      kdays = 0
      knf = 0
      kdd = NINT(REAL(nn_secave,wp)/rdt) 
      kd1 = NINT(86400._wp / REAL(nn_secave,wp)) 
      ktt = kd1 * ndays
      zdv1 = 1._wp/REAL(kdd,wp)
      zdvd = 1._wp/REAL(ktt,wp)

      ALLOCATE ( zifld(jpi,jpj,jpfea+2) )
      ALLOCATE ( zimea(jpi,jpj,jpfea) )
      ALLOCATE ( ziuse(jpi,jpj,jpfea,ktt) )
      ALLOCATE ( zofld(jpi,jpj       ) )
      ALLOCATE ( mxld (jpi,jpj       ) )

      zofld(:,:) = 0._wp

      call ini_zifld

      np = NINT ( SUM ( tmask(:,:,1) ) )
      ALLOCATE( zi1d(np,jpfea) )
      ALLOCATE( zo1d(np,    1) )

      IF(lwp) WRITE(numout,*) '       Number of timesteps for averaging                           kdd :', kdd
      IF(lwp) WRITE(numout,*) '       Number of timesteps per day                                 kd1 :', kd1
      IF(lwp) WRITE(numout,*) '       Number of total fields to store                             ktt :', ktt
      IF(lwp) WRITE(numout,*) '       Number of sea points in this domain                          np :', np
      IF(lwp) WRITE(numout,*) '       Divider (kdd)                                              zdv1 :', zdv1,1./zdv1
      IF(lwp) WRITE(numout,*) '       Divider (ktt)                                              zdvd :', zdvd,1./zdvd

      ! We will have only one model in this test
      ALLOCATE ( NMOD(1) )

      ! Initialization
      kvrb = 0
      IF(lwp) kvrb = 10
      CALL ANNIF_INIT ( NMOD, kvrb, numout, numout)

      ! Reading
      CALL ANNIF_RDMOD( NMOD(1), TRIM(cfnmod) )

#if defined key_infero
      ! Initialize Infero as needed
      IF(ln_infero) THEN
            ! Assuming cfnmod is a NetCDF file
            ic = LEN_TRIM ( cfnmod ) - 3
            IF(lwp) WRITE(numout,*) ' Initializing Infero model :',TRIM(cfnmod(1:ic))
            CALL nemo2infero_init(dnninfmod,1,TRIM(cfnmod(1:ic)),0,TRIM(input_layer_name),&
            & TRIM(output_layer_name),numout,lwp)
      ENDIF
#endif

   END SUBROUTINE dnn_qcorr_init

   SUBROUTINE ini_zifld

      IMPLICIT NONE
      ! (Re-)Initialize the flds
      zifld(:,:,:)  = 0._wp
      zifld(:,:,6)  = -1.E+20_wp
      zifld(:,:,18) = -1.E+20_wp
      zifld(:,:,25) = +1.E+20_wp
      zifld(:,:,26) = +1.E+20_wp

   END SUBROUTINE ini_zifld

   SUBROUTINE dnn_qcorr(kt)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: kt
      REAL(wp) :: ztmp

      IF( kt .eq. nit000 ) CALL dnn_qcorr_init
      IF(.NOT. ln_dnnqcorr ) RETURN

      IF( ln_climcorr ) THEN
              
         CALL fld_read( kt, nn_fsbc, sf_qcl )   ! Read data
         qns(:,:) = qns(:,:) + sf_qcl(1)%fnow(:,:,1)*tmask(:,:,1)*rn_qampl

      ELSE

         IF( MOD(kt-1,kdd) .eq. 0 ) THEN

              IF (kt .gt. nit000 ) THEN

                      IF(lwp) WRITE(numout,*) '  dnn_qcorr : cycle at (kt,kd):', kt, kdays
                      call dtfix(kdays)
                      call dtmea
                      IF ( kdays .eq. ktt ) ln_spinup = .TRUE.
                      IF ( ln_spinup ) call dnn_qrp(kt)

              ENDIF

              kdays=kdays+1
              IF(kdays .gt. ktt) THEN
                      kdays = 1
              ENDIF

         ENDIF

         CALL collect(kt)

         IF( ln_spinup .AND. ( MOD( kt-1, nn_fdnq ) == 0 .OR. nn_fdnq .le. 1 ) ) THEN
           qns(:,:) = qns(:,:) + zofld(:,:)*tmask(:,:,1)*rn_qampl
           call iom_put("dnn_qrp",zofld(:,:)*tmask(:,:,1)*rn_qampl)
         ENDIF

      ENDIF

   END SUBROUTINE dnn_qcorr

   SUBROUTINE dtfix(kd)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: kd
      INTEGER             :: ji,jj,jk

      DO jk= 3, jpfea
         ziuse(:,:,jk,kd) = zifld(:,:,jk)*zdv1
      ENDDO

      IF(.NOT. ln_subdaily) THEN
        DO jj=1,jpj
           DO ji=1,jpi
              ziuse(ji,jj,6, kd) = zifld(ji,jj,6) -zifld(ji,jj,25)
              ziuse(ji,jj,18,kd) = zifld(ji,jj,18)-zifld(ji,jj,26)
           ENDDO
        ENDDO
      ENDIF

      call ini_zifld

   END SUBROUTINE dtfix

   SUBROUTINE dtmea

      IMPLICIT NONE
      INTEGER  :: jd, jk, ji, jj, ntot
      REAL(wp) :: zc

      IF( ln_subdaily ) THEN
              zc = 1._wp
              ntot = kd1
      ELSE
              zc = 1._wp
              ntot = ndays
      ENDIF

      zimea(:,:,:) = 0._wp
      DO jk= 3, jpfea
         DO jd=1,ntot
            zimea(:,:,jk) = zimea(:,:,jk) + ziuse(:,:,jk,jd) * zdvd * zc
         ENDDO
      ENDDO
      zimea(:,:,1) = glamt   (:,:)
      zimea(:,:,2) = gphit   (:,:)

      IF(ln_subdaily) THEN
        DO jj=1,jpj
           DO ji=1,jpi
              zimea(ji,jj,6 ) = MAXVAL ( ziuse(ji,jj, 4,:) ) - MINVAL( ziuse(ji,jj, 4,:) )
              zimea(ji,jj,18) = MAXVAL ( ziuse(ji,jj,20,:) ) - MINVAL( ziuse(ji,jj,20,:) )
           ENDDO
        ENDDO
      ENDIF

   END SUBROUTINE dtmea

   SUBROUTINE collect(kt)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: kt

      INTEGER  :: ji,jj,jk
      REAL(wp), DIMENSION(jpi,jpj) :: z2d, z2f

      z2d(:,:)  = 0._wp
      z2f(:,:)  = 0._wp
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               z2d(ji,jj) = z2d(ji,jj) + e3t_n(ji,jj,jk) * tsn(ji,jj,jk,jp_tem) * tmask(ji,jj,jk)
               z2f(ji,jj) = z2f(ji,jj) + e3t_n(ji,jj,jk) * tsn(ji,jj,jk,jp_sal) * tmask(ji,jj,jk)
            END DO
         END DO
      END DO

      CALL mx_ld(kt)

      zifld(:,:,3) = zifld(:,:,3) + rau0_rcp * z2d(:,:)
      zifld(:,:,4) = zifld(:,:,4) + mxld    (:,:)
      zifld(:,:,5) = zifld(:,:,5) + tprecip (:,:)
      WHERE ( mxld (:,:) > zifld(:,:,6) ) zifld(:,:,6) = mxld (:,:)
      zifld(:,:,7) = zifld(:,:,7) + qemp_oce(:,:)
      zifld(:,:,8) = zifld(:,:,8) + blk_aux(:,:,3)
      zifld(:,:,9) = zifld(:,:,9) + blk_aux(:,:,1)
      zifld(:,:,10)= zifld(:,:,10) + qns_oce (:,:)
      zifld(:,:,11)= zifld(:,:,11) + blk_aux(:,:,2)
      zifld(:,:,12)= zifld(:,:,12) + qsr_oce * ( 1._wp - at_i_b )
      zifld(:,:,13)= zifld(:,:,13) + ( qsr_oce + qns_oce ) * ( 1._wp - at_i_b ) + qemp_oce
      zifld(:,:,14)= zifld(:,:,14) + rnf     (:,:)
      zifld(:,:,15)= zifld(:,:,15) + rau0 * z2f(:,:)
      zifld(:,:,16)= zifld(:,:,16) + sfx     (:,:)
      zifld(:,:,17)= zifld(:,:,17) + tsn(:,:,1,jp_sal)
      WHERE ( tsn(:,:,1,jp_tem) > zifld(:,:,18) ) zifld(:,:,18) = tsn(:,:,1,jp_tem)
      zifld(:,:,19)= zifld(:,:,19) + taum    (:,:)
      zifld(:,:,20)= zifld(:,:,20) + tsn(:,:,1,jp_tem)
      zifld(:,:,21)= zifld(:,:,21) + emp(:,:) - rnf(:,:) 
      zifld(:,:,22)= zifld(:,:,22) + wndm    (:,:)
      zifld(:,:,23)= zifld(:,:,23) + sshn    (:,:)
      zifld(:,:,24)= zifld(:,:,24) + REAL( nday_year, wp )
      WHERE ( mxld (:,:) < zifld(:,:,25) ) zifld(:,:,25) = mxld (:,:)
      WHERE ( tsn(:,:,1,jp_tem) < zifld(:,:,26) ) zifld(:,:,26) = tsn(:,:,1,jp_tem)

   END SUBROUTINE collect

   SUBROUTINE dnn_qrp(kt)

      IMPLICIT NONE
      INTEGER,INTENT(IN)  :: kt
      INTEGER  :: ji,jj,kp,jf,sav
      LOGICAL  :: llv
      REAL(wp) :: ztmp1(1),ztmp2(1)

      knf=knf+1
      IF(lwp) WRITE(numout,*) '  dnn_qcorr : inference (kt, knf):',kt,knf

      llv=.false.
      kp=0
      DO jj=1,jpj
         DO ji=1,jpi
            IF( tmask(ji,jj,1) .GT. 0.5_wp ) THEN
                kp=kp+1
                zi1d(kp,:) = zimea(ji,jj,:)
                CALL ANNIF_NORMALIZ(NMOD(1),zi1d(kp,:),1,1)
            ENDIF
         ENDDO
      ENDDO

      IF(ln_infero) THEN
#if defined key_infero
         CALL nemo2infero_predict(dnninfmod, np, jpfea, zi1d, np, 1, zo1d)
#endif
      ELSE IF (ln_eophis) THEN
         CALL inferences( kt, zi1d, zo1d )
      ELSE
         CALL ANNIF_FWD(NMOD(1), np, zi1d(:,:), zo1d(:,:))
      ENDIF

      llv=.false.
      kp=0
      zofld(:,:) = 0._wp
      DO jj=1,jpj
         DO ji=1,jpi
            IF( tmask(ji,jj,1) .GT. 0.5_wp ) THEN
                kp=kp+1
                CALL ANNIF_NORMALIZ(NMOD(1),zo1d(kp,:),2,2)
                zofld(ji,jj) = zo1d(kp,1)
            ENDIF
         ENDDO
      ENDDO

      ztmp1(1) = minval( zo1d(:,1))
      ztmp2(1) = maxval( zo1d(:,1))
      IF(lwp) WRITE(numout,*) '  LOC Min/Max : ',ztmp1(1),ztmp2(1)
      CALL mpp_min( "dnn_qrp", ztmp1(1) )
      CALL mpp_max( "dnn_qrp", ztmp2(1) )
      IF(lwp) WRITE(numout,*) '  GLO Min/Max : ',ztmp1(1),ztmp2(1)

   END SUBROUTINE dnn_qrp

   SUBROUTINE mx_ld(kt)
      !
      INTEGER, INTENT(IN) :: kt
      INTEGER  ::   ji, jj, jk      ! dummy loop indices
      INTEGER  ::   iikn, iiki, ikt ! local integer
      REAL(wp) ::   zN2_c           ! local scalar
      INTEGER, DIMENSION(jpi,jpj) ::   nmln   ! 2D workspace
      REAL(wp) ::   rho_c = 0.01_wp    !: density criterion for mixed layer depth
      !!----------------------------------------------------------------------
      !
      ! w-level of the mixing and mixed layers
      IF( kt .eq. nit000) THEN
          CALL eos_rab( tsb, rab_b )       ! before local thermal/haline expension ratio at T-points
          CALL bn2    ( tsb, rab_b, rn2b ) ! before Brunt-Vaisala frequency
      ENDIF
      nmln(:,:)  = nlb10               ! Initialization to the number of w ocean point
      mxld(:,:)  = 0._wp               ! here hmlp used as a dummy variable, integrating vertically N^2
      zN2_c = grav * rho_c * r1_rau0   ! convert density criteria into N^2 criteria
      DO jk = nlb10, jpkm1
         DO jj = 1, jpj                ! Mixed layer level: w-level 
            DO ji = 1, jpi
               ikt = mbkt(ji,jj)
               mxld(ji,jj) = mxld(ji,jj) + MAX( rn2b(ji,jj,jk) , 0._wp ) * e3w_n(ji,jj,jk)
               IF( mxld(ji,jj) < zN2_c )   nmln(ji,jj) = MIN( jk , ikt ) + 1   ! Mixed layer level
            END DO
         END DO
      END DO
      !
      ! depth of the mixing and mixed layers
      DO jj = 1, jpj
         DO ji = 1, jpi
            iikn = nmln(ji,jj)
            mxld (ji,jj) = gdepw_n(ji,jj,iikn  ) * ssmask(ji,jj)    ! Mixed layer depth
         END DO
      END DO
      !
   END SUBROUTINE mx_ld

#else
   CONTAINS
   SUBROUTINE dnn_qcorr(kt)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: kt

   END SUBROUTINE dnn_qcorr
#endif

END MODULE dnnqcorr
