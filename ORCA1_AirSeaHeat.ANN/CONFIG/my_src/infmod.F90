MODULE infmod
   !!======================================================================
   !!                       ***  MODULE  infmod  ***
   !! Machine Learning Inferences : manage connexion with external ML codes 
   !!======================================================================
   !! History :  4.2.1  ! 2023-09  (A. Barge)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   naminf          : machine learning models formulation namelist
   !!   inferences_init : initialization of Machine Learning based models
   !!   inferences      : ML based models
   !!   inf_snd         : send data to external trained model
   !!   inf_rcv         : receive inferences from external trained model
   !!----------------------------------------------------------------------
   USE oce             ! ocean fields
   USE dom_oce         ! ocean domain fields
   USE cpl_oasis3      ! OASIS3 coupling
   USE timing
   USE iom
   USE in_out_manager
   USE lib_mpp

   IMPLICIT NONE
   PRIVATE

   PUBLIC inf_alloc          ! function called in inferences_init 
   PUBLIC inf_dealloc        ! function called in inferences_final
   PUBLIC inferences_init    ! routine called in nemogcm.F90
   PUBLIC inferences         ! routine called in stpmlf.F90
   PUBLIC inferences_final   ! routine called in nemogcm.F90

   INTEGER, PARAMETER ::   jps_zi1d = 1    ! sea temperature
   INTEGER, PARAMETER ::   jps_tmsk = 2    ! sea salinity
   INTEGER, PARAMETER ::   jps_inf = 2   ! total number of sendings for inferences

   INTEGER, PARAMETER ::   jpr_zo1d = 1   ! density inferences-computed
   INTEGER, PARAMETER ::   jpr_inf = 1   ! total number of inference receptions

   INTEGER, PARAMETER ::   jpinf = MAX(jps_inf,jpr_inf) ! Maximum number of exchanges

   TYPE( DYNARR ), SAVE, DIMENSION(jpinf) ::  infsnd, infrcv  ! sent/received inferences

   !
   !!-------------------------------------------------------------------------
   !!                    Namelist for the Inference Models
   !!-------------------------------------------------------------------------
   !                           !!** naminf namelist **
   !TYPE ::   FLD_INF              !: Field informations ...  
   !   CHARACTER(len = 32) ::         ! 
   !END TYPE FLD_INF
   !
   LOGICAL , PUBLIC ::   ln_inf    !: activate module for inference models
   
   !!-------------------------------------------------------------------------

CONTAINS

   INTEGER FUNCTION inf_alloc()
      !!----------------------------------------------------------------------
      !!             ***  FUNCTION inf_alloc  ***
      !!----------------------------------------------------------------------
      INTEGER :: ierr
      INTEGER :: jn
      !!----------------------------------------------------------------------
      ierr = 0
      !
      DO jn = 1, jpinf
         IF( srcv(ntypinf,jn)%laction ) ALLOCATE( infrcv(jn)%z3(jpi,jpj,srcv(ntypinf,jn)%nlvl), STAT=ierr )
         IF( ssnd(ntypinf,jn)%laction ) ALLOCATE( infsnd(jn)%z3(jpi,jpj,ssnd(ntypinf,jn)%nlvl), STAT=ierr )
         inf_alloc = MAX(ierr,0)
      END DO
      !
   END FUNCTION inf_alloc

   
   INTEGER FUNCTION inf_dealloc()
      !!----------------------------------------------------------------------
      !!             ***  FUNCTION inf_dealloc  ***
      !!----------------------------------------------------------------------
      INTEGER :: ierr
      INTEGER :: jn
      !!----------------------------------------------------------------------
      ierr = 0
      !
      DO jn = 1, jpinf
         IF( srcv(ntypinf,jn)%laction ) DEALLOCATE( infrcv(jn)%z3, STAT=ierr )
         IF( ssnd(ntypinf,jn)%laction ) DEALLOCATE( infsnd(jn)%z3, STAT=ierr )
         inf_dealloc = MAX(ierr,0)
      END DO
      !
   END FUNCTION inf_dealloc


   SUBROUTINE inferences_init 
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE inferences_init  ***
      !!
      !! ** Purpose :   Initialisation of the models that rely on external inferences
      !!
      !! ** Method  :   * Read naminf namelist
      !!                * create data for models
      !!----------------------------------------------------------------------
      !
      INTEGER ::   ios   ! Local Integer
      !!----------------------------------------------------------------------
      !
      ! ======================================== !
      !     Define exchange needs for Models     !
      ! ======================================== !
      !
      ! default definitions of ssnd snd srcv
      srcv(ntypinf,:)%laction = .FALSE.  ;  srcv(ntypinf,:)%clgrid = 'T'  ;  srcv(ntypinf,:)%nsgn = 1.
      srcv(ntypinf,:)%nct = 1  ;  srcv(ntypinf,:)%nlvl = 1
      !
      ssnd(ntypinf,:)%laction = .FALSE.  ;  ssnd(ntypinf,:)%clgrid = 'T'  ;  ssnd(ntypinf,:)%nsgn = 1.
      ssnd(ntypinf,:)%nct = 1  ;  ssnd(ntypinf,:)%nlvl = 1
       
      ! ---------------------------------- !
      !  Python ANN, Storto et al. (2024)  !
      ! ---------------------------------- !

      ! sending of inputs
      ssnd(ntypinf,jps_zi1d)%clname = 'E_OUT_0'
      ssnd(ntypinf,jps_zi1d)%laction = .TRUE.
      ssnd(ntypinf,jps_zi1d)%nlvl = 24

      ssnd(ntypinf,jps_tmsk)%clname = 'E_OUT_1'
      ssnd(ntypinf,jps_tmsk)%laction = .TRUE.

      ! reception of ANN
      srcv(ntypinf,jpr_zo1d)%clname = 'E_IN_0'
      srcv(ntypinf,jpr_zo1d)%laction = .TRUE.

      ! ------------------------------ !

      ! ================================= !
      !   Define variables for coupling
      ! ================================= !
      CALL cpl_var(jpinf, jpinf, 1, ntypinf)
      !
      IF( inf_alloc() /= 0 )     CALL ctl_stop( 'STOP', 'inf_alloc : unable to allocate arrays' )
      !
   END SUBROUTINE inferences_init


   SUBROUTINE inferences( kt, i1d, o1d )
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE inferences  ***
      !!
      !! ** Purpose :   update the ocean data with the ML based models
      !!
      !! ** Method  :   *  
      !!                * 
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt               ! ocean time step
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: i1d
      REAL(wp), DIMENSION(:,:), INTENT(out) :: o1d
      !
      INTEGER :: jj, ji, kp                           ! loop indexes
      INTEGER :: isec, info, jn                       ! local integer
      REAL(wp), DIMENSION(jpi,jpj,jpk)   ::  zdata    ! sending buffer
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('inferences')
      !
      isec = ( kt - nit000 ) * NINT( rdt ) - 86400  ! Date of exchange - remove shift for synchronicity
      IF (lwp) print*, 'Time: ', isec
      info = OASIS_idle
      !
      ! ------  Prepare data to send ------
      !
      ! Inputs
      kp = 0
      zdata = 0.0_wp
      IF( ssnd(ntypinf,jps_zi1d)%laction ) THEN
         DO jj = 1, jpj
            DO ji = 1, jpj
               IF( tmask(ji,jj,1) .GT. 0.5_wp ) THEN
                  kp = kp + 1
                  zdata(ji,jj,1:ssnd(ntypinf,jps_zi1d)%nlvl) = i1d(kp,:)
               ENDIF
            ENDDO
         ENDDO
         infsnd(jps_zi1d)%z3(:,:,1:ssnd(ntypinf,jps_zi1d)%nlvl) = zdata(:,:,1:ssnd(ntypinf,jps_zi1d)%nlvl)
      ENDIF  
      !
      ! t-grid mask
      IF( ssnd(ntypinf,jps_tmsk)%laction ) THEN
         infsnd(jps_tmsk)%z3(:,:,1:ssnd(ntypinf,jps_tmsk)%nlvl) = tmask(:,:,1:ssnd(ntypinf,jps_tmsk)%nlvl)
      ENDIF
      !
      ! ========================
      !   Proceed all sendings
      ! ========================
      !
      DO jn = 1, jpinf
         IF ( ssnd(ntypinf,jn)%laction ) THEN
            CALL cpl_snd( jn, isec, ntypinf, infsnd(jn)%z3, info)
         ENDIF
      END DO
      !
      ! .... some external operations ....
      !
      ! ==========================
      !   Proceed all receptions
      ! ==========================
      !
      DO jn = 1, jpinf
         IF( srcv(ntypinf,jn)%laction ) THEN
            CALL cpl_rcv( jn, isec, ntypinf, infrcv(jn)%z3, info)
         ENDIF
      END DO
      !
      ! ------ Distribute receptions  ------
      !
      ! Outputs
      kp = 0
      o1d = 0.0_wp
      IF( srcv(ntypinf,jpr_zo1d)%laction ) THEN
         DO jj = 1, jpj
            DO ji = 1, jpj
               IF( tmask(ji,jj,1) .GT. 0.5_wp ) THEN
                  kp = kp + 1
                  o1d(kp,1) = infrcv(jpr_zo1d)%z3(jj,ji,1)
               ENDIF
            ENDDO
         ENDDO
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('inferences')
      !
   END SUBROUTINE inferences


   SUBROUTINE inferences_final
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE inferences_final  ***
      !!
      !! ** Purpose :   Free memory used for inferences modules
      !!
      !! ** Method  :   * Deallocate arrays
      !!----------------------------------------------------------------------
      !
      IF( inf_dealloc() /= 0 )     CALL ctl_stop( 'STOP', 'inf_dealloc : unable to free memory' )
      !
   END SUBROUTINE inferences_final 
   !!=======================================================================
END MODULE infmod
