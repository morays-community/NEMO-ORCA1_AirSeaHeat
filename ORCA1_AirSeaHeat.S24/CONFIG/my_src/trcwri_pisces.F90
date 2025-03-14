MODULE trcwri_pisces
   !!======================================================================
   !!                       *** MODULE trcwri ***
   !!    PISCES :   Output of PISCES tracers
   !!======================================================================
   !! History :   1.0  !  2009-05 (C. Ethe)  Original code
   !!----------------------------------------------------------------------
#if defined key_top

   !!----------------------------------------------------------------------
   !! trc_wri_pisces   :  outputs of concentration fields
   !!----------------------------------------------------------------------
   USE trc         ! passive tracers common variables 
   USE sms_pisces  ! PISCES variables
   USE iom         ! I/O manager
   USE ioipsl         !
   USE lib_mpp         ! MPP library
   USE dianam         ! build name of file (routine)
   USE in_out_manager ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_wri_pisces 

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcwri_pisces.F90 10069 2018-08-28 14:12:24Z nicolasmartin $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
#if defined key_iomput 
CONTAINS
   SUBROUTINE trc_wri_pisces
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_wri_trc  ***
      !!
      !! ** Purpose :   output passive tracers fields 
      !!---------------------------------------------------------------------
      CHARACTER (len=20)           :: cltra
      REAL(wp)                     :: zfact
      INTEGER                      :: ji, jj, jk, jn
      REAL(wp), DIMENSION(jpi,jpj) :: zdic, zo2min, zdepo2min
      !!---------------------------------------------------------------------
 
      ! write the tracer concentrations in the file
      ! ---------------------------------------
      IF( ln_p2z ) THEN
         DO jn = jp_pcs0, jp_pcs1
            cltra = TRIM( ctrcnm(jn) )                  ! short title for tracer
            CALL iom_put( cltra, trn(:,:,:,jn) )
         END DO
      ELSE
         DO jn = jp_pcs0, jp_pcs1
            zfact = 1.0e+6 
            IF( jn == jpno3 .OR. jn == jpnh4 ) zfact = rno3 * 1.0e+6 
            IF( jn == jppo4  )                 zfact = po4r * 1.0e+6
            cltra = TRIM( ctrcnm(jn) )                  ! short title for tracer
            IF( iom_use( cltra ) )  CALL iom_put( cltra, trn(:,:,:,jn) * zfact )
         END DO

         IF( iom_use( "INTDIC" ) ) THEN                     !   DIC content in kg/m2
            zdic(:,:) = 0.
            DO jk = 1, jpkm1
               zdic(:,:) = zdic(:,:) + trn(:,:,jk,jpdic) * e3t_n(:,:,jk) * tmask(:,:,jk) * 12.
            ENDDO
            CALL iom_put( 'INTDIC', zdic )     
         ENDIF
         !
         IF( iom_use( "O2MIN" ) .OR. iom_use ( "ZO2MIN" ) ) THEN  ! Oxygen minimum concentration and depth 
            zo2min   (:,:) = trn(:,:,1,jpoxy) * tmask(:,:,1)
            zdepo2min(:,:) = gdepw_n(:,:,1)   * tmask(:,:,1)
            DO jk = 2, jpkm1
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     IF( tmask(ji,jj,jk) == 1 ) then
                        IF( trn(ji,jj,jk,jpoxy) < zo2min(ji,jj) ) then
                           zo2min   (ji,jj) = trn(ji,jj,jk,jpoxy)
                           zdepo2min(ji,jj) = gdepw_n(ji,jj,jk)
                        ENDIF
                     ENDIF
                  END DO
               END DO
            END DO
            !
            CALL iom_put('O2MIN' , zo2min     )                              ! oxygen minimum concentration
            CALL iom_put('ZO2MIN', zdepo2min  )                              ! depth of oxygen minimum concentration
             !
         ENDIF
     ENDIF
      !
   END SUBROUTINE trc_wri_pisces

#else

   INTEGER ::   nid_T, nz_T, nh_T, ndim_T, ndim_hT   ! grid_T file
   INTEGER ::          nb_T              , ndim_bT   ! grid_T file
   INTEGER ::   ndex(1)                              ! ???
   INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: ndex_hT
   INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: ndex_T
   INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: ndex_bT

CONTAINS
   !!----------------------------------------------------------------------
   !!  Dummy module :                                     No passive tracer
   !!----------------------------------------------------------------------
   SUBROUTINE trc_wri_pisces (kt) 
      !
      INTEGER, INTENT(IN) :: kt
      !
      LOGICAL ::   ll_print = .FALSE.                        ! =T print and flush numout
      CHARACTER (len=40) ::   clhstnam, clop, clmx           ! local names
      INTEGER  ::   inum = 11                                ! temporary logical unit
      INTEGER  ::   ji, jj, jk                               ! dummy loop indices
      INTEGER  ::   iimi, iima, ipk, it, itmod, ijmi, ijma   ! local integers
      INTEGER  ::   jn, ierr
      REAL(wp) ::   zsto, zout, zmax, zjulian                ! local scalars
      CHARACTER (len=20)           :: cltra
      REAL(wp)                     :: zfact
      REAL(wp), DIMENSION(jpi,jpj) :: zdic, zo2min, zdepo2min

      clop = "x"
      zsto=rdt
      clop = "ave("//TRIM(clop)//")"
      zout = nn_write * rdt
      zmax = ( nitend - nit000 + 1 ) * rdt

      ! Define indices of the horizontal output zoom and vertical limit storage
      iimi = 1      ;      iima = jpi
      ijmi = 1      ;      ijma = jpj
      ipk = jpk

      ! define time axis
      it = kt
      itmod = kt - nit000 + 1

      IF( kt == nit000 ) THEN

         ierr = 0
         ALLOCATE( ndex_hT(jpi*jpj) , ndex_T(jpi*jpj*jpk) , STAT=ierr )
         !
         CALL ymds2ju( nyear, nmonth, nday, rdt, zjulian )
         zjulian = zjulian - adatrj   !   set calendar origin to the beginning of the experiment

         CALL dia_nam( clhstnam, nn_write, 'pisces' )
         CALL histbeg( clhstnam, jpi, glamt, jpj, gphit,           & 
            &          iimi, iima-iimi+1, ijmi, ijma-ijmi+1,       &
            &          nit000-1, zjulian, rdt, nh_T, nid_T, domain_id=nidom, snc4chunks=snc4set )
         CALL histvert( nid_T, "deptht", "Vertical T levels",      &
            &           "m", ipk, gdept_1d, nz_T, "down" )
         CALL wheneq( jpi*jpj*ipk, tmask, 1, 1., ndex_T , ndim_T  )      ! volume
         CALL wheneq( jpi*jpj    , tmask, 1, 1., ndex_hT, ndim_hT )      ! surface

         IF( ln_p2z ) THEN
             DO jn = jp_pcs0, jp_pcs1
                cltra = TRIM( ctrcnm(jn) )
                CALL histdef( nid_T, cltra, TRIM(cltra)//" Concentration", "mol/L",&
                & jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
             ENDDO
         ELSE
             DO jn = jp_pcs0, jp_pcs1
                cltra = TRIM( ctrcnm(jn) )
                CALL histdef( nid_T, cltra, TRIM(cltra)//" Concentration", "mol/L",&
                & jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
             END DO
         ENDIF

         CALL histdef( nid_T, "INTDIC", "Integrated DIC" , "mol/L"  ,   &
         &    jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "O2MIN", "Minimum O2" , "mol/L"  ,   &
         &    jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "ZO2MIN", "Depth of Minimum O2" , "m"  ,   &
         &    jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )

         CALL histend( nid_T, snc4chunks=snc4set )

      ENDIF


      IF( ln_p2z ) THEN
         DO jn = jp_pcs0, jp_pcs1
            cltra = TRIM( ctrcnm(jn) )                  ! short title for tracer
            CALL histwrite( nid_T, cltra, it, trn(:,:,:,jn), ndim_T , ndex_T  )
         END DO
      ELSE
         DO jn = jp_pcs0, jp_pcs1
            zfact = 1.0e+6 
            IF( jn == jpno3 .OR. jn == jpnh4 ) zfact = rno3 * 1.0e+6 
            IF( jn == jppo4  )                 zfact = po4r * 1.0e+6
            cltra = TRIM( ctrcnm(jn) )                  ! short title for tracer
            CALL histwrite( nid_T, cltra, it, trn(:,:,:,jn)*zfact, ndim_T , ndex_T  )
         END DO
      ENDIF

      zdic(:,:) = 0.
      DO jk = 1, jpkm1
           zdic(:,:) = zdic(:,:) + trn(:,:,jk,jpdic) * e3t_n(:,:,jk) * tmask(:,:,jk) * 12.
      ENDDO
      CALL histwrite( nid_T, 'INTDIC', it, zdic, ndim_hT , ndex_hT  )
      !
      zo2min   (:,:) = trn(:,:,1,jpoxy) * tmask(:,:,1)
      zdepo2min(:,:) = gdepw_n(:,:,1)   * tmask(:,:,1)
      DO jk = 2, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( tmask(ji,jj,jk) == 1 ) then
                  IF( trn(ji,jj,jk,jpoxy) < zo2min(ji,jj) ) then
                     zo2min   (ji,jj) = trn(ji,jj,jk,jpoxy)
                     zdepo2min(ji,jj) = gdepw_n(ji,jj,jk)
                  ENDIF
               ENDIF
            END DO
         END DO
      END DO
      CALL histwrite( nid_T, 'O2MIN', it, zo2min, ndim_hT , ndex_hT  )
      CALL histwrite( nid_T, 'ZO2MIN', it, zdepo2min, ndim_hT , ndex_hT  )
      !
      IF( kt == nitend ) THEN
         CALL histclo( nid_T )
      ENDIF
    END SUBROUTINE trc_wri_pisces
#endif
#endif
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcwri_pisces.F90 10069 2018-08-28 14:12:24Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!======================================================================
END MODULE trcwri_pisces
