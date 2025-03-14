MODULE sbcssr
   !!======================================================================
   !!                       ***  MODULE  sbcssr  ***
   !! Surface module :  heat and fresh water fluxes a restoring term toward observed SST/SSS
   !!======================================================================
   !! History :  3.0  !  2006-06  (G. Madec)  Original code
   !!            3.2  !  2009-04  (B. Lemaire)  Introduce iom_put
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_ssr       : add to sbc a restoring term toward SST/SSS climatology
   !!   sbc_ssr_init  : initialisation of surface restoring
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE sbc_oce        ! surface boundary condition
   USE phycst         ! physical constants
   USE sbcrnf         ! surface boundary condition : runoffs
   !
   USE fldread        ! read input fields
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager
   USE lib_mpp        ! distribued memory computing library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)  
   USE stopack        ! Stochastic physics package

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_ssr        ! routine called in sbcmod
   PUBLIC   sbc_ssr_init   ! routine called in sbcmod
   PUBLIC   sbc_ssr_alloc  ! routine called in sbcmod

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   erp   !: evaporation damping   [kg/m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   qrp   !: heat flux damping        [w/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   coefice   !: under ice relaxation coefficient

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   wei   !: heat flux damping        [w/m2]
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   zdc   !: heat flux damping        [w/m2]

   !                                   !!* Namelist namsbc_ssr *
   INTEGER, PUBLIC ::   nn_sstr         ! SST/SSS restoring indicator
   INTEGER, PUBLIC ::   nn_sssr         ! SST/SSS restoring indicator
   INTEGER, PUBLIC ::   nn_qtor         ! SST/SSS restoring indicator
   LOGICAL,SAVE    ::   ln_qtdata=.false. ! restoring from file
   REAL(wp)        ::   rn_dqdt         ! restoring factor on SST and SSS
   REAL(wp)        ::   rn_deds         ! restoring factor on SST and SSS
   REAL(wp)        ::   rn_qtor         ! restoring factor on SST and SSS
   LOGICAL         ::   ln_sssr_bnd     ! flag to bound erp term 
   REAL(wp)        ::   rn_sssr_bnd     ! ABS(Max./Min.) value of erp term [mm/day]
   INTEGER         ::   nn_sssr_ice     ! Control of restoring under ice
   INTEGER         ::   nn_sstr_ice     ! Control of restoring under ice

   REAL(wp) , ALLOCATABLE, DIMENSION(:) ::   buffer   ! Temporary buffer for exchange
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_sst   ! structure of input SST (file informations, fields read)
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_sss   ! structure of input SSS (file informations, fields read)
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_qto   ! structure of input SSS (file informations, fields read)
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_qtd   ! structure of input SSS (file informations, fields read)

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: sbcssr.F90 12276 2019-12-20 11:14:26Z cetlod $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sbc_ssr( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE sbc_ssr  ***
      !!
      !! ** Purpose :   Add to heat and/or freshwater fluxes a damping term
      !!                toward observed SST and/or SSS.
      !!
      !! ** Method  : - Read namelist namsbc_ssr
      !!              - Read observed SST and/or SSS
      !!              - at each nscb time step
      !!                   add a retroaction term on qns    (nn_sstr = 1)
      !!                   add a damping term on sfx        (nn_sssr = 1)
      !!                   add a damping term on emp        (nn_sssr = 2)
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in   ) ::   kt   ! ocean time step
      !!
      INTEGER  ::   ji, jj   ! dummy loop indices
      REAL(wp) ::   zerp     ! local scalar for evaporation damping
      REAL(wp) ::   zqrp     ! local scalar for heat flux damping
      REAL(wp) ::   zsrp     ! local scalar for unit conversion of rn_deds factor
      REAL(wp) ::   zerp_bnd ! local scalar for unit conversion of rn_epr_max factor
      REAL(wp) ::   zqtor(jpi,jpj) ! timescale
      INTEGER  ::   ierror   ! return error code
      !!
      CHARACTER(len=100) ::  cn_dir          ! Root directory for location of ssr files
      !!----------------------------------------------------------------------
      !
      IF( nn_sstr + nn_sssr + nn_qtor /= 0 ) THEN
         !
         IF( nn_sstr == 1)   CALL fld_read( kt, nn_fsbc, sf_sst )   ! Read SST data and provides it at kt
         IF( nn_sssr >= 1)   CALL fld_read( kt, nn_fsbc, sf_sss )   ! Read SSS data and provides it at kt
         IF( nn_qtor >= 1)   THEN
                 CALL fld_read( kt, nn_fsbc, sf_qto )
                 IF ( ln_qtdata ) THEN
                      CALL fld_read( kt, nn_fsbc, sf_qtd )   ! Read SST data and provides it at kt
                      zqtor(:,:) = sf_qtd(1)%fnow(:,:,1)
                 ELSE
                      zqtor(:,:) = rn_qtor
                 ENDIF
         ENDIF
         !
         !                                         ! ========================= !
         IF( MOD( kt-1, nn_fsbc ) == 0 ) THEN      !    Add restoring term     !
            !                                      ! ========================= !
            !
            IF( nn_sstr == 1 ) THEN                                   !* Temperature restoring term

               qrp = rn_dqdt
               IF( ln_stopack .AND. nn_spp_dqdt > 0 ) CALL spp_gen(kt, qrp,nn_spp_dqdt,rn_dqdt_sd,jk_spp_dqdt )

               IF( nn_sstr_ice == 1 ) THEN
                 DO jj = 1, jpj
                    DO ji = 1, jpi
                       zqrp = qrp(ji,jj)  * ( sst_m(ji,jj) - sf_sst(1)%fnow(ji,jj,1) ) * tmask(ji,jj,1) * wei(ji,jj) * &
                            & exp ( -fr_i(ji,jj)*fr_i(ji,jj)/ 0.16_wp )
                       qns(ji,jj) = qns(ji,jj) + zqrp
                       qrp(ji,jj) = zqrp
                    END DO
                 END DO
               ELSE
                 DO jj = 1, jpj
                    DO ji = 1, jpi
                       zqrp = qrp(ji,jj) * ( sst_m(ji,jj) - sf_sst(1)%fnow(ji,jj,1) ) * tmask(ji,jj,1) * wei(ji,jj)
                       qns(ji,jj) = qns(ji,jj) + zqrp
                       qrp(ji,jj) = zqrp
                    END DO
                 END DO
               ENDIF
            ENDIF
            !
            IF( nn_sstr .ne. 1 ) qrp(:,:) =0._wp
            !
            IF( nn_qtor == 1 ) THEN                                   !* Temperature restoring term
               IF( nn_sstr_ice == 1 ) THEN
                 DO jj = 1, jpj
                    DO ji = 1, jpi
                       zqrp = zqtor(ji,jj) * ( qns(ji,jj) + qsr(ji,jj) - sf_qto(1)%fnow(ji,jj,1) ) * tmask(ji,jj,1) * wei(ji,jj) * &
                            & exp ( -fr_i(ji,jj)*fr_i(ji,jj)/ 0.16_wp )
                       qns(ji,jj) = qns(ji,jj) + zqrp
                       qrp(ji,jj) = qrp(ji,jj) + zqrp
                    END DO
                 END DO
               ELSE
                 DO jj = 1, jpj
                    DO ji = 1, jpi
                       zqrp = zqtor(ji,jj) * ( qns(ji,jj) + qsr(ji,jj) - sf_qto(1)%fnow(ji,jj,1) ) * tmask(ji,jj,1) * wei(ji,jj)
                       qns(ji,jj) = qns(ji,jj) + zqrp
                       qrp(ji,jj) = qrp(ji,jj) + zqrp
                    END DO
                 END DO
               ENDIF
            ENDIF
            !
            IF( nn_qtor == 2 ) THEN                                   !* Temperature restoring term
               IF( nn_sstr_ice == 1 ) THEN
                 DO jj = 1, jpj
                    DO ji = 1, jpi
                       zqrp = zqtor(ji,jj) * ( qns(ji,jj) + qsr(ji,jj) - sf_qto(1)%fnow(ji,jj,1) ) * tmask(ji,jj,1) * wei(ji,jj) * &
                            & exp ( -fr_i(ji,jj)*fr_i(ji,jj)/ 0.16_wp )
                       qsr(ji,jj) = qns(ji,jj) + zqrp
                       qrp(ji,jj) = qrp(ji,jj) + zqrp
                    END DO
                 END DO
               ELSE
                 DO jj = 1, jpj
                    DO ji = 1, jpi
                       zqrp = zqtor(ji,jj) * ( qns(ji,jj) + qsr(ji,jj) - sf_qto(1)%fnow(ji,jj,1) ) * tmask(ji,jj,1) * wei(ji,jj)
                       qsr(ji,jj) = qns(ji,jj) + zqrp
                       qrp(ji,jj) = qrp(ji,jj) + zqrp
                    END DO
                 END DO
               ENDIF
            ENDIF
            !
            !
            IF( nn_sssr /= 0 .AND. nn_sssr_ice /= 1 ) THEN
              ! use fraction of ice ( fr_i ) to adjust relaxation under ice if nn_sssr_ice .ne. 1
              ! n.b. coefice is initialised and fixed to 1._wp if nn_sssr_ice = 1

               DO jj = 1, jpj
                  DO ji = 1, jpi
                     SELECT CASE ( nn_sssr_ice )
                       CASE ( 0 )    ;  coefice(ji,jj) = 1._wp - fr_i(ji,jj)              ! no/reduced damping under ice
                       CASE  DEFAULT ;  coefice(ji,jj) = 1._wp + ( nn_sssr_ice - 1 ) * fr_i(ji,jj) ! reinforced damping (x nn_sssr_ice) under ice )
                     END SELECT
                  END DO
               END DO
            ENDIF
            !
            IF( nn_sssr == 1 ) THEN                                   !* Salinity damping term (salt flux only (sfx))

               erp = rn_deds / rday                                  ! from [mm/day] to [kg/m2/s]
               IF( ln_stopack .AND. nn_spp_dedt > 0 ) CALL spp_gen(kt, erp, nn_spp_dedt, rn_dedt_sd, jk_spp_deds )

               DO jj = 1, jpj
                  DO ji = 1, jpi
                     zerp = erp(ji,jj) * ( 1. - 2.*rnfmsk(ji,jj) )   &      ! No damping in vicinity of river mouths
                        &        *   coefice(ji,jj)            &      ! Optional control of damping under sea-ice
                        &        * ( sss_m(ji,jj) - sf_sss(1)%fnow(ji,jj,1) ) * tmask(ji,jj,1) * wei(ji,jj)
                     sfx(ji,jj) = sfx(ji,jj) + zerp                 ! salt flux
                     erp(ji,jj) = zerp / MAX( sss_m(ji,jj), 1.e-20 ) ! converted into an equivalent volume flux (diagnostic only)
                  END DO
               END DO
               !
            ELSEIF( nn_sssr == 2 ) THEN                               !* Salinity damping term (volume flux (emp) and associated heat flux (qns)
               erp = rn_deds / rday                                  ! from [mm/day] to [kg/m2/s]
               IF( ln_stopack .AND. nn_spp_dedt > 0 ) CALL spp_gen(kt, erp, nn_spp_dedt, rn_dedt_sd, jk_spp_deds )

               zerp_bnd = rn_sssr_bnd / rday                          !       -              -    
               DO jj = 1, jpj
                  DO ji = 1, jpi                            
                     zerp = erp(ji,jj) * ( 1. - 2.*rnfmsk(ji,jj) )   &      ! No damping in vicinity of river mouths
                        &        *   coefice(ji,jj)            &      ! Optional control of damping under sea-ice
                        &        * ( sss_m(ji,jj) - sf_sss(1)%fnow(ji,jj,1) )   &
                        &        / MAX(  sss_m(ji,jj), 1.e-20   ) * tmask(ji,jj,1) * wei(ji,jj)
                     IF( ln_sssr_bnd )   zerp = SIGN( 1., zerp ) * MIN( zerp_bnd, ABS(zerp) )
                     emp(ji,jj) = emp (ji,jj) + zerp
                     qns(ji,jj) = qns(ji,jj) - zerp * rcp * sst_m(ji,jj)
                     erp(ji,jj) = zerp
                  END DO
               END DO
            ENDIF
            !
         ENDIF
         !
      ENDIF
      !
   END SUBROUTINE sbc_ssr

 
   SUBROUTINE sbc_ssr_init
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_ssr_init  ***
      !!
      !! ** Purpose :   initialisation of surface damping term
      !!
      !! ** Method  : - Read namelist namsbc_ssr
      !!              - Read observed SST and/or SSS if required
      !!---------------------------------------------------------------------
      INTEGER  ::   ji, jj   ! dummy loop indices
      REAL(wp) ::   zerp     ! local scalar for evaporation damping
      REAL(wp) ::   zqrp     ! local scalar for heat flux damping
      REAL(wp) ::   zsrp     ! local scalar for unit conversion of rn_deds factor
      REAL(wp) ::   zerp_bnd ! local scalar for unit conversion of rn_epr_max factor
      INTEGER  ::   ierror   ! return error code
      REAL(wp) ::   rn_dcoast = 0._wp ! distance to coast weigthing
      !!
      CHARACTER(len=100) ::  cn_dir          ! Root directory for location of ssr files
      TYPE(FLD_N) ::   sn_sst, sn_sss,sn_qto,sn_qtd        ! informations about the fields to be read
      NAMELIST/namsbc_ssr/ cn_dir, nn_sstr, nn_sssr, rn_dqdt, rn_deds, sn_sst, &
              & sn_sss, ln_sssr_bnd, rn_sssr_bnd, nn_sssr_ice, rn_dcoast, nn_sstr_ice,&
              & sn_qto, nn_qtor, rn_qtor, ln_qtdata, sn_qtd
      INTEGER     ::  ios, inum
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'sbc_ssr : SST and/or SSS damping term '
         WRITE(numout,*) '~~~~~~~ '
      ENDIF
      ! 
      REWIND( numnam_ref )              ! Namelist namsbc_ssr in reference namelist : 
      READ  ( numnam_ref, namsbc_ssr, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namsbc_ssr in reference namelist' )

      REWIND( numnam_cfg )              ! Namelist namsbc_ssr in configuration namelist :
      READ  ( numnam_cfg, namsbc_ssr, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namsbc_ssr in configuration namelist' )
      IF(lwm) WRITE ( numond, namsbc_ssr )

      IF(lwp) THEN                 !* control print
         WRITE(numout,*) '   Namelist namsbc_ssr :'
         WRITE(numout,*) '      SST restoring term (Yes=1)             nn_sstr        = ', nn_sstr
         WRITE(numout,*) '         dQ/dT (restoring magnitude on SST)     rn_dqdt     = ', rn_dqdt, ' W/m2/K'
         WRITE(numout,*) '      SSS damping term (Yes=1, salt   flux)  nn_sssr        = ', nn_sssr
         WRITE(numout,*) '                       (Yes=2, volume flux) '
         WRITE(numout,*) '         dE/dS (restoring magnitude on SST)     rn_deds     = ', rn_deds, ' mm/day'
         WRITE(numout,*) '         flag to bound erp term                 ln_sssr_bnd = ', ln_sssr_bnd
         WRITE(numout,*) '         ABS(Max./Min.) erp threshold           rn_sssr_bnd = ', rn_sssr_bnd, ' mm/day'
         WRITE(numout,*) '      Cntrl of surface restoration under ice nn_sssr_ice    = ', nn_sssr_ice
         WRITE(numout,*) '          ( 0 = no restoration under ice)'
         WRITE(numout,*) '          ( 1 = restoration everywhere  )'
         WRITE(numout,*) '          (>1 = enhanced restoration under ice  )'
         WRITE(numout,*) '      Cntrl of SST     restoration under ice nn_sssr_ice    = ', nn_sstr_ice
         WRITE(numout,*) '          ( 0 = restoration everywhere  )'
         WRITE(numout,*) '          ( 1 = partial restoration under ice  )'
         WRITE(numout,*) '      Use distance to coast weigthing (if >0) (m) rn_dcoast = ', rn_dcoast
         WRITE(numout,*) '      QTOT damping term (Yes=1)              nn_qtor        = ', nn_qtor
         WRITE(numout,*) '          restoring magnitude                rn_qtor        = ', rn_qtor
         WRITE(numout,*) '          restoring magnitude from file      ln_qtdata      = ', ln_qtdata
      ENDIF
      !
      IF( nn_sstr == 1 ) THEN      !* set sf_sst structure & allocate arrays
         !
         ALLOCATE( sf_sst(1), STAT=ierror )
         IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_ssr: unable to allocate sf_sst structure' )
         ALLOCATE( sf_sst(1)%fnow(jpi,jpj,1), STAT=ierror )
         IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_ssr: unable to allocate sf_sst now array' )
         !
         ! fill sf_sst with sn_sst and control print
         CALL fld_fill( sf_sst, (/ sn_sst /), cn_dir, 'sbc_ssr', 'SST restoring term toward SST data', 'namsbc_ssr', no_print )
         IF( sf_sst(1)%ln_tint )   ALLOCATE( sf_sst(1)%fdta(jpi,jpj,1,2), STAT=ierror )
         IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_ssr: unable to allocate sf_sst data array' )
         !
      ENDIF
      !
      IF( nn_sssr >= 1 ) THEN      !* set sf_sss structure & allocate arrays
         !
         ALLOCATE( sf_sss(1), STAT=ierror )
         IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_ssr: unable to allocate sf_sss structure' )
         ALLOCATE( sf_sss(1)%fnow(jpi,jpj,1), STAT=ierror )
         IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_ssr: unable to allocate sf_sss now array' )
         !
         ! fill sf_sss with sn_sss and control print
         CALL fld_fill( sf_sss, (/ sn_sss /), cn_dir, 'sbc_ssr', 'SSS restoring term toward SSS data', 'namsbc_ssr', no_print )
         IF( sf_sss(1)%ln_tint )   ALLOCATE( sf_sss(1)%fdta(jpi,jpj,1,2), STAT=ierror )
         IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_ssr: unable to allocate sf_sss data array' )
         !
      ENDIF
      !
      IF( nn_qtor >= 1 ) THEN      !* set sf_sss structure & allocate arrays
         !
         ALLOCATE( sf_qto(1), STAT=ierror )
         IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_ssr: unable to allocate sf_qto structure' )
         ALLOCATE( sf_qto(1)%fnow(jpi,jpj,1), STAT=ierror )
         IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_ssr: unable to allocate sf_qto now array' )
         !
         ! fill sf_sss with sn_sss and control print
         CALL fld_fill( sf_qto, (/ sn_qto /), cn_dir, 'sbc_ssr', 'QTO restoring term toward QTO data', 'namsbc_ssr', no_print )
         IF( sf_qto(1)%ln_tint )   ALLOCATE( sf_qto(1)%fdta(jpi,jpj,1,2), STAT=ierror )
         IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_ssr: unable to allocate sf_qto data array' )
         !
         IF( ln_qtdata ) THEN
            ALLOCATE( sf_qtd(1), STAT=ierror )
            IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_ssr: unable to allocate sf_qtd structure' )
            ALLOCATE( sf_qtd(1)%fnow(jpi,jpj,1), STAT=ierror )
            IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_ssr: unable to allocate sf_qtd now array' )
            CALL fld_fill( sf_qtd, (/ sn_qtd /), cn_dir, 'sbc_ssr', 'QTO restoring weigth', 'namsbc_ssr', no_print )
            IF( sf_qtd(1)%ln_tint )   ALLOCATE( sf_qtd(1)%fdta(jpi,jpj,1,2), STAT=ierror )
            IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_ssr: unable to allocate sf_qtd data array' )
         ENDIF
         !
      ENDIF
      !
      coefice(:,:) = 1._wp         !  Initialise coefice to 1._wp ; will not need to be changed if nn_sssr_ice=1
      !                            !* Initialize qrp and erp if no restoring 
      IF( nn_sstr /= 1                   )   qrp(:,:) = 0._wp
      IF( nn_sssr /= 1 .OR. nn_sssr /= 2 )   erp(:,:) = 0._wp

      IF( rn_dcoast > 0._wp ) THEN
            CALL iom_open( 'dist.coast.nc', inum )
            CALL iom_get ( inum, jpdom_data, 'Tcoast2d', zdc )
            CALL iom_close( inum )
            DO jj=1,jpj
               DO ji=1,jpi
                  wei(ji,jj) = 1._wp - exp ( - zdc(ji,jj)**5 / ( rn_dcoast**5 ) )
               ENDDO
            ENDDO
      ELSE
            wei(:,:) = 1._wp
      ENDIF
      !
   END SUBROUTINE sbc_ssr_init
         
   INTEGER FUNCTION sbc_ssr_alloc()
      !!----------------------------------------------------------------------
      !!               ***  FUNCTION sbc_ssr_alloc  ***
      !!----------------------------------------------------------------------
      sbc_ssr_alloc = 0       ! set to zero if no array to be allocated
      IF( .NOT. ALLOCATED( erp ) ) THEN
         ALLOCATE( qrp(jpi,jpj), erp(jpi,jpj), coefice(jpi,jpj), &
                   wei(jpi,jpj), zdc(jpi,jpj), STAT= sbc_ssr_alloc )
         !
         IF( lk_mpp                  )   CALL mpp_sum ( 'sbcssr', sbc_ssr_alloc )
         IF( sbc_ssr_alloc /= 0 )   CALL ctl_warn('sbc_ssr_alloc: failed to allocate arrays.')
         !
      ENDIF
   END FUNCTION
      
   !!======================================================================
END MODULE sbcssr
