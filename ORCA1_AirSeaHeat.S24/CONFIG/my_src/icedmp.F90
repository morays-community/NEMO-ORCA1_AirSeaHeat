MODULE icedmp
   !!======================================================================
   !!                       ***  MODULE  icedmp  ***
   !! Ice damping module :  relax sea-ice concentraion and thickness to obs 
   !!======================================================================
   !! History :  4.0  !  2021-07  (A. Storto)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   ice_dmp       : sea-ice restoring
   !!   ice_dmp_init  : initialisation of sea-ice restoring
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE ice            ! ice variables
   USE dom_oce        ! ocean space and time domain
   USE sbc_oce        ! surface module
   USE phycst         ! physical constants
   !
   USE fldread        ! read input fields
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager
   USE lib_mpp        ! distribued memory computing library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)  

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_dmp        ! routine called in ice_stp
   PUBLIC   ice_dmp_init   ! routine called in ice_stp

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   wei
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   zmd
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   zbs
   !
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   weic
   !                                   !!* Namelist namice_dmp *
   INTEGER, PUBLIC ::   nn_sicr         ! SIC restoring indicator
   INTEGER, PUBLIC ::   nn_sitr         ! SIT restoring indicator
   REAL(wp)        ::   zn_dsic         ! restoring factor on SIC
   REAL(wp)        ::   zn_dsit         ! restoring factor on SIT
   INTEGER, PUBLIC ::   nn_sicm         ! SIC restoring method
   INTEGER, PUBLIC ::   nn_sicmv        ! SIC restoring method for thickness
   REAL(wp)        ::   rn_dcoast       ! Distance to coast filtering
   REAL(wp)        ::   rn_minicediff, zhicifmin
   REAL(wp)        ::   rn_amin, rn_amax

   REAL(wp) , ALLOCATABLE, DIMENSION(:) ::   buffer   ! Temporary buffer for exchange
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_sic   ! structure of input SIC (file informations, fields read)
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_sit   ! structure of input SIT (file informations, fields read)
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_cwe   ! structure of input SIC weigths (file informations, fields read)

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: sbcssr.F90 12276 2019-12-20 11:14:26Z cetlod $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_dmp( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE ice_dmp  ***
      !!
      !! ** Purpose :   Damping to observed sea-ice variables.
      !!
      !!              - Read observed SIC and/or SIT
      !!                   add a retroaction term on SIC and/or SIT
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in   ) ::   kt   ! ocean time step
      !!
      INTEGER  ::   ji, jj   ! dummy loop indices
      INTEGER  ::   ierror   ! return error code
      INTEGER  ::   ki, kj, kc
      !!
      CHARACTER(len=100) ::  cn_dir          ! Root directory for location of ssr files
      TYPE(FLD_N) ::   sn_sic, sn_sit        ! informations about the fields to be read
      REAL(wp) ::  zmd(jpi,jpj), zbs(jpi,jpj), zdiff, zs1, zdiffm(jpl)
      REAL(wp) ::  zs2, zzh
      !!------------------------------------------------------------------------
      !
      IF( nn_sicr + nn_sitr /= 0 ) THEN
         !
         IF( nn_sicr == 1)   CALL fld_read( kt, nn_fsbc, sf_sic )   ! Read SIC data and provides it at kt
         IF( nn_sicr == 1 .AND. nn_sicm == 2 )   CALL fld_read( kt, nn_fsbc, sf_cwe )   ! Read SIC weights data
         IF( nn_sitr == 1)   CALL fld_read( kt, nn_fsbc, sf_sit )   ! Read SIT data and provides it at kt
         !
         IF( nn_sicr == 1 ) THEN                                   !* SIC restoring term
                 !
                 zmd(:,:) = at_i_b(:,:)
                 zbs(:,:) = sf_sic(1)%fnow(:,:,1)
                 DO jj = 1, jpj
                    DO ji = 1, jpi
                       !
                       IF( ABS ( zmd(ji,jj) - zbs(ji,jj) ) .gt. rn_minicediff ) THEN
                             !
                             zdiff = ( zbs(ji,jj) - zmd(ji,jj) ) * wei(ji,jj) * zn_dsic
                             zs1 = SUM( a_i_b(ji,jj,: ) )
                             !
                             IF( nn_sicm .eq. 1 ) THEN
                                 IF ( zs1 .GT. epsi20 ) THEN
                                     zdiffm (:) = a_i_b(ji,jj,:) / zs1
                                 ELSE
                                     zdiffm (:) = 1._wp / REAL( jpl, wp)
                                 ENDIF
                             ELSEIF( nn_sicm .eq. 2 ) THEN
                                 ki  = INT( ( glamt(ji,jj) + 180._wp ) / 10._wp ) +1
                                 kj  = INT( ( gphit(ji,jj) +  90._wp ) / 10._wp ) +1
                                 kc  = MIN( 5,INT( zs1 / 0.2_wp ) + 1 )
                                 zs2 = SUM ( weic(ki,kj,kc,:) )
                                 IF( zs2 .GT. epsi20 ) THEN 
                                     zdiffm (:) = weic(ki,kj,kc,:) / zs2
                                 ELSE
                                     IF ( zs1 .GT. epsi20 ) THEN
                                          zdiffm (:) = a_i_b(ji,jj,:) / zs1
                                     ELSE
                                          zdiffm (:) = 1._wp / REAL( jpl, wp)
                                     ENDIF
                                 ENDIF
                             ELSEIF( nn_sicm .eq. 3 ) THEN
!                                zdiffm (:) = 
                             ENDIF
                             ! Adjust
                             IF( at_i(ji,jj) + zdiff > rn_amax ) THEN
                                 zdiff = rn_amax - at_i(ji,jj)
                             ENDIF
                             IF( at_i(ji,jj) + zdiff < rn_amin ) THEN
                                 zdiff = rn_amin - at_i(ji,jj)
                             ENDIF
                             ! Apply to now and back fields
                             at_i(ji,jj)     = at_i(ji,jj)   + zdiff
                             a_i(ji,jj,:)    = a_i(ji,jj,:)  + zdiffm(:) * zdiff
                             at_i_b(ji,jj)   = at_i(ji,jj)   + zdiff
                             a_i_b(ji,jj,:)  = a_i(ji,jj,:)  + zdiffm(:) * zdiff
                             IF( zs1 .le. epsi20 .and. hm_i(ji,jj) .lt. zhicifmin .AND. at_i(ji,jj) .gt. epsi20 ) THEN
                                 hm_i(ji,jj)  = zhicifmin
                                 h_i(ji,jj,1) = zhicifmin
                                 h_i_b(ji,jj,1) = zhicifmin
                                 v_i(ji,jj,1) = h_i (ji,jj,1) * a_i  (ji,jj,1)
                                 v_i_b(ji,jj,1) = h_i_b (ji,jj,1) * a_i_b  (ji,jj,1)
                             ENDIF
                             !  Take care of thickness
                             IF( nn_sicmv .EQ. 1 ) THEN ! Change thickness (useless??)
                                h_i_b(ji,jj,:) = v_i_b(ji,jj,:) / a_i_b(ji,jj,:) 
                                h_i  (ji,jj,:) = v_i  (ji,jj,:) / a_i  (ji,jj,:) 
                                IF( at_i(ji,jj) > epsi20 ) THEN
                                      zzh = 1._wp / at_i(ji,jj)
                                ELSE
                                      zzh = 0._wp
                                ENDIF
                                hm_i(ji,jj) = vt_i(ji,jj) * zzh
                             ELSEIF( nn_sicmv .EQ. 2 ) THEN
                                v_i_b(ji,jj,:) = v_i_b(ji,jj,:) + zdiffm(:)*zdiff*h_i_b(ji,jj,:)
                                v_i  (ji,jj,:) = v_i  (ji,jj,:) + zdiffm(:)*zdiff*h_i_b(ji,jj,:)
                                vt_i(ji,jj)    = SUM( v_i (ji,jj,:) )
                             ELSEIF( nn_sicmv .EQ. 3 ) THEN
                                v_i_b(ji,jj,:) = h_i_b(ji,jj,:) * a_i_b(ji,jj,:)
                                v_i  (ji,jj,:) = h_i  (ji,jj,:) * a_i  (ji,jj,:)
                                vt_i(ji,jj)    = SUM( v_i (ji,jj,:) )
                             ENDIF
                             !
                       ENDIF
                    END DO
                 END DO
            !
         ENDIF
            !
      ENDIF
         !
   END SUBROUTINE ice_dmp

   SUBROUTINE ice_dmp_init
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE ice_dmp_init  ***
      !!
      !! ** Purpose :   initialisation of surface damping term
      !!
      !! ** Method  : - Read namelist namice_dmp
      !!              - Read observed SST and/or SSS if required
      !!---------------------------------------------------------------------
      INTEGER  ::   ji, jj   ! dummy loop indices
      INTEGER  ::   ierror   ! return error code
      REAL(wp) ::   rn_dcoast = 0._wp ! distance to coast weigthing
      REAL(wp)        ::   rn_dsic         ! restoring factor on SST and SSS
      REAL(wp)        ::   rn_dsit         ! restoring factor on SST and SSS
      !!
      CHARACTER(len=100) ::  cn_dir           ! Root directory for location of ssr files
      TYPE(FLD_N) ::   sn_sic, sn_sit, sn_cwe ! informations about the fields to be read
      INTEGER     ::  ios, inum
      INTEGER     ::  jc,ia(3)
      REAL(wp)    ::  zs
      NAMELIST/namice_dmp/ cn_dir, nn_sicr, nn_sitr, rn_dsic, rn_dsit, sn_sic, &
                         & sn_sit, rn_dcoast, nn_sicm, nn_sicmv, rn_minicediff,&
                         & zhicifmin, rn_amin, rn_amax
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'ice_dmp : SST and/or SSS damping term '
         WRITE(numout,*) '~~~~~~~ '
      ENDIF
      ! 
      IF( ice_dmp_alloc() .ne. 0 ) CALL ctl_stop('STOP', 'ice_dmp: unable to allocate arrays' )
      ! 
      REWIND( numnam_ice_ref )              ! Namelist namice_dmp
      READ  ( numnam_ice_ref, namice_dmp, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namice_dmp in reference namelist' )
      REWIND( numnam_ice_cfg )              ! Namelist namice_dmp in configuration namelist
      READ  ( numnam_ice_cfg, namice_dmp, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namice_dmp in configuration namelist' )
      IF(lwm) WRITE ( numoni, namice_dmp )
      ! 
      IF(lwp) THEN                 !* control print
         WRITE(numout,*) '   Namelist namice_dmp :'
         WRITE(numout,*) '      SIC restoring term (Yes=1)             nn_sicr        = ', nn_sicr
         WRITE(numout,*) '      SIT restoring term (Yes=1)             nn_sitr        = ', nn_sitr
         WRITE(numout,*) '      SIC method for spreading               nn_sicm        = ', nn_sicm
         WRITE(numout,*) '      SIC method for spreading on thickness  nn_sicmv       = ', nn_sicmv
         WRITE(numout,*) '      Restoring magnitude on SIC             rn_dsic        = ', rn_dsic
         WRITE(numout,*) '      Restoring magnitude on SIT             rn_dsit        = ', rn_dsit
         WRITE(numout,*) '      Use distance to coast weigth (m,if >0) rn_dcoast      = ', rn_dcoast
         WRITE(numout,*) '      Minumum mod-obs diff in SIC for nudg)  rn_minicediff  = ', rn_minicediff
         WRITE(numout,*) '      Minumum sea-ice thickness for SIC>0    zhicifmin      = ', zhicifmin
         WRITE(numout,*) '      Min/Max sea-ice concentration rn_amin, rn_amax        = ', rn_amin, rn_amax
      ENDIF
      !
      IF( nn_sicr == 1 ) THEN      !* set sf_sst structure & allocate arrays
         !
         IF( nn_sicm  .lt. 1 .or. nn_sicm  .gt. 2 ) CALL ctl_stop( 'STOP', 'ice_dmp: nn_sicm not valid' )
         IF( nn_sicmv .lt. 0 .or. nn_sicmv .gt. 3 ) CALL ctl_stop( 'STOP', 'ice_dmp: nn_sicm not valid' )
         !
         ALLOCATE( sf_sic(1), STAT=ierror )
         IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'ice_dmp: unable to allocate sf_sic structure' )
         ALLOCATE( sf_sic(1)%fnow(jpi,jpj,1), STAT=ierror )
         IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'ice_dmp: unable to allocate sf_sic now array' )
         !
         CALL fld_fill( sf_sic, (/ sn_sic /), cn_dir, 'ice_dmp', 'SIC restoring data', 'namice_dmp', no_print )
         IF( sf_sic(1)%ln_tint )   ALLOCATE( sf_sic(1)%fdta(jpi,jpj,1,2), STAT=ierror )
         IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'ice_dmp: unable to allocate sf_sst data array' )
         !
         IF( nn_sicm .eq. 2 ) THEN
            ALLOCATE( weic(36,18,5,jpl) )
            OPEN(9018,FILE='weigc.txt',STATUS='OLD')
            DO jc=1,5
               DO jj=1,18
                  DO ji=1,36
                     READ(9018,*) ia(1),ia(2),ia(3),&
                   & weic(ji,jj,jc,1),weic(ji,jj,jc,2),weic(ji,jj,jc,3),&
                   & weic(ji,jj,jc,4),weic(ji,jj,jc,5)
                     zs = SUM( weic(ji,jj,jc,:) )
                     IF( zs .gt. epsi20 ) weic(ji,jj,jc,:) = weic(ji,jj,jc,:) / zs
                  ENDDO
               ENDDO
            ENDDO
            CLOSE(9018)
         ENDIF
         !
         zn_dsic = REAL(nn_fsbc,wp) * rn_rdt / (rn_dsic*86400._wp)
         !
      ENDIF
      !
      IF( nn_sitr == 1 ) THEN      !* set sf_sst structure & allocate arrays
         !
         ALLOCATE( sf_sit(1), STAT=ierror )
         IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'ice_dmp: unable to allocate sf_sit structure' )
         ALLOCATE( sf_sit(1)%fnow(jpi,jpj,1), STAT=ierror )
         IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'ice_dmp: unable to allocate sf_sit now array' )
         !
         CALL fld_fill( sf_sit, (/ sn_sit /), cn_dir, 'ice_dmp', 'SIT restoring data', 'namice_dmp', no_print )
         IF( sf_sit(1)%ln_tint )   ALLOCATE( sf_sit(1)%fdta(jpi,jpj,1,2), STAT=ierror )
         IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'ice_dmp: unable to allocate sf_sst data array' )
         !
         zn_dsit = REAL(nn_fsbc,wp) * rn_rdt / (rn_dsit**86400._wp)
         !
      ENDIF

      !
      IF( rn_dcoast > 0._wp ) THEN
            CALL iom_open( 'dist.coast.nc', inum )
            CALL iom_get ( inum, jpdom_data, 'Tcoast2d', wei )
            CALL iom_close( inum )
            DO jj=1,jpj
               DO ji=1,jpi
                  wei(ji,jj) = 1._wp - exp ( - wei(ji,jj) / ( rn_dcoast ) )
               ENDDO
            ENDDO
      ELSE
            wei(:,:) = 1._wp
      ENDIF
      !
   END SUBROUTINE ice_dmp_init
         
   INTEGER FUNCTION ice_dmp_alloc()
      !!----------------------------------------------------------------------
      !!               ***  FUNCTION ice_dmp_alloc  ***
      !!----------------------------------------------------------------------
      ice_dmp_alloc = 0       ! set to zero if no array to be allocated
      IF( .NOT. ALLOCATED( zmd ) ) THEN
         ALLOCATE( zmd(jpi,jpj), zbs(jpi,jpj), wei(jpi,jpj), STAT= ice_dmp_alloc )
         !
         IF( lk_mpp                  )   CALL mpp_sum ( 'icedmp', ice_dmp_alloc )
         IF( ice_dmp_alloc /= 0 )   CALL ctl_warn('ice_dmp_alloc: failed to allocate arrays.')
         !
      ENDIF
   END FUNCTION
      
   !!======================================================================
END MODULE icedmp
