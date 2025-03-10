MODULE sbcblk_corr
   !
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE sbc_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   USE fldread        ! read input fields
   USE lib_fortran    ! to use key_nosignedzero
   !
   USE iom            ! I/O manager library
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! distribued memory computing library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE prtctl         ! Print control
   !
   IMPLICIT NONE
   PRIVATE
   !
   PUBLIC   sbc_blk_corr
   !
   INTEGER , PARAMETER ::   jpfld   = 8           ! maximum number of files to read
   INTEGER , PARAMETER ::   jp_wnds = 1           ! index of 10m wind velocity (i-component) (m/s)    at T-point
   INTEGER , PARAMETER ::   jp_wndd = 2           ! index of 10m wind velocity (j-component) (m/s)    at T-point
   INTEGER , PARAMETER ::   jp_tair = 3           ! index of 10m air temperature             (Kelvin)
   INTEGER , PARAMETER ::   jp_humi = 4           ! index of specific humidity               ( % )
   INTEGER , PARAMETER ::   jp_qsr  = 5           ! index of solar heat                      (W/m2)
   INTEGER , PARAMETER ::   jp_qlw  = 6           ! index of Long wave                       (W/m2)
   INTEGER , PARAMETER ::   jp_prec = 7           ! index of total precipitation (rain+snow) (Kg/m2/s)
   INTEGER , PARAMETER ::   jp_slp  = 8           ! index of sea level pressure              (Pa)
   !
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf   ! structure of input fields (file informations, fields read)
   CHARACTER(LEN=4), PARAMETER :: cvlst(jpfld) = &
   & (/'wnds','wndd','tair','humi','qsr ','qlw ','prec','slp '/)
   !
   LOGICAL, SAVE       :: ln_init = .false.
   LOGICAL, SAVE       :: ln_vcorr ( jpfld ) = .true. 
   !
   REAL(wp), DIMENSION(:,:), ALLOCATABLE :: wsp, wdr
   !
   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: sbcblk.F90 15613 2021-12-22 09:35:54Z cetlod $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sbc_blk_corr(kt,u10,v10,slp,t2m,hum,prc,snw,qsr,qlw)
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_blk_corr  ***
      !!
      !! ** Purpose :   correct the bulk formulae input fields
      !!
      !! ** Method  : 
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(IN)  ::  kt
      REAL(wp),INTENT(INOUT), DIMENSION(jpi,jpj) :: u10,v10,slp,t2m,hum,prc,snw,qsr,qlw
      !!
      !!
      INTEGER  ::   jfpr, jfld            ! dummy loop indice and argument
      INTEGER  ::   ios, ierror, ioptio   ! Local integer
      INTEGER  ::   nn_cgstoch
      INTEGER  ::   ji,jj
      !!
      CHARACTER(len=100)            ::   cn_dir                ! Root directory for location of atmospheric forcing files
      TYPE(FLD_N), DIMENSION(jpfld) ::   slf_i                 ! array of namelist informations on the fields to read
      TYPE(FLD_N) ::   sn_wnds, sn_wndd, sn_humi, sn_qsr       ! informations about the fields to be read
      TYPE(FLD_N) ::   sn_qlw , sn_tair, sn_prec, sn_slp
      !
      REAL(wp), PARAMETER :: r2d = 57.29578_wp
      REAL(wp) :: z1,z2
      NAMELIST/namsbc_blk_corr/ sn_wnds, sn_wndd, sn_humi, sn_qsr,               &
      &                         sn_qlw , sn_tair, sn_prec, sn_slp, cn_dir
      !!---------------------------------------------------------------------
      IF(.not. ln_init) THEN
              !
              IF(lwp) THEN
                  WRITE(numout,*)
                  WRITE(numout,*) 'sbc_blk_corr : correction of bulk formulae inputs'
                  WRITE(numout,*) '~~~~~~~~~~~~'
              ENDIF
              !
              ln_init = .TRUE.
              !                                   !** read bulk namelist  
              REWIND( numnam_ref )                !* Namelist namsbc_blk in reference namelist : bulk parameters
              READ  ( numnam_ref, namsbc_blk_corr, IOSTAT = ios, ERR = 901)
901           IF( ios /= 0 )   CALL ctl_nam ( ios , 'namsbc_blk_corr in reference namelist' )
              !
              REWIND( numnam_cfg )                !* Namelist namsbc_blk in configuration namelist : bulk parameters
              READ  ( numnam_cfg, namsbc_blk_corr, IOSTAT = ios, ERR = 902 )
902           IF( ios >  0 )   CALL ctl_nam ( ios , 'namsbc_blk_corr in configuration namelist' )
              !
              IF(lwm) WRITE( numond, namsbc_blk_corr )
              !
              !                                      !- store namelist information in an array
              slf_i(jp_wnds) = sn_wnds   ;   slf_i(jp_wndd) = sn_wndd
              slf_i(jp_qsr ) = sn_qsr    ;   slf_i(jp_qlw ) = sn_qlw
              slf_i(jp_tair) = sn_tair   ;   slf_i(jp_humi) = sn_humi
              slf_i(jp_prec) = sn_prec   ;   slf_i(jp_slp)  = sn_slp
              !
              jfld = jpfld
              !                                      !- allocate the bulk structure
              ALLOCATE( sf(jfld), STAT=ierror )
              IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_blk_corr: unable to allocate sf structure' )
              !                                      !- fill the bulk structure with namelist informations
              CALL fld_fill( sf, slf_i, cn_dir, 'sbc_blk_corr', 'surface boundary condition -- bulk formulae correction', 'namsbc_blk_corr' )
              !
              DO jfpr = 1, jfld
                 !
                 IF( TRIM(sf(jfpr)%clrootname) == 'NOT USED' ) THEN    !--  not used field  --!   (only now allocated and set to zero)
                    ln_vcorr(jfpr) = .false.
                 ELSE                                                  !-- used field --!
                    ALLOCATE( sf(jfpr)%fnow(jpi,jpj,1) )
                    IF( slf_i(jfpr)%ln_tint )   ALLOCATE( sf(jfpr)%fdta(jpi,jpj,1,2) )
                 ENDIF
                 !
                 IF(lwp) WRITE(numout,*) 'Correction of ',trim(cvlst(jfpr)),' : ',ln_vcorr(jfpr)
              ENDDO
              !
              ALLOCATE ( wsp(jpi,jpj), wdr(jpi,jpj) )
              wsp(:,:) = 0._wp ;       wdr(:,:) = 0._wp
      ENDIF
      ! 
      IF( COUNT( (/ln_vcorr(jp_wnds),ln_vcorr(jp_wndd)/) ) .EQ. 1 ) THEN
              CALL ctl_stop( 'STOP', 'sbc_blk_corr : inconsistent jp_wnds / jp_wndd' )
      ENDIF
      ! 
      CALL fld_read( kt, nn_fsbc, sf )             ! input fields provided at the current time-step
      ! 
      IF ( ln_vcorr( jp_tair ) ) t2m(:,:) = t2m(:,:) + sf(jp_tair)%fnow(:,:,1) * tmask(:,:,1)
      IF ( ln_vcorr( jp_slp  ) ) slp(:,:) = slp(:,:) + sf(jp_slp )%fnow(:,:,1) * tmask(:,:,1)
      IF ( ln_vcorr( jp_qsr  ) ) THEN
!              z1=minval( qsr, mask=(tmask(:,:,1) > 0.5_wp) )
!              z2=maxval( qsr, mask=(tmask(:,:,1) > 0.5_wp) )
!              CALL mpp_min ( 'blk', z1) ; CALL mpp_max ( 'blk', z2)
!              IF(lwp) WRITE(numout,*) 'Inner Minmax for qsr :',z1, z2

!              z1=minval( sf(jp_qsr )%fnow(:,:,1), mask=(tmask(:,:,1) > 0.5_wp) )
!              z2=maxval( sf(jp_qsr )%fnow(:,:,1), mask=(tmask(:,:,1) > 0.5_wp) )
!              CALL mpp_min ( 'blk', z1) ; CALL mpp_max ( 'blk', z2)
!              IF(lwp) WRITE(numout,*) 'Inner Minmax for qsr corr :',z1, z2

              qsr(:,:) = qsr(:,:) * sf(jp_qsr )%fnow(:,:,1) 

!              z1=minval( qsr, mask=(tmask(:,:,1) > 0.5_wp) )
!              z2=maxval( qsr, mask=(tmask(:,:,1) > 0.5_wp) )
!              CALL mpp_min ( 'blk', z1) ; CALL mpp_max ( 'blk', z2)
!              IF(lwp) WRITE(numout,*) 'Inner Minmax after for qsr :',z1, z2
      ENDIF
      IF ( ln_vcorr( jp_qlw  ) ) qlw(:,:) = qlw(:,:) * sf(jp_qlw )%fnow(:,:,1) 
      IF ( ln_vcorr( jp_humi ) ) hum(:,:) = hum(:,:) * sf(jp_humi)%fnow(:,:,1) 
      IF ( ln_vcorr( jp_wnds ) .AND. ln_vcorr( jp_wndd ) ) THEN
           wsp(:,:) = sqrt ( u10(:,:)*u10(:,:) + v10(:,:)*v10(:,:) )
           wdr(:,:) = atan2( v10(:,:)          , u10(:,:)           ) * r2d
           wsp(:,:) = wsp(:,:) + sf(jp_wnds)%fnow(:,:,1) * tmask(:,:,1)
           wdr(:,:) = wdr(:,:) + sf(jp_wndd)%fnow(:,:,1) * tmask(:,:,1)
           u10(:,:) = wsp(:,:) * cos( wdr(:,:) / r2d )  
           v10(:,:) = wsp(:,:) * sin( wdr(:,:) / r2d )
      ENDIF
      IF ( ln_vcorr( jp_prec ) ) THEN
           wsp(:,:) = 0._wp
           WHERE( prc(:,:) .GT. 1.e-9_wp ) wsp(:,:) = snw(:,:) / prc(:,:)
           WHERE( wsp(:,:) .LT. 0._wp ) wsp(:,:) = 0._wp
           WHERE( wsp(:,:) .GT. 1._wp ) wsp(:,:) = 1._wp
           prc(:,:) = prc(:,:) * sf(jp_prec)%fnow(:,:,1) 
           snw(:,:) = wsp(:,:) * prc(:,:)
      ENDIF
      ! 
   END SUBROUTINE sbc_blk_corr
   !!======================================================================
END MODULE sbcblk_corr
