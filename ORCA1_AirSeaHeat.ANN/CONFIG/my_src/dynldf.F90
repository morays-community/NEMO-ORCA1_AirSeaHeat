MODULE dynldf
   !!======================================================================
   !!                       ***  MODULE  dynldf  ***
   !! Ocean physics:  lateral diffusivity trends 
   !!=====================================================================
   !! History :  2.0  ! 2005-11  (G. Madec)  Original code (new step architecture)
   !!            3.7  ! 2014-01  (F. Lemarie, G. Madec)  restructuration/simplification of ahm specification,
   !!                 !                                  add velocity dependent coefficient and optional read in file
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_ldf      : update the dynamics trend with the lateral diffusion
   !!   dyn_ldf_init : initialization, namelist read, and parameters control
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   USE ldfdyn         ! lateral diffusion: eddy viscosity coef.
   USE dynldf_lap_blp ! lateral mixing   (dyn_ldf_lap & dyn_ldf_blp routines)
   USE dynldf_iso     ! lateral mixing                 (dyn_ldf_iso routine )
   USE trd_oce        ! trends: ocean variables
   USE trddyn         ! trend manager: dynamics   (trd_dyn      routine)
   !
   USE prtctl         ! Print control
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! distribued memory computing library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE timing         ! Timing
   USE stopack        ! Stochastic physics package

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_ldf       ! called by step module 
   PUBLIC   dyn_ldf_init  ! called by opa  module 

   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: ahmt0, ahmf0

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: dynldf.F90 10068 2018-08-28 14:09:04Z nicolasmartin $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_ldf( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_ldf  ***
      !! 
      !! ** Purpose :   compute the lateral ocean dynamics physics.
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   ztrdu, ztrdv
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('dyn_ldf')
      !
      IF( l_trddyn )   THEN                      ! temporary save of momentum trends
         ALLOCATE( ztrdu(jpi,jpj,jpk) , ztrdv(jpi,jpj,jpk) )
         ztrdu(:,:,:) = ua(:,:,:) 
         ztrdv(:,:,:) = va(:,:,:) 
      ENDIF
      !
      IF( ln_stopack .AND. nn_spp_ahm1 .GT. 0) THEN
           IF( kt .EQ. nit000 ) ahmt0 = ahmt
           ahmt = ahmt0
           CALL spp_ahm(kt,ahmt,nn_spp_ahm1,rn_ahm1_sd,jk_spp_ahm1)
      ENDIF
      IF( ln_stopack .AND. nn_spp_ahm2 .GT. 0) THEN
           IF( kt .EQ. nit000 ) ahmf0 = ahmf
           ahmf = ahmf0
           CALL spp_ahm(kt,ahmf,nn_spp_ahm2,rn_ahm2_sd,jk_spp_ahm2)
      ENDIF
      !
      SELECT CASE ( nldf_dyn )                   ! compute lateral mixing trend and add it to the general trend
      !
      CASE ( np_lap   )    ;   CALL dyn_ldf_lap( kt, ub, vb, ua, va, 1 )      ! iso-level    laplacian
      CASE ( np_lap_i )    ;   CALL dyn_ldf_iso( kt )                         ! rotated      laplacian
      CASE ( np_blp   )    ;   CALL dyn_ldf_blp( kt, ub, vb, ua, va    )      ! iso-level bi-laplacian
      !
      END SELECT

      IF( l_trddyn ) THEN                        ! save the horizontal diffusive trends for further diagnostics
         ztrdu(:,:,:) = ua(:,:,:) - ztrdu(:,:,:)
         ztrdv(:,:,:) = va(:,:,:) - ztrdv(:,:,:)
         CALL trd_dyn( ztrdu, ztrdv, jpdyn_ldf, kt )
         DEALLOCATE ( ztrdu , ztrdv )
      ENDIF
      !                                          ! print sum trends (used for debugging)
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' ldf  - Ua: ', mask1=umask,   &
         &                       tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
      !
      IF( ln_timing )   CALL timing_stop('dyn_ldf')
      !
   END SUBROUTINE dyn_ldf


   SUBROUTINE dyn_ldf_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_ldf_init  ***
      !! 
      !! ** Purpose :   initializations of the horizontal ocean dynamics physics
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN                     !==  Namelist print  ==!
         WRITE(numout,*)
         WRITE(numout,*) 'dyn_ldf_init : Choice of the lateral diffusive operator on dynamics'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namdyn_ldf: already read in ldfdyn module'
         WRITE(numout,*) '      see ldf_dyn_init report for lateral mixing parameters'
         WRITE(numout,*)
         !
         SELECT CASE( nldf_dyn )             ! print the choice of operator
         CASE( np_no_ldf )   ;   WRITE(numout,*) '   ==>>>   NO lateral viscosity'
         CASE( np_lap    )   ;   WRITE(numout,*) '   ==>>>   iso-level laplacian operator'
         CASE( np_lap_i  )   ;   WRITE(numout,*) '   ==>>>   rotated laplacian operator with iso-level background'
         CASE( np_blp    )   ;   WRITE(numout,*) '   ==>>>   iso-level bi-laplacian operator'
         END SELECT
      ENDIF
      IF( ln_stopack .AND. nn_spp_ahm1 .GT. 0) ALLOCATE( ahmt0(jpi,jpj,jpk) )
      IF( ln_stopack .AND. nn_spp_ahm2 .GT. 0) ALLOCATE( ahmf0(jpi,jpj,jpk) )
      !
   END SUBROUTINE dyn_ldf_init

   !!======================================================================
END MODULE dynldf
