MODULE CGSTOCH

   USE par_kind
   USE oce, ONLY : sshn ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE sbc_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   USE iom            ! I/O routines
   USE lib_mpp

IMPLICIT NONE

PRIVATE

PUBLIC :: CGSTOCH_INIT, SCATTER1, SCATTER2, SCATTER3, GATHER1, SCATTER2b, SCATTER2c, &
        & q_sat0, rho_air0, L_vap0, cp_air0

LOGICAL, PUBLIC :: ln_cgstoch = .false.
LOGICAL, PUBLIC :: ln_cgstbia = .true. 
LOGICAL, PUBLIC :: ln_cgsst = .true. ,ln_cgssq = .true.,ln_cgwnd= .true.,&
                 & ln_cgrho= .true.,ln_cgtau= .true.,ln_cgflx= .true.
REAL(wp), ALLOCATABLE, DIMENSION(:,:), PUBLIC :: &
        & zst0, ztpot0, zsq0, sf_humi0, wndm0, Cd_atm0, &
        & Ch_atm0, Ce_atm0, t_zu0, q_zu0, sf_slp0, &
        & zU_zu0, cdn_oce0, chn_oce0, cen_oce0, zrhoa0, &
        & sst_sgst, ssc_sgst, cg_tmask, zps, zps2, &
        & taum0, zwnd_i0, zwnd_j0, &
        & zqla0,zevap0,zqsb0,sf_tair0

INTEGER :: ni, nj, ncg 
INTEGER, PUBLIC :: ncg_jpi, ncg_jpj
LOGICAL, SAVE :: llinit_rnd = .false.
LOGICAL, SAVE :: llinit_rnd2= .false.

!... Ziggurat Section

! integer,  parameter       ::  dp=wp
real(kind=dp), parameter  ::  m1=2147483648.0_dp,   m2=2147483648.0_dp, &
                              half=0.5_dp, ve=0.003949659822581572_dp, &
                              vn=0.00991256303526217_dp
real(kind=dp)             ::  dn=3.442619855899_dp, tn=3.442619855899_dp, &
                              de=7.697117470131487_dp, &
                              te=7.697117470131487_dp
real(kind=dp)             ::  q
integer,  save            ::  iz, jz, jsr=123456789, kn(0:127), ke(0:255), hz
real(kind=dp), save       ::  wn(0:127), fn(0:127), we(0:255), fe(0:255)
logical,  save            ::  initialized=.false.
public                    ::  zigset, shr3, uni, rnor, rexp

interface uni
   module procedure uni_scalar, uni_vec, uni_vec_given_seed
end interface uni
interface rnor
   module procedure rnor_scalar,rnor_vec
end interface rnor

!... END Ziggurat Section

contains
!
SUBROUTINE CGSTOCH_INIT( iaux )

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: iaux

        INTEGER :: inum
        REAL(wp) :: zs
        REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: ztst, ztst0, ztst1

        ncg = iaux

        IF( lwp ) THEN
                write(numout,*) ' CGSTOCH '
                write(numout,*) ' ~~~~~~~'
                write(numout,*) '     Coarse-grained stochastic bulk formulas'
                write(numout,*) '     Initialization'
        ENDIF

        ! Read the standard deviations
        ALLOCATE ( sst_sgst(jpi,jpj), ssc_sgst(jpi,jpj) )
        IF( lwp ) write(numout,*) '     Reading cgstoch.nc input file'
        CALL iom_open( 'cgstoch', inum )
        CALL iom_get( inum, jpdom_autoglo, 'sst', sst_sgst)
        CALL iom_get( inum, jpdom_autoglo, 'ssc', ssc_sgst)
        CALL iom_close( inum )

        ni = ncg * jpi
        nj = ncg * jpj
        ln_cgstoch = .TRUE.
        IF( lwp ) write(numout,*) '     CG domain fact: ',ncg
        IF( lwp ) write(numout,*) '     CG domain size: ',ni,nj

        ! Allocating ancillary variables
        ALLOCATE ( zst0(ni,nj), ztpot0(ni,nj), zsq0(ni,nj), sf_humi0(ni,nj), wndm0(ni,nj), Cd_atm0(ni,nj), &
        & Ch_atm0(ni,nj), Ce_atm0(ni,nj), t_zu0(ni,nj), q_zu0(ni,nj), sf_slp0(ni,nj), &
        & zU_zu0(ni,nj), cdn_oce0(ni,nj), chn_oce0(ni,nj), cen_oce0(ni,nj), zrhoa0(ni,nj), &
        & taum0(ni,nj), zwnd_i0(ni,nj), zwnd_j0(ni,nj), &
        & zqla0(ni,nj), zevap0(ni,nj), zqsb0(ni,nj), sf_tair0(ni,nj) )

        CALL SYSTEM_CLOCK(inum)
        IF(nmember .ne. 0) THEN
                inum = inum+(nmember*100)
        ENDIF
        CALL zigset(inum)
        IF( lwp ) write(numout,*) '     Ziggurat seed : ',inum

        ncg_jpi = ni
        ncg_jpj = nj

        ALLOCATE( cg_tmask(ni,nj) )
        CALL SCATTER1(tmask,cg_tmask)

        ! TEST the routines
        ALLOCATE (ztst(jpi,jpj), ztst0(ni,nj),ztst1(jpi,jpj) )
        ztst = sshn * tmask(:,:,1)
        ztst0= 0._wp
        ztst1= 0._wp
        IF( lwp ) write(numout,*) '     CH1:',SUM(ABS(ztst)),SUM(ABS(ztst0)),SUM(ABS(ztst1))
        CALL SCATTER1t(ztst,ztst0)
        IF( lwp ) write(numout,*) '     CH2:',SUM(ABS(ztst)),SUM(ABS(ztst0)),SUM(ABS(ztst1))
        CALL GATHER1t(ztst0,ztst1)
        IF( lwp ) write(numout,*) '     CH3:',SUM(ABS(ztst)),SUM(ABS(ztst0)),SUM(ABS(ztst1))
        zs=SUM( ABS( ztst-ztst1) )
        IF( lwp ) write(numout,*) '     TEST: sum of dfs', zs
        DEALLOCATE ( ztst, ztst0, ztst1)

END SUBROUTINE CGSTOCH_INIT

SUBROUTINE SCATTER1(fldin,fldou)
        IMPLICIT NONE
        REAL(wp), INTENT(IN) :: fldin(jpi,jpj)
        REAL(wp), INTENT(OUT):: fldou( ni, nj)
        INTEGER  :: ji,jj,ki,kj
        DO jj=1,jpj*ncg
           kj=(jj-1)/ncg+1
           DO ji=1,jpi*ncg
              ki=(ji-1)/ncg+1
              fldou(ji,jj) = fldin(ki,kj)
           ENDDO
        ENDDO
END SUBROUTINE SCATTER1

SUBROUTINE SCATTER1t(fldin,fldou)
        IMPLICIT NONE
        REAL(wp), INTENT(IN) :: fldin(jpi,jpj)
        REAL(wp), INTENT(OUT):: fldou( ni, nj)
        INTEGER  :: ji,jj,ki,kj
!        write(5000+narea,*) 'ST',jpi,jpj,ni,nj,ncg
        DO jj=1,jpj*ncg
           kj=(jj-1)/ncg+1
           DO ji=1,jpi*ncg
              ki=(ji-1)/ncg+1
              fldou(ji,jj) = fldin(ki,kj)
 !             write(5000+narea,'(4I4,2F12.3)') ji,jj,ki,kj,fldin(ki,kj),fldou(ji,jj)
           ENDDO
        ENDDO
  !      write(5000+narea,*) 'EN'
END SUBROUTINE SCATTER1t

SUBROUTINE SCATTER2b(fldin,fldou,zstdv)
        IMPLICIT NONE
        REAL(wp), INTENT(IN) :: fldin(jpi,jpj)
        REAL(wp), INTENT(IN) :: zstdv(jpi,jpj)
        REAL(wp), INTENT(OUT):: fldou( ni, nj)

        INTEGER  :: ji,jj,ki,kj,kk,ki2,kj2
        REAL(wp)             :: zz(ni*nj)
        REAL(wp)             :: zs(jpi,jpj)
        REAL(wp)             :: zcg2, zr, zp
        REAL(wp)             :: zaa(ncg,ncg), zbb(ncg,ncg)
        REAL(wp)             :: taus, taut, zpo

        taus = REAL( ncg*ncg, wp )
        taut = REAL( nn_fsbc*rn_rdt/86400._wp, wp ) 
        zcg2 = REAL( ncg*ncg, wp )
        CALL SCATTER1(fldin,fldou)
        zaa(:,:) = exp(-1._wp/taus )
        zaa(1,1) = exp(-1._wp/taut )
        zbb(:,:) = sqrt ( 1._wp - zaa*zaa )
        if(.not.llinit_rnd) then
                ALLOCATE ( zps(jpi,jpj) )
                IF(lwp) write(numout,*) ' init random field'
                ZZ(1:jpi*jpj)=RNOR(jpi*jpj)
                kk=0
                zpo=zz(1)
                DO jj=1,jpj
                 DO ji=1,jpi
                   kk=kk+1
                   zps(ji,jj) = zaa(ncg,ncg)*zpo+zbb(ncg,ncg)*zz(kk)
                 ENDDO
                ENDDO
                llinit_rnd = .true.
        endif

        kk = 0
        zs = 0._wp
        ZZ=RNOR(ni*nj)
        DO kj=1,jpj
           DO ki=1,jpi
              zpo = zps(ki,kj)
              DO kj2=1,ncg 
               DO ki2=1,ncg 
                 kk=kk+1
                 ji=(ki-1)*ncg+ki2
                 jj=(kj-1)*ncg+kj2
                 zp = zaa(ki2,kj2) * zpo + zbb(ki2,kj2) * zz(kk)
                 fldou(ji,jj) = fldou(ji,jj) + zp
                 zpo = zp
               ENDDO
              ENDDO
              zps(ki,kj) = zpo
           ENDDO
        ENDDO

END SUBROUTINE SCATTER2b

SUBROUTINE SCATTER2c(fldin,fldou,zstdv)
        IMPLICIT NONE
        REAL(wp), INTENT(IN) :: fldin(jpi,jpj)
        REAL(wp), INTENT(IN) :: zstdv(jpi,jpj)
        REAL(wp), INTENT(OUT):: fldou( ni, nj)

        INTEGER  :: ji,jj,ki,kj,kk,ki2,kj2
        REAL(wp)             :: zz(ni*nj)
        REAL(wp)             :: zs(jpi,jpj)
        REAL(wp)             :: zcg2, zr, zp
        REAL(wp)             :: zaa(ncg,ncg), zbb(ncg,ncg)
        REAL(wp)             :: taus, taut, zpo

        taus = REAL( ncg*ncg, wp )
        taut = REAL( nn_fsbc*rn_rdt/86400._wp, wp ) 
        zcg2 = REAL( ncg*ncg, wp )
        CALL SCATTER1(fldin,fldou)
        zaa(:,:) = exp(-1._wp/taus )
        zaa(1,1) = exp(-1._wp/taut )
        zbb(:,:) = sqrt ( 1._wp - zaa*zaa )
        if(.not.llinit_rnd2) then
                ALLOCATE ( zps2(jpi,jpj) )
                IF(lwp) write(numout,*) ' init random field'
                ZZ(1:jpi*jpj)=RNOR(jpi*jpj)
                kk=0
                zpo=zz(1)
                DO jj=1,jpj
                 DO ji=1,jpi
                   kk=kk+1
                   zps2(ji,jj) = zaa(ncg,ncg)*zpo+zbb(ncg,ncg)*zz(kk)
                 ENDDO
                ENDDO
                llinit_rnd2 = .true.
        endif

        kk = 0
        zs = 0._wp
        ZZ=RNOR(ni*nj)
        DO kj=1,jpj
           DO ki=1,jpi
              zpo = zps2(ki,kj)
              DO kj2=1,ncg 
               DO ki2=1,ncg 
                 kk=kk+1
                 ji=(ki-1)*ncg+ki2
                 jj=(kj-1)*ncg+kj2
                 zp = zaa(ki2,kj2) * zpo + zbb(ki2,kj2) * zz(kk)
                 fldou(ji,jj) = fldou(ji,jj) + zp
                 zpo = zp
               ENDDO
              ENDDO
              zps2(ki,kj) = zpo
           ENDDO
        ENDDO

END SUBROUTINE SCATTER2c

SUBROUTINE SCATTER2(fldin,fldou,zstdv)
        IMPLICIT NONE
        REAL(wp), INTENT(IN) :: fldin(jpi,jpj)
        REAL(wp), INTENT(IN) :: zstdv(jpi,jpj)
        REAL(wp), INTENT(OUT):: fldou( ni, nj)

        INTEGER  :: ji,jj,ki,kj,kk
        REAL(wp)             :: zz(ni*nj)
        REAL(wp)             :: zs(jpi,jpj)
        REAL(wp)             :: zp( ni, nj)
        REAL(wp)             :: zcg2, zr

        zcg2 = REAL( ncg*ncg, wp )
        CALL SCATTER1(fldin,fldou)
        ZZ=RNOR(ni*nj)
        kk = 0
        zs = 0._wp
        DO jj=1,jpj*ncg
           kj=(jj-1)/ncg+1
           DO ji=1,jpi*ncg
              ki=(ji-1)/ncg+1
              kk=kk+1
              zp(ji,jj) = tmask(ki,kj,1)*zz(kk)*zstdv(ki,kj)
              zs(ki,kj) = zs(ki,kj) + zp(ji,jj) / zcg2
           ENDDO
        ENDDO
        IF( .NOT. ln_cgstbia ) zs(:,:) = 0._wp
        !
        DO jj=1,jpj*ncg
           kj=(jj-1)/ncg+1
           DO ji=1,jpi*ncg
              ki=(ji-1)/ncg+1
              zr = zp(ji,jj) - zs(ki,kj)
              fldou(ji,jj) = fldou(ji,jj) + zr
           ENDDO
        ENDDO

END SUBROUTINE SCATTER2

SUBROUTINE SCATTER3(fldin,fldou,zstdv)
        IMPLICIT NONE
        REAL(wp), INTENT(IN) :: fldin(jpi,jpj)
        REAL(wp), INTENT(IN) :: zstdv(jpi,jpj)
        REAL(wp), INTENT(OUT):: fldou( ni, nj)

        INTEGER  :: ji,jj,ki,kj,kk
        REAL(wp)             :: zz(ni*nj)
        REAL(wp)             :: zs(jpi,jpj)
        REAL(wp)             :: zp( ni, nj)
        REAL(wp)             :: zcg2, zr

        zcg2 = REAL( ncg*ncg, wp )
        CALL SCATTER1(fldin,fldou)
        ZZ=RNOR(ni*nj)
        kk = 0
        zs = 0._wp
        DO jj=1,jpj*ncg
           kj=(jj-1)/ncg+1
           DO ji=1,jpi*ncg
              ki=(ji-1)/ncg+1
              kk=kk+1
              zp(ji,jj) = tmask(ki,kj,1)*zz(kk)*zstdv(ki,kj)
              zs(ki,kj) = zs(ki,kj) + zp(ji,jj)*zp(ji,jj) / zcg2
           ENDDO
        ENDDO
        !
        IF( .NOT. ln_cgstbia ) zs(:,:) = 0._wp
        !
        DO jj=1,jpj*ncg
           kj=(jj-1)/ncg+1
           DO ji=1,jpi*ncg
              ki=(ji-1)/ncg+1
              zr = zp(ji,jj) - SQRT ( zs(ki,kj) )
              fldou(ji,jj) = SQRT ( fldou(ji,jj)*fldou(ji,jj) - zr*zr )
           ENDDO
        ENDDO

END SUBROUTINE SCATTER3

SUBROUTINE GATHER1t(fldin,fldou)
        IMPLICIT NONE
        REAL(wp), INTENT(IN) :: fldin( ni, nj)
        REAL(wp), INTENT(OUT):: fldou(jpi,jpj)
        INTEGER  :: ji,jj,ki,kj
        REAL(wp)             :: zcg2
        fldou = 0._wp
        zcg2 = REAL( ncg*ncg, wp )
        !write(6000+narea,*) 'ST',ni,nj,jpi,jpj,ncg,zcg2
        DO jj=1,jpj*ncg
           kj=(jj-1)/ncg+1
           DO ji=1,jpi*ncg
              ki=(ji-1)/ncg+1
              fldou(ki,kj) = fldou(ki,kj) + fldin(ji,jj) / zcg2
        !      write(6000+narea,'(4I4,2F12.3)') ji,jj,ki,kj,fldin(ji,jj),fldou(ki,kj)
           ENDDO
        ENDDO
        !write(6000+narea,*) 'EN'
        !fldou = fldou * tmask(:,:,1)
END SUBROUTINE GATHER1t

SUBROUTINE GATHER1(fldin,fldou)
        IMPLICIT NONE
        REAL(wp), INTENT(IN) :: fldin( ni, nj)
        REAL(wp), INTENT(OUT):: fldou(jpi,jpj)
        INTEGER  :: ji,jj,ki,kj
        REAL(wp)             :: zcg2
        fldou = 0._wp
        zcg2 = REAL( ncg*ncg, wp )
        DO jj=1,jpj*ncg
           kj=(jj-1)/ncg+1
           DO ji=1,jpi*ncg
              ki=(ji-1)/ncg+1
              fldou(ki,kj) = fldou(ki,kj) + fldin(ji,jj) / zcg2
           ENDDO
        ENDDO
        !fldou = fldou * tmask(:,:,1)
END SUBROUTINE GATHER1

   FUNCTION q_sat0( ptak, pslp, reps )
      !!----------------------------------------------------------------------------------
      !!                           ***  FUNCTION q_sat  ***
      !!
      !! ** Purpose : Specific humidity at saturation in [kg/kg]
      !!              Based on accurate estimate of "e_sat"
      !!              aka saturation water vapor (Goff, 1957)
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://sourceforge.net/p/aerobulk)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION( ni, nj), INTENT(in) ::   ptak    ! air temperature                       [K]
      REAL(wp), DIMENSION( ni, nj), INTENT(in) ::   pslp    ! sea level atmospheric pressure       [Pa]
      REAL(wp), DIMENSION( ni, nj)             ::   q_sat0  ! Specific humidity at saturation   [kg/kg]
      REAL(wp)                    , INTENT(in) ::   reps
      !
      INTEGER  ::   ji, jj         ! dummy loop indices
      REAL(wp) ::   ze_sat, ztmp   ! local scalar
      !!----------------------------------------------------------------------------------
      !
      DO jj = 1, nj
         DO ji = 1, ni
            !
            ztmp = rt0 / ptak(ji,jj)
            !
            ! Vapour pressure at saturation [hPa] : WMO, (Goff, 1957)
            ze_sat = 10.**( 10.79574*(1. - ztmp) - 5.028*LOG10(ptak(ji,jj)/rt0)        &
               &    + 1.50475*10.**(-4)*(1. - 10.**(-8.2969*(ptak(ji,jj)/rt0 - 1.)) )  &
               &    + 0.42873*10.**(-3)*(10.**(4.76955*(1. - ztmp)) - 1.) + 0.78614  )
               !
            q_sat0(ji,jj) = reps * ze_sat/( 0.01_wp*pslp(ji,jj) - (1._wp - reps)*ze_sat )   ! 0.01
            !
         END DO
      END DO
      !
   END FUNCTION q_sat0

      FUNCTION rho_air0( ptak, pqa, pslp, R_dry, rctv0 )
      !!-------------------------------------------------------------------------------
      !!                           ***  FUNCTION rho_air  ***
      !!
      !! ** Purpose : compute density of (moist) air using the eq. of state of the atmosphere
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://sourceforge.net/p/aerobulk)
      !!-------------------------------------------------------------------------------
      REAL(wp), DIMENSION( ni, nj), INTENT(in) ::   ptak      ! air temperature             [K]
      REAL(wp), DIMENSION( ni, nj), INTENT(in) ::   pqa       ! air specific humidity   [kg/kg]
      REAL(wp), DIMENSION( ni, nj), INTENT(in) ::   pslp      ! pressure in                [Pa]
      REAL(wp), DIMENSION( ni, nj)             ::   rho_air0  ! density of moist air   [kg/m^3]
      REAL(wp)                    , INTENT(in) ::   R_dry, rctv0
      !!-------------------------------------------------------------------------------
      !
      rho_air0 = pslp / (  R_dry*ptak * ( 1._wp + rctv0*pqa )  )
      !
   END FUNCTION rho_air0

      FUNCTION L_vap0( psst )
      !!---------------------------------------------------------------------------------
      !!                           ***  FUNCTION L_vap  ***
      !!
      !! ** Purpose : Compute the latent heat of vaporization of water from temperature
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://sourceforge.net/p/aerobulk)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(ni,nj)             ::   L_vap0   ! latent heat of vaporization   [J/kg]
      REAL(wp), DIMENSION(ni,nj), INTENT(in) ::   psst   ! water temperature                [K]
      !!----------------------------------------------------------------------------------
      !
      L_vap0 = (  2.501 - 0.00237 * ( psst(:,:) - rt0)  ) * 1.e6
      !
   END FUNCTION L_vap0

   FUNCTION cp_air0( pqa )
      !!-------------------------------------------------------------------------------
      !!                           ***  FUNCTION cp_air  ***
      !!
      !! ** Purpose : Compute specific heat (Cp) of moist air
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://sourceforge.net/p/aerobulk)
      !!-------------------------------------------------------------------------------
      REAL(wp), DIMENSION( ni, nj), INTENT(in) ::   pqa      ! air specific humidity         [kg/kg]
      REAL(wp), DIMENSION( ni, nj)             ::   cp_air0  ! specific heat of moist air   [J/K/kg]
      REAL(wp), PARAMETER ::   Cp_dry = 1005.0       !: Specic heat of dry air, constant pressure      [J/K/kg]
      REAL(wp), PARAMETER ::   Cp_vap = 1860.0       !: Specic heat of water vapor, constant pressure  [J/K/kg]
      !!-------------------------------------------------------------------------------
      !
      Cp_air0 = Cp_dry + Cp_vap * pqa
      !
   END FUNCTION cp_air0


! Marsaglia & Tsang generator for random normals & random exponentials.
! Translated from C by Alan Miller (amiller@bigpond.net.au)

! Marsaglia, G. & Tsang, W.W. (2000) `The ziggurat method for generating
! random variables', J. Statist. Software, v5(8).

! This is an electronic journal which can be downloaded from:
! http://www.jstatsoft.org/v05/i08

! N.B. It is assumed that all integers are 32-bit.
! N.B. The value of M2 has been halved to compensate for the lack of
!      unsigned integers in Fortran.

! Latest version - 1 January 2001

subroutine zigset(jsrseed)
integer, intent(in) :: jsrseed
integer  :: i
!  set the seed
jsr = jsrseed
!  tables for rnor
q = vn*exp(half*dn*dn)
kn(0) = (dn/q)*m1
kn(1) = 0
wn(0) = q/m1
wn(127) = dn/m1
fn(0) = 1.0_dp
fn(127) = exp(-half*dn*dn)
do  i = 126, 1, -1
   dn = sqrt(-2.0_dp * log(vn/dn + exp(-half*dn*dn)))
   kn(i+1) = (dn/tn)*m1
   tn = dn
   fn(i) = exp(-half*dn*dn)
   wn(i) = dn/m1
end do
!  tables for rexp
q = ve*exp(de)
ke(0) = (de/q)*m2
ke(1) = 0
we(0) = q/m2
we(255) = de/m2
fe(0) = 1.0_dp
fe(255) = exp(-de)
do  i = 254, 1, -1
   de = -log(ve/de + exp(-de))
   ke(i+1) = m2 * (de/te)
   te = de
   fe(i) = exp(-de)
   we(i) = de/m2
end do
initialized = .true.
end subroutine zigset

!  generate random 32-bit integers
function shr3() result(ival)
integer  ::  ival
jz = jsr
jsr = ieor(jsr, ishft(jsr,  13))
jsr = ieor(jsr, ishft(jsr, -17))
jsr = ieor(jsr, ishft(jsr,   5))
ival = jz + jsr
end function shr3

!  generate uniformly distributed random numbers
function uni_scalar() result(ran)
real(kind=dp)  ::  ran
ran = half + 0.2328306e-9_dp * shr3()
end function uni_scalar
!
function uni_vec(n) result(ran)
integer, intent(in) :: n
real(kind=dp)       :: ran(n)
integer             :: i
do i=1,n
   ran(i) = half + 0.2328306e-9_dp * shr3()   
end do
end function uni_vec
!
function uni_vec_given_seed(n,seed) result(ran)
! return n uniform variates
integer, intent(in) :: n
integer, intent(in) :: seed
real(kind=dp)       :: ran(n)
integer             :: i
call zigset(seed)
do i=1,n
   ran(i) = half + 0.2328306e-9_dp * shr3()   
end do
end function uni_vec_given_seed
!
function rnor_scalar() result(fn_val)
!  generate random normal variate 
real(kind=dp)             ::  fn_val
real(kind=dp), parameter  ::  r = 3.442620_dp
real(kind=dp)             ::  x, y
if (.not. initialized) call zigset(jsr)
hz = shr3()
iz = iand(hz, 127)
if (abs(hz) < kn(iz)) then
   fn_val = hz * wn(iz)
else
   do
      if (iz == 0) then
         do
            x = -0.2904764_dp* log(uni())
            y = -log(uni())
            if (y+y >= x*x) exit
         end do
         fn_val = r+x
         if (hz <= 0) fn_val = -fn_val
         return
      end if
      x = hz * wn(iz)
      if (fn(iz) + uni()*(fn(iz-1)-fn(iz)) < exp(-half*x*x)) then
         fn_val = x
         return
      end if
      hz = shr3()
      iz = iand(hz, 127)
      if (abs(hz) < kn(iz)) then
         fn_val = hz * wn(iz)
         return
      end if
   end do
end if
end function rnor_scalar
!
function rnor_vec(n) result(ran)
! generate n random normal variates
integer, intent(in) :: n
real(kind=dp)       :: ran(n)
integer             :: i
do i=1,n
   ran(i) = rnor_scalar()
end do
end function rnor_vec
!
function rexp() result(fn_val)
!  generate random exponential variate
real(kind=dp)  ::  fn_val
real(kind=dp)  ::  x
if (.not. initialized) call zigset(jsr)
jz = shr3()
iz = iand(jz, 255)
if (abs(jz) < ke(iz)) then
   fn_val = abs(jz) * we(iz)
   return
end if
do
   if (iz == 0) then
      fn_val = 7.69711 - log(uni())
      return
   end if
   x = abs(jz) * we(iz)
   if (fe(iz) + uni()*(fe(iz-1) - fe(iz)) < exp(-x)) then
      fn_val = x
      return
   end if
   jz = shr3()
   iz = iand(jz, 255)
   if (abs(jz) < ke(iz)) then
      fn_val = abs(jz) * we(iz)
      return
   end if
end do
end function rexp

END MODULE CGSTOCH
