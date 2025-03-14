!-- SETUP --
#define NEMO
#undef  MAIN
!-----------
#ifdef MAIN
#undef NEMO
#endif
!-----------

MODULE PERLIN

!!
!! ------
!!
!! PERLIN NOISE FORTRAN IMPLEMENTATION
!! 
!! This module provides function and suborutines
!! to generate multiple Perlin noise fields,
!! for use inside and outside NEMO.
!! It is specifically conceived for stochastic physics
!! packages, with an independent perturbation steup/field
!! for each parameter/tendency.
!!
!! Andrea Storto @ CNR.it  -  Aug 2023
!!
!! For NEMO, the module needs the NEMO cpp macro to be defined
!!
!! Typical Usage within NEMO:
!!
!!       ! Initialization
!!       IF( kt .eq. nit000 ) THEN
!!           ! PNOISE_START requires the number of random fields to generate, 
!!           ! and the local domain size, and optionally the random generator
!!           ! seed for ziggurat (to be different for each ensmember)
!!           ! and the logical to abort if restart fields are not available
!!           CALL PNOISE_START(3,jpi,jpj,ksd,.false.)
!!           !  PNOISE_INIT(K1,K2,N,M,LP,NOC,NLA,RPE,RTA,INMD,RDC)
!!           !  Where K1,K2 are the first to last random field with the same
!!           !  parameters and the other refer to the Perlin noise parameters 
!!           !  (see below for details)
!!           !  Here we group two fields (with same parameters)
!!           !  and a third one with different ones
!!           CALL PNOISE_INIT(1,2,20,20,(/.TRUE.,.FALSE./),2,2,0.5_wp,12._wp,1,50000._wp)
!!           CALL PNOISE_INIT(3,3,16,16,(/.TRUE.,.FALSE./),1,2,0.5_wp,12._wp,1,250000._wp)
!!       ENDIF
!!
!!       ! Noise generation at each timestep / update
!!       CALL PNOISE_ADVANCE()
!!       z2d(:,:) = PNOISE(1)%NOISE(:,:) * tmask(:,:,1)
!!
!!       ! End of NEMO example
!!
!!       PARAMETERS OF THE PERLIN NOISE TO PASS TO PNOISE_INIT
!!
!!       N,M  : Size of the Perlin grid
!!       LP   : Logical(2) for cyclic conditions in x,y
!!       NOC  : Number of octaves
!!       NLA  : Lacunarity
!!       RPE  : Persistency
!!       RTA  : Tau (time decorrelation scale for AR(1)
!!       INMD : Initialization (1=from scracth; 2=from restart file)
!!       RDC  : Distance to coast decay (m) if greater than 1km
!!
!!       Further to PNOISE_STARTPNOISE_START, PNOISE_INIT, PNOISE_ADVANCE,
!!       other useful subroutines are
!!       PERLIN_NOISE : low-level noise generation routine
!!       PNOISE_RESTART_WRITE : to write noise fields into NEMO restarts
!!
!! ------
!!

#ifdef NEMO
USE par_kind
USE DOM_OCE
USE LIB_MPP
USE in_out_manager
USE iom
USE ioipsl
#endif

USE ziggurat

IMPLICIT NONE

PRIVATE

PUBLIC PERLIN_NOISE, PNOISE_START, PNOISE_INIT, PNOISE_ADVANCE, PNOISE_RESTART_WRITE, PERR

#ifndef NEMO
PUBLIC WRITE_PN, WRITE_PN3
#endif

! For use outside NEMO
#ifndef NEMO
INTEGER, PARAMETER :: WP=8
LOGICAL, PARAMETER :: lwp=.TRUE.
INTEGER, PARAMETER :: mpprank = 0
INTEGER, PARAMETER :: numout  = 6
#endif

! Compound/Derived type for Perlin Noise
TYPE PNOISE_TYPE
        INTEGER  :: SZ_X, SZ_Y, NX, NY
        INTEGER  :: NN_OCTAVES,NN_LACUNARITY
        INTEGER  :: NN_INITMODE
        REAL(wp) :: RN_PERSISTENCY
        LOGICAL  :: LN_PER(2)
        REAL(wp) :: RN_TAU
        REAL(wp) :: RN_WA
        REAL(wp) :: RN_WB
        REAL(wp) :: RN_DCOAST
        REAL(wp), ALLOCATABLE  :: PNOUT(:,:)
        REAL(wp), ALLOCATABLE  :: NOISE(:,:)
        REAL(wp), ALLOCATABLE  :: MASK (:,:)
END TYPE
PUBLIC PNOISE_TYPE

! Compound/Derived type for Perlin Noise setup
TYPE PNOISE_SETUP
        INTEGER  :: SZ_X
        INTEGER  :: SZ_Y
        LOGICAL  :: LN_PERX
        LOGICAL  :: LN_PERY
        INTEGER  :: NN_OCTAVES
        INTEGER  :: NN_LACUNARITY
        REAL(wp) :: RN_PERSISTENCY
END TYPE
PUBLIC PNOISE_SETUP

! Compound/Derived variable
TYPE(PNOISE_TYPE), ALLOCATABLE, PUBLIC :: PNOISE(:)

INTEGER, SAVE :: PNSIZE= 0
INTEGER, SAVE :: PNTCNT=-1
INTEGER, SAVE :: KNOUT=-1

LOGICAL, PARAMETER :: ln_debug = .FALSE.
LOGICAL, SAVE      :: ln_strict_restart = .FALSE.
LOGICAL, PUBLIC    :: ln_pnoise_adv_init =.FALSE.

! Error procedures
INTERFACE PERR
  MODULE PROCEDURE PERR1, PERR2, PERR3, PERR2B
END INTERFACE

CONTAINS

SUBROUTINE PNOISE_START(NI,KX,KY,KSEED,LNSTRRS,KKOUT)
IMPLICIT NONE
INTEGER, INTENT(IN) :: NI,KX,KY
INTEGER, OPTIONAL   :: KSEED
LOGICAL, OPTIONAL   :: LNSTRRS
INTEGER :: JI, NSEED
INTEGER, OPTIONAL :: KKOUT
ALLOCATE( PNOISE(NI) )
PNSIZE = NI
IF(PRESENT(KKOUT)) THEN
        KNOUT=KKOUT
ELSE
        KNOUT=-1
ENDIF
IF(KNOUT>0) WRITE(KNOUT,*) ' PNOISE :: Start with fields=',PNSIZE

DO JI=1,NI
   PNOISE(JI)%NX=KX
   PNOISE(JI)%NY=KY
   PNOISE(JI)%RN_TAU=-1._wp
   ALLOCATE ( PNOISE(JI)%PNOUT( KX, KY ) )
   ALLOCATE ( PNOISE(JI)%NOISE( KX, KY ) )
   ALLOCATE ( PNOISE(JI)%MASK ( KX, KY ) )
   PNOISE(JI)%PNOUT(:,:) = 0._wp
   PNOISE(JI)%NOISE(:,:) = 0._wp
   PNOISE(JI)%MASK (:,:) = 0._wp
ENDDO
IF(PRESENT(KSEED)) THEN
    NSEED = KSEED
ELSE
    CALL SYSTEM_CLOCK(NSEED)
ENDIF
CALL ZIGSET(NSEED)
IF(PRESENT(LNSTRRS)) ln_strict_restart = LNSTRRS
END SUBROUTINE PNOISE_START

SUBROUTINE PNOISE_INIT(K1,K2,N,M,LP,NOC,NLA,RPE,RTA,INMD,RNDC)
IMPLICIT NONE
INTEGER, INTENT(IN) :: K1,K2,N,M,NOC,NLA,INMD
LOGICAL, INTENT(IN) :: LP(2)
REAL(WP),INTENT(IN) :: RPE,RTA,RNDC
INTEGER :: JI
CHARACTER(LEN=11)   :: CRVAR
INTEGER             :: id1
!
IF( PNSIZE .EQ. 0) CALL PERR('PNOISE_INIT: PNOISE_START has not been called correctly')
PNTCNT = 0
!
IF( K1 .GT. PNSIZE .OR. K2 .GT. PNSIZE) THEN
        CALL PERR('PNOISE_INIT, K1 or K2 exceed the size of fields',K1,K2)
ENDIF
!
IF( NOC .LT. 1 .OR. NOC .GT. 5) THEN
        CALL PERR('PNOISE_INIT, Wrong number of octaves',NOC)
ENDIF
!
IF( NLA .LT. 1 .OR. NLA .GT. 5) THEN
        CALL PERR('PNOISE_INIT, Wrong lacunarity',NLA)
ENDIF
!
IF( RTA .LT. 1._wp .OR. RTA .GT. 1.E+6_wp ) THEN
        CALL PERR('PNOISE_INIT, Incorrect TAU',RTA)
ENDIF
!
IF( INMD .NE.1 .AND. INMD .NE. 2 ) THEN
        CALL PERR('PNOISE_INIT, INMD not supported',INMD)
ENDIF
!
IF( RPE .LT. 1.e-9_wp .OR. RPE .GE. 1._wp ) THEN
        CALL PERR('PNOISE_INIT, Incorrect RPE',RPE)
ENDIF
!
IF( RNDC .LT. 0._wp .OR. RNDC .GT. 1.E+6_wp ) THEN
        CALL PERR('PNOISE_INIT, Incorrect RNDC',RNDC)
ENDIF
!
IF(KNOUT>0) WRITE(KNOUT,*) ' PNOISE :: Init fields=',K1,K2
!
DO JI=K1,K2
   PNOISE(JI)%SZ_X = N
   PNOISE(JI)%SZ_Y = M
   PNOISE(JI)%NN_OCTAVES = NOC
   PNOISE(JI)%NN_LACUNARITY = NLA
   PNOISE(JI)%NN_INITMODE = INMD
   PNOISE(JI)%LN_PER = LP
   PNOISE(JI)%RN_PERSISTENCY = RPE
   PNOISE(JI)%RN_TAU = RTA
   PNOISE(JI)%RN_DCOAST = RNDC
   PNOISE(JI)%RN_WA  = EXP( -1._wp / RTA )
   PNOISE(JI)%RN_WB  = SQRT( 1._wp - PNOISE(JI)%RN_WA*PNOISE(JI)%RN_WA )

   IF ( PNOISE(JI)%NN_INITMODE .EQ. 2 ) THEN
#ifdef NEMO
        WRITE(CRVAR,'(A,I4.4)') 'pnoise_',JI
        id1 = iom_varid( numror, TRIM(CRVAR), ldstop = .FALSE. )
        IF( id1 > 0 ) THEN
            CALL iom_get( numror, jpdom_autoglo, TRIM(CRVAR), PNOISE(JI)%NOISE(:,:))
        ELSE
                IF( ln_strict_restart ) CALL PERR('Restart field for var '//TRIM(CRVAR)//' not found')
        ENDIF
#else
        CALL PERR('Restart mode only in NEMO')
#endif
        PNOISE(JI)%NN_INITMODE = 0
   ENDIF
ENDDO
END SUBROUTINE PNOISE_INIT

SUBROUTINE PNOISE_RESTART_WRITE( kt, cdrw )
  INTEGER         , INTENT(in) ::   kt         ! ocean time-step
  CHARACTER(len=*), INTENT(in) ::   cdrw       ! "READ"/"WRITE" flag
#ifdef NEMO
  CHARACTER(LEN=11)   :: CRVAR
  INTEGER             :: JI
  IF( TRIM(cdrw) == 'WRITE' ) THEN   ! Create restart file
     DO JI=1,PNSIZE
        WRITE(CRVAR,'(A,I4.4)') 'pnoise_',JI
        CALL iom_rstput( kt, nitrst, numrow, CRVAR, PNOISE(JI)%NOISE(:,:) )
     ENDDO
  ENDIF
#endif
END SUBROUTINE PNOISE_RESTART_WRITE

SUBROUTINE PNOISE_ADVANCE(nn_write)
IMPLICIT NONE
INTEGER :: JI, INUM
REAL(WP):: ZSUM, ZCNT
REAL(WP), ALLOCATABLE :: ZDC(:,:,:)
REAL(WP), ALLOCATABLE :: PN3(:,:,:)
CHARACTER(LEN=99) :: CFOUT
INTEGER, OPTIONAL :: nn_write
LOGICAL :: ln_write=.false., ln_iomput=.false.

IF( PNTCNT .LT. 0) CALL PERR('PNOISE_ADVANCE: PNOISE_INIT has not been called correctly')

IF(.NOT.ln_pnoise_adv_init) ln_pnoise_adv_init = .TRUE.

IF( PNTCNT .EQ. 0) THEN
#ifdef NEMO
    IF( ANY(PNOISE(1:PNSIZE)%RN_DCOAST .GT. 1000._wp) ) THEN
            ALLOCATE( ZDC(jpi,jpj,jpk) ); ZDC=1._wp
            CALL iom_open('dist.coast', inum )
            CALL iom_get(inum,jpdom_data,'Tcoast',zdc)
            CALL iom_close( inum )
            DO JI=1,PNSIZE
                 IF( PNOISE(JI)%RN_DCOAST .GT. 1000._wp ) THEN
                     PNOISE(JI)%MASK(:,:) = 1._wp - exp(-zdc(:,:,1)/PNOISE(JI)%RN_DCOAST )
                 ELSE
                     PNOISE(JI)%MASK(:,:) = TMASK(:,:,1)
                 ENDIF
            ENDDO
            DEALLOCATE ( ZDC )
    ELSE
            DO JI=1,PNSIZE
                 PNOISE(JI)%MASK(:,:) = TMASK(:,:,1)
            ENDDO
    ENDIF
#else
            DO JI=1,PNSIZE
                 PNOISE(JI)%MASK(:,:) = 1._wp
            ENDDO
#endif
ENDIF

PNTCNT = PNTCNT +1

DO JI=1,PNSIZE

  CALL PERLIN_NOISE(PNOISE(JI)%SZ_X,PNOISE(JI)%SZ_Y,PNOISE(JI)%NX,PNOISE(JI)%NY,&
     & PNOISE(JI)%LN_PER,PNOISE(JI)%PNOUT,PNOISE(JI)%NN_OCTAVES,&
     & PNOISE(JI)%NN_LACUNARITY,PNOISE(JI)%RN_PERSISTENCY )

  ZSUM =  SUM( PNOISE(JI)%PNOUT, MASK=( PNOISE(JI)%MASK(:,:) .GT. 0._wp ))
  ZCNT =  REAL( COUNT( (PNOISE(JI)%MASK(:,:) .GT. 0._wp )), WP)
#ifdef key_mpp_mpi
  CALL MPP_SUM ( 'pnoise', ZSUM )
  CALL MPP_SUM ( 'pnoise', ZCNT )
#endif  
  ZSUM = ZSUM / ZCNT

  ZSUM =  SUM( (PNOISE(JI)%PNOUT-ZSUM)*(PNOISE(JI)%PNOUT-ZSUM), MASK=( PNOISE(JI)%MASK(:,:) .GT. 0._wp ))
#ifdef key_mpp_mpi
  CALL MPP_SUM ( 'pnoise', ZSUM )
#endif 
  ZSUM = ZSUM / ZCNT

  PNOISE(JI)%PNOUT = PNOISE(JI)%MASK(:,:)*PNOISE(JI)%PNOUT
  PNOISE(JI)%PNOUT = PNOISE(JI)%PNOUT / SQRT( ZSUM )

  IF( ln_debug .AND. lwp .AND. JI .EQ. 1 ) WRITE(numout,*) ' >1> ',PNTCNT,PNOISE(JI)%NN_INITMODE

  IF( PNOISE(JI)%NN_INITMODE .EQ. 0 ) THEN ! Already initialized

      IF( ln_debug .AND. lwp .AND. JI .EQ. 1 ) WRITE(numout,*) ' >2> ', PNOISE(JI)%RN_TAU,PNOISE(JI)%RN_WA,PNOISE(JI)%RN_WB

      IF( PNOISE(JI)%RN_TAU .GT. 0._WP ) THEN
          PNOISE(JI)%NOISE(:,:) = PNOISE(JI)%NOISE(:,:)*PNOISE(JI)%RN_WA + PNOISE(JI)%PNOUT*PNOISE(JI)%RN_WB
      ELSE
          PNOISE(JI)%NOISE(:,:) = PNOISE(JI)%PNOUT
      ENDIF 

      IF( ln_debug .AND. lwp .AND. JI .EQ. 1 ) WRITE(numout,*) '>3> ',SUM( PNOISE(JI)%PNOUT*PNOISE(JI)%PNOUT),&
                                                        SUM( PNOISE(JI)%NOISE*PNOISE(JI)%NOISE)
  ELSEIF( PNOISE(JI)%NN_INITMODE .EQ. 1 ) THEN ! From scratch

      PNOISE(JI)%NOISE(:,:) = PNOISE(JI)%PNOUT
      PNOISE(JI)%NN_INITMODE = 0

  ENDIF 

ENDDO

IF( PRESENT(nn_write) ) THEN
 IF(nn_write>0) ln_iomput = .TRUE.
 IF(nn_write>1) ln_write = .TRUE.
ENDIF

 IF( ln_iomput ) THEN
   IF(ln_write) ALLOCATE(PN3(PNOISE(1)%NX,PNOISE(1)%NY,PNSIZE)) 
   DO JI=1,PNSIZE
#ifdef NEMO
      WRITE(CFOUT,'(A,I2.2,A)') 'pnoise_',JI
      CALL iom_put(TRIM(CFOUT),PNOISE(JI)%NOISE(:,:))
#endif
      IF(ln_write) PN3(:,:,JI) = PNOISE(JI)%NOISE(:,:)
      IF(ln_write.AND.KNOUT>0) WRITE(KNOUT,*) ' PNOISE :: Field ',JI,'=',&
              & MINVAL(PN3(:,:,JI)),MAXVAL(PN3(:,:,JI))
   ENDDO
   IF(ln_write) THEN
#ifdef key_mpp_mpi
      WRITE(CFOUT,'(A,I8.8,A,I4.4,A)') 'pnoise_out_',PNTCNT,'_',narea,'.nc'
#else
      WRITE(CFOUT,'(A,I8.8,A)') 'pnoise_out_',PNTCNT,'.nc'
#endif
      CALL WRITE_PN3(TRIM(CFOUT),PNOISE(1)%NX,PNOISE(1)%NY,PNSIZE,PN3)
      DEALLOCATE(PN3)
   ENDIF
 ENDIF

END SUBROUTINE PNOISE_ADVANCE

SUBROUTINE PERLIN_NOISE_INT(N, M, NN, MM, LPER, PN)
IMPLICIT NONE

INTEGER, INTENT(IN) :: N, M, NN, MM
LOGICAL, INTENT(IN) :: LPER(2)
REAL(WP),   INTENT(OUT) :: PN(NN,MM)

INTEGER :: I, J, K, KK, IERR
INTEGER, PARAMETER :: KROOT=0
REAL(WP) :: AR(2,N,M), ZZ(2*N*M), DX, DY
#ifdef NEMO
REAL(WP) ::  XS(JPIGLO), YS(JPJGLO)
#ifdef key_mpp_mpi
INCLUDE 'mpif.h'
#endif
#else
REAL(WP) ::  XS(NN), YS(MM)
#endif

IF( mpprank .EQ. 0 ) THEN
#ifdef FORTRAN_RANDOM
    CALL random_stdnormalv( ZZ(:) ) ! faster so far
#else
    DO J=1,2*n*m
      ZZ(J)= rnor()
    ENDDO
#endif
ENDIF
#ifdef key_mpp_mpi
      CALL mpi_bcast( ZZ, 2*N*M , mpi_double,  &
      &               kroot, mpi_comm_oce, ierr )
#endif

KK=0
DO K=1,M
  DO J=1,N
    DO I=1,2
       KK=KK+1
       AR(I,J,K) = ZZ(KK)
    ENDDO
  ENDDO
ENDDO

AR = AR / SQRT( SUM( AR*AR ) )

IF(LPER(1)) AR(:,N,:) = AR(:,1,:)
IF(LPER(2)) AR(:,:,M) = AR(:,:,1)

#ifdef NEMO

DX = REAL(N-1,WP) / REAL(JPIGLO,WP)
DY = REAL(M-1,WP) / REAL(JPJGLO,WP)

DO K=1,JPIGLO
   XS(K)=1._WP+(K-1)*DX
ENDDO
DO K=1,JPJGLO
   YS(K)=1._WP+(K-1)*DY
ENDDO

DO J=1,JPJ
   DO I=1,JPI
       PN(I,J) = FF(N,M,XS(MIG(I)),YS(MJG(J)),AR)
   ENDDO
ENDDO

#else

DX = REAL(N-1) / REAL(NN)
DY = REAL(M-1) / REAL(MM)

DO K=1,NN
   XS(K)=1._WP+(K-1)*DX
ENDDO
DO K=1,MM
   YS(K)=1._WP+(K-1)*DY
ENDDO

DO J=1,MM
   DO I=1,NN
       PN(I,J) = FF(N,M,XS(I),YS(J),AR)
   ENDDO
ENDDO

#endif

END SUBROUTINE PERLIN_NOISE_INT

REAL FUNCTION FF(K,L,X,Y,VF)
IMPLICIT NONE

INTEGER,INTENT(IN) :: K,L
REAL(WP), INTENT(IN)   :: X,Y,VF(2,K,L)
INTEGER            :: I,J
REAL(WP), DIMENSION(2) :: UU, V1, V2, V3, V4
REAL(WP), DIMENSION(2) :: U1, U2, U3, U4
REAL(WP)           :: A1, A2, A3, A4, P, Q, B1, B2

I = FLOOR(X)
J = FLOOR(Y)

IF(I<1 .OR. J<1 .OR. I>K .OR. J>L) THEN
   IF(lwp) WRITE(numout,*) 'WRONG INDEX',I,J,K,L
   CALL PERR("WRONG INDEXING")
ENDIF

UU = (/X, Y/)

V1(:) = VF(:,I  ,J  )
V2(:) = VF(:,I+1,J  )
V3(:) = VF(:,I  ,J+1)
V4(:) = VF(:,I+1,J+1)

U1(:) = UU(:) - (/I  , J  /)
U2(:) = UU(:) - (/I+1, J  /)
U3(:) = UU(:) - (/I  , J+1/)
U4(:) = UU(:) - (/I+1, J+1/)

A1    = SUM( V1 * U1 )
A2    = SUM( V2 * U2 )
A3    = SUM( V3 * U3 )
A4    = SUM( V4 * U4 )

P     = SS( X - I)
Q     = SS( Y - J)

B1    = (1._WP-P)*A1+P*A2
B2    = (1._WP-P)*A3+P*A4

FF    = (1._WP-Q)*B1 + Q*B2

END FUNCTION FF

REAL FUNCTION SS ( PP )
  IMPLICIT NONE
  REAL(WP) :: PP
  SS = 3*PP*PP - 2*PP*PP*PP
END FUNCTION SS

SUBROUTINE PERLIN_NOISE(KX, KY, NX, NY, LL_PER, PNOISE, NOCT, LACUNARITY, PERSISTENCE)
IMPLICIT NONE

INTEGER, INTENT(IN) :: KX, KY, NX, NY, NOCT
LOGICAL, INTENT(IN) :: LL_PER(2)
INTEGER, INTENT(IN) :: LACUNARITY
REAL(WP),INTENT(IN) :: PERSISTENCE
REAL(WP),INTENT(OUT):: PNOISE(NX, NY)
!
REAL(WP)            :: PNOIST(NX, NY)
REAL(WP)            :: AM
INTEGER             :: JK, FR

! Default to persistence=0.5, lacunarity=2

IF( NOCT .LT.1 .OR. NOCT .GT.8 )  CALL PERR('NUMBER OF OCTAVES NOT SUPPORTED',NOCT)
IF( LACUNARITY.LT. 1. )           CALL PERR('LACUNARITY        NOT SUPPORTED',LACUNARITY)
IF( PERSISTENCE .GT.1.)           CALL PERR('PERSISTENCE       NOT SUPPORTED',PERSISTENCE)

PNOISE = 0.

FR=1
AM=1.
DO JK=1, NOCT
       CALL PERLIN_NOISE_INT(KX*FR, KY*FR, NX, NY, LL_PER, PNOIST)
       PNOISE = PNOISE + AM*PNOIST
       FR = FR * LACUNARITY
       AM = AM * PERSISTENCE
ENDDO

END SUBROUTINE PERLIN_NOISE

SUBROUTINE PERR1(CC)
IMPLICIT NONE
CHARACTER(LEN=*) :: CC
WRITE(0,'(X,A)') TRIM(CC)
CALL MYABORT()
END SUBROUTINE PERR1

SUBROUTINE PERR2(CC,II)
IMPLICIT NONE
CHARACTER(LEN=*) :: CC
INTEGER          :: II
WRITE(0,'(X,A,X,I9)') TRIM(CC),II
CALL MYABORT()
END SUBROUTINE PERR2

SUBROUTINE PERR2B(CC,II,JJ)
IMPLICIT NONE
CHARACTER(LEN=*) :: CC
INTEGER          :: II,JJ
WRITE(0,'(X,A,X,I9,X,I9)') TRIM(CC),II,JJ
CALL MYABORT()
END SUBROUTINE PERR2B

SUBROUTINE PERR3(CC,ZZ)
IMPLICIT NONE
CHARACTER(LEN=*) :: CC
REAL(WP)         :: ZZ
WRITE(0,'(X,A,X,F12.3)') TRIM(CC),ZZ
CALL MYABORT()
END SUBROUTINE PERR3

SUBROUTINE MYABORT
#ifdef NEMO
CALL ctl_stop( 'STOP', 'Perlin module error' )
#else
CALL ABORT()
#endif
END SUBROUTINE MYABORT

SUBROUTINE WRITE_PN(CF,IM,JM,PN)
USE NETCDF
IMPLICIT NONE

CHARACTER(LEN=*) :: CF
INTEGER, INTENT(IN) :: IM, JM
REAL(WP), INTENT(IN) :: PN(IM,JM)

INTEGER :: NCID, X_DIMID, Y_DIMID, DIMIDS1(2),VARID

  CALL CHECK( NF90_CREATE(TRIM(CF), NF90_CLOBBER, NCID) )
  CALL CHECK( NF90_DEF_DIM(NCID, 'x', IM, X_DIMID) )
  CALL CHECK( NF90_DEF_DIM(NCID, 'y', JM, Y_DIMID) )
  DIMIDS1=  (/X_DIMID, Y_DIMID/)
  CALL CHECK( NF90_DEF_VAR(NCID, 'noise',NF90_REAL,DIMIDS1,VARID) )
  CALL CHECK( NF90_ENDDEF(NCID) )
  CALL CHECK( NF90_PUT_VAR(NCID, VARID, PN ))
  CALL CHECK( NF90_CLOSE(NCID) )

END SUBROUTINE WRITE_PN

SUBROUTINE WRITE_PN3(CF,IM,JM,KM,PN)
USE NETCDF
IMPLICIT NONE

CHARACTER(LEN=*) :: CF
INTEGER, INTENT(IN) :: IM, JM, KM
REAL(WP), INTENT(IN) :: PN(IM,JM,KM)

INTEGER :: NCID, X_DIMID, Y_DIMID, K_DIMID, DIMIDS1(3),VARID

  CALL CHECK( NF90_CREATE(TRIM(CF), NF90_CLOBBER, NCID) )
  CALL CHECK( NF90_DEF_DIM(NCID, 'x', IM, X_DIMID) )
  CALL CHECK( NF90_DEF_DIM(NCID, 'y', JM, Y_DIMID) )
  CALL CHECK( NF90_DEF_DIM(NCID, 'p', KM, K_DIMID) )
  DIMIDS1=  (/X_DIMID, Y_DIMID, K_DIMID/)
  CALL CHECK( NF90_DEF_VAR(NCID, 'noise',NF90_REAL,DIMIDS1,VARID) )
  CALL CHECK( NF90_ENDDEF(NCID) )
  CALL CHECK( NF90_PUT_VAR(NCID, VARID, PN ))
  CALL CHECK( NF90_CLOSE(NCID) )

END SUBROUTINE WRITE_PN3

SUBROUTINE CHECK(STATUS)
USE NETCDF
IMPLICIT NONE
INTEGER, INTENT ( IN) :: STATUS

  IF(STATUS /= NF90_NOERR) THEN
     CALL PERR(TRIM(NF90_STRERROR(STATUS)))
  END IF

END SUBROUTINE CHECK

END MODULE PERLIN

!!!!!!
!
! Standalone test is enabled with the MAIN cpp macro
! and without the NEMO cpp macro
!
!!!!!!

#ifdef MAIN

PROGRAM TEST

USE PERLIN
USE NETCDF
IMPLICIT NONE
INTEGER, PARAMETER :: WP=8
INTEGER,PARAMETER :: JPI=362, JPJ=292
INTEGER ::  KI=10 ,  KJ=10
LOGICAL :: LP(2)=(/.TRUE.,.TRUE./)
INTEGER :: NOCT=2, LAC=2
REAL(WP):: PERS=0.5
REAL(WP):: PN(JPI,JPJ)

! Test low-level function
CALL PERLIN_NOISE(KI,KJ,JPI,JPJ,LP,PN,NOCT,LAC,PERS)
CALL WRITE_PN('test_ll.nc',JPI,JPJ,PN)

! Test APIs
CALL PNOISE_START(1,jpi,jpj,123456,.false.,6)
CALL PNOISE_INIT(1,1,KI,KJ,LP,NOCT,LAC,PERS,21.0_wp,1,200000._wp)
CALL PNOISE_ADVANCE()
CALL WRITE_PN('test_api.nc',JPI,JPJ,PNOISE(1)%NOISE(:,:))

END PROGRAM TEST
#endif
