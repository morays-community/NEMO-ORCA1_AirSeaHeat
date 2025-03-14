MODULE ANNIF

#define wp 8

! Artifical Neural Networks Inference in Fortran (ANNIF)
! ------------------------------------------------------
!
! Bridging module for online inference in Fortran90
! Including TL and AD operators for data assimilation.
!
! Version 1:    05/JUL/2023
! ---------
!               Support for dense neural networks and
!               "relu"/"tanh"/"sigm" activation (only).
!               Read models in intermediate NetCDF
!               format
!
! Author(s): Andrea.Storto (@ CNR .it)
!

USE NETCDF

IMPLICIT NONE

    PRIVATE

    TYPE ANNIF_TYPE
        CHARACTER(LEN=199) :: cnnmodf
        INTEGER :: numfea, numneu, numout, numwrk, numlay, numhdl
        REAL(wp), ALLOCATABLE  :: wei(:,:,:), bia(:,:)
        REAL(wp), ALLOCATABLE  :: rin(:,:), rou(:,:)
        CHARACTER(LEN=4)   :: cnnact
        INTEGER :: innact
        INTEGER :: iverb  = 0
        INTEGER :: inorma = 0
    END TYPE
    PUBLIC ANNIF_TYPE

    CHARACTER(LEN=10), PARAMETER :: CINST = ' ANNIF >> ', &
                                  & CINSE = ' ANNIF ERR'

    INTEGER, SAVE   :: IOUNO = 6
    INTEGER, SAVE   :: IOUNE = 0

    PUBLIC ANNIF_RDMOD, ANNIF_INIT, ANNIF_FWD, ANNIF_TL, ANNIF_AD, ANNIF_END, &
           ANNIF_TESTAD, ANNIF_TESTTL, ANNIF_NORMALIZ

    LOGICAL, SAVE :: LLIRS = .FALSE.

    INTERFACE POPREAL
        MODULE PROCEDURE POPREAL_4, POPREAL_8
    END INTERFACE
    INTERFACE ANNIF_FWD
        MODULE PROCEDURE ANNIF_FWD_1, ANNIF_FWD_2
    END INTERFACE
    INTERFACE PUSHREAL
        MODULE PROCEDURE PUSHREAL_4, PUSHREAL_8
    END INTERFACE

CONTAINS

SUBROUTINE ANNIF_RDMOD( IMOD, CFILE )

!  Reading NN model from NetCDF
!  ---

IMPLICIT NONE

        TYPE(ANNIF_TYPE), INTENT(INOUT) :: IMOD
        CHARACTER(LEN=*) :: CFILE
        CHARACTER(LEN=4) :: CNORMA

        IMOD%cnnmodf = TRIM(CFILE)

        IF( IMOD%iverb .gt. 0 ) THEN
          WRITE(IOUNO,*) 
          WRITE(IOUNO,*) CINST//' Reading file '//TRIM(IMOD%cnnmodf)
          WRITE(IOUNO,*) 
        ENDIF

        CALL GETGATT(IMOD%cnnmodf, "num_features" , IMOD%numfea )
        CALL GETGATT(IMOD%cnnmodf, "num_neurons"  , IMOD%numneu )
        CALL GETGATT(IMOD%cnnmodf, "num_outputs"  , IMOD%numout )
        CALL GETGATT(IMOD%cnnmodf, "num_workdim"  , IMOD%numwrk )
        CALL GETGATT(IMOD%cnnmodf, "num_layers"   , IMOD%numlay )
        CALL GETGATT(IMOD%cnnmodf, "num_hiddenlay", IMOD%numhdl )
        CALL GETGATC(IMOD%cnnmodf, "activation"   , IMOD%cnnact )

        IF( IMOD%numhdl .NE. IMOD%numlay-2 ) THEN
                WRITE(IOUNE,*) CINSE,' num_layers, num_hiddenlay =',IMOD%numhdl,IMOD%numlay
                WRITE(IOUNE,*) CINSE//' num_layers -2 != num_hiddenlay'
                WRITE(IOUNE,*) CINSE//' not supported, aborting'
                CALL ABORT()
        ENDIF

        IF( IMOD%cnnact .NE. 'relu' .AND. IMOD%cnnact .NE. 'none' &
              & .AND. IMOD%cnnact .NE. 'tanh' .AND. IMOD%cnnact .NE. 'sigm' ) THEN
                WRITE(IOUNE,*) CINSE,' activation =',IMOD%cnnact
                WRITE(IOUNE,*) CINSE//' not supported, aborting'
                CALL ABORT()
        ENDIF
        IMOD%innact = 0
        IF( IMOD%cnnact .EQ. 'relu' ) IMOD%innact = 1
        IF( IMOD%cnnact .EQ. 'tanh' ) IMOD%innact = 2
        IF( IMOD%cnnact .EQ. 'sigm' ) IMOD%innact = 3

        CALL GETGATC(IMOD%cnnmodf, "normalization", CNORMA )
        SELECT CASE(CNORMA)
                CASE('none') 
                        IMOD%inorma = 0 ;;
                CASE('unif') 
                        IMOD%inorma = 1 ;;
                CASE('uni0') 
                        IMOD%inorma = 2 ;;
                CASE('norm')
                        IMOD%inorma = 3 ;;
                CASE DEFAULT
                        WRITE(IOUNE,*) CINSE,' normalizaz =',CNORMA
                        WRITE(IOUNE,*) CINSE//' not supported, aborting'
                        CALL ABORT()
        END SELECT

        ALLOCATE( IMOD%wei (IMOD%numwrk, IMOD%numwrk, IMOD%numlay), &
                & IMOD%bia (IMOD%numwrk, IMOD%numlay), &
                & IMOD%rin (IMOD%numfea, 2) , &
                & IMOD%rou (IMOD%numout, 2) )

        CALL GETVAR3(IMOD%cnnmodf,"weights",IMOD%wei )
        CALL GETVAR2(IMOD%cnnmodf,"bias",IMOD%bia )

        CALL GETVAR2(IMOD%cnnmodf,"range_inputs" ,IMOD%rin )
        CALL GETVAR2(IMOD%cnnmodf,"range_outputs",IMOD%rou )

        IF( IMOD%iverb .gt. 1 ) CALL ANNIF_PRINT(IMOD)

END SUBROUTINE ANNIF_RDMOD

SUBROUTINE ANNIF_PRINT(IMOD)

!  Print info 
!  ---

        IMPLICIT NONE

        TYPE(ANNIF_TYPE) :: IMOD

        WRITE(IOUNO,*)
        WRITE(IOUNO,*) CINST//''
        WRITE(IOUNO,*) CINST//' ANN Setup'
        WRITE(IOUNO,*) CINST//' Filename        : ', TRIM(IMOD%cnnmodf)
        WRITE(IOUNO,*) CINST//' No of features  : ', IMOD%numfea
        WRITE(IOUNO,*) CINST//' No of neurons   : ', IMOD%numneu
        WRITE(IOUNO,*) CINST//' No of outputs   : ', IMOD%numout
        WRITE(IOUNO,*) CINST//' No of workdim   : ', IMOD%numwrk
        WRITE(IOUNO,*) CINST//' No of layers    : ', IMOD%numlay
        WRITE(IOUNO,*) CINST//' No of hidden lay: ', IMOD%numhdl
        WRITE(IOUNO,*) CINST//' Activation type : ', IMOD%cnnact
        WRITE(IOUNO,*) CINST//' Activation numb : ', IMOD%innact
        WRITE(IOUNO,*) CINST//''
        WRITE(IOUNO,*)

END SUBROUTINE ANNIF_PRINT

SUBROUTINE ANNIF_END(IMOD)

!  Finalizing NN model
!  ---

        IMPLICIT NONE

        TYPE(ANNIF_TYPE), INTENT(INOUT) :: IMOD(:)
        INTEGER :: NM, JM

        NM = SIZE ( IMOD )
        DO JM=1,NM
           DEALLOCATE( IMOD(JM)%wei, IMOD(JM)%bia )
        ENDDO

END SUBROUTINE ANNIF_END

SUBROUTINE ANNIF_INIT( IMOD, IVRB, IOU, IEU )

!  Initializing NN model
!  ---

        IMPLICIT NONE

        TYPE(ANNIF_TYPE), INTENT(INOUT) :: IMOD(:)
        INTEGER, INTENT(IN) :: IVRB
        INTEGER, OPTIONAL, INTENT(IN) :: IOU, IEU
        INTEGER :: JM, NM

        NM = SIZE ( IMOD )
        DO JM=1,NM
          IMOD(JM)%iverb = MAX(0,MIN(IVRB, 10))
        ENDDO

        IF(PRESENT(IOU )) IOUNO = IOU
        IF(PRESENT(IEU )) IOUNE = IEU

END SUBROUTINE ANNIF_INIT

SUBROUTINE ANNIF_NORMALIZ(IMOD,ZV,IDATA,IMODE)

!  Data normalization
!  ---

        IMPLICIT NONE
        TYPE(ANNIF_TYPE), INTENT(INOUT) :: IMOD
        REAL(wp), INTENT(INOUT) :: ZV(:)
        INTEGER, INTENT(IN)    :: IDATA, IMODE
        INTEGER                :: JI

        IF( IMOD%inorma .EQ. 0 ) RETURN

        IF( IDATA .EQ. 1 ) THEN   ! FEATURES

                IF ( IMODE.EQ.1 ) THEN
                        IF( IMOD%inorma .EQ. 1 ) THEN
                          DO JI=1, IMOD%numfea
                             ZV(JI) = -1. + 2.*(ZV(JI)- IMOD%rin(JI,1))/(IMOD%rin(JI,2)-IMOD%rin(JI,1))
                          ENDDO
                        ELSEIF ( IMOD%inorma .EQ. 2 ) THEN
                          DO JI=1, IMOD%numfea
                             ZV(JI) = (ZV(JI)- IMOD%rin(JI,1))/(IMOD%rin(JI,2)-IMOD%rin(JI,1))
                          ENDDO
                        ELSEIF ( IMOD%inorma .EQ. 3 ) THEN
                          DO JI=1, IMOD%numfea
                             ZV(JI) = (ZV(JI)- IMOD%rin(JI,1))/IMOD%rin(JI,2)
                          ENDDO
                        ENDIF
                ELSE IF ( IMODE.EQ.2 ) THEN
                        IF( IMOD%inorma .EQ. 1 ) THEN
                          DO JI=1, IMOD%numfea
                             ZV(JI) = (ZV(JI)+1.)*(IMOD%rin(JI,2)-IMOD%rin(JI,1))/2. + IMOD%rin(JI,1)
                          ENDDO
                        ELSEIF ( IMOD%inorma .EQ. 2 ) THEN
                          DO JI=1, IMOD%numfea
                             ZV(JI) = ZV(JI)*(IMOD%rin(JI,2)-IMOD%rin(JI,1)) + IMOD%rin(JI,1)
                          ENDDO
                        ELSEIF ( IMOD%inorma .EQ. 3 ) THEN
                          DO JI=1, IMOD%numfea
                             ZV(JI) = IMOD%rin(JI,2)*ZV(JI) + IMOD%rin(JI,1)
                          ENDDO
                        ENDIF
                ELSE
                        WRITE(IOUNE,*) CINSE//' IMODE neither 1 or 2'
                        WRITE(IOUNE,*) CINSE//' not supported, aborting'
                        CALL ABORT()
                ENDIF

         ELSEIF ( IDATA .EQ. 2 ) THEN

                IF ( IMODE.EQ.1 ) THEN
                        IF( IMOD%inorma .EQ. 1 ) THEN
                          DO JI=1, IMOD%numout
                             ZV(JI) = -1. + 2.*(ZV(JI)- IMOD%rou(JI,1))/(IMOD%rou(JI,2)-IMOD%rou(JI,1))
                          ENDDO
                        ELSEIF ( IMOD%inorma .EQ. 2 ) THEN
                          DO JI=1, IMOD%numout
                             ZV(JI) = (ZV(JI)- IMOD%rou(JI,1))/(IMOD%rou(JI,2)-IMOD%rou(JI,1))
                          ENDDO
                        ELSEIF ( IMOD%inorma .EQ. 3 ) THEN
                          DO JI=1, IMOD%numout
                             ZV(JI) = (ZV(JI)- IMOD%rou(JI,1))/IMOD%rou(JI,2)
                          ENDDO
                        ENDIF
                ELSE IF ( IMODE.EQ.2 ) THEN
                        IF( IMOD%inorma .EQ. 1 ) THEN
                          DO JI=1, IMOD%numout
                             ZV(JI) = (ZV(JI)+1.)*(IMOD%rou(JI,2)-IMOD%rou(JI,1))/2. + IMOD%rou(JI,1)
                          ENDDO
                        ELSEIF ( IMOD%inorma .EQ. 2 ) THEN
                          DO JI=1, IMOD%numout
                             ZV(JI) = ZV(JI)*(IMOD%rou(JI,2)-IMOD%rou(JI,1)) + IMOD%rou(JI,1)
                          ENDDO
                        ELSEIF ( IMOD%inorma .EQ. 3 ) THEN
                          DO JI=1, IMOD%numout
                             ZV(JI) = IMOD%rou(JI,2)*ZV(JI) + IMOD%rou(JI,1)
                          ENDDO
                        ENDIF
                ELSE

                        WRITE(IOUNE,*) CINSE//' IMODE neither 1 or 2'
                        WRITE(IOUNE,*) CINSE//' not supported, aborting'
                        CALL ABORT()

                ENDIF

         ELSE

                WRITE(IOUNE,*) CINSE//' IDATA neither 1 or 2'
                WRITE(IOUNE,*) CINSE//' not supported, aborting'
                CALL ABORT()

         ENDIF


END SUBROUTINE ANNIF_NORMALIZ

SUBROUTINE ANNIF_FWD_1(IMOD, INPU, OUTP)

!  NN forward model Wrapper
!  ---

        IMPLICIT NONE
        TYPE(ANNIF_TYPE), INTENT(IN) :: IMOD
        REAL(wp), INTENT(IN)       :: INPU(IMOD%numfea)
        REAL(wp), INTENT(OUT)      :: OUTP(IMOD%numout)
        INTEGER                :: IAC

        IAC=IMOD%innact

        CALL NNINFERW(IMOD%numlay,IMOD%numfea,IMOD%numout,IMOD%numneu,&
        & INPU, OUTP, IMOD%wei, IMOD%bia, IMOD%numwrk, IAC)

END SUBROUTINE ANNIF_FWD_1

SUBROUTINE ANNIF_FWD_2(IMOD, NC, INPU, OUTP)

!  NN forward model Wrapper
!  ---

        IMPLICIT NONE
        TYPE(ANNIF_TYPE), INTENT(IN) :: IMOD
        INTEGER, INTENT(IN) :: NC
        REAL(wp), INTENT(IN)       :: INPU(NC,IMOD%numfea)
        REAL(wp), INTENT(OUT)      :: OUTP(NC,IMOD%numout)
        INTEGER                :: IAC, JC

        IAC=IMOD%innact

        DO JC=1,NC
           CALL NNINFERW(IMOD%numlay,IMOD%numfea,IMOD%numout,IMOD%numneu,&
           & INPU(JC,:), OUTP(JC,:), IMOD%wei, IMOD%bia, IMOD%numwrk, IAC)
        ENDDO

END SUBROUTINE ANNIF_FWD_2

SUBROUTINE ANNIF_TL(IMOD, INPU, INPU_TL, OUTP, OUTP_TL)

!  NN tangent-linear model Wrapper
!  ---

        IMPLICIT NONE
        TYPE(ANNIF_TYPE), INTENT(IN) :: IMOD
        REAL(wp), INTENT(IN)       :: INPU(IMOD%numfea)
        REAL(wp), INTENT(OUT)      :: OUTP(IMOD%numout)
        REAL(wp), INTENT(IN)       :: INPU_TL(IMOD%numfea)
        REAL(wp), INTENT(OUT)      :: OUTP_TL(IMOD%numout)
        INTEGER                :: IAC

        IAC=IMOD%innact

        CALL NNINFER_WD(IMOD%numlay,IMOD%numfea,IMOD%numout,IMOD%numneu,&
        & INPU, INPU_TL, OUTP, OUTP_TL, IMOD%wei, IMOD%bia, IMOD%numwrk, IAC)

END SUBROUTINE ANNIF_TL

SUBROUTINE ANNIF_AD(IMOD, INPU, INPU_AD, OUTP, OUTP_AD)

!  NN adjoint model Wrapper
!  ---

        IMPLICIT NONE
        TYPE(ANNIF_TYPE), INTENT(IN) :: IMOD
        REAL(wp), INTENT(IN)       :: INPU(IMOD%numfea)
        REAL(wp), INTENT(OUT)      :: OUTP(IMOD%numout)
        REAL(wp), INTENT(INOUT)    :: INPU_AD(IMOD%numfea)
        REAL(wp), INTENT(INOUT)    :: OUTP_AD(IMOD%numout)
        INTEGER                :: IAC

        IAC=IMOD%innact
        CALL NNINFER_WB(IMOD%numlay,IMOD%numfea,IMOD%numout,IMOD%numneu,&
        & INPU, INPU_AD, OUTP, OUTP_AD, IMOD%wei, IMOD%bia, IMOD%numwrk, IAC)

END SUBROUTINE ANNIF_AD

!
! Internal routines
!

SUBROUTINE NNINFERW(NUML, NFEA, NOUT, NNEU, INPU, OUTP, WEI, BIA, NWRK, IACT)

!  Forward ANN operator
!  ---

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NUML, NFEA, NOUT, NNEU, NWRK, IACT
      REAL(wp), INTENT(IN)    :: INPU(NFEA)
      REAL(wp), INTENT(OUT)   :: OUTP(NOUT)
      REAL(wp) :: XXW (NWRK)
      REAL(wp) :: XXW2(NWRK)
      REAL(wp) :: WEI(NWRK, NWRK, NUML)
      REAL(wp) :: BIA(      NWRK, NUML)
      INTEGER :: JI, JJ, JK, JL

      ! First layer
      XXW=0.
      DO JJ=1,NNEU
         DO JK=1,NFEA
            XXW(JJ)=XXW(JJ)+INPU(JK)*WEI(JK,JJ,1)
         ENDDO
      ENDDO
      XXW(1:NNEU) =  XXW(1:NNEU) + BIA(1:NNEU,1)
      
      ! RELU activation
      IF(IACT.EQ.1) THEN
        DO JJ=1,NNEU
                 XXW(JJ) = MAX(0., XXW(JJ) )
        ENDDO
      ! TANH activation
      ELSEIF(IACT.EQ.2) THEN
        DO JJ=1,NNEU
                 XXW(JJ) = TANH( XXW(JJ) )
        ENDDO
      ! SIGMOID activation
      ELSEIF(IACT.EQ.3) THEN
        DO JJ=1,NNEU
                 XXW(JJ) = 1./(1.+EXP(-XXW(JJ) ))
        ENDDO
      ENDIF

      ! Hidden layers
      DO JL=2,NUML-1
         XXW2=0.
         DO JJ=1,NNEU
           DO JK=1,NNEU
              XXW2(JJ)=XXW2(JJ)+XXW(JK)*WEI(JK,JJ,JL)
           ENDDO
         ENDDO
         XXW2(1:NNEU) =  XXW2(1:NNEU) + BIA(1:NNEU,JL)
      
         ! RELU activation
         IF(IACT.EQ.1) THEN
           DO JJ=1,NNEU
                 XXW2(JJ) = MAX(0., XXW2(JJ) )
           ENDDO
         ! TANH activation
         ELSEIF(IACT.EQ.2) THEN
           DO JJ=1,NNEU
                 XXW2(JJ) = TANH( XXW2(JJ) )
           ENDDO
         ! SIGMOID activation
         ELSEIF(IACT.EQ.3) THEN
           DO JJ=1,NNEU
                 XXW2(JJ) = 1./(1.+EXP(-XXW2(JJ) ))
           ENDDO
         ENDIF
         XXW(1:NNEU) = XXW2(1:NNEU)
      ENDDO
      
      ! Last layer
      OUTP=0.
      DO JJ=1,NOUT
         DO JK=1,NNEU
            OUTP(JJ) = OUTP(JJ) + XXW(JK)*WEI(JK,JJ,NUML)
         ENDDO
      ENDDO
      OUTP(1:NOUT) =  OUTP(1:NOUT) + BIA(1:NOUT,NUML)

END SUBROUTINE NNINFERW

SUBROUTINE NNINFER_WD(numl, nfea, nout, nneu, inpu, inpud, outp, outpd, wei, bia, nwrk,iact)

!  Tangent-linear ANN operator
!  ---

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: numl, nfea, nout, nneu, nwrk,iact
  REAL(wp), INTENT(IN) :: inpu(nfea)
  REAL(wp), INTENT(IN) :: inpud(nfea)
  REAL(wp), INTENT(OUT) :: outp(nout)
  REAL(wp), INTENT(OUT) :: outpd(nout)
  REAL(wp) :: xxw(nwrk)
  REAL(wp) :: xxwd(nwrk)
  REAL(wp) :: xxw2(nwrk)
  REAL(wp) :: xxw2d(nwrk)
  REAL(wp) :: wei(nwrk, nwrk, numl)
  REAL(wp) :: bia(nwrk, numl)
  REAL(wp) :: zt
  INTEGER :: ji, jj, jk, jl
  INTRINSIC MAX
! First layer
  xxw = 0.
  xxwd = 0.0
  DO jj=1,nneu
    DO jk=1,nfea
      xxwd(jj) = xxwd(jj) + wei(jk, jj, 1)*inpud(jk)
      xxw(jj) = xxw(jj) + inpu(jk)*wei(jk, jj, 1)
    END DO
  END DO
  xxw(1:nneu) = xxw(1:nneu) + bia(1:nneu, 1)
! RELU activation
  IF(IACT.EQ.1) THEN
   DO jj=1,nneu
    IF (0. .LT. xxw(jj)) THEN
      xxw(jj) = xxw(jj)
    ELSE
      xxwd(jj) = 0.
      xxw(jj) = 0.
    END IF
   END DO
! TANH activation
  ELSEIF(IACT.EQ.2) THEN
   DO jj=1,nneu
      xxwd(jj) = (1.0-TANH(xxw(jj))**2)*xxwd(jj)
      xxw(jj) = TANH(xxw(jj))
   END DO
! SIGMOID activation
  ELSEIF(IACT.EQ.3) THEN
   DO jj=1,nneu
      zt = 1./(EXP(-xxw(jj))+1.)
      xxwd(jj) = zt*EXP(-xxw(jj))*xxwd(jj)/(EXP(-xxw(jj))+1.)
      xxw(jj) = zt
   END DO
  ENDIF
! Hidden layers
  DO jl=2,numl-1
    xxw2 = 0.
    xxw2d = 0.
    DO jj=1,nneu
      DO jk=1,nneu
        xxw2d(jj) = xxw2d(jj) + wei(jk, jj, jl)*xxwd(jk)
        xxw2(jj) = xxw2(jj) + xxw(jk)*wei(jk, jj, jl)
      END DO
    END DO
    xxw2(1:nneu) = xxw2(1:nneu) + bia(1:nneu, jl)
! RELU activation
    IF(IACT.EQ.1) THEN
     DO jj=1,nneu
      IF (0. .LT. xxw2(jj)) THEN
        xxw2(jj) = xxw2(jj)
      ELSE
        xxw2d(jj) = 0.
        xxw2(jj) = 0.
      END IF
     END DO
! TANH activation
    ELSEIF(IACT.EQ.2) THEN
     DO jj=1,nneu
        xxw2d(jj) = (1.0-TANH(xxw2(jj))**2)*xxw2d(jj)
        xxw2(jj) = TANH(xxw2(jj))
     END DO
! SIGMOID activation
    ELSEIF(IACT.EQ.3) THEN
     DO jj=1,nneu
        zt = 1./(EXP(-xxw2(jj))+1.)
        xxw2d(jj) = zt*EXP(-xxw2(jj))*xxw2d(jj)/(EXP(-xxw2(jj))+1.)
        xxw2(jj) = zt
     END DO
    ENDIF
    xxwd(1:nneu) = xxw2d(1:nneu)
    xxw(1:nneu) = xxw2(1:nneu)
  END DO
! Last layer
  outp = 0.
  outpd = 0.
  DO jj=1,nout
    DO jk=1,nneu
      outpd(jj) = outpd(jj) + wei(jk, jj, numl)*xxwd(jk)
      outp(jj) = outp(jj) + xxw(jk)*wei(jk, jj, numl)
    END DO
  END DO
  outp(1:nout) = outp(1:nout) + bia(1:nout, numl)
END SUBROUTINE NNINFER_WD

SUBROUTINE NNINFER_WB(numl, nfea, nout, nneu, inpu, inpub, &
                & outp, outpb, wei, bia, nwrk, iact)

!  Adjoint ANN operator
!  ---

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: numl, nfea, nout, nneu, nwrk,iact
  REAL(wp), INTENT(IN) :: inpu(nfea)
  REAL(wp) :: inpub(nfea)
  REAL(wp) :: outp(nout)
  REAL(wp) :: outpb(nout)
  REAL(wp) :: xxw(nwrk)
  REAL(wp) :: xxwb(nwrk)
  REAL(wp) :: xxw2(nwrk)
  REAL(wp) :: xxw2b(nwrk)
  REAL(wp) :: wei(nwrk, nwrk, numl)
  REAL(wp) :: bia(nwrk, numl)
  REAL(wp) :: zt
  INTEGER :: ji, jj, jk, jl
  INTRINSIC MAX
  INTEGER*4 :: branch

! First layer
  xxw = 0.
  DO jj=1,nneu
    DO jk=1,nfea
      xxw(jj) = xxw(jj) + inpu(jk)*wei(jk, jj, 1)
    END DO
  END DO
  xxw(1:nneu) = xxw(1:nneu) + bia(1:nneu, 1)
! RELU activation
  IF (iact .EQ. 1) THEN
    DO jj=1,nneu
      IF (0. .LT. xxw(jj)) THEN
        CALL PUSHCONTROL1B(0)
        xxw(jj) = xxw(jj)
      ELSE
        xxw(jj) = 0.
        CALL PUSHCONTROL1B(1)
      END IF
    END DO
    CALL PUSHCONTROL2B(3)
  ELSE IF (iact .EQ. 2) THEN
! TANH activation
    DO jj=1,nneu
      CALL PUSHREAL(xxw(jj))
      xxw(jj) = TANH(xxw(jj))
    END DO
    CALL PUSHCONTROL2B(2)
  ELSE IF (iact .EQ. 3) THEN
! SIGMOID activation
    DO jj=1,nneu
      CALL PUSHREAL(xxw(jj))
      xxw(jj) = 1./(1.+EXP(-xxw(jj)))
    END DO
    CALL PUSHCONTROL2B(1)
  ELSE
    CALL PUSHCONTROL2B(0)
  END IF
! Hidden layers
  DO jl=2,numl-1
    xxw2 = 0.
    DO jj=1,nneu
      DO jk=1,nneu
        xxw2(jj) = xxw2(jj) + xxw(jk)*wei(jk, jj, jl)
      END DO
    END DO
    xxw2(1:nneu) = xxw2(1:nneu) + bia(1:nneu, jl)
! RELU activation
    IF (iact .EQ. 1) THEN
      DO jj=1,nneu
        IF (0. .LT. xxw2(jj)) THEN
          CALL PUSHCONTROL1B(0)
          xxw2(jj) = xxw2(jj)
        ELSE
          xxw2(jj) = 0.
          CALL PUSHCONTROL1B(1)
        END IF
      END DO
      CALL PUSHCONTROL2B(0)
    ELSE IF (iact .EQ. 2) THEN
! TANH activation
      DO jj=1,nneu
        CALL PUSHREAL(xxw2(jj))
        xxw2(jj) = TANH(xxw2(jj))
      END DO
      CALL PUSHCONTROL2B(1)
    ELSE IF (iact .EQ. 3) THEN
! SIGMOID activation
      DO jj=1,nneu
        CALL PUSHREAL(xxw2(jj))
        xxw2(jj) = 1./(1.+EXP(-xxw2(jj)))
      END DO
      CALL PUSHCONTROL2B(2)
    ELSE
      CALL PUSHCONTROL2B(3)
    END IF
    xxw(1:nneu) = xxw2(1:nneu)
  END DO
  xxwb = 0.0
  DO jj=nout,1,-1
    DO jk=nneu,1,-1
      xxwb(jk) = xxwb(jk) + wei(jk, jj, numl)*outpb(jj)
    END DO
  END DO
  DO jl=numl-1,2,-1
    xxw2b = 0.0
    xxw2b(1:nneu) = xxw2b(1:nneu) + xxwb(1:nneu)
    xxwb(1:nneu) = 0.0
    CALL POPCONTROL2B(branch)
    IF (branch .LT. 2) THEN
      IF (branch .EQ. 0) THEN
        DO jj=nneu,1,-1
          CALL POPCONTROL1B(branch)
          IF (branch .NE. 0) xxw2b(jj) = 0.0
        END DO
      ELSE
        DO jj=nneu,1,-1
          CALL POPREAL(xxw2(jj))
          xxw2b(jj) = (1.0-TANH(xxw2(jj))**2)*xxw2b(jj)
        END DO
      END IF
    ELSE IF (branch .EQ. 2) THEN
      DO jj=nneu,1,-1
        CALL POPREAL(xxw2(jj))
        zt = EXP(-xxw2(jj)) + 1.
        xxw2b(jj) = EXP(-xxw2(jj))*xxw2b(jj)/zt**2
      END DO
    END IF
    DO jj=nneu,1,-1
      DO jk=nneu,1,-1
        xxwb(jk) = xxwb(jk) + wei(jk, jj, jl)*xxw2b(jj)
      END DO
    END DO
  END DO
  CALL POPCONTROL2B(branch)
  IF (branch .LT. 2) THEN
    IF (branch .NE. 0) THEN
      DO jj=nneu,1,-1
        CALL POPREAL(xxw(jj))
        zt = EXP(-xxw(jj)) + 1.
        xxwb(jj) = EXP(-xxw(jj))*xxwb(jj)/zt**2
      END DO
    END IF
  ELSE IF (branch .EQ. 2) THEN
    DO jj=nneu,1,-1
      CALL POPREAL(xxw(jj))
      xxwb(jj) = (1.0-TANH(xxw(jj))**2)*xxwb(jj)
    END DO
  ELSE
    DO jj=nneu,1,-1
      CALL POPCONTROL1B(branch)
      IF (branch .NE. 0) xxwb(jj) = 0.0
    END DO
  END IF
  inpub = 0.
  DO jj=nneu,1,-1
    DO jk=nfea,1,-1
      inpub(jk) = inpub(jk) + wei(jk, jj, 1)*xxwb(jj)
    END DO
  END DO
  outpb = 0.

END SUBROUTINE NNINFER_WB

SUBROUTINE ANNIF_TESTAD(IMOD,NTRIALS,INPU)

!  Adjoint test wrapper
!  ---

IMPLICIT NONE

        TYPE(ANNIF_TYPE), INTENT(IN) :: IMOD
        INTEGER, INTENT(IN)    :: NTRIALS
        REAL(wp), INTENT(IN)       :: INPU(IMOD%numfea)
        REAL(wp)                   :: PIN(IMOD%numfea),PI0(IMOD%numfea)
        REAL(wp)                   :: ROU(IMOD%numout),ROUD(IMOD%numout)
        REAL(wp)                   :: ROUB(IMOD%numout),ROUB0(IMOD%numout)
        REAL(wp)                   :: ZPREC(NTRIALS), ZA, ZB
        INTEGER                :: JT

        ! Random generator init
        IF(.NOT.LLIRS) CALL INIT_SEED

        ! Loop over number of tests
        DO JT=1,NTRIALS

          CALL random_stdnormalv(IMOD%numfea,PIN)
          CALL random_stdnormalv(IMOD%numout,ROUB)
          ROUD=0.
          PI0=0.
          CALL ANNIF_TL(IMOD,INPU,PIN,ROU,ROUD)
          ROUB0=ROUB
          CALL ANNIF_AD(IMOD,INPU,PI0,ROU,ROUB)
          ZA=DOT_PRODUCT( PIN, PI0 )
          ZB=DOT_PRODUCT( ROUD, ROUB0 )
          ZPREC(JT) = 100*ABS(ZA-ZB)/ABS(ZA)

        ENDDO

        ! Print results
        WRITE(IOUNO,*) 
        WRITE(IOUNO,*) CINST//''
        WRITE(IOUNO,*) CINST//' ADJOINT TEST'
        WRITE(IOUNO,*) CINST//' No of trials      : ', NTRIALS
        WRITE(IOUNO,*) CINST//' Minimum error  (%): ', MINVAL(ZPREC)
        WRITE(IOUNO,*) CINST//' Average error  (%): ', SUM(ZPREC)/NTRIALS
        WRITE(IOUNO,*) CINST//' Maximum error  (%): ', MAXVAL(ZPREC)
        WRITE(IOUNO,*) CINST//''
        WRITE(IOUNO,*) 

END SUBROUTINE ANNIF_TESTAD

SUBROUTINE ANNIF_TESTTL(IMOD,NTRIALS,INPU)

!  Tangent linear test wrapper
!  ---

IMPLICIT NONE

        TYPE(ANNIF_TYPE), INTENT(IN) :: IMOD
        INTEGER, INTENT(IN)    :: NTRIALS
        REAL(wp), INTENT(IN)       :: INPU(IMOD%numfea)
        REAL(wp)                   :: PIN(IMOD%numfea)
        REAL(wp)                   :: ROU(IMOD%numout),ROUD(IMOD%numout)
        REAL(wp)                   :: ROU2(IMOD%numout)
        REAL(wp)                   :: ZPREC(NTRIALS)
        REAL(wp)                   :: ZA(IMOD%numout), ZB(IMOD%numout)
        INTEGER                :: JT

        ! Random generator init
        IF(.NOT.LLIRS) CALL INIT_SEED

        CALL ANNIF_FWD(IMOD,INPU,ROU)

        ! Loop over number of tests
        DO JT=1,NTRIALS

          CALL random_stdnormalv(IMOD%numfea,PIN)
          PIN=PIN*0.01
          !!WRITE(*,*) '*--*',JT
          !!WRITE(*,*) MINVAL(INPU),MAXVAL(INPU)
          !!WRITE(*,*) MINVAL(PIN),MAXVAL(PIN)
          !!WRITE(*,*) MINVAL(INPU+PIN),MAXVAL(INPU+PIN)
          !!WRITE(*,*) MAXVAL(ABS(PIN))
          !WHERE (INPU+PIN > 1.) PIN=0.
          !WHERE (INPU+PIN < 1.) PIN=0.
          CALL ANNIF_FWD(IMOD,INPU+PIN,ROU2)
          CALL ANNIF_TL(IMOD,INPU,PIN,ROU,ROUD)
          ZA=ROU2-ROU
          ZB=ROUD
          !WRITE(*,*) SUM(ABS(ZA)),SUM(ABS(ZB))
          ZPREC(JT) = 100*ABS(SUM(ABS(ZA))-SUM(ABS(ZB)))/SUM(ABS(ZA))

        ENDDO

        ! Print results
        WRITE(IOUNO,*) 
        WRITE(IOUNO,*) CINST//''
        WRITE(IOUNO,*) CINST//' TANGENT-LINEAR TEST'
        WRITE(IOUNO,*) CINST//' No of trials      : ', NTRIALS
        WRITE(IOUNO,*) CINST//' Minimum error  (%): ', MINVAL(ZPREC)
        WRITE(IOUNO,*) CINST//' Average error  (%): ', SUM(ZPREC)/NTRIALS
        WRITE(IOUNO,*) CINST//' Maximum error  (%): ', MAXVAL(ZPREC)
        WRITE(IOUNO,*) CINST//''
        WRITE(IOUNO,*) 

END SUBROUTINE ANNIF_TESTTL

!
! NetCDF reading routines
!

SUBROUTINE GETGATT(CF,CA,AT)
USE NETCDF
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: CF, CA
  INTEGER , INTENT(OUT) :: AT
  INTEGER :: NCID, VID
  CALL CHECK( NF90_OPEN(CF,NF90_NOWRITE,NCID) )
  CALL CHECK( NF90_GET_ATT(NCID, NF90_GLOBAL, CA, AT) )
  CALL CHECK( NF90_CLOSE(NCID) )
END SUBROUTINE GETGATT

SUBROUTINE GETGATC(CF,CA,CAT)
USE NETCDF
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN)  :: CF, CA
  CHARACTER(LEN=*), INTENT(OUT) :: CAT
  INTEGER :: NCID, VID
  CALL CHECK( NF90_OPEN(CF,NF90_NOWRITE,NCID) )
  CALL CHECK( NF90_GET_ATT(NCID, NF90_GLOBAL, CA, CAT) )
  CALL CHECK( NF90_CLOSE(NCID) )
END SUBROUTINE GETGATC

SUBROUTINE GETVAR3(CF,CV,VV)
USE NETCDF
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: CF, CV
  REAL(wp), DIMENSION(:,:,:), INTENT(OUT) :: VV
  INTEGER :: NCID, VID
  CALL CHECK( NF90_OPEN(CF,NF90_NOWRITE,NCID) )
  CALL CHECK( NF90_INQ_VARID(NCID,CV,VID) )
  CALL CHECK( NF90_GET_VAR(NCID,VID,VV) )
  CALL CHECK( NF90_CLOSE(NCID) )
END SUBROUTINE GETVAR3

SUBROUTINE GETVAR2(CF,CV,VV)
USE NETCDF
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: CF, CV
  REAL(wp), DIMENSION(:,:), INTENT(OUT) :: VV
  INTEGER :: NCID, VID
  CALL CHECK( NF90_OPEN(CF,NF90_NOWRITE,NCID) )
  CALL CHECK( NF90_INQ_VARID(NCID,CV,VID) )
  CALL CHECK( NF90_GET_VAR(NCID,VID,VV) )
  CALL CHECK( NF90_CLOSE(NCID) )
END SUBROUTINE GETVAR2

SUBROUTINE CHECK(STATUS)
    USE NETCDF
    IMPLICIT NONE
    INTEGER, INTENT ( IN) :: STATUS
    IF(STATUS /= NF90_NOERR) THEN
      PRINT *, TRIM(NF90_STRERROR(STATUS))
      CALL ABORT()
    END IF
END SUBROUTINE CHECK

!
! Random generation routines
!

subroutine random_stduniform(u)
   implicit none
   real(wp),intent(out) :: u
   real(wp) :: r
   call random_number(r)
   u = 1 - r
end subroutine random_stduniform

subroutine random_stduniformv(n,u)
   implicit none
   real(wp),intent(out) :: u(n)
   real(wp) :: r
   integer :: i, n
   do i = 1 ,n
     call random_number(r)
     u(i) = 1 - r
   enddo
end subroutine random_stduniformv

subroutine random_stdnormal(x)
   implicit none
   real(wp),intent(out) :: x
   real(wp),parameter :: pi=3.14159265
   real(wp) :: u1,u2
   call random_stduniform(u1)
   call random_stduniform(u2)
   x = sqrt(-2*log(u1))*cos(2*pi*u2)
end subroutine random_stdnormal

subroutine random_stdnormalv(n,x)
   implicit none
   real(wp),intent(out) :: x(n)
   real(wp),parameter :: pi=3.14159265
   real(wp) :: u1,u2
   integer :: i, n
   do i = 1 ,n
     call random_stduniform(u1)
     call random_stduniform(u2)
     x(i) = sqrt(-2*log(u1))*cos(2*pi*u2)
   enddo
end subroutine random_stdnormalv

subroutine init_seed
   implicit none
        integer       :: n
        call random_seed(size=n)
        call random_seed(put=urandom_seed(n))  ! Put seed array into PRNG.
        LLIRS = .TRUE.
end subroutine init_seed

function urandom_seed(n, stat) result(seed)
   implicit none
        !! Returns a seed array filled with random values from `/dev/urandom`.
        integer, intent(in)            :: n
        integer, intent(out), optional :: stat
        integer                        :: seed(n)
        integer                        :: fu, rc

        open (access='stream', action='read', file='/dev/urandom', &
              form='unformatted', iostat=rc, newunit=fu)
        if (present(stat)) stat = rc
        if (rc == 0) read (fu) seed
        close (fu)
end function urandom_seed

!... Routines for TAPENADE portability

SUBROUTINE PUSHREAL_4(ZZ)
        REAL(4), INTENT(INOUT) :: ZZ
        CALL PUSHREAL4(ZZ)
END SUBROUTINE PUSHREAL_4
SUBROUTINE PUSHREAL_8(ZZ)
        REAL(8), INTENT(INOUT) :: ZZ
        CALL PUSHREAL8(ZZ)
END SUBROUTINE PUSHREAL_8
SUBROUTINE POPREAL_4(ZZ)
        REAL(4), INTENT(INOUT) :: ZZ
        CALL POPREAL4(ZZ)
END SUBROUTINE POPREAL_4
SUBROUTINE POPREAL_8(ZZ)
        REAL(8), INTENT(INOUT) :: ZZ
        CALL POPREAL8(ZZ)
END SUBROUTINE POPREAL_8

END MODULE ANNIF
