PROGRAM TEST

!
! Testing the use of the ANNIF module
! for online inference and TL/AD support
!
! A.S. @ cnr.it

        USE ANNIF
        IMPLICIT NONE

        REAL, ALLOCATABLE, DIMENSION(:) :: RIN, ROU, PIN, ROUD ,ROU2, ROUB, PIB, PI0, ROUB0
        CHARACTER(LEN=90) :: CFNN
        INTEGER :: JI
        TYPE(ANNIF_TYPE), ALLOCATABLE :: NNMOD(:)

        ! This is the NN model in NetCDF format
        CFNN='test/my_model_p4.nc'

        ! We will have only one model in this test
        ALLOCATE ( NNMOD(1) )

        ! Initialization
        CALL ANNIF_INIT ( NNMOD, 10 )

        ! Reading
        CALL ANNIF_RDMOD( NNMOD(1), CFNN )

        ! Allocat some arrays
        ALLOCATE( RIN(NNMOD(1)%numfea),  ROU(NNMOD(1)%numout), PIN(NNMOD(1)%numfea) )
        ALLOCATE( ROUD(NNMOD(1)%numout), ROU2(NNMOD(1)%numout) )
        ALLOCATE( ROUB(NNMOD(1)%numout), ROUB0(NNMOD(1)%numout), PIB(NNMOD(1)%numfea), PI0(NNMOD(1)%numfea) )

        ! Read some one sample
        OPEN(UNIT=21,FILE='test/sample1d.out')
        DO JI=1,NNMOD(1)%numfea
                READ(21,*) RIN(JI)
        ENDDO

        ! Forward modelling
        CALL ANNIF_FWD(NNMOD(1),RIN,ROU)
        WRITE(*,*) ' FWD MODEL, first ten output variables:'
        WRITE(*,*) ROU(1:10)

        PIN=0.
        PIN(72:75)=0.1

        CALL ANNIF_TL(NNMOD(1),RIN,PIN,ROU,ROUD)
        WRITE(*,*) ' TL MODEL, first ten output variables:'
        WRITE(*,*) ROU(1:10)
        WRITE(*,*) ROUD(1:10)

        PIB=0. 
        ROUD=0.
        ROUD(1:3)=0.1
        CALL ANNIF_AD(NNMOD(1),RIN,PIB,ROU,ROUD)
        WRITE(*,*) ' AD MODEL, first ten output variables:'
        WRITE(*,*) ROU(1:10)
        WRITE(*,*) PIB(1:10)

        STOP

!!SUBROUTINE ANNIF_AD(IMOD, INPU, INPU_AD, OUTP, OUTP_AD)
!!
!!!  NN adjoint model Wrapper
!!!  ---
!!
!!        IMPLICIT NONE
!!        TYPE(ANNIF_TYPE), INTENT(IN) :: IMOD
!!        REAL, INTENT(IN)       :: INPU(IMOD%numfea)
!!        REAL, INTENT(OUT)      :: OUTP(IMOD%numout)
!!        REAL, INTENT(INOUT)    :: INPU_AD(IMOD%numfea)
!!        REAL, INTENT(INOUT)    :: OUTP_AD(IMOD%numout)
!!        INTEGER                :: IAC
!!
!!        IAC=IMOD%innact
!!        CALL NNINFER_WB(IMOD%numlay,IMOD%numfea,IMOD%numout,IMOD%numneu,&
!!        & INPU, INPU_AD, OUTP, OUTP_AD, IMOD%wei, IMOD%bia, IMOD%numwrk, IAC)
!!
!!##


        ! Perform TL/AD tests
        CALL ANNIF_TESTTL(NNMOD(1),1000,RIN)
        CALL ANNIF_TESTAD(NNMOD(1),1000,RIN)

        WRITE(*,*) 
        WRITE(*,*) ' Changing activation to RELU'
        WRITE(*,*) ' to test TL/AD'
        WRITE(*,*)

        NNMOD(1)%cnnact = 'relu'
        NNMOD(1)%innact = 1

        CALL ANNIF_TESTTL(NNMOD(1),1000,RIN)
        CALL ANNIF_TESTAD(NNMOD(1),1000,RIN)

        WRITE(*,*)
        WRITE(*,*) ' Changing activation to LINEAR'
        WRITE(*,*) ' to test TL/AD'
        WRITE(*,*)

        NNMOD(1)%cnnact = 'none'
        NNMOD(1)%innact = 0

        CALL ANNIF_TESTTL(NNMOD(1),1000,RIN)
        CALL ANNIF_TESTAD(NNMOD(1),1000,RIN)

        WRITE(*,*)
        WRITE(*,*) ' Changing activation to SIGMOID'
        WRITE(*,*) ' to test TL/AD'
        WRITE(*,*)

        NNMOD(1)%cnnact = 'sigm'
        NNMOD(1)%innact = 3

        CALL ANNIF_TESTTL(NNMOD(1),1000,RIN)
        CALL ANNIF_TESTAD(NNMOD(1),1000,RIN)


END PROGRAM TEST
