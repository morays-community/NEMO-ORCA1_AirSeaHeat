module nemo2infero
#ifdef key_infero

!!
!! Simple wrapper for INFERO for online inference in the NEMO model
!! Assume only one ML model is needed each time
!! A.Storto@CNR - 2023-06
!!
!! Use it as follows:
!! 1. Initialization:
!!    CALL nemo2infero_init(instance, ktype, cpath, kntype, laynam1, laynam2, kloguni, lwpi )
!!      instance: infero_instance type, defining data and model
!!      ktype:   Type of model (1: TF;   2: ONNX;   3: TFLite)
!!      cpath:   Directory where is the model
!!      kntype:  Normalization type (0: NO; 1: (-1,1); 2: (0,1); 3: mean/sd)
!!      laynam1, laynam2 : layer names in the NN model 
!!   (to be found eg with $ saved_model_cli show --dir TF_model_path --tag_set serve --signature_def serving_default )
!!      kloguni: Output unit for verbose diagnostics (e.g. numout in NEMO)
!!      lwpi   : Verbose diagnostics enabled
!!
!! 2. Inference:
!!    CALL nemo2infero_predict(instance,l,m,ina,n,p,oua,r1i,r2i,r1o,r2o)  (for 2D arrays, e.g. dense layers)
!!    CALL nemo2infero_predict(instance,l,m,n,p,ina,q,r,s,t,oua,r1i,r2i,r1o,r2o)  (for 4D arrays, e.g. convolutional layers)
!!      instance: infero_instance type, defining data and model
!!      ina :    Input array of dim  (l,m) or (l,m,n,p)
!!      oua :    Output array of dim (n,p) or (q,r,s,t)
!!      r1i,r2i :  For basic normalization (with dimensions m (2d) or m,n,p (4d)
!!      r1o,r2o :  For basic normalization (with dimensions p (2d) or r,s,t (4d)
!!      Note the dimensions must be consistent with the ML model read 
!!      during initialization
!!
!! 3. Finalize:
!!    CALL nemo2infero_finalize()
!!      with    optional <instance> sets free all arrays / model / tensors
!!      without optional <instance> finalizes infero
!! 

use lib_mpp
use inferof
use par_kind, only : wp
use iso_c_binding, only : c_double, c_int, c_float, c_char, c_null_char, c_ptr

implicit none
private

logical, save  :: ln_n2i_init = .false.
logical, save  :: lwi = .false.
integer, parameter :: sp = SELECTED_REAL_KIND( 6, 37)
integer, parameter :: dp = wp
integer, parameter :: clen = 1024
integer, save      :: iloguni = 6, ierruni = 0

INTERFACE nemo2infero_predict
MODULE PROCEDURE nemo2infero_predict_4dp, &
                 nemo2infero_predict_2dp
END INTERFACE

public :: nemo2infero_init, nemo2infero_predict, nemo2infero_finalize

! Compound data type to allow multiple infero instances
! Normalization values could be in
type infero_instance
        type(infero_model) :: model
        type(infero_tensor_set) :: iset, oset
        real(c_float), allocatable :: it2f(:,:), ot2f(:,:), it4f(:,:,:,:), ot4f(:,:,:,:)
        character(len=128) :: layer_names(2)
        logical :: ll_first_inference = .true. 
        logical :: ll_init = .false.
        integer :: knorm
end type
public infero_instance

contains

 subroutine nemo2infero_init(instance, ktype, cpath, kntype, laynam1, laynam2, kloguni, lwpi )

  implicit none
  type(infero_instance), intent(inout) :: instance
  integer, intent(in) :: ktype, kntype
  character(len=*), intent(in) :: cpath
  character(len=clen) :: model_type
  character(len=clen) :: yaml_config
  character(len=clen) :: cnstr
  character(len=*), intent(in) :: laynam1
  character(len=*), intent(in) :: laynam2
  integer, intent(in), optional :: kloguni
  logical, intent(in), optional :: lwpi

  if(present(kloguni)) iloguni = kloguni
  if(present(lwpi)) lwi = lwpi

  if(lwi) then
     write(iloguni,*) ''
     write(iloguni,*) ' INITIALIZING MODULE NEMO2INFERO '
     write(iloguni,*) ''
  endif

  if(.not. ln_n2i_init ) then 
      call infero_check(infero_initialise())
      ln_n2i_init = .true.
  endif

  select case (ktype)
        case(1)
           model_type = 'tf_c'
        case(2)
           model_type = 'onnx'
        case(3)
           model_type = 'tflite'
        case(4)
           model_type = 'trt'
        case default
           call abor1('Model type not recognized')
   end select

  ! YAML config string
  yaml_config = "---"//NEW_LINE('A') &
  //"  path: "//TRIM(cpath)//NEW_LINE('A') &
  //"  type: "//TRIM(model_type)//c_null_char

  if(lwi) &
  & write(iloguni,*) ' N2I :: Initiliazing model '//trim(model_type)//' from path '//trim(cpath)

  ! Get the model
  call infero_check(instance%model%initialise_from_yaml_string(yaml_config))

  if(lwi) then
     call infero_check(instance%model%print_statistics())
     call infero_check(instance%model%print_config())
  endif

  instance%layer_names(1) = laynam1
  instance%layer_names(2) = laynam2

  instance%knorm = kntype
  select case (kntype )
        case (0) 
           cnstr = 'None'
        case (1) 
           cnstr = 'Range (-1,1)'
        case (2) 
           cnstr = 'Range (0,1)'
        case (3) 
           cnstr = 'Normal'
        case default
           call abor1('Normalization not recognized')
  end select

  if(lwi) write(iloguni,*) ' N2I :: Normalization : '//trim(cnstr)

  instance%ll_init = .true.

 end subroutine nemo2infero_init

  subroutine nemo2infero_predict_4dp(instance,id1,id2,id3,id4,ina,od1,od2,od3,od4,oua,r1i,r2i,r1o,r2o)

  implicit none
  type(infero_instance), intent(inout) :: instance
  integer, intent(in) :: id1,id2,id3,id4
  integer, intent(in) :: od1,od2,od3,od4
  real(dp), intent(in)  :: ina(id1,id2,id3,id4)
  real(dp), intent(out) :: oua(od1,od2,od3,od4)
  real(dp), intent(in), optional  :: r1i(id2,id3,id4),r2i(id2,id3,id4),r1o(od2,od3,od4),r2o(od2,od3,od4)

  ! input and output tensors
  real(c_float)              :: cr1i(id2,id3,id4),cr2i(id2,id3,id4)
  integer :: jc

  if( .not. ln_n2i_init )     call abor1('INFERO is not initialized, aborting')
  if( .not. instance%ll_init) call abor1('INFERO Model is not initialized, aborting')
  if( id1 .ne. od1 ) call abor1('INFERO Model dimension mismatch, aborting')
  if( instance%knorm .gt. 0 .AND. (.not. present(r1i) .or. .not. present(r2i) &
          .or. .not. present(r1o) .or. .not. present(r2o) ) ) &
          call abor1('Mandatory range for normalization not provided, aborting')

  ! Setup tensors
  IF( instance%ll_first_inference ) THEN
     if(allocated(instance%it4f) .or. allocated(instance%ot4f)) &
     & call abor1('nemo2infero_predict_4dp : it4f/ot4f already allocated')
     allocate( instance%it4f(id1,id2,id3,id4) )
     allocate( instance%ot4f(od1,od2,od3,od4) )
     ! Input layer
     call infero_check(instance%iset%initialise())
     call infero_check(instance%iset%push_tensor(instance%it4f, TRIM(instance%layer_names(1))))
     call infero_check(instance%iset%print())
     ! Output layer
     call infero_check(instance%oset%initialise())
     call infero_check(instance%oset%push_tensor(instance%ot4f, TRIM(instance%layer_names(2))))
     call infero_check(instance%oset%print())
     instance%ll_first_inference = .false.
  ENDIF

  ! Cast & normalize
  instance%it4f = REAL ( ina, KIND=c_float )
  if ( instance%knorm .gt. 0 ) THEN
          cr1i = REAL( r1i, KIND=c_float ) 
          cr2i = REAL( r2i, KIND=c_float )
          if ( instance%knorm .eq. 1 ) THEN
                  DO jc=1,id1
                     instance%it4f(jc,:,:,:) = -1. + (2*(instance%it4f(jc,:,:,:) - cr1i(:,:,:))/(cr2i(:,:,:)-cr1i(:,:,:)))
                  ENDDO
          elseif ( instance%knorm .eq. 2 ) THEN
                  DO jc=1,id1
                     instance%it4f(jc,:,:,:) = (instance%it4f(jc,:,:,:)-cr1i(:,:,:))/(cr2i(:,:,:)-cr1i(:,:,:))
                  ENDDO
          elseif ( instance%knorm .eq. 3 ) THEN
                  DO jc=1,id1
                     instance%it4f(jc,:,:,:) = (instance%it4f(jc,:,:,:)-cr1i(:,:,:))/cr2i(:,:,:)
                  ENDDO
          endif
  ENDIF

  ! Inference
  call infero_check(instance%model%infer(instance%iset, instance%oset))

  ! Cast back
  oua = REAL( instance%ot4f, KIND=dp )


  if ( instance%knorm .eq. 1 ) then
          DO jc=1,id1
             oua(jc,:,:,:) = (( oua(jc,:,:,:) + 1. ) * (r2o(:,:,:)-r1o(:,:,:)))/2. + r1o(:,:,:)
          ENDDO
  elseif ( instance%knorm .eq. 2 ) then
          DO jc=1,id1
             oua(jc,:,:,:) = oua(jc,:,:,:)*(r2o(:,:,:)-r1o(:,:,:)) + r1o(:,:,:)
          ENDDO
  elseif ( instance%knorm .eq. 3 ) then
          DO jc=1,id1
             oua(jc,:,:,:) = r2o(:,:,:)*oua(jc,:,:,:)+r1o(:,:,:)
          ENDDO
  endif

 end subroutine nemo2infero_predict_4dp

 subroutine nemo2infero_predict_2dp(instance,id1,id2,ina,od1,od2,oua,r1i,r2i,r1o,r2o)

  implicit none
  type(infero_instance), intent(inout) :: instance
  integer, intent(in) :: id1,id2
  integer, intent(in) :: od1,od2
  real(dp), intent(in)  :: ina(id1,id2)
  real(dp), intent(out) :: oua(od1,od2)
  real(dp), intent(in), optional  :: r1i(id2),r2i(id2),r1o(od2),r2o(od2)

  ! input and output tensors
  real(c_float)              :: cr1i(id2), cr2i(id2)
  integer :: jc

  if( .not. ln_n2i_init )     call abor1('INFERO is not initialized, aborting')
  if( .not. instance%ll_init) call abor1('INFERO Model is not initialized, aborting')
  if( id1 .ne. od1 ) call abor1('INFERO Model dimension mismatch, aborting')
  if( instance%knorm .gt. 0 .AND. (.not. present(r1i) .or. .not. present(r2i) &
          .or. .not. present(r1o) .or. .not. present(r2o) ) ) &
          call abor1('Mandatory range for normalization not provided, aborting')

  ! Setup tensors
  IF( instance%ll_first_inference ) THEN
     if(allocated(instance%it2f) .or. allocated(instance%ot2f)) &
     & call abor1('nemo2infero_predict_2dp : it2f/ot2f already allocated')
     allocate( instance%it2f(id1,id2) )
     allocate( instance%ot2f(od1,od2) )
     ! Input layer
     call infero_check(instance%iset%initialise())
     call infero_check(instance%iset%push_tensor(instance%it2f, TRIM(instance%layer_names(1))))
     call infero_check(instance%iset%print())
     ! Output layer
     call infero_check(instance%oset%initialise())
     call infero_check(instance%oset%push_tensor(instance%ot2f, TRIM(instance%layer_names(2))))
     call infero_check(instance%oset%print())
     instance%ll_first_inference = .false.
  ENDIF

  ! Cast & normalize
  instance%it2f = REAL ( ina, KIND=c_float )
  if ( instance%knorm .gt. 0 ) THEN
          cr1i = REAL( r1i, KIND=c_float ) 
          cr2i = REAL( r2i, KIND=c_float )
          if ( instance%knorm .eq. 1 ) THEN
                  DO jc=1,id1
                     instance%it2f(jc,:) = -1. + (2*(instance%it2f(jc,:) - cr1i(:))/(cr2i(:)-cr1i(:)))
                  ENDDO
          elseif ( instance%knorm .eq. 2 ) THEN
                  DO jc=1,id1
                     instance%it2f(jc,:) = (instance%it2f(jc,:)-cr1i(:))/(cr2i(:)-cr1i(:))
                  ENDDO
          elseif ( instance%knorm .eq. 3 ) THEN
                  DO jc=1,id1
                     instance%it2f(jc,:) = (instance%it2f(jc,:)-cr1i(:))/cr2i(:)
                  ENDDO
          endif
  ENDIF

  ! Inference
  call infero_check(instance%model%infer(instance%iset, instance%oset))

  ! Cast back
  oua = REAL( instance%ot2f, KIND=dp )

  if ( instance%knorm .eq. 1 ) then
          DO jc=1,id1
             oua(jc,:) = (( oua(jc,:) + 1. ) * (r2o(:)-r1o(:)))/2. + r1o(:)
          ENDDO
  elseif ( instance%knorm .eq. 2 ) then
          DO jc=1,id1
             oua(jc,:) = oua(jc,:)*(r2o(:)-r1o(:)) + r1o(:)
          ENDDO
  elseif ( instance%knorm .eq. 3 ) then
          DO jc=1,id1
             oua(jc,:) = r2o(:)*oua(jc,:)+r1o(:)
          ENDDO
  endif

 end subroutine nemo2infero_predict_2dp

 subroutine nemo2infero_finalize(instance)

  implicit none
  type(infero_instance), intent(inout), optional :: instance

  ! free the tensors and the model
  if(present(instance)) then
     write(iloguni,*) ' Finalizing the INFERO instance'
     call infero_check(instance%iset%free())
     call infero_check(instance%oset%free())
     call infero_check(instance%model%free())
     if(allocated(instance%it2f)) deallocate(instance%it2f)
     if(allocated(instance%it4f)) deallocate(instance%it4f)
     if(allocated(instance%ot2f)) deallocate(instance%ot2f)
     if(allocated(instance%ot4f)) deallocate(instance%ot4f)
  else
     write(iloguni,*) ' Finalizing INFERO'
     call infero_check(infero_finalise())
  endif

 end subroutine nemo2infero_finalize

 subroutine abor1(ctxt)

  implicit none
  character(len=*) :: ctxt
  write(ierruni,*) ' Error : '//trim(ctxt)
  call ctl_stop('STOP',trim(ctxt))

 end subroutine abor1

#else

public :: nemo2infero_init, nemo2infero_predict, nemo2infero_finalize
character(len=128), public :: layer_names(2)
type infero_instance
        logical :: ll_init = .false.
end type
public infero_instance

contains
 subroutine nemo2infero_init
         CALL ctl_stop( 'STOP', 'nemo2infero : NEMO compiled without key_infero' )
 end subroutine nemo2infero_init
 subroutine nemo2infero_predict
         CALL ctl_stop( 'STOP', 'nemo2infero : NEMO compiled without key_infero' )
 end subroutine nemo2infero_predict
 subroutine nemo2infero_finalize
         CALL ctl_stop( 'STOP', 'nemo2infero : NEMO compiled without key_infero' )
 end subroutine nemo2infero_finalize

#endif
end module nemo2infero
