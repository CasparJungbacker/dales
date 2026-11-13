!> cuFFT backend for FFT operations.
module modfft_cufft_t

#ifdef USE_CUDA
  use cufft
#endif
  use fortran_support, only: finish
  use modfft_t,        only: fft_t, FORWARD_FFT, BACKWARD_FFT 
  use modprecision,    only: pois_r

  implicit none

  private

  character(len=*), parameter :: modname = 'modfft_cufft_t'

  public :: fft_cufft_t

  !> Implementation of the FFT type using the cuFFT backend.
  type, extends(fft_t) :: fft_cufft_t
    private
    integer :: direction
    integer :: nbatch
    integer :: size
    integer :: plan = 0 !< cuFFT plan handle.
  contains
    procedure :: init => cufft_init
    procedure :: get_worksize => cufft_get_worksize
    procedure :: set_workspace => cufft_set_workspace
    procedure :: exec => cufft_exec
    procedure :: finalize => cufft_finalize
  end type fft_cufft_t

contains

  !> Initialize the cuFFT FFT plan.
  subroutine cufft_init(this, direction, size, nbatch, in, out)

    class(fft_cufft_t), intent(inout) :: this

    integer,               intent(in) :: direction !< Direction of the FFT (FORWARD_FFT or BACKWARD_FFT).
    integer,               intent(in) :: size      !< Size of the FFT.
    integer,               intent(in) :: nbatch    !< Number of FFTs to perform in a batch.
    real(pois_r), pointer, intent(in) :: in(:)     !< Pointer to the input data.
    real(pois_r), pointer, intent(in) :: out(:)    !< Pointer to the output data.

    character(len=*), parameter :: routine = modname//'/cufft_init'

    integer           :: istat
    integer           :: fft_kind
    integer, pointer  :: null_ptr => null()
    character(len=64) :: msg

#ifdef USE_CUDA
    select case (direction)
      case (FORWARD_FFT)
#if POIS_PRECISION == 32
        fft_kind = CUFFT_R2C
#else
        fft_kind = CUFFT_C2R
#endif
      case (BACKWARD_FFT)
#if POIS_PRECISION == 32
        fft_kind = CUFFT_D2Z
#else
        fft_kind = CUFFT_Z2D
#endif
      case default
        call finish(routine, 'Invalid FFT direction')
    end select

    this%size = size
    this%nbatch = nbatch

    istat = cufftPlanMany( &
      this%plan, &
      1, &
      size, &
      null_ptr, &
      0, &
      0, &
      null_ptr, &
      0, &
      0, &
      fft_kind, &
      nbatch &
    )

    if (istat /= 0) then
      write(msg, *) "cuFFT returned non-zero exit code: ", istat
      call finish(routine, msg)
    end if
#else
    call finish(routine, "DALES was not compiled with CUDA support")
#endif

  end subroutine cufft_init   

  !> Get the workspace size required for a cuFFT FFT.
  function cufft_get_worksize(this) result(worksize)

    class(fft_cufft_t), intent(in) :: this

    character(len=*), parameter :: routine = modname//'/cufft_get_worksize'

    integer                 :: worksize       !< Required workspace size.
    integer(int_ptr_kind()) :: worksize_bytes !< Workspace size in bytes.
    integer                 :: istat          !< cuFFT return code.
    character(len=64)       :: msg

#ifdef USE_CUDA
    istat = cufftGetSize(this%plan, worksize_bytes)

    if (istat /= 0) then
      write(msg, *) "cuFFT returned non-zero exit code: ", istat
      call finish(routine, msg)
    end if

    ! Convert bytes => number of elements
    worksize = worksize_bytes / (storage_size(1.0_pois_r) / 8)
#endif

  end function cufft_get_worksize

  subroutine cufft_set_workspace(this, workspace)

    class(fft_cufft_t), intent(inout) :: this

    real(pois_r), pointer, intent(in) :: workspace(:) !< Pointer to the workspace. The underlying array is assumed to be allocated on GPU already.

    character(len=*), parameter :: routine = modname//'/cufft_set_workspace' 

    integer           :: istat !< cuFFT return code.
    character(len=64) :: msg

#ifdef USE_CUDA
    !$acc host_data use_device(workspace)
    !$omp target data use_device_addr(workspace)
    istat = cufftSetWorkArea(this%plan, workspace)
    !$acc end host_data
    !$omp end target data

    if (istat /= 0) then
      write(msg, *) "cuFFT returned non-zero exit code: ", istat
      call finish(routine, msg)
    end if
#endif

  end subroutine cufft_set_workspace

  !> Execute the FFT using the cuFFT backend.
  subroutine cufft_exec(this, in)
    
    class(fft_cufft_t), intent(in) :: this
    
    real(pois_r), intent(inout) :: in(:) !< Pointer to the input/output data. The underlying array is assumed to be allocated on GPU already.

    character(len=*), parameter :: routine = modname//'/cufft_exec'

    integer           :: istat        !< cuFFT return code.
    integer           :: ibatch       !< Batch index.
    integer           :: batch_offset !< Offset for the current batch in the input/output array.
    character(len=64) :: msg

#ifdef USE_CUDA
    if (this%direction == FORWARD_FFT) then
      !$acc host_data use_device(in)
      !$omp target data use_device_addr(in)
#if POIS_PRECISION == 32
      istat = cufftExecR2C(this%plan, in, in)
#else
      istat = cufftExecD2Z(this%plan, in, in)
#endif
      !$acc end host_data
      !$omp end target data

      ! Post-process signal
      !$acc parallel loop default(present)
      !$omp target teams loop
      do ibatch = 0, this%nbatch - 1
        batch_offset = ibatch * (this%size + 2)
        in(2 + batch_offset) = in(this%size + 1 + batch_offset)
      end do
    else
      ! Pre-process signal
      !$acc parallel loop default(present)
      !$omp target teams loop
      do ibatch = 0, this%nbatch - 1
        batch_offset = ibatch * (this%size + 2)
        in(this%size + 1 + batch_offset) = in(2 + batch_offset)
        in(2 + batch_offset) = 0
      end do

      !$acc host_data use_device(in)
      !$omp target data use_device_addr(in)
#if POIS_PRECISION == 32
      istat = cufftExecC2R(this%plan, in, in)
#else
      istat = cufftExecZ2D(this%plan, in, in)
#endif
      !$acc end host_data
      !$omp end target data
    end if

    if (istat /= 0) then
      write(msg, *) "cuFFT returned non-zero exit code: ", istat
      call finish(routine, msg)
    end if
#endif

  end subroutine cufft_exec

  !> Destroy the cuFFT FFT plan.
  subroutine cufft_finalize(this)

    class(fft_cufft_t), intent(inout) :: this 

    character(len=*), parameter :: routine = modname//'/cufft_finalize'

    integer           :: istat !< cuFFT return code.
    character(len=64) :: msg

#ifdef USE_CUDA
    istat = cufftDestroy(this%plan)

    if (istat /= 0) then
      write(msg, *) "cuFFT returned non-zero exit code: ", istat
      call finish(routine, msg)
    end if
#endif

  end subroutine cufft_finalize

end module modfft_cufft_t
