module modfft_fftw_t

  use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_associated, c_int
#ifdef USE_FFTW
  use modfftw_interface, only: fftw_plan_many_r2r_if, fftw_execute_r2r_if, &
                               fftw_destroy_plan, fftwf_destroy_plan, C_FFTW_R2R_KIND, &
                               FFTW_R2HC, FFTW_HC2R, FFTW_MEASURE
#endif
  use modfft_t, only: fft_t, FORWARD_FFT, BACKWARD_FFT
  use modprecision, only: pois_r
  use fortran_support, only: finish

  implicit none

  private

  character(len=*), parameter :: modname = 'modfft_fftw_t'

  public :: fft_fftw_t

  type, extends(fft_t) :: fft_fftw_t
    private
    type(c_ptr) :: plan = c_null_ptr
  contains
    procedure :: init => fftw_init
    procedure :: get_worksize => fftw_get_worksize
    procedure :: set_workspace => fftw_set_workspace
    procedure :: exec => fftw_exec
    procedure :: finalize => fftw_finalize
  end type fft_fftw_t

contains

  subroutine fftw_init(this, direction, size, nbatch, in, out)

    class(fft_fftw_t), intent(inout) :: this

    integer,      intent(in) :: direction
    integer,      intent(in) :: size
    integer,      intent(in) :: nbatch
    real(pois_r), pointer, intent(in) :: in(:)
    real(pois_r), pointer, intent(in) :: out(:)

    character(len=*), parameter :: routine = modname//'/fftw_init'

    integer :: embed(1)
    integer :: kind(1)

#ifdef USE_FFTW
    if (size <= 0 .or. nbatch <= 0) call finish(routine, 'FFT size and batch count must be positive')

    ! TODO: implement discrete cosine transform
    select case (direction)
      case (FORWARD_FFT)
        kind(1) = FFTW_R2HC
      case (BACKWARD_FFT)
        kind(1) = FFTW_HC2R
      case default
        call finish(routine, 'Invalid FFT direction')
    end select

    embed(1) = size

    this%plan = fftw_plan_many_r2r_if( &
      rank=1, &
      n=embed, &
      howmany=nbatch, &
      in=in, &
      inembed=embed, &
      istride=1, &
      idist=size, &
      out=out, &
      onembed=embed, &
      ostride=1, &
      odist=size, &
      kind=kind, &
      flags=FFTW_MEASURE &
    )

    if (.not. c_associated(this%plan)) call finish(routine, 'Failed to create FFTW plan')

#else
    call finish(routine, "DALES was not compiled with FFTW support. Please recompile with FFTW enabled.")
#endif

  end subroutine fftw_init

  function fftw_get_worksize(this) result(worksize)

    class(fft_fftw_t), intent(in) :: this
    integer :: worksize

    worksize = -1

  end function fftw_get_worksize

  subroutine fftw_set_workspace(this, workspace)

    class(fft_fftw_t), intent(inout) :: this
    real(pois_r), target, intent(in) :: workspace(:)

  end subroutine fftw_set_workspace

  subroutine fftw_exec(this, in)

    class(fft_fftw_t), intent(in) :: this

    real(pois_r), intent(inout) :: in(:)

    character(len=*), parameter :: routine = modname//'/fftw_exec'

    if (.not. c_associated(this%plan)) call finish(routine, 'FFTW plan is not initialized')

#ifdef USE_FFTW
    call fftw_execute_r2r_if(this%plan, in, in)
#endif

  end subroutine fftw_exec

  subroutine fftw_finalize(this)

    class(fft_fftw_t), intent(inout) :: this
#ifdef USE_FFTW
    if (c_associated(this%plan)) then
#if POIS_PRECISION == 32
      call fftwf_destroy_plan(this%plan)
#else
      call fftw_destroy_plan(this%plan)
#endif
      this%plan = c_null_ptr
    end if
#endif

  end subroutine fftw_finalize

end module modfft_fftw_t