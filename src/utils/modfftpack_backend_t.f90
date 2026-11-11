module modfftpack_backend_t

  use fortran_support, only: finish
  use fftnew, only: rfftf, rfftb, rffti
  use modfft_backend_t, only: fft_backend_t, FORWARD_FFT, BACKWARD_FFT
  use modprecision, only: pois_r

  implicit none

  private

  character(len=*), parameter :: modname = 'modfftpack_backend_t'

  public :: fftpack_backend_t

  type, extends(fft_backend_t) :: fftpack_backend_t
    private
    integer :: direction
    integer :: nbatch
    integer :: size
    !> Internal workspace. FFTPACK actually stores values in the work array, so we don't use the externally supplide one.
    real(pois_r), allocatable :: work(:)
  contains
    procedure :: init => fftpack_init
    procedure :: get_worksize => fftpack_get_worksize
    procedure :: set_workspace => fftpack_set_workspace
    procedure :: exec => fftpack_exec
    procedure :: finalize => fftpack_finalize
  end type fftpack_backend_t

contains

  !! Initialize an FFT with the FFTPACK backend.
  !!
  !! @param[in] direction Direction of the FFT (FORWARD_FFT or BACKWARD_FFT).
  !! @param[in] size Size of the FFT.
  !! @param[in] nbatch Number of FFTs to perform in a batch.
  !! @param[in] in Input data.
  !! @param[out] out Output data.
  subroutine fftpack_init(this, direction, size, nbatch, in, out)

    class(fftpack_backend_t), intent(inout) :: this

    integer,      intent(in) :: direction
    integer,      intent(in) :: size
    integer,      intent(in) :: nbatch
    real(pois_r), intent(in) :: in(:)

    real(pois_r), intent(out) :: out(:)

    character(len=*), parameter :: routine = modname//'//fftpack_init'

    this%size = size
    this%nbatch = nbatch

    select case(direction)
      case(FORWARD_FFT, BACKWARD_FFT)
        this%direction = direction
      case default
        call finish(routine, 'Invalid FFT direction')
    end select

    allocate(this%work(this%get_worksize()))
    call rffti(this%size, this%work)

  end subroutine fftpack_init

  !> Get the required size of the internal workspace for the FFTPACK backend.
  !!
  !! @return Required size of the internal workspace.
  function fftpack_get_worksize(this) result(worksize)

    class(fftpack_backend_t), intent(in) :: this

    character(len=*), parameter :: routine = modname//'//fftpack_get_worksize'

    integer :: worksize

    worksize = 2 * this%size + 15 

  end function fftpack_get_worksize

  !> Set the internal workspace pointer to a user-allocated buffer.
  !!
  !! The supplied workspace is actually not used by the FFTPACK backend, 
  !! but this routine is provided to satisfy the interface of the fft_backend_t type.
  !!
  !! @param[in] workspace Buffer to use as workspace.
  subroutine fftpack_set_workspace(this, workspace)

    class(fftpack_backend_t), intent(inout) :: this

    real(pois_r), target, intent(in) :: workspace(:)

    character(len=*), parameter :: routine = modname//'//fftpack_set_workspace'

  end subroutine fftpack_set_workspace

  !> Execute FFT with the FFTPACK backend.
  !!
  !! @param[inout] in Input data.
  subroutine fftpack_exec(this, in)

    class(fftpack_backend_t), intent(in) :: this

    real(pois_r), intent(inout) :: in(:)

    character(len=*), parameter :: routine = modname//'//fftpack_exec'

    integer :: ibatch

    if (this%direction == FORWARD_FFT) then
      do ibatch = 0, this%nbatch - 1
        call rfftf(this%size, in(ibatch * this%size + 1:), this%work)
      end do
    else
      do ibatch = 0, this%nbatch - 1
        call rfftb(this%size, in(ibatch * this%size + 1:), this%work)
      end do
    end if

  end subroutine fftpack_exec

  !> Finalize the FFTPACK backend.
  subroutine fftpack_finalize(this)

    class(fftpack_backend_t), intent(inout) :: this

    character(len=*), parameter :: routine = modname//'//fftpack_finalize'

    deallocate(this%work)

  end subroutine fftpack_finalize

end module modfftpack_backend_t