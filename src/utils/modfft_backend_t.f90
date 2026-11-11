module modfft_backend_t

  use modprecision, only: pois_r

  implicit none

  private

  public :: FORWARD_FFT
  public :: BACKWARD_FFT
  public :: fft_backend_t

  integer, parameter :: FORWARD_FFT = 1
  integer, parameter :: BACKWARD_FFT = -1

  !> Abstract type for FFT implementations. Executes batched 1D FFTs.
  type, abstract :: fft_backend_t
  contains
    !> Constructor.
    procedure(fft_init), deferred :: init
    !> Get the work size required for the FFT.
    procedure(fft_get_worksize), deferred :: get_worksize
    !> Set the workspace for the FFT.
    procedure(fft_set_workspace), deferred :: set_workspace
    !> Execute the FFT.
    procedure(fft_exec), deferred :: exec
    !> Destructor.
    procedure(fft_finalize), deferred :: finalize
  end type fft_backend_t

  abstract interface
    !> Initialize an FFT.
    !!
    !! @param[in] direction Direction of the FFT (FORWARD_FFT or BACKWARD_FFT).
    !! @param[in] size Size of the FFT.
    !! @param[in] nbatch Number of FFTs to perform in a batch.
    !! @param[in] in Input data.
    !! @param[out] out Output data.
    subroutine fft_init(this, direction, size, nbatch, in, out)
      import :: fft_backend_t, pois_r
      class(fft_backend_t), intent(inout) :: this
      integer, intent(in) :: direction
      integer, intent(in) :: size
      integer, intent(in) :: nbatch
      real(pois_r), intent(in) :: in(:)
      real(pois_r), intent(out) :: out(:)
    end subroutine fft_init

    !> Get the work size required for the FFT.
    !!
    !! @return Required work size.
    function fft_get_worksize(this) result(worksize)
      import :: fft_backend_t
      class(fft_backend_t), intent(in) :: this
      integer :: worksize
    end function fft_get_worksize

    !> Set the workspace for the FFT.
    !!
    !! @param[in] workspace Workspace to use for the FFT.
    subroutine fft_set_workspace(this, workspace)
      import :: fft_backend_t, pois_r
      class(fft_backend_t), intent(inout) :: this
      real(pois_r), target, intent(in) :: workspace(:)
    end subroutine fft_set_workspace

    !> Execute the FFT.
    !!
    !! @param[inout] in Input data to transform. Output data is written here.
    subroutine fft_exec(this, in)
      import :: fft_backend_t, pois_r
      class(fft_backend_t), intent(in) :: this
      real(pois_r), intent(inout) :: in(:)
    end subroutine fft_exec

    !> Finalize the FFT.
    subroutine fft_finalize(this)
      import :: fft_backend_t
      class(fft_backend_t), intent(inout) :: this
    end subroutine fft_finalize

  end interface

end module modfft_backend_t