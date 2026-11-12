!> Provides a factory method for creating FFT plans.
submodule(modfft_t) submodfft_factory

  use fortran_support,  only: finish
  use modfft_fftpack_t, only: fft_fftpack_t
  use modfft_fftw_t,    only: fft_fftw_t

  implicit none

  character(len=*), parameter :: modname = 'submodfft_factory'

contains

  !> Create and initialize an FFT plan.
  !!
  !! @param[inout] object Allocatable object of type fft_t to be created.
  !! @param[in] backend_id ID of the FFT backend to use.
  !! @param[in] direction Direction of the FFT (FORWARD_FFT or BACKWARD_FFT).
  !! @param[in] size Size of the FFT.
  !! @param[in] nbatch Number of FFTs to perform in a batch.
  !! @param[in] in Input data.
  !! @param[in] out Output data.
  module subroutine fft_make_plan(object, backend_id, direction, size, nbatch, in, out)

    class(fft_t), allocatable, intent(inout) :: object

    integer,      intent(in) :: backend_id
    integer,      intent(in) :: direction
    integer,      intent(in) :: size
    integer,      intent(in) :: nbatch
    real(pois_r), pointer, intent(in) :: in(:)
    real(pois_r), pointer, intent(in) :: out(:)

    character(len=*), parameter :: routine = modname//'/fft_make_plan'

    if (allocated(object)) deallocate(object)

    select case(backend_id)
      case (0)
        allocate(fft_fftpack_t :: object)
      case (100)
        allocate(fft_fftw_t :: object)
      case default
        call finish(routine, 'Unknown FFT backend ID')
    end select

    call object%init(direction, size, nbatch, in, out)

  end subroutine fft_make_plan

end submodule submodfft_factory