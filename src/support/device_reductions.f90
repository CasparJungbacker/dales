!> Wrapper module for CUDA reductions.
module device_reductions
  
  use iso_c_binding
  use iso_fortran_env, only: real32, real64
  use openacc

  implicit none

  private

  interface
    subroutine reduction_profile_halo_float(input, output, nh, itot, jtot, &
                                            ktot, stream) &
      bind(C, name='reduction_profile_halo_float')
      import :: c_ptr, acc_handle_kind
      type(c_ptr), value :: input, output
      integer, value :: nh, itot, jtot, ktot
      integer(acc_handle_kind), value :: stream
    end subroutine reduction_profile_halo_float

    subroutine reduction_profile_halo_double(input, output, nh, itot, jtot, &
                                             ktot, stream) &
      bind(C, name='reduction_profile_halo_double')
      import :: c_ptr, acc_handle_kind
      type(c_ptr), value :: input, output
      integer, value :: nh, itot, jtot, ktot
      integer(acc_handle_kind), value :: stream
    end subroutine reduction_profile_halo_double
  end interface

  interface reduction_profile
    module procedure :: reduction_profile_real32
    module procedure :: reduction_profile_real64
  end interface reduction_profile

contains

  subroutine reduction_profile_real32(input, nh, itot, jtot, ktot, output, &
                                      opt_stream)

    real(real32), target, intent(in) :: input(:,:,:)
    integer,              intent(in) :: nh
    integer,              intent(in) :: itot
    integer,              intent(in) :: jtot
    integer,              intent(in) :: ktot

    real(real32), target, intent(out) :: output(:)

    integer, optional, intent(in) :: opt_stream

    integer(acc_handle_kind) :: stream

    if (present(opt_stream)) then
      stream = acc_get_cuda_stream(opt_stream)
    else
      stream = acc_get_cuda_stream(1)
    end if

    !$acc host_data use_device(input, output)
    call reduction_profile_halo_float(c_loc(input), c_loc(output), nh, itot, &
                                      jtot, ktot, stream)
    !$acc end host_data

  end subroutine reduction_profile_real32

  subroutine reduction_profile_real64(input, nh, itot, jtot, ktot, output, &
                                      opt_stream)

    real(real64), target, intent(in) :: input(:,:,:)
    integer,              intent(in) :: nh
    integer,              intent(in) :: itot
    integer,              intent(in) :: jtot
    integer,              intent(in) :: ktot

    real(real64), target, intent(out) :: output(:)

    integer, optional, intent(in) :: opt_stream

    integer(acc_handle_kind) :: stream

    if (present(opt_stream)) then
      stream = acc_get_cuda_stream(opt_stream)
    else
      stream = acc_get_cuda_stream(1)
    end if

    !$acc host_data use_device(input, output)
    call reduction_profile_halo_double(c_loc(input), c_loc(output), nh, itot, &
                                       jtot, ktot, stream)
    !$acc end host_data

  end subroutine reduction_profile_real64

end module device_reductions