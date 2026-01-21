module modreduction_cuda

  use iso_fortran_env
  use iso_c_binding
  use openacc
 
  implicit none

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

end module modreduction_cuda
