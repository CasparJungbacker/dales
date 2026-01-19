program test_bind
  use iso_c_binding
  use openacc
  
  implicit none
  
  interface
    subroutine reduction_2d_float(input, output, nh, itot, jtot) &
      bind(C, name="reduction_2d_float")
      import :: c_ptr, c_int
      type(c_ptr), value :: input, output
      integer(c_int), value :: nh, itot, jtot
    end subroutine
  end interface
  
  real(c_float), allocatable, target :: a(:,:)
  real(c_float), target :: res
  integer(c_int) :: itot, jtot, i, j, nh

  itot = 64
  jtot = 64
  nh = 1
    
  allocate(a(itot+2*nh,jtot+2*nh))
    
  do j = 1, jtot + 2 * nh
    do i = 1, itot + 2 * nh
      if (i > nh .and. i <= itot + nh .and. j > nh .and. j <= jtot + nh) then
        a(i,j) = 1
      else
        a(i,j) = -5000
      end if
    end do
  end do

  res = 0
  
  !$acc enter data copyin(a(:,:)) create(res)

  !$acc update device(res)
  
  !$acc host_data use_device(a, res)
  call reduction_2d_float(c_loc(a(1,1)), c_loc(res), nh, itot, jtot)
  !$acc end host_data
  
  
  print *, "the result is: ", res

end program
