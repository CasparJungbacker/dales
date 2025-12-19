module modpluto

  use, intrinsic :: iso_fortran_env, only: real32, real64

  use pluto_module, only: pluto, pluto_memory_resource, pluto_allocator

  implicit none

  private

  type(pluto_allocator) :: allocator

  public :: init_pluto
  public :: get_field
  public :: release_field

  interface get_field
    procedure :: get_field_r8
  end interface

contains

  !> Setup the Pluto allocator.
  subroutine init_pluto()

    allocator = pluto%make_allocator(pluto%host_resource())

  end subroutine init_pluto

  subroutine get_field_r8(field, size)

    real(real64), pointer, intent(out) :: field(:)
    integer, intent(in) :: size

    call allocator%allocate(field, [size])

  end subroutine get_field_r8

  subroutine release_field(field)

    real(real64), pointer, intent(inout) :: field(:)

    call allocator%deallocate(field)

  end subroutine release_field

end module modpluto