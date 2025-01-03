module modtracer_type

  use modprecision, only: field_r
 
  implicit none
  
  type tracer_t
  ! Fixed tracer properties
      ! Tracer name
      character(len=16) :: tracname
      ! Tracer long name
      character(len=64) :: traclong="dummy long name"
      ! Tracer unit
      character(len=16) :: unit="dummy unit"
      ! Moleculare mass of tracer (g mol-1)
      real(field_r)     :: molar_mass=-999.
      ! Tracer index in sv0, svm, svp
      integer           :: trac_idx=-1
      ! Boolean if tracer is emitted 
      logical           :: lemis=.false.
      ! Boolean if tracer is reactive
      logical           :: lreact=.false.
      ! Boolean if tracer is deposited
      logical           :: ldep=.false.
      ! Boolean if in A-gs
      logical           :: lags=.false.
      ! Boolean if in cloud microphysics
      logical           :: lmicro=.false.
      ! Boolean if tracer is nudged
      logical           :: lnudge=.false.
      ! Boolean if in aerosol microphysics
      logical           :: laero=.false.
      ! ! Static tracer properties:
      logical           :: lsurfsource = .false.
      real(field_r)     :: surface_source = 0
      ! real :: diffusivity

  contains
    procedure :: print_properties => tracer_print_properties
  end type tracer_t

  type, public :: tracer_ptr_t
    type(tracer_t), pointer :: ptr => null()
  end type tracer_ptr_t

contains

  subroutine tracer_print_properties(self)

    class(tracer_t), intent(in) :: self

    write(*,*) "Tracer: ", self%tracname
    write(*,*) "  long name  : ", trim(self%traclong)
    write(*,*) "  unit       : ", trim(self%unit)
    write(*,*) "  molar mass : ", self%molar_mass
    write(*,*) "  index      : ", self%trac_idx
    write(*,*) "  lemis      : ", self%lemis
    write(*,*) "  lreact     : ", self%lreact
    write(*,*) "  ldep       : ", self%ldep
    write(*,*) "  lags       : ", self%lags
    write(*,*) "  lmicro     : ", self%lmicro
    write(*,*) "  laero      : ", self%laero

  end subroutine tracer_print_properties

end module modtracer_type
