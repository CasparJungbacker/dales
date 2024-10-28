!> \file modaeroso
module modaerosol
  use modglobal,      only: ifnamopt, fname_options, checknamelisterror, &
                            cexpnr, i1, j1, k1
  use modmpi,         only: myid
  use modprecision,   only: field_r
  use modtracer_type, only: tracer_t, tracer_ptr_t
  use modtracers,     only: add_tracer, allocate_tracers, tracer_prop, &
                            get_tracer
  use modstat_nc
  use go,             only: goSplitString_s ! TODO: move this to utils
  use utils,          only: dales_error

  implicit none

  private

  ! Parameters
  integer,       parameter :: maxmodes = 9
  character(3),  parameter :: modenames(9) = (/ 'nus', 'ais', 'acs', 'cos', &
                                                'aii', 'aci', 'coi', 'inr', &
                                                'inc' /)
  character(22), parameter :: longnames(9) = (/ 'soluble nucleation    ', &
                                                'soluble Aitken        ', &
                                                'soluble accumulation  ', &
                                                'soluble coarse        ', &
                                                'insoluble Aitken      ', &
                                                'insoluble accumulation', &
                                                'insoluble coarse      ', &
                                                'in-rain               ', &
                                                'in-cloud              ' /)
  real(field_r), parameter :: sigma_g(maxmodes) = (/ 1.59, 1.59, 1.59, 2.00, &
                                                     1.59, 1.59, 2.00, -999., &
                                                     -999. /)
  integer,       parameter :: iNUS = 1, iAIS = 2, iACS = 3, iCOS = 4, &
                              iAII = 5, iACI = 6, iCOI = 7, iINR = 8, &
                              iINC = 9

  ! Types
  !> Aerosol type
  type aerosol_t
    character(3)  :: name       !< Short name.
    character(64) :: long_name  !< Long name.
    real(field_r) :: rho        !< Density.
    real(field_r) :: kappa      !< Hygroscopicity.
    integer       :: nmodes = 0 !< Number of modes this aerosol participates in.

    type(tracer_ptr_t), allocatable :: tracers(:) !< List of tracers for mass
                                                  !! concentrations.
    character(3),       allocatable :: modes(:)   !< List of modes this aerosol
                                                  !! participates in.
   contains
    procedure, pass(self) :: construct => aerosol_construct !< Constructor
  end type aerosol_t

  !> Pointer to an aerosol
  type aerosol_ptr_t
    type(aerosol_t), pointer :: ptr => null()
  end type aerosol_ptr_t

  !> M7 mode type
  type mode_t
    ! General properties
    character(3)       :: name          !< Short name of the mode.
    character(64)      :: long_name     !< Full name of the mode.
    integer            :: nspecies = 0  !< Number of species in the mode.
    real(field_r)      :: sigma_g       !< Geometric standard deviation.
    logical            :: lactivation   !< Aerosols can be activated.
    type(tracer_ptr_t) :: N             !< Tracer for number concentration.

    type(aerosol_ptr_t), allocatable :: aerosols(:) !< List of aerosols.

    ! Fields
    real(field_r), allocatable :: conc(:,:,:,:) !< Mass/number concentrations.
    real(field_r), allocatable :: tend(:,:,:,:) !< Tendencies.
   contains
    procedure, pass(self) :: construct => mode_construct     !< Constructor.
    procedure, pass(self) :: add_aerosol => mode_add_aerosol !< Add an aerosol to the mode.
    procedure, pass(self) :: allocate => mode_allocate       !< Allocate workspace.
  end type mode_t

  ! Variables
  logical      :: laerosol = .false. !< Switch for enabling/disabling interactive aerosols.
  type(mode_t) :: modes(maxmodes)    !< List of modes.

  type(aerosol_t), allocatable :: aerosols(:)

  ! Procedures
  public :: initaerosol
  public :: exitaerosol

contains
  !> Read input files and setup aerosols and M7 modes.
  subroutine initaerosol
    integer       :: imod, ierr, ncid, nvars, iaer, mode_loc
    character(3)  :: name
    character(64) :: long_name
    character(27) :: modes_str
    real(field_r) :: rho, kappa

    integer, allocatable :: varids(:)

    namelist /NAMAEROSOL/ laerosol

    ! Read input
    if (myid == 0) then
      ! Namelist
      open(ifnamopt, file=fname_options, status='old', iostat=ierr)
      read(ifnamopt, NAMAEROSOL, iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMAEROSOL')
      close(ifnamopt)
    end if

    if (.not. laerosol) return

    ! Setup the modes
    do imod = 1, maxmodes
      call modes(imod) % construct(name=modenames(imod), long_name=longnames(imod), sigma_g=sigma_g(imod))
    end do

    if (myid == 0) then
      call nchandle_error(nf90_open("aerosol."//cexpnr//".nc", NF90_NOWRITE, ncid))
      call nchandle_error(nf90_inquire(ncid, nVariables=nvars))
      
      allocate(aerosols(nvars), varids(nvars))

      call nchandle_error(nf90_inq_varids(ncid, nvars, varids))

      do iaer = 1, nvars
        call nchandle_error( &
          nf90_inquire_variable(ncid, varids(iaer), name=name) &
        )
        call nchandle_error(&
          nf90_get_att(ncid, varids(iaer), "long_name", long_name) &
        )
        call nchandle_error(&
          nf90_get_att(ncid, varids(iaer), "rho", rho) &
        )
        call nchandle_error(&
          nf90_get_att(ncid, varids(iaer), "kappa", kappa) &
        )
        call nchandle_error(&
          nf90_get_att(ncid, varids(iaer), "modes", modes_str) &
        )

        modes_str = trim(modes_str)//",inr,inc" ! Manually enable in-rain and in-cloud modes

        call add_tracer(trim(name)//"_inc", &
                        long_name=trim(long_name)//", in-cloud mode", &
                        laero=.true.)
        call add_tracer(trim(name)//"_inr", &
                        long_name=trim(long_name)//", in-rain mode", &
                        laero=.true.)

        call aerosols(iaer) % construct(name, long_name, rho, kappa, modes_str)
      end do
      call nchandle_error(nf90_close(ncid))
    end if

    ! Ok, now add the aerosols to the corresponding modes
    do iaer = 1, nvars
      do imod = 1, size(aerosols(iaer) % modes)
        mode_loc = findloc(modenames, aerosols(iaer) % modes(imod), dim=1) 

        if (.not. mode_loc > 0) then
          call dales_error("Mode '"//aerosols(iaer) % modes(imod)// &
                           "' enabled for aerosol "// &
                           trim(aerosols(iaer) % long_name)// &
                           " does not exist!")
        end if
        call modes(mode_loc) % add_aerosol(aerosols(iaer)) 
      end do
    end do

    ! Finally, allocate memory
    do imod = 1, maxmodes
      call modes(imod) % allocate
    end do

    do imod = 1, maxmodes
      if (modes(imod) % nspecies < 1) cycle
      write(6,*) "Mode: ", modes(imod) % name, " ("//trim(modes(imod) % long_name)//")"
      write(6,'(13A)',advance='no') "    Species: "
      do iaer = 1, modes(imod) % nspecies
        write(6,'(4A)',advance='no') modes(imod) % aerosols(iaer) % ptr % name, " "
      end do
      write(6,*)
    end do

    deallocate(varids)

  end subroutine initaerosol

  !> Deallocates memory
  subroutine exitaerosol
    deallocate(aerosols)
  end subroutine exitaerosol

  !> Construct a mode
  !!
  !! \param name Short name.
  !! \param long_name Long name.
  !! \param sigma_g Geometric standard deviation.
  subroutine mode_construct(self, name, long_name, sigma_g)
    class(mode_t), intent(inout) :: self
    character(3),  intent(in)    :: name
    character(*),  intent(in)    :: long_name
    real(field_r), intent(in)    :: sigma_g

    self % name = trim(name)
    self % long_name = trim(long_name)
    self % sigma_g = sigma_g

    ! In-rain and in-cloud modes do not have activation
    if (name == 'inr' .or. name == 'inc') then
      self % lactivation = .false.
    else
      self % lactivation = .true.
    end if
  
  end subroutine mode_construct

  !> \brief Add an aerosol to a mode.
  !! On the first call, also sets up a tracer for the number concentration.
  !!
  !! \param aerosol Aerosol to add.
  subroutine mode_add_aerosol(self, aerosol)
    class(mode_t),           intent(inout) :: self
    type(aerosol_t), target, intent(in)    :: aerosol

    integer :: itrac
    type(aerosol_ptr_t), allocatable :: tmp(:)

    self % nspecies = self % nspecies + 1

    if (.not. allocated(self % aerosols)) then
      allocate(self % aerosols(1))
    else ! Expand the list of aerosols
      allocate(tmp(self % nspecies))
      tmp(1:self % nspecies - 1) = self % aerosols(1:self % nspecies - 1)
      call move_alloc(tmp, self % aerosols)
    end if

    self % aerosols(self % nspecies) % ptr => aerosol

    ! Check if we have a number concentration setup
    if (.not. associated(self % N % ptr)) then
      call add_tracer("N_"//self % name, long_name="number concentration, "// &
                      trim(self % long_name)//" mode", laero=.true., isv=itrac)
      self % N % ptr => get_tracer(itrac)
    end if                      
  end subroutine mode_add_aerosol

  !> Allocate memory for mass/number concentrations and tendencies.
  subroutine mode_allocate(self)
    class(mode_t), intent(inout) :: self 

    if (self % nspecies < 1) return

    allocate(self % conc(2:i1,2:j1,k1,self % nspecies), &
             self % tend(2:i1,2:j1,k1,self % nspecies))
  end subroutine mode_allocate

  !> \brief Construct an aerosol
  !! 
  !! \param name Short name
  !! \param long_name Long name
  !! \param rho Density
  !! \param kappa Hygroscopicity
  !! \param modes Comma-separated list of mode names the aerosol participates
  !! in. For example, to enable the soluble accumulation and coarse modes, provide
  !! `modes = "acs,cos".
  subroutine aerosol_construct(self, name, long_name, rho, kappa, modes)
    class(aerosol_t), intent(inout) :: self
    character(3),     intent(in)    :: name
    character(64),    intent(in)    :: long_name
    real(field_r),    intent(in)    :: rho, kappa
    character(*),     intent(in)    :: modes

    integer :: nmodes, ierr, imod
    character(3) :: modes_sep(maxmodes)

    self % name = name
    self % long_name = long_name
    self % rho = rho
    self % kappa = kappa

    ! Separate the comma-separated values into parts
    call goSplitString_s(modes, nmodes, modes_sep, ierr, ',')

    self % nmodes = nmodes
    allocate(self % modes(nmodes), self % tracers(nmodes))

    do imod = 1, nmodes
      self % modes(imod) = modes_sep(imod)
      ! Look for the tracers
      self % tracers(imod) % ptr => get_tracer(trim(self % name)//"_"// &
                                               trim(self % modes(imod)))
    end do
  end subroutine aerosol_construct
end module modaerosol