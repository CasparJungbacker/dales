!  This file is part of DALES.
!
! DALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! DALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!  Copyright 1993-2024 Delft University of Technology, Wageningen
!  University, Utrecht University, KNMI, TNO
!
!> Definitions and functions for M7 aerosols microphysics.
!!
!! \author Marco de Bruine
!! \author Caspar Jungbacker, TU Delft
!!
!! \see https://gmd.copernicus.org/articles/12/5177/2019/
module modaerosol
  use modglobal,      only: ifnamopt, fname_options, checknamelisterror, &
                            cexpnr, i1, j1, k1, ih, jh, pi, nsv
  use modmath,        only: erfcinv, inv_sqrt_two
  use modmpi,         only: myid, D_MPI_BCAST, commwrld, mpierr
  use modprecision,   only: field_r
  use modtracer_type, only: tracer_t, tracer_ptr_t
  use modtracers,     only: add_tracer, allocate_tracers, tracer_prop, &
                            get_tracer
  use modstat_nc
  use go,             only: goSplitString_s ! TODO: move this to utils
  use utils,          only: dales_error

  implicit none

  private
   
  save

  ! Parameters
  integer, public,      parameter :: maxmodes = 9
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
  integer,       public, parameter :: iNUS = 1, iAIS = 2, iACS = 3, iCOS = 4, &
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

    character(3),       allocatable :: modes(:)   !< List of modes this aerosol
                                                  !! participates in.
   contains
    procedure, pass(self) :: construct => aerosol_construct !< Constructor
  end type aerosol_t

  !> Pointer to an aerosol
  !type aerosol_ptr_t
  !  type(aerosol_t), pointer :: ptr => null()
  !end type aerosol_ptr_t

  !> M7 mode type
  type, public :: mode_t
    ! General properties
    character(3)       :: name          !< Short name of the mode.
    character(64)      :: long_name     !< Full name of the mode.
    integer            :: nspecies = 0  !< Number of species in the mode.
    real(field_r)      :: sigma_g       !< Geometric standard deviation.
    logical            :: lactivation   !< Aerosols can be activated.
    type(tracer_ptr_t) :: N             !< Tracer for number concentration.

    type(tracer_ptr_t), allocatable :: species(:) !< List of tracers for species
    real(field_r),      allocatable :: rho(:)     !< Densities
    real(field_r),      allocatable :: kappa(:)   !< Hygroscopicities

    ! Fields
    real(field_r), allocatable :: conc(:,:,:,:) !< Mass/number concentrations.
    real(field_r), allocatable :: tend(:,:,:,:) !< Tendencies.
   contains
    procedure, pass(self) :: construct => mode_construct     !< Constructor.
    procedure, pass(self) :: add_aerosol => mode_add_aerosol !< Add an aerosol to the mode.
    procedure, pass(self) :: allocate => mode_allocate       !< Allocate workspace.
    procedure, pass(self) :: copy_in => mode_copy_in
    procedure, pass(self) :: copy_out => mode_copy_out
  end type mode_t

  ! Variables
  logical, public, protected :: laerosol = .false. !< Switch for enabling/disabling interactive aerosols.
  type(mode_t), public :: modes(maxmodes)    !< List of modes.

  type(aerosol_t), public, allocatable :: aerosols(:)
  integer, public, protected :: naero_trac = 0

  ! Procedures
  public :: initaerosol
  public :: exitaerosol
  public :: activation

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

    call D_MPI_BCAST(laerosol, 1, 0, commwrld, mpierr)

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
        write(6,'(4A)',advance='no') modes(imod) % species(iaer) % ptr % tracname, " "
      end do
      write(6,*)
    end do

    do imod = 1, maxmodes
      naero_trac = naero_trac + modes(imod) % nspecies + 1
    end do

    deallocate(varids)

  end subroutine initaerosol

  !> Deallocates memory
  subroutine exitaerosol
     if (.not. laerosol) return
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

    call add_tracer(trim(name)//"_n", laero=.true.)

    self % N % ptr => get_tracer(trim(name) // "_n")
  
  end subroutine mode_construct

  !> \brief Add an aerosol to a mode.
  !! On the first call, also sets up a tracer for the number concentration.
  !!
  !! \param aerosol Aerosol to add.
  subroutine mode_add_aerosol(self, aerosol)
    class(mode_t),           intent(inout) :: self
    type(aerosol_t), target, intent(in)    :: aerosol

    integer :: itrac
    type(tracer_ptr_t), allocatable :: tmp_spec(:)
    real(field_r),      allocatable :: tmp_rho(:)
    real(field_r),      allocatable :: tmp_kappa(:)

    self % nspecies = self % nspecies + 1

    if (.not. allocated(self % species)) then
      allocate(self % species(1), self % rho(1), self % kappa(1))
    else ! Expand the list of aerosols
      allocate(tmp_spec(self % nspecies), tmp_rho(self % nspecies), &
               tmp_kappa(self % nspecies))

      tmp_spec(1:self % nspecies - 1) = self % species(1:self % nspecies - 1)
      tmp_rho(1:self % nspecies - 1) = self % rho(1:self % nspecies - 1)
      tmp_kappa(1:self % nspecies - 1) = self % kappa(1:self % nspecies - 1)

      call move_alloc(tmp_spec, self % species)
      call move_alloc(tmp_rho, self % rho)
      call move_alloc(tmp_kappa, self % kappa)
    end if

    self % species(self % nspecies) % ptr => &
          get_tracer(trim(aerosol % name) // "_" // trim(self % name))

    ! For efficiency
    self % rho(self % nspecies) = aerosol % rho
    self % kappa(self % nspecies) = aerosol % kappa

  end subroutine mode_add_aerosol

  !> Allocate memory for mass/number concentrations and tendencies.
  subroutine mode_allocate(self)
    class(mode_t), intent(inout) :: self 

    if (self % nspecies < 1) return

    allocate(self % conc(2:i1,2:j1,k1,self % nspecies + 1), &
             self % tend(2:i1,2:j1,k1,self % nspecies + 1))

  end subroutine mode_allocate

  !> Copies aerosol fields to work space.
  !!
  !! \param sv Tracer fields.
  subroutine mode_copy_in(self, sv)
    class(mode_t), intent(inout) :: self    
    real(field_r), intent(in)    :: sv(2-ih:i1+ih,2-jh:j1+jh,1:k1,1:nsv)

    integer :: iaer, sv_idx
    integer :: i, j, k

    if (self % nspecies < 1) return

    ! First, the number concentration
    sv_idx = self % N % ptr % trac_idx
    do k = 1, k1
      do j = 2, j1
        do i = 2, i1
          self % conc(i,j,k,1) = sv(i,j,k,sv_idx)
        end do
      end do
    end do

    ! And the mass concentrations
    do iaer = 1, self % nspecies
      sv_idx = self % species(iaer) % ptr % trac_idx
      do k = 1, k1
        do j = 2, j1
          do i = 2, i1
            self % conc(i,j,k,iaer+1) = sv(i,j,k,sv_idx)
          end do
        end do
      end do
    end do

  end subroutine mode_copy_in

  !> Copies computed tendencies to svp fields.
  !!
  !! \param svp Tracer tendency fields.
  subroutine mode_copy_out(self, svp, svm, delt)
    class(mode_t), intent(in)    :: self
    real(field_r), intent(inout) :: svp(2-ih:i1+ih,2-jh:j1+jh,1:k1,1:nsv)
    real(field_r), intent(in)    :: svm(2-ih:i1+ih,2-jh:j1+jh,1:k1,1:nsv)
    real(field_r), intent(in)    :: delt

    integer :: iaer, sv_idx
    integer :: i, j, k
    
    if (self % nspecies < 1) return

    sv_idx = self % N % ptr % trac_idx
    print *, self % name, "  ", sv_idx
    do k = 1, k1
      do j = 2, j1
        do i = 2, i1
          svp(i,j,k,sv_idx) = svp(i,j,k,sv_idx) + self % tend(i,j,k,1)
          svp(i,j,k,sv_idx) = max(svp(i,j,k,sv_idx), -svm(i,j,k,sv_idx)/delt)
        end do
      end do
    end do

    do iaer = 1, self % nspecies
      sv_idx = self % species(iaer) % ptr % trac_idx
      do k = 1, k1
        do j = 2, j1
          do i = 2, i1
            svp(i,j,k,sv_idx) = svp(i,j,k,sv_idx) + self % tend(i,j,k,iaer+1)
            svp(i,j,k,sv_idx) = max(svp(i,j,k,sv_idx), -svm(i,j,k,sv_idx)/delt)
          end do
        end do
      end do
    end do

  end subroutine mode_copy_out

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
    allocate(self % modes(nmodes))

    do imod = 1, nmodes
      self % modes(imod) = modes_sep(imod)
    end do
  end subroutine aerosol_construct

  subroutine activation
    use modfields,    only: w0
    use modmicrodata, only: delt
    call activation_pn15(modes(iAIS), modes(iACS), modes(iCOS), modes(iINC), w0, delt)
  end subroutine activation

  subroutine activation_pn15(m_ais, m_acs, m_cos, m_inc, w, delt)
    type(mode_t),  intent(inout) :: m_ais
    type(mode_t),  intent(inout) :: m_acs
    type(mode_t),  intent(inout) :: m_cos
    type(mode_t),  intent(inout) :: m_inc
    real(field_r), intent(in)    :: w(2-ih:i1+ih,2-jh:j1+jh,1:k1)
    real(field_r), intent(in)    :: delt

    integer :: i, j, k, s
    integer :: imod, iaer

    real(field_r), parameter :: r_crit = 35E-9

    real(field_r) :: &
      mode_total_mass, &
      mode_mean_rho, &
      mode_median_diameter, &
      f_activated, &
      N_activated, &
      dNcdt
    real(field_r) :: fn, fm, tend_n, tend_m

    associate(qa_ais => m_ais % conc(:,:,:,2:), N_ais => m_ais % conc(:,:,:,1), &
              qap_ais => m_ais % tend(:,:,:,2:), Np_ais => m_ais % tend(:,:,:,1), & 
              qa_acs => m_acs % conc(:,:,:,2:),  N_acs => m_acs % conc(:,:,:,1), &
              qap_acs => m_acs % tend(:,:,:,2:), Np_acs => m_acs % tend(:,:,:,1), &
              qa_cos => m_cos % conc(:,:,:,2:), N_cos => m_cos % conc(:,:,:,1), &
              qap_cos => m_cos % tend(:,:,:,2:), Np_cos => m_cos % tend(:,:,:,1), &
              qa_inc => m_inc % conc(:,:,:,2:), Nc => m_inc % conc(:,:,:,1), &
              qap_inc => m_inc % tend(:,:,:,2:), Ncp => m_inc % conc(:,:,:,1))

    do k = 1, k1
      do j = 2, j1
        do i = 2, i1
          ! Determine the fraction of AIS aerosol with r > 35 nm
          mode_total_mass = sum(qa_ais(i,j,k,:))
          mode_mean_rho = sum(qa_ais(i,j,k,:) / m_ais % rho(:))

          mode_median_diameter = ((6 * mode_total_mass) / (pi * N_ais(i,j,k) * mode_mean_rho * 1E9))**(1.0_field_r/3) &
                                * exp((-3 * log(m_ais % sigma_g)**2) / 2)
          f_activated = 1 - 0.5_field_r * erfc(-log(2 * r_crit / mode_median_diameter) * inv_sqrt_two) * log(m_ais % sigma_g)
          N_activated = 1E-6 * (f_activated * N_ais(i,j,k) + N_acs(i,j,k) + N_cos(i,j,k))
          dNcdt = 1E6 / delt * (0.1 * (w(i,j,k) * 100 * N_activated / (w(i,j,k) * 100 + 0.023_field_r * N_activated)))**1.27_field_r &
                  - 1E-6 * Nc(i,j,k)
          dNcdt = max(dNcdt, 0.0_field_r)

          ! COS mode
          fn = dNcdt * delt / N_cos(i,j,k)
          fm = 1 - 0.5_field_r * erfc(erfcinv(2 * fn) - 3 * log(m_cos % sigma_g) * inv_sqrt_two) ! Overflow for fn -> 1?
          fn = merge(1.0_field_r, fn, fn > 1.0_field_r)
          fm = merge(1.0_field_r, fm, fn > 1.0_field_r)

          tend_n = fn * N_cos(i,j,k) / delt
          tend_n = max(0.0_field_r, tend_n) ! Make sure that we don't have negative activation
          Np_cos(i,j,k) = Np_cos(i,j,k) - tend_n
          Ncp(i,j,k) = Ncp(i,j,k) + tend_n
          
          do s = 1, m_cos % nspecies
            tend_m = fm * qa_cos(i,j,k,s) / delt
            tend_m = max(0.0_field_r, tend_m)
            qap_cos(i,j,k,s) = qap_cos(i,j,k,s) - tend_m
            qap_inc(i,j,k,s) = qap_inc(i,j,k,s) + tend_m
          end do

          dNcdt = merge(dNcdt - N_cos(i,j,k) / delt, 0.0_field_r, dNcdt * delt > N_cos(i,j,k))

          ! ACS mode
          fn = dNcdt * delt / N_acs(i,j,k)
          fm = 1 - 0.5_field_r * erfc(erfcinv(2 * fn) - 3 * log(m_acs % sigma_g) * inv_sqrt_two)
          fn = merge(1.0_field_r, fn, fn > 1.0_field_r)
          fm = merge(1.0_field_r, fm, fn > 1.0_field_r)

          tend_n = fn * N_acs(i,j,k) / delt
          tend_n = max(0.0_field_r, tend_n)
          Np_acs(i,j,k) = Np_acs(i,j,k) - tend_n
          Ncp(i,j,k) = Ncp(i,j,k) + tend_n
          
          do s = 1, m_acs % nspecies
            tend_m = fm * qa_acs(i,j,k,s) / delt
            tend_m = max(0.0_field_r, tend_m)
            qap_acs(i,j,k,s) = qap_acs(i,j,k,s) - tend_m
            qap_inc(i,j,k,s) = qap_inc(i,j,k,s) + tend_m
          end do

          dNcdt = merge(dNcdt - N_acs(i,j,k) / delt, 0.0_field_r, dNcdt * delt > N_acs(i,j,k))

          ! AIS mode
          fn = dNcdt * delt / N_ais(i,j,k)
          fm = 1 - 0.5_field_r * erfc(erfcinv(2 * fn) - 3 * log(m_ais % sigma_g) * inv_sqrt_two)
          fn = merge(1.0_field_r, fn, fn > 1.0_field_r)
          fm = merge(1.0_field_r, fm, fn > 1.0_field_r)

          tend_n = fn * N_ais(i,j,k) / delt
          tend_n = max(0.0_field_r, tend_n)
          Np_ais(i,j,k) = Np_ais(i,j,k) - tend_n
          Ncp(i,j,k) = Ncp(i,j,k) + tend_n
          
          do s = 1, m_ais % nspecies
            tend_m = fm * qa_ais(i,j,k,s) / delt
            tend_m = max(0.0_field_r, tend_m)
            qap_ais(i,j,k,s) = qap_ais(i,j,k,s) - tend_m
            qap_inc(i,j,k,s) = qap_inc(i,j,k,s) + tend_m
          end do
        end do
      end do
    end do

    end associate
            
  end subroutine activation_pn15

end module modaerosol
