! This file is part of DALES.
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
! Copyright 1993-2024 The DALES team.
!
!> Definitions and functions for M7 aerosols microphysics.
!!
!! \author Marco de Bruine
!! \author Caspar Jungbacker, TU Delft
!!
!! \see https://gmd.copernicus.org/articles/12/5177/2019/
module modaerosol
  use modglobal,      only: ifnamopt, fname_options, checknamelisterror, &
                            cexpnr, i1, j1, k1, ih, jh, pi, nsv, rhow, kmax
  use modmath,        only: erfcinv, inv_sqrt_two
  use modmpi,         only: myid, D_MPI_BCAST, commwrld, mpierr
  use modprecision,   only: field_r
  use modtracer_type, only: tracer_t, tracer_ptr_t
  use modtracers,     only: add_tracer, allocate_tracers, tracer_prop, &
                            get_tracer
  use modstat_nc
  use go,             only: goSplitString_s ! TODO: move this to utils
  use utils,          only: dales_error
  use modtimer,       only: timer_tic, timer_toc
  use modlookuptable, only: LT2_t, LT2_create, LT2_set_col, LT2_get_col

  implicit none

  private
   
  save

  ! Parameters
  character(*),  parameter :: modname = "modaerosol"

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

  integer, parameter, public :: maxmodes = 9
  integer, parameter, public :: iNUS = 1, iAIS = 2, iACS = 3, iCOS = 4, &
                                iAII = 5, iACI = 6, iCOI = 7, iINR = 8, &
                                iINC = 9

  real(field_r), parameter :: sigma_g(maxmodes) = (/ 1.59, 1.59, 1.59, 2.00, &
                                                     1.59, 1.59, 2.00, -999., &
                                                     -999. /)
  real(field_r), parameter :: cldrad(10) = log([5., 10., 15., 20., 25., 30., &
                                                35., 40., 45., 50.])
  real(field_r), parameter :: rainrate(5) = log([0.01, 0.1, 1., 10., 100.])
  real(field_r), parameter :: eps0 = 1E-20
  
  !> Aerosol type, used in initialization
  type aerosol_t
    character(3)  :: name
    character(64) :: long_name
    real(field_r) :: rho
    real(field_r) :: kappa
    integer       :: nmodes
    character(3)  :: modes(maxmodes)
  end type aerosol_t

  !> M7 mode type
  type, public :: mode_t
    ! General properties
    logical            :: enabled       !< Mode is enabled.
    character(3)       :: name          !< Short name of the mode.
    character(64)      :: long_name     !< Full name of the mode.
    integer            :: nspecies = 0  !< Number of species in the mode.
    real(field_r)      :: sigma_g       !< Geometric standard deviation.
    real(field_r)      :: log_sigma_g   !< Log of the GSD
    logical            :: lactivation   !< Aerosols can be activated.

    integer,       allocatable :: trac_idx(:)
    integer,       allocatable :: aero_idx(:)
    real(field_r), allocatable :: rho(:)     !< Densities
    real(field_r), allocatable :: kappa(:)   !< Hygroscopicities

    ! Fields
    ! Eventually, we want to get rid of these
    real(field_r), allocatable :: conc(:,:,:,:) !< Mass/number concentrations.
    real(field_r), allocatable :: tend(:,:,:,:) !< Tendencies.
   contains
    procedure, pass(self) :: construct => mode_construct     !< Constructor.
    procedure, pass(self) :: add_aerosol => mode_add_aerosol !< Add an aerosol to the mode.
    procedure, pass(self) :: allocate => mode_allocate       !< Allocate workspace.
    procedure, pass(self) :: copy_in => mode_copy_in         !< Populate workspace.
    procedure, pass(self) :: copy_out => mode_copy_out       !< Copy out tendencies.
  end type mode_t

  ! Variables
  logical,      public, protected :: laerosol = .false. !< Switch for enabling/disabling interactive aerosols.
  type(mode_t), public            :: modes(maxmodes)   !< List of modes.

  type(aerosol_t), allocatable                    :: aerosols(:)
  integer,         allocatable, public, protected :: idx_tab(:,:)

  ! Scavenging lookup tables
  type(LT2_t) :: inc_tab_m
  type(LT2_t) :: inc_tab_n
  type(LT2_t) :: blc_tab_m
  type(LT2_t) :: blc_tab_n

  ! Procedures
  public :: initaerosol
  public :: exitaerosol
  public :: activation
  public :: scavenging

contains
  !> Read input files and setup aerosols and M7 modes.
  subroutine initaerosol
    integer       :: imod, ierr, ncid, nvars, iaer, mode_loc
    character(3)  :: name
    character(64) :: long_name
    character(27) :: modes_str
    real(field_r) :: rho, kappa
    character(3)  :: modes_list(maxmodes)
    integer       :: aero_idx_in_mode

    integer, allocatable :: varids(:)

    ! Values for lookup tables
    include "scavenging.inc"

    namelist /NAMAEROSOL/ laerosol

    ! Read input
    if (myid == 0) then
      ! Namelist
      open(ifnamopt, file=fname_options, status='old', iostat=ierr)
      read(ifnamopt, NAMAEROSOL, iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMAEROSOL')
      close(ifnamopt)
    end if

    call d_mpi_bcast(laerosol, 1, 0, commwrld, mpierr)

    if (.not. laerosol) return

#ifdef DALES_GPU
    call dales_error("Aerosol microphysics are not supported on GPU!")
#endif

    ! Setup the modes
    do imod = 1, maxmodes
      call modes(imod) % construct(name=modenames(imod), &
                                   long_name=longnames(imod), &
                                   sigma_g=sigma_g(imod))
    end do

    if (myid == 0) then
      call nchandle_error(nf90_open("aerosol."//cexpnr//".nc", NF90_NOWRITE, &
                          ncid))
      call nchandle_error(nf90_inquire(ncid, nVariables=nvars))
    end if

    call d_mpi_bcast(nvars, 1, 0, commwrld, mpierr)
    
    allocate(aerosols(nvars), varids(nvars), idx_tab(maxmodes, nvars))
    idx_tab(:,:) = -1

    if (myid == 0) then
      call nchandle_error(nf90_inq_varids(ncid, nvars, varids))
    end if

    do iaer = 1, nvars
      if (myid == 0) then
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

        ! Manually enable in-rain and in-cloud modes
        modes_str = trim(modes_str)//",inr,inc"
      end if

      call d_mpi_bcast(name, 3, 0, commwrld, mpierr)
      call d_mpi_bcast(long_name, 64, 0, commwrld, mpierr)
      call d_mpi_bcast(rho, 1, 0, commwrld, mpierr)
      call d_mpi_bcast(kappa, 1, 0, commwrld, mpierr)
      call d_mpi_bcast(modes_str, 27, 0, commwrld, mpierr)

      call add_tracer(trim(name)//"_inc", &
                      long_name=trim(long_name)//", in-cloud mode", &
                      laero=.true.)
      call add_tracer(trim(name)//"_inr", &
                      long_name=trim(long_name)//", in-rain mode", &
                      laero=.true.)
      
      aerosols(iaer) % name = name
      aerosols(iaer) % long_name = long_name
      aerosols(iaer) % rho = rho
      aerosols(iaer) % kappa = kappa

      call goSplitString_s(modes_str, aerosols(iaer) % nmodes, &
                           aerosols(iaer) % modes, ierr, ',')
    end do
    
    if (myid == 0) then
      call nchandle_error(nf90_close(ncid))
    end if

    ! Ok, now add the aerosols to the corresponding modes
    do iaer = 1, nvars
      do imod = 1, aerosols(iaer) % nmodes
        mode_loc = findloc(modenames, aerosols(iaer) % modes(imod), dim=1)

        if (.not. mode_loc > 0) then
          call dales_error("Mode '"//aerosols(iaer) % modes(imod)// &
                           "' enabled for aerosol "// &
                           trim(aerosols(iaer) % long_name)// &
                           " does not exist!")
        end if
        call modes(mode_loc) % add_aerosol(aerosols(iaer), iaer, &
                                           aero_idx_in_mode) 

        idx_tab(mode_loc,iaer) = aero_idx_in_mode
      end do
    end do


    ! Finally, allocate memory
    do imod = 1, maxmodes
      call modes(imod) % allocate
    end do

    do imod = 1, maxmodes
      if (modes(imod) % nspecies < 1) cycle
      write(6,*) "Mode: ", modes(imod) % name, &
                 " ("//trim(modes(imod) % long_name)//")"
      write(6,*)
    end do

    deallocate(varids)

    ! Setup the lookup tables for scavenging routines
    inc_tab_m = LT2_create([cldrad(1), aerrad(1)], [cldrad(10), aerrad(60)], &
                           [10, 60], 1)
    call LT2_set_col(inc_tab_m, 1, cldrad, aerrad, scavenging_eff_incloud_m)

    inc_tab_n = LT2_create([cldrad(1), aerrad(1)], [cldrad(10), aerrad(60)], &
                           [10, 60], 1)
    call LT2_set_col(inc_tab_n, 1, cldrad, aerrad, scavenging_eff_incloud_n)

    blc_tab_m = LT2_create([rainrate(1), aerrad_blc(1)], &
                           [rainrate(5), aerrad_blc(100)], &
                           [5, 100], 1) 
    call LT2_set_col(blc_tab_m, 1, rainrate, aerrad_blc, &
                     scavenging_eff_belowcloud_m)

    blc_tab_n = LT2_create([rainrate(1), aerrad_blc(1)], &
                           [rainrate(5), aerrad_blc(100)], &
                           [5, 100], 1) 
    call LT2_set_col(blc_tab_n, 1, rainrate, aerrad_blc, &
                     scavenging_eff_belowcloud_n)


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
    if (sigma_g > 0) then
      self % log_sigma_g = log(sigma_g)
    end if

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
  subroutine mode_add_aerosol(self, aerosol, aero_idx, aero_idx_in_mode)
    class(mode_t),   intent(inout) :: self
    type(aerosol_t), intent(in)    :: aerosol
    integer,         intent(in)    :: aero_idx
    integer,         intent(out)   :: aero_idx_in_mode

    integer :: itrac, isv

    integer,       allocatable :: tmp_idx(:)
    integer,       allocatable :: tmp_aero_idx(:)
    real(field_r), allocatable :: tmp_rho(:)
    real(field_r), allocatable :: tmp_kappa(:)

    self % enabled  = .true. 

    ! First species, also add number concentration
    if (self % nspecies == 0) then
      self % nspecies = self % nspecies + 1

      allocate(self % trac_idx(2), self % rho(1), &
               self % kappa(1), self % aero_idx(1))

      ! Add a tracer for the number concentration
      if (self % name == 'inc') then
        call add_tracer('Nc', laero=.true., isv=isv)
      else if (self % name == 'inr') then
        call add_tracer('Nr', laero=.true., isv=isv)
      else
        call add_tracer(trim(self % name) // "_n", laero=.true., isv=isv)
      end if

      ! Set trac index for number concentration
      self % trac_idx(1) = isv
    else ! Expand the list of aerosols
      self % nspecies = self % nspecies + 1

      ! Some Fortran magic, expands the allocatable arrays
      allocate(tmp_idx(self % nspecies + 1), tmp_rho(self % nspecies), &
               tmp_kappa(self % nspecies), tmp_aero_idx(self % nspecies))
      
      tmp_idx(1:self % nspecies) = self % trac_idx(1:self % nspecies)
      tmp_aero_idx(1:self % nspecies - 1) = &
        self % aero_idx(1:self % nspecies - 1)
      tmp_rho(1:self % nspecies - 1) = self % rho(1:self % nspecies - 1)
      tmp_kappa(1:self % nspecies - 1) = self % kappa(1:self % nspecies - 1)

      call move_alloc(tmp_idx, self % trac_idx)
      call move_alloc(tmp_aero_idx, self % aero_idx)
      call move_alloc(tmp_rho, self % rho)
      call move_alloc(tmp_kappa, self % kappa)
    end if

    ! Add tracer for the new aerosol. This probably already exists, but
    ! add_tracer will give us the trac_idx in that case.
    call add_tracer(trim(aerosol % name) // "_" // trim(self % name), &
                    laero=.true., isv=isv)
 
    self % trac_idx(self % nspecies + 1) = isv
    self % aero_idx(self % nspecies) = aero_idx
    self % rho(self % nspecies) = aerosol % rho
    self % kappa(self % nspecies) = aerosol % kappa

    aero_idx_in_mode = self % nspecies

  end subroutine mode_add_aerosol

  !> Allocate memory for mass/number concentrations and tendencies.
  subroutine mode_allocate(self)
    class(mode_t), intent(inout) :: self 

    if (self % nspecies < 1) return

    allocate(self % conc(2:i1,2:j1,1:k1,self % nspecies + 1), &
             self % tend(2:i1,2:j1,1:k1,self % nspecies + 1))

    self % conc(:,:,:,:) = 0.0_field_r
    self % tend(:,:,:,:) = 0.0_field_r

    !$acc enter data copyin(self%conc(2:i1,2:j1,1:k1,self%nspecies + 1), &
    !$acc&                  self%tend(2:i1,2:j1,1:k1,self%nspecies + 1))

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

    ! And the mass concentrations
    !$acc parallel loop gang vector default(present) private(sv_idx)
    do iaer = 1, self % nspecies + 1
      sv_idx = self % trac_idx(iaer)
      !$acc loop collapse(3)
      do k = 1, k1
        do j = 2, j1
          do i = 2, i1
            self % conc(i,j,k,iaer) = max(sv(i,j,k,sv_idx), 0.0_field_r)
          end do
        end do
      end do
    end do

  end subroutine mode_copy_in

  !> Copies computed tendencies to svp fields.
  !!
  !! \param svp Tracer tendency fields.
  subroutine mode_copy_out(self, svp, svm, delt)
    class(mode_t), intent(inout) :: self
    real(field_r), intent(inout) :: svp(2-ih:i1+ih,2-jh:j1+jh,1:k1,1:nsv)
    real(field_r), intent(in)    :: svm(2-ih:i1+ih,2-jh:j1+jh,1:k1,1:nsv)
    real(field_r), intent(in)    :: delt

    integer :: iaer, sv_idx
    integer :: i, j, k
    character(3) :: name
    
    if (self % nspecies < 1) return

    !$acc parallel loop gang vector default(present) private(sv_idx)
    do iaer = 1, self % nspecies + 1
      sv_idx = self % trac_idx(iaer)
      !$acc loop collapse(3)
      do k = 1, k1
        do j = 2, j1
          do i = 2, i1
            svp(i,j,k,sv_idx) = svp(i,j,k,sv_idx) + self % tend(i,j,k,iaer)
            svp(i,j,k,sv_idx) = max(svp(i,j,k,sv_idx), -svm(i,j,k,sv_idx)/delt)
            self % tend(i,j,k,iaer) = 0.0_field_r
          end do
        end do
      end do
    end do

  end subroutine mode_copy_out

  subroutine activation
    use modfields,    only: w0
    use modmicrodata, only: delt, qcmask

    call activation_pn15(modes(iAIS), modes(iACS), modes(iCOS), modes(iINC), &
                         w0, delt)
  end subroutine activation

  !> \brief Aerosol activation based on updraft velocity
  !!
  !! \see https://doi.org/10.5194/acp-15-9217-2015
  !! 
  !! \param m_ais Aitken soluble mode.
  !! \param m_acs Accumulation soluble mode.
  !! \param m_cos Coarse soluble mode.
  !! \param m_inc In-cloud mode.
  !! \param w Vertical velocity.
  !! \param delt Time step size.
  subroutine activation_pn15(m_ais, m_acs, m_cos, m_inc, w, delt)
    use modmicrodata, only: qcmask
    type(mode_t),  intent(inout) :: m_ais
    type(mode_t),  intent(inout) :: m_acs
    type(mode_t),  intent(inout) :: m_cos
    type(mode_t),  intent(inout) :: m_inc
    real(field_r), intent(in)    :: w(2-ih:i1+ih,2-jh:j1+jh,1:k1)
    real(field_r), intent(in)    :: delt

    integer :: i, j, k, s
    integer :: imod, iaer
    integer :: my_number, my_target

    character(*),  parameter :: routine = modname//"::activation_pn15"
    real(field_r), parameter :: r_crit = 35E-9

    real(field_r) :: &
      mode_total_mass, &
      mode_mean_rho, &
      mode_median_diameter, &
      f_activated, &
      N_activated, &
      dNcdt
    real(field_r) :: fn, fm, tend_n, tend_m
    real(field_r) :: w0

    call timer_tic(routine, 1)

    do k = 1, k1
      do j = 2, j1
        do i = 2, i1
          if (qcmask(i,j,k)) then
            N_activated = 0

            ! Step 1: compute how much particles can activate in the AIS mode
            if (m_ais%enabled) then
              mode_total_mass = 0
              mode_mean_rho = 0

              do s = 1, m_ais % nspecies
                mode_total_mass = mode_total_mass + m_ais%conc(i,j,k,s)
                mode_mean_rho = mode_mean_rho &
                                + m_ais%conc(i,j,k,s) / m_ais % rho(s)
              end do

              mode_mean_rho = mode_total_mass / (mode_mean_rho + eps0)

              mode_median_diameter = ((6 * mode_total_mass) / (pi * &
                                     m_ais%conc(i,j,k,1) * mode_mean_rho * 1E9 + &
                                     eps0)) **(1.0_field_r/3) &
                                     * exp((-3 * m_ais % log_sigma_g**2) / 2)
              f_activated = 1 - 0.5_field_r * erfc(-log(2 * r_crit / &
                            mode_median_diameter) * inv_sqrt_two) * m_ais % sigma_g
              N_activated = 1E-6 * f_activated * m_ais%conc(i,j,k,1)
            end if

            ! Step 2: compute tendency of CCN
            if (m_acs%enabled) then
              N_activated = N_activated + 1E-6 * m_acs%conc(i,j,k,1)
            end if

            if (m_cos%enabled) then
              N_activated = N_activated + 1E-6 * m_cos%conc(i,j,k,1)
            end if

            w0 = max(0.0_field_r, w(i+1,j+1,k))
            N_activated = max(0.0_field_r, N_activated)
            dNcdt = 1E6 / delt * (0.1 * (w0 * 100 * N_activated / (w0 * 100 + &
                    0.023_field_r * N_activated + eps0)))**1.27_field_r &
                    - 1E-6 * m_inc%conc(i,j,k,1)
            dNcdt = max(dNcdt, 0.0_field_r)

            ! Step 3: move aerosol mass + number from free modes to the in-cloud
            !         mode, starting from the largest mode
            ! COS mode
            if (m_cos % enabled) then
              fn = dNcdt * delt / (m_cos%conc(i,j,k,1) + eps0)
              fn = max(min(fn, 1.0_field_r), 0.0_field_r)
              fm = 1 - 0.5_field_r * erfc(erfcinv(2 * fn) &
                   - 3 * m_cos % log_sigma_g * inv_sqrt_two)
              fm = merge(1.0_field_r, fm, fn > 1.0_field_r)

              tend_n = fn * m_cos%conc(i,j,k,1) / delt
              tend_n = max(0.0_field_r, tend_n)
              m_cos%tend(i,j,k,1) = m_cos%tend(i,j,k,1) - tend_n
              m_inc%tend(i,j,k,1) = m_inc%tend(i,j,k,1) + tend_n
              
              do s = 1, m_cos % nspecies
                tend_m = fm * m_cos%conc(i,j,k,s+1) / delt
                tend_m = max(0.0_field_r, tend_m)
                my_number = m_cos % aero_idx(s)
                my_target = idx_tab(iINC, my_number)
                m_cos%tend(i,j,k,s+1) = m_cos%tend(i,j,k,s+1) - tend_m
                m_inc%tend(i,j,k,my_target) = m_inc%tend(i,j,k,my_target) + tend_m
              end do

              dNcdt = merge(dNcdt - m_cos%conc(i,j,k,1) / delt, 0.0_field_r, &
                            dNcdt * delt > m_cos%conc(i,j,k,1))
            end if

            ! ACS mode
            if (m_acs % enabled) then
              fn = dNcdt * delt / (m_acs%conc(i,j,k,1) + eps0)
              fn = max(min(fn, 1.0_field_r), 0.0_field_r)
              fm = 1 - 0.5_field_r * &
                erfc(erfcinv(2 * fn) - 3 * m_acs % log_sigma_g * inv_sqrt_two)
              fm = merge(1.0_field_r, fm, fn > 1.0_field_r)

              tend_n = fn * m_acs%conc(i,j,k,1) / delt
              tend_n = max(0.0_field_r, tend_n)
              m_acs%tend(i,j,k,1) = m_acs%tend(i,j,k,1) - tend_n
              m_inc%tend(i,j,k,1) = m_inc%tend(i,j,k,1) + tend_n
              
              do s = 1, m_acs % nspecies
                tend_m = fm * m_acs%conc(i,j,k,s+1) / delt
                tend_m = max(0.0_field_r, tend_m)
                my_number = m_acs % aero_idx(s)
                my_target = idx_tab(iINC, my_number)
                m_acs%tend(i,j,k,s) = m_acs%tend(i,j,k,s) - tend_m
                m_inc%tend(i,j,k,my_target) = m_inc%tend(i,j,k,my_target) + tend_m
              end do

              dNcdt = merge(dNcdt - m_acs%conc(i,j,k,1) / delt, 0.0_field_r, &
                            dNcdt * delt > m_acs%conc(i,j,k,1))
            end if

            ! AIS mode
            if (m_ais%enabled) then
              fn = dNcdt * delt / (m_ais%conc(i,j,k,1) + eps0)
              fn = max(min(fn, 1.0_field_r), 0.0_field_r)
              fm = 1 - 0.5_field_r * &
                erfc(erfcinv(2 * fn) - 3 * m_ais % log_sigma_g * inv_sqrt_two)
              fm = merge(1.0_field_r, fm, fn > 1.0_field_r)

              tend_n = fn * m_ais%conc(i,j,k,1) / delt
              tend_n = max(0.0_field_r, tend_n)
              m_ais%tend(i,j,k,1) = m_ais%tend(i,j,k,1) - tend_n
              m_inc%tend(i,j,k,1) = m_inc%tend(i,j,k,1) + tend_n
              
              do s = 1, m_ais % nspecies - 1
                tend_m = fm * m_ais%conc(i,j,k,s) / delt
                tend_m = max(0.0_field_r, tend_m)
                my_number = m_ais % aero_idx(s)
                my_target = idx_tab(iINC, my_number)
                m_ais%tend(i,j,k,s) = m_ais%tend(i,j,k,s) - tend_m
                m_inc%tend(i,j,k,my_target) = m_inc%tend(i,j,k,my_target) + tend_m
              end do
            end if
          end if
        end do
      end do
    end do

    call timer_toc(routine)
            
  end subroutine activation_pn15


  !> \brief Computes aerosol scavenging by rain and cloud droplets.
  !!
  !! Based on the approach by Croft et al. (2009) and Croft et al. (2010).
  !! \see https://doi.org/10.5194/acp-10-1511-2010
  !! \see https://doi.org/10.5194/acp-9-4653-2009
  !!
  !! \param ql Liquid water mixing ratio.
  !! \param sed_qr Sedimentation rate of rain droplets.
  !! \param Nc Cloud droplet number concentration.
  !! \param rhof Density of full levels.
  !! \param delt Time step size.
  !! \param modes List of aerosol modes.
  subroutine scavenging(ql, sed_qr, Nc, qrmask, rhof, delt, modes)
    use modmicrodata, only: qcmask
    real(field_r), intent(in)    :: ql(2-ih:i1+ih,2-jh:j1+jh,1:k1)
    real(field_r), intent(in)    :: sed_qr(2:i1,2:j1,1:k1)
    real(field_r), intent(in)    :: Nc(2:i1,2:j1,1:k1)
    logical,       intent(in)    :: qrmask(2:i1,2:j1,1:k1)
    real(field_r), intent(in)    :: rhof(1:k1)
    real(field_r), intent(in)    :: delt
    type(mode_t),  intent(inout) :: modes(maxmodes)

    character(*), parameter :: routine = modname//"::scavenging"

    integer       :: i, j, k, s, imod
    type(mode_t)  :: mode
    integer       :: my_numb, target_idx
    real(field_r) :: mode_mean_mass, mode_mean_rho
    real(field_r) :: mean_cloud_droplet_size, rain_rate
    real(field_r) :: mean_aerosol_radius
    real(field_r) :: f_scav_inc_m, f_scav_inc_n
    real(field_r) :: f_scav_blc_m, f_scav_blc_n
    real(field_r) :: tend_n, tend_m
    real(field_r) :: rainrate

    call timer_tic(routine, 1)

    do imod = 1, maxmodes - 2 ! Exclude in-cloud and in-rain modes
      mode = modes(imod)
      if (.not. mode%enabled) cycle
      do k = 1, kmax
        do j = 2, j1
          do i = 2, i1 
            if (qcmask(i,j,k) .and. Nc(i,j,k) > 1E3) then
            ! Mode mean properties
            mode_mean_mass = 0
            mode_mean_rho = 0

            do s = 1, mode%nspecies
              mode_mean_mass = mode_mean_mass + mode%conc(i,j,k,s+1)
              mode_mean_rho = mode_mean_rho + &
                (mode%conc(i,j,k,s+1) / mode%rho(s))
            end do

            mode_mean_rho = mode_mean_mass / (mode_mean_rho + eps0)

            if (mode_mean_mass > 0 .and. mode_mean_rho > 0) then
              ! Compute mean cloud droplet size and rain rate.
              ! Make sure both stay within the bounds of the lookup table
              mean_cloud_droplet_size = max( &
                1E6 * (3 * ql(i,j,k) * rhof(k) / &
                       (4 * pi * Nc(i,j,k) * rhow + eps0)), &
                5.001 &
              )
              mean_cloud_droplet_size = min(mean_cloud_droplet_size, 49.999)

              ! Compute mean aerosol radius in this mode
              mean_aerosol_radius = 0.5 * (6 * mode_mean_mass &
                / (pi * mode%conc(i,j,k,1) * mode_mean_rho + eps0)) &
                **(1.0_field_r / 3) &
                * exp((-3 * mode%log_sigma_g * mode%log_sigma_g) / 2)

              mean_aerosol_radius = min(100 * mean_aerosol_radius, 8E-3)
              mean_aerosol_radius = max(mean_aerosol_radius, 1E-8)

              ! Compute how much aerosol is washed out (number and mass)
              f_scav_inc_m = LT2_get_col(inc_tab_m, 1, &
                                         log(mean_cloud_droplet_size), &
                                         log(mean_aerosol_radius))
              f_scav_inc_m = 1E-6 * Nc(i,j,k) * f_scav_inc_m

              f_scav_inc_n = LT2_get_col(inc_tab_n, 1, &
                                         log(mean_cloud_droplet_size), &
                                         log(mean_aerosol_radius))
              f_scav_inc_n = 1E-6 * Nc(i1,j1,k) * f_scav_inc_n

              f_scav_inc_m = merge(1 / delt, f_scav_inc_m, &
                f_scav_inc_m * delt > 1 .or. f_scav_inc_n * delt > 1)
              f_scav_inc_n = merge(1 / delt, f_scav_inc_n, &
                f_scav_inc_m * delt > 1 .or. f_scav_inc_n * delt > 1)

              ! Remove aerosol from the free modes
              tend_n = f_scav_inc_n * max(0.0_field_r, mode%conc(i,j,k,1))
              mode%conc(i,j,k,1) = mode%conc(i,j,k,1) - tend_n
              do s = 1, mode%nspecies
                tend_m  = f_scav_inc_m * max(0.0_field_r, mode%conc(i,j,k,s+1))
                mode%tend(i,j,k,s+1) = mode%tend(i,j,k,s+1) - tend_m
                my_numb = mode%aero_idx(s) 
                target_idx = idx_tab(iINC,my_numb)
                modes(iINC)%tend(i,j,k,target_idx) = &
                  modes(iINC)%tend(i,j,k,target_idx) + tend_m
              end do
            end if
          end if
          end do
        end do
      end do
      ! Below-cloud
      do k = 1, kmax
        do j = 2, j1
          do i = 2, i1 
            if (qrmask(i,j,k) .and. sed_qr(i,j,k)*3600 > 0.01_field_r) then
              ! Mode mean properties
              mode_mean_mass = 0
              mode_mean_rho = 0

              do s = 1, mode%nspecies
                mode_mean_mass = mode_mean_mass + mode%conc(i,j,k,s+1)
                mode_mean_rho = mode_mean_rho + (mode%conc(i,j,k,s+1) / &
                                                 mode%rho(s))
              end do

              mode_mean_rho = mode_mean_mass / (mode_mean_rho + eps0)

              if (mode_mean_mass > 0 .and. mode_mean_rho > 0) then
                ! Compute mean cloud droplet size and rain rate.
                ! Make sure both stay within the bounds of the lookup table
                mean_cloud_droplet_size = (3 * ql(i,j,k) * rhof(k) / &
                  (4 * pi * Nc(i,j,k) * rhow)) * 1E6
                mean_cloud_droplet_size = min(max(mean_cloud_droplet_size, 0.01), 99.99)

                rainrate = min(max(sed_qr(i,j,k) * 3600, 0.01001), 99.999)

                ! Compute mean aerosol radius in this mode
                mean_aerosol_radius = 0.5 * (6 * mode_mean_mass / &
                  (pi * mode%conc(i,j,k,1) * mode_mean_rho + eps0)) &
                  **(1.0_field_r / 3) &
                  * exp((-3 * mode%log_sigma_g * mode%log_sigma_g) / 2)

                mean_aerosol_radius = min(0.9999E3_field_r, &
                                          mean_aerosol_radius * 1E6)
                mean_aerosol_radius = max(mean_aerosol_radius, 1.001E-3_field_r)

                ! Compute how much aerosol is washed out (number and mass)
                f_scav_blc_m = LT2_get_col(blc_tab_m, 1, &
                                           log(rainrate), &
                                           log(mean_aerosol_radius))

                f_scav_blc_n = LT2_get_col(blc_tab_n, 1, &
                                           log(rainrate), &
                                           log(mean_aerosol_radius))

                f_scav_blc_m = merge(1 / delt, f_scav_blc_m, &
                  f_scav_blc_m * delt > 1 .or. f_scav_blc_n * delt > 1)
                f_scav_blc_n = merge(1 / delt, f_scav_blc_n, &
                  f_scav_blc_m * delt > 1 .or. f_scav_blc_n * delt > 1)

                ! Remove aerosol from the free modes
                tend_n = f_scav_blc_n * max(0.0_field_r, mode%conc(i,j,k,1))
                mode%tend(i,j,k,1) = mode%tend(i,j,k,1) - tend_n
                do s = 1, mode%nspecies
                  tend_m  = f_scav_blc_m * max(0.0_field_r, mode%conc(i,j,k,s+1))
                  mode%tend(i,j,k,s+1) = mode%tend(i,j,k,s+1) - tend_m
                  my_numb = mode%aero_idx(s) 
                  target_idx = idx_tab(iINR,my_numb)
                  modes(iINR) % tend(i,j,k,target_idx) = &
                    modes(iINR) % tend(i,j,k,target_idx) + tend_m
                end do
              end if
            end if
          end do
        end do
      end do
    end do

    call timer_toc(routine)

  end subroutine scavenging

end module modaerosol
