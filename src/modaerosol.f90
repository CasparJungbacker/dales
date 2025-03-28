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
  use modmath,        only: inv_sqrt_two
  use modmicrodata,   only: qcmin, qa_inc, qa_inr, qap_inc, qap_inr
  use modfields,      only: sv0, svp, svm
  use modmicrodata,   only: qcmin, qa_inc, qa_inr, qap_inc, qap_inr, delt
  use modmpi,         only: myid, D_MPI_BCAST, commwrld, mpierr
  use modprecision,   only: field_r
  use modtracers,     only: add_tracer, allocate_tracers, tracer_prop, &
                            get_tracer_index
  use modstat_nc
  use go,             only: goSplitString_s ! TODO: move this to utils
  use utils,          only: dales_error
  use modtimer,       only: timer_tic, timer_toc
  use modlookuptable, only: LT2_t, LT2_create, LT2_set_col, LT2_get_col, LT2_get_col_inline

  implicit none

  private
   
  save

  public :: aerosol_prepare
  public :: aerosol_finalize

  public :: aerosol_get_index_in_mode
  public :: aerosol_get_index_in_cloud
  public :: aerosol_get_type_in_cloud

  !public :: aero_redistribute
  public :: aerosol_cloud_to_rain

  public :: aerosol_names
  public :: n_species_active
  public :: rho_a
  public :: inc_idx
  public :: species_active

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
  character(len=3),  parameter :: aerosol_names(5) = [ &
    'so4', 'ss ', 'pom', 'bc ', 'du ' & 
  ]
  character(len=26), parameter :: aerosol_longnames(5) = [ &
    'sulfate                   ', &
    'sea salt                  ', &
    'particulate organic matter', &
    'black carbon              ', &
    'dust                      ' &
  ]

  integer, parameter, public :: maxmodes = 9
  integer, parameter :: maxspecies = 5
  integer, parameter, public :: iNUS = 1, iAIS = 2, iACS = 3, iCOS = 4, &
                                iAII = 5, iACI = 6, iCOI = 7, iINR = 8, &
                                iINC = 9
  integer, parameter :: iso4 = 1, iss = 2, ipom = 3, ibc = 4, idu = 5

  real(field_r), parameter :: eps0 = 1E-20

  ! Aerosol chemical properties
  real(field_r), parameter :: rho_a(maxspecies) = &
    [1841, 2165, 1800, 1300, 2650]
  !$acc declare copyin(rho_a)
  
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
    integer :: aerosol_type(5)
    integer :: tracer_idx(5)

    real(field_r), allocatable :: n(:,:,:)
    real(field_r), allocatable :: q(:,:,:,:)
    real(field_r), allocatable :: np(:,:,:)
    real(field_r), allocatable :: qp(:,:,:,:)

    integer,       allocatable :: trac_idx(:)
    integer,       allocatable :: aero_idx(:)
    real(field_r), allocatable :: rho(:)     !< Densities
    real(field_r), allocatable :: kappa(:)   !< Hygroscopicities

    ! Fields
    ! Eventually, we want to get rid of these
    !real(field_r), allocatable :: conc(:,:,:,:) !< Mass/number concentrations.
    !real(field_r), allocatable :: tend(:,:,:,:) !< Tendencies.
   contains
    procedure, pass(self) :: construct => mode_construct     !< Constructor.
    procedure, pass(self) :: add_aerosol => mode_add_aerosol !< Add an aerosol to the mode.
    procedure, pass(self) :: allocate => mode_allocate       !< Allocate workspace.
    procedure, pass(self) :: copy_in => mode_copy_in         !< Populate workspace.
    procedure, pass(self) :: copy_out => mode_copy_out       !< Copy out tendencies.
  end type mode_t

  ! Variables
  logical,      public, protected :: laerosol = .false. !< Switch for enabling/disabling interactive aerosols.
  logical                         :: lscavenging = .true.
  type(mode_t), public, target    :: modes(maxmodes)   !< List of modes.

  type(aerosol_t), allocatable                    :: aerosols(:)
  integer,         allocatable, public, protected :: idx_tab(:,:)
  integer, protected :: inc_idx(maxspecies)
  integer            :: inc_type(maxspecies)
  logical, protected :: species_active(maxspecies)
  integer, protected :: n_species_active

  ! Switches for enabling/disabling aerosol species
  logical :: lso4, lss, lpom, lbc, ldu

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

  function aerosol_get_index_in_mode(itype, mode) result(idx)

    integer,      intent(in) :: itype
    type(mode_t), intent(in) :: mode
    !$acc routine seq

    integer :: idx

    ! Linear search
    do idx = 1, mode%nspecies
      if (itype == mode%aerosol_type(idx)) return
    end do

    idx = -1

  end function aerosol_get_index_in_mode

  !> Get array index of aerosol in the in-cloud and in-rain categories
  function aerosol_get_index_in_cloud(itype) result(idx)

    integer, intent(in) :: itype
    !$acc routine seq

    integer :: idx

    idx = inc_idx(itype)

  end function aerosol_get_index_in_cloud

  function aerosol_get_type_in_cloud(idx) result(itype)

    integer, intent(in) :: idx
    !$acc routine seq

    integer :: itype

    itype = inc_type(idx)

  end function aerosol_get_type_in_cloud

  !> Read input files and setup aerosols and M7 modes.
  subroutine initaerosol
    character(len=*), parameter :: routine = modname//"::initaerosol"
    integer       :: imod, ierr, ncid, nvars, iaer, mode_loc
    character(3)  :: name
    character(64) :: long_name
    character(27) :: modes_str
    real(field_r) :: rho, kappa
    character(3)  :: modes_list(maxmodes)
    integer       :: aero_idx_in_mode
    real(field_r), parameter :: sigma_g(maxmodes) = (/ 1.59, 1.59, 1.59, 2.00, &
                                                       1.59, 1.59, 2.00, -999., &
                                                       -999. /)
    real(field_r), parameter :: cldrad(10) = log([5., 10., 15., 20., 25., 30., &
                                                  35., 40., 45., 50.])
    real(field_r), parameter :: rainrate(5) = log([0.01, 0.1, 1., 10., 100.])

    integer, allocatable :: varids(:)

    ! Values for lookup tables
    include "scavenging.inc"

    namelist /NAMAEROSOL/ laerosol, lscavenging, lso4, lss, lpom, lbc, ldu

    ! Read input
    if (myid == 0) then
      ! Namelist
      open(ifnamopt, file=fname_options, status='old', iostat=ierr)
      read(ifnamopt, NAMAEROSOL, iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMAEROSOL')
      close(ifnamopt)
    end if

    call d_mpi_bcast(laerosol, 1, 0, commwrld, mpierr)
    call d_mpi_bcast(lscavenging, 1, 0, commwrld, mpierr)
    call d_mpi_bcast(lso4, 1, 0, commwrld, mpierr)
    call d_mpi_bcast(lss, 1, 0, commwrld, mpierr)
    call d_mpi_bcast(lpom, 1, 0, commwrld, mpierr)
    call d_mpi_bcast(lbc, 1, 0, commwrld, mpierr)
    call d_mpi_bcast(ldu, 1, 0, commwrld, mpierr)

    ! Setup the modes
    do imod = 1, maxmodes
      call mode_construct(modes(imod), name=modenames(imod), &
                          long_name=longnames(imod), sigma_g=sigma_g(imod))
    end do

    if (.not. laerosol) return

    ! Add aerosols to modes

    ! active_matrix represents the following table, where x=1
    !                | NUS | AIS | ACS | COS | AII | ACI | COI
    ! SO4            |  x  |  x  |  x  |  x  |     |     |
    ! Sea Salt       |     |     |  x  |  x  |     |     |
    ! Organic Matter |     |  x  |  x  |  x  |  x  |     |
    ! Black Carbon   |     |  x  |  x  |  x  |  x  |     |
    ! Dust           |     |     |  x  |  x  |     |  x  |  x

    ! Sulphuric acid
    if (lso4) then
      block
        integer :: my_modes(4) = [iNUS, iAIS, iACS, iCOS]
        do imod = 1, size(my_modes)
          call mode_add_aerosol(modes(my_modes(imod)), itype=iso4)
        end do
        call add_tracer('so4_c', long_name='so4 in-cloud mass concentration', unit='kg/kg')
        call add_tracer('so4_r', long_name='so4 in-rain mass concentration', unit='kg/kg')
      end block
    end if

    ! Sea salt
    if (lss) then
      block
        integer :: my_modes(2) = [iACS, iCOS]
        do imod = 1, size(my_modes)
          call mode_add_aerosol(modes(my_modes(imod)), itype=iss)
        end do
        call add_tracer('ss_c', long_name='sea salt in-cloud mass concentration', unit='kg/kg')
        call add_tracer('ss_r', long_name='sea salt in-rain mass concentration', unit='kg/kg')
      end block
    end if

    ! Primary organic matter
    if (lpom) then
      block
        integer :: my_modes(4) = [iAIS, iACS, iCOS, iAII]
        do imod = 1, size(my_modes)
          call mode_add_aerosol(modes(my_modes(imod)), itype=ipom)
        end do
        call add_tracer('pom_c', long_name='organic matter in-cloud mass concentration', unit='kg/kg')
        call add_tracer('pom_r', long_name='organic matter in-rain mass concentration', unit='kg/kg')
      end block
    end if

    ! Black carbon
    if (lbc) then
      block
        integer :: my_modes(4) = [iAIS, iACS, iCOS, iAII]
        do imod = 1, size(my_modes)
          call mode_add_aerosol(modes(my_modes(imod)), itype=ibc)
        end do
        call add_tracer('bc_c', long_name='black carbon in-cloud mass concentration', unit='kg/kg')
        call add_tracer('bc_r', long_name='black carbon in-rain mass concentration', unit='kg/kg')
      end block
    end if

    ! Mineral dust
    if (ldu) then
      block
        integer :: my_modes(4) = [iACS, iCOS, iACI, iCOI]
        do imod = 1, size(my_modes)
          call mode_add_aerosol(modes(my_modes(imod)), itype=idu)
        end do
        call add_tracer('du_c', long_name='mineral dust in-cloud mass concentration', unit='kg/kg')
        call add_tracer('du_r', long_name='mineral dust in-rain mass concentration', unit='kg/kg')
      end block
    end if

    ! Compute indices of aerosols in the array of in-cloud mass
    block

      integer :: i, j
      integer :: index

      species_active(1) = lso4
      species_active(2) = lss
      species_active(3) = lpom
      species_active(4) = lbc
      species_active(5) = ldu

      n_species_active = count(species_active, dim=1)
      
      do i = 1, maxspecies
        index = 0
        if (.not. species_active(i)) cycle
        do j = 1, i
          index = index + merge(1, 0, species_active(j))
        end do
        inc_idx(i) = index
        if (index > 0) inc_type(index) = i
      end do

    end block

    ! Finally, allocate memory
    do imod = 1, maxmodes
      call mode_allocate(modes(imod))
    end do

    ! "Temporary" arrays for in-cloud and in-rain categories of aerosol
    allocate(qa_inc(2:i1,2:j1,k1,n_species_active), &
             qa_inr(2:i1,2:j1,k1,n_species_active), &
             qap_inc(2:i1,2:j1,k1,n_species_active), &
             qap_inr(2:i1,2:j1,k1,n_species_active))

    qa_inc = 0
    qa_inr = 0
    qap_inc = 0
    qap_inr = 0

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

    !$acc enter data copyin(inc_tab_m, inc_tab_m%x1(1:10), inc_tab_m%x2(1:60), &
    !$acc                   inc_tab_m%rows_cols(1:10,1:60,1), &
    !$acc                   inc_tab_n, inc_tab_n%x1(1:10), inc_tab_n%x2(1:60), &
    !$acc                   inc_tab_n%rows_cols(1:10,1:60,1), &
    !$acc                   blc_tab_m, blc_tab_m%x1(1:5), blc_tab_m%x2(1:100), &
    !$acc                   blc_tab_m%rows_cols(1:5,1:100,1), &
    !$acc                   blc_tab_n, blc_tab_n%x1(1:5), blc_tab_n%x2(1:100), &
    !$acc                   blc_tab_n%rows_cols(1:5,1:100,1))
    !$acc enter data copyin(inc_tab_m%n_points, inc_tab_m%x_min, inc_tab_m%inv_fac)
    !$acc enter data copyin(idx_tab)

  end subroutine initaerosol

  !> Deallocates memory
  subroutine exitaerosol
    if (.not. laerosol) return
  end subroutine exitaerosol

  ! Prepares aerosol fields for microphysics calculations
  subroutine aerosol_prepare

    character(len=*), parameter :: routine = modname//'/aerosol_prepare'

    integer :: i, j, k, s, imod
    integer :: itype, idx_c, idx_r

    if (.not. laerosol) return

    call timer_tic(routine, 1)

    ! TODO: Possible optimization: replace this with pointers
    ! need to make sure that the in-cloud species are contiguous in sv array
    ! or: copy them while transposing to (s,k,j,i)

    ! Copy ambient mass and number concentrations to temp fields
    do imod = 1, maxmodes
      call modes(imod) % copy_in(sv0)
    end do

    do s = 1, n_species_active
      itype = aerosol_get_type_in_cloud(s)
      idx_c = get_tracer_index(trim(aerosol_names(itype))//"_c")
      idx_r = get_tracer_index(trim(aerosol_names(itype))//"_r")
      do k = 1, kmax
        do j = 2, j1
          do i = 2, i1
            ! Copy mass concentrations
            qa_inc(i,j,k,s) = max(sv0(i,j,k,idx_c), 0.0_field_r)
            qa_inr(i,j,k,s) = max(sv0(i,j,k,idx_r), 0.0_field_r)
            ! Reset tendency fields
            qap_inc(i,j,k,s) = 0
            qap_inr(i,j,k,s) = 0
          end do
        end do
      end do
    end do

    call timer_toc(routine)

  end subroutine aerosol_prepare

  subroutine aerosol_finalize

    character(len=*), parameter :: routine = modname//'/aerosol_finalize'

    integer :: i, j, k, s, imod
    integer :: itype, idx_c, idx_r

    if (.not. laerosol) return

    call timer_tic(routine, 1)

    if (laerosol) then
      do imod = 1, maxmodes
        call modes(imod) % copy_out(svp, svm, delt)
      end do
    end if

    do s = 1, n_species_active
      itype = aerosol_get_type_in_cloud(s)
      idx_c = get_tracer_index(trim(aerosol_names(itype))//"_c")
      idx_r = get_tracer_index(trim(aerosol_names(itype))//"_r")
      do k = 1, k1
        do j = 2, j1
          do i = 2, i1
            svp(i,j,k,idx_c) = svp(i,j,k,idx_c) + qap_inc(i,j,k,s)
            svp(i,j,k,idx_r) = svp(i,j,k,idx_r) + qap_inr(i,j,k,s)
          end do
        end do
      end do
    end do

    call timer_toc(routine)

  end subroutine aerosol_finalize

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
  
  end subroutine mode_construct

  subroutine mode_add_aerosol(self, itype)
    class(mode_t), intent(inout) :: self
    integer, intent(in) :: itype
    
    integer :: isv

    self%enabled = .true.
    self%nspecies = self%nspecies + 1
    self%aerosol_type(self%nspecies) = itype

    ! Setup a tracer for the mass concentration
    call add_tracer(trim(aerosol_names(itype))//"_"//self%name, &
                    long_name=aerosol_longnames(itype), isv=isv)

    self%tracer_idx(self%nspecies) = isv

  end subroutine mode_add_aerosol

  !> Allocate memory for mass/number concentrations and tendencies.
  subroutine mode_allocate(self)
    class(mode_t), intent(inout) :: self 

    ! Make sure static data is available on GPU, even if we don't use this mode
    !$acc enter data copyin(self, self%enabled)

    if (.not. self%enabled) return

    allocate(self%n(2:i1,2:j1,1:k1), self%np(2:i1,2:j1,1:k1), &
             self%q(2:i1,2:j1,1:k1,self%nspecies), &
             self%qp(2:i1,2:j1,1:k1,self%nspecies))

    self%n(:,:,:) = 0
    self%np(:,:,:) = 0
    self%q(:,:,:,:) = 0
    self%qp(:,:,:,:) = 0

    !$acc enter data copyin(self%n, self%np, self%q, self%qp)

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

    sv_idx = get_tracer_index(self%name//'_n')

    !$acc parallel loop gang vector collapse(3) default(present) async
    do k = 1, k1
      do j = 2, j1
        do i = 2, i1
          self%n(i,j,k) = max(sv(i,j,k,sv_idx), 0.0_field_r)
        end do
      end do
    end do

    ! And the mass concentrations
    do iaer = 1, self%nspecies
      sv_idx = get_tracer_index(aerosol_names(self%aerosol_type(iaer))//'_'//self%name)
      !$acc parallel loop gang vector collapse(3) default(present) async
      do k = 1, k1
        do j = 2, j1
          do i = 2, i1
            self%q(i,j,k,iaer) = max(sv(i,j,k,sv_idx), 0.0_field_r)
          end do
        end do
      end do
    end do

    !$acc wait

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

    sv_idx = get_tracer_index(self%name//'_n')

    !$acc parallel loop gang vector collapse(3) default(present) async
    do k = 1, k1
      do j = 2, j1
        do i = 2, i1
          svp(i,j,k,sv_idx) = svp(i,j,k,sv_idx) + &
                              max(self%np(i,j,k), -svm(i,j,k,sv_idx) / delt)
          self%np(i,j,k) = 0
        end do
      end do
    end do

    do iaer = 1, self%nspecies
      sv_idx = get_tracer_index(aerosol_names(self%aerosol_type(iaer))//'_'//self%name)
      !$acc parallel loop gang vector collapse(3) default(present) async
      do k = 1, k1
        do j = 2, j1
          do i = 2, i1
            svp(i,j,k,sv_idx) = svp(i,j,k,sv_idx) + &
                                max(self%qp(i,j,k,iaer), -svm(i,j,k,sv_idx) / delt)
            self%qp(i,j,k,iaer) = 0
          end do
        end do
      end do
    end do

    !$acc wait

  end subroutine mode_copy_out

  subroutine activation
    use modfields,    only: w0, sv0
    use modmicrodata, only: delt, qcmask, iNc, Ncp, qa_inc, qap_inc

    call activation_pn15(w0, sv0(:,:,:,iNc), delt, modes(iAIS), modes(iACS), &
                         modes(iCOS), Ncp, qap_inc)
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
  subroutine activation_pn15(w, nc, delt, m_ais, m_acs, m_cos, ncp, qap_inc)
    use modmicrodata, only: qcmask
    real(field_r), intent(in)    :: w(2:,2:,:)
    real(field_r), intent(in)    :: nc(2:,2:,:)
    real(field_r), intent(in)    :: delt
    type(mode_t),  intent(inout) :: m_ais
    type(mode_t),  intent(inout) :: m_acs
    type(mode_t),  intent(inout) :: m_cos
    real(field_r), intent(inout) :: ncp(2:,2:,:)
    real(field_r), intent(inout) :: qap_inc(2:,2:,:,:)

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
    real(field_r) :: mass, num, rho
    real(field_r) :: w0

    call timer_tic(routine, 1)

    !$acc parallel loop collapse(3) default(present) &
    !$acc private(N_activated, mode_total_mass, mode_mean_rho, &
    !$acc         mode_median_diameter, f_activated, N_activated, dNcdt, fn, &
    !$acc         fm, tend_n, tend_m, mass, num, rho, w0, my_number, my_target)
    do k = 1, k1
      do j = 2, j1
        do i = 2, i1
          if (qcmask(i,j,k)) then
            N_activated = 0

            ! Step 1: compute how much particles can activate in the AIS mode
            if (m_ais%enabled) then
              mode_total_mass = 0
              mode_mean_rho = 0

              ! Inner loop, perhaps better to reorder arrays to (isv,k,j,i)
              !$acc loop seq
              do s = 1, m_ais % nspecies
                mass = m_ais%q(i,j,k,s)
                rho = rho_a(m_ais%aerosol_type(s))
                mode_total_mass = mode_total_mass + mass
                mode_mean_rho = mode_mean_rho + mass / rho
              end do

              mode_mean_rho = mode_total_mass / (mode_mean_rho + eps0)

              num = m_ais%n(i,j,k)

              mode_median_diameter = ((6 * mode_total_mass) / (pi * &
                                     num * mode_mean_rho + &
                                     eps0))**(1.0_field_r/3) &
                                     * exp((-3 * m_ais % log_sigma_g**2) / 2)
              f_activated = 1 - 0.5_field_r * erfc(-log(2 * r_crit / &
                            mode_median_diameter) * inv_sqrt_two) * m_ais % sigma_g
              N_activated = 1E-6 * f_activated * num
            end if

            ! Step 2: compute tendency of CCN
            if (m_acs%enabled) then
              N_activated = N_activated + 1E-6 * m_acs%n(i,j,k)
            end if

            if (m_cos%enabled) then
              N_activated = N_activated + 1E-6 * m_cos%n(i,j,k)
            end if

            w0 = max(0.0_field_r, w(i,j,k))
            N_activated = max(0.0_field_r, N_activated)
            dNcdt = 1E6 / delt * (0.1 * (w0 * 100 * N_activated / (w0 * 100 + &
                    0.023_field_r * N_activated + eps0)))**1.27_field_r &
                    - 1E-6 * Nc(i,j,k)
            dNcdt = max(dNcdt, 0.0_field_r)

            ! Step 3: move aerosol mass + number from free modes to the in-cloud
            !         mode, starting from the largest mode
            ! COS mode
            if (m_cos % enabled) then
              num = m_cos%n(i,j,k)
              fn = dNcdt * delt / (num + eps0)
              fn = max(min(fn, 1.0_field_r), 0.0_field_r)
              fm = 1 - 0.5_field_r * erfc(erfinv(1 - (2 * fn)) &
                   - 3 * m_cos % log_sigma_g * inv_sqrt_two)
              fm = merge(1.0_field_r, fm, fn > 1.0_field_r)

              tend_n = fn * num / delt
              tend_n = max(0.0_field_r, tend_n)
              m_cos%np(i,j,k) = m_cos%np(i,j,k) - tend_n
              Ncp(i,j,k) = Ncp(i,j,k) + tend_n
              
              !$acc loop seq
              do s = 1, m_cos % nspecies
                mass = m_cos%q(i,j,k,s)
                tend_m = max(0.0_field_r, fm * mass / delt)
                my_target = inc_idx(m_cos%aerosol_type(s)) 
                m_cos%qp(i,j,k,s) = m_cos%qp(i,j,k,s) - tend_m
                qap_inc(i,j,k,my_target) = qap_inc(i,j,k,my_target) + tend_m
              end do

              ! dNcdt = dNcdt - tend_n?
              dNcdt = merge(dNcdt - num / delt, 0.0_field_r, &
                            dNcdt * delt > num)
            end if

            ! ACS mode
            if (m_acs % enabled) then
              num = m_acs%n(i,j,k)
              fn = dNcdt * delt / (num + eps0)
              fn = max(min(fn, 1.0_field_r), 0.0_field_r)
              fm = 1 - 0.5_field_r * &
                erfc(erfinv(1- (2 * fn)) - 3 * m_acs % log_sigma_g * inv_sqrt_two)
              fm = merge(1.0_field_r, fm, fn > 1.0_field_r)

              tend_n = fn * num / delt
              tend_n = max(0.0_field_r, tend_n)
              m_acs%np(i,j,k) = m_acs%np(i,j,k) - tend_n
              Ncp(i,j,k) = Ncp(i,j,k) + tend_n
              
              do s = 1, m_acs % nspecies
                mass = m_acs%q(i,j,k,s)
                tend_m = fm * mass / delt
                tend_m = max(0.0_field_r, tend_m)
                my_target = inc_idx(m_acs%aerosol_type(s))
                m_acs%qp(i,j,k,s) = m_acs%qp(i,j,k,s) - tend_m
                qap_inc(i,j,k,my_target) = qap_inc(i,j,k,my_target) + tend_m
              end do

              dNcdt = merge(dNcdt - num / delt, 0.0_field_r, &
                            dNcdt * delt > num)
            end if

            ! AIS mode
            if (m_ais%enabled) then
              num = m_ais%n(i,j,k)
              fn = dNcdt * delt / (num + eps0)
              fn = max(min(fn, 1.0_field_r), 0.0_field_r)
              fm = 1 - 0.5_field_r * &
                erfc(erfinv(1 - (2 * fn)) - 3 * m_ais % log_sigma_g * inv_sqrt_two)
              fm = merge(1.0_field_r, fm, fn > 1.0_field_r)

              tend_n = fn * num / delt
              tend_n = max(0.0_field_r, tend_n)
              m_ais%np(i,j,k) = m_ais%np(i,j,k) - tend_n
              Ncp(i,j,k) = Ncp(i,j,k) + tend_n
              
              do s = 1, m_ais%nspecies
                mass = m_ais%q(i,j,k,s)
                tend_m = fm * mass / delt
                tend_m = max(0.0_field_r, tend_m)
                my_target = inc_idx(m_ais%aerosol_type(s))
                m_ais%qp(i,j,k,s) = m_ais%qp(i,j,k,s) - tend_m
                qap_inc(i,j,k,my_target) = qap_inc(i,j,k,my_target) + tend_m
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
  subroutine scavenging(ql, sed_qr, Nc, qrmask, rhof, delt, qap_inc, qap_inr)
    real(field_r), intent(in)    :: ql(2:,2:,:)
    real(field_r), intent(in)    :: sed_qr(2:i1,2:j1,1:k1)
    real(field_r), intent(in)    :: Nc(2:i1,2:j1,1:k1)
    logical,       intent(in)    :: qrmask(2:i1,2:j1,1:k1)
    real(field_r), intent(in)    :: rhof(1:k1)
    real(field_r), intent(in)    :: delt
    real(field_r), intent(inout) :: qap_inc(2:,2:,:,:)
    real(field_r), intent(inout) :: qap_inr(2:,2:,:,:)

    character(*), parameter :: routine = modname//"::scavenging"

    integer       :: i, j, k, s, imod
    type(mode_t), pointer  :: mode
    integer       :: my_numb, target_idx
    real(field_r) :: mass, num
    real(field_r) :: mode_mean_mass, mode_mean_rho
    real(field_r) :: mean_cloud_droplet_size, rain_rate
    real(field_r) :: mean_aerosol_radius
    real(field_r) :: f_scav_inc_m, f_scav_inc_n
    real(field_r) :: f_scav_blc_m, f_scav_blc_n
    real(field_r) :: tend_n, tend_m
    real(field_r) :: rainrate
    real(field_r) :: rho
    real(field_r) :: ql_, nc_

    if (.not. lscavenging) return

    call timer_tic(routine, 1)

    do imod = 1, maxmodes ! Exclude in-cloud and in-rain modes
      mode => modes(imod)
      if (.not. mode%enabled) cycle
      !$acc parallel loop gang vector collapse(3) default(present) &
      !$acc private(mode_mean_mass, mode_mean_rho, mass, num, &
      !$acc         mean_cloud_droplet_size, mean_aerosol_radius, &
      !$acc         f_scav_inc_m, f_scav_inc_n, tend_m, tend_n, my_numb, &
      !$acc         target_idx) &
      !$acc async(1)
      do k = 1, kmax
        do j = 2, j1
          do i = 2, i1 
            ql_ = ql(i,j,k)
            nc_ = Nc(i,j,k)
            if (ql_ > qcmin .and. nc_ > 1E3) then
            ! Mode mean properties
            mode_mean_mass = 0
            mode_mean_rho = 0

            !$acc loop seq
            do s = 1, mode%nspecies
              mass = mode%q(i,j,k,s)
              rho = rho_a(mode%aerosol_type(s))
              mode_mean_mass = mode_mean_mass + mass
              mode_mean_rho = mode_mean_rho + (mass / rho)
            end do

            mode_mean_rho = mode_mean_mass / (mode_mean_rho + eps0)

            if (mode_mean_mass > 0 .and. mode_mean_rho > 0) then
              ! Compute mean cloud droplet size and rain rate.
              ! Make sure both stay within the bounds of the lookup table
              mean_cloud_droplet_size = max( &
                1E6_field_r * (3 * ql_ * rhof(k) / &
                       (4 * pi * nc_ * rhow + eps0)), &
                5.001_field_r &
              )
              mean_cloud_droplet_size = min(mean_cloud_droplet_size, &
                                           49.999_field_r)

              ! Compute mean aerosol radius in this mode
              num = mode%n(i,j,k)
              mean_aerosol_radius = 0.5 * (6 * mode_mean_mass &
                / (pi * num * mode_mean_rho + eps0)) &
                **(1.0_field_r / 3) &
                * exp((-3 * mode%log_sigma_g * mode%log_sigma_g) / 2)

              mean_aerosol_radius = min(100 * mean_aerosol_radius, 8E-3_field_r)
              mean_aerosol_radius = max(mean_aerosol_radius, 1E-8_field_r)

              ! Compute how much aerosol is washed out (number and mass)
              f_scav_inc_m = LT2_get_col_inline(inc_tab_m, 1, &
                                         log(mean_cloud_droplet_size), &
                                         log(mean_aerosol_radius))
              f_scav_inc_m = 1E-6 * Nc(i,j,k) * f_scav_inc_m

              f_scav_inc_n = LT2_get_col_inline(inc_tab_n, 1, &
                                         log(mean_cloud_droplet_size), &
                                         log(mean_aerosol_radius))
              f_scav_inc_n = 1E-6 * Nc(i,j,k) * f_scav_inc_n

              f_scav_inc_m = merge(1 / delt, f_scav_inc_m, &
                f_scav_inc_m * delt > 1 .or. f_scav_inc_n * delt > 1)
              f_scav_inc_n = merge(1 / delt, f_scav_inc_n, &
                f_scav_inc_m * delt > 1 .or. f_scav_inc_n * delt > 1)

              ! Remove aerosol from the free modes
              tend_n = f_scav_inc_n * max(0.0_field_r, num)
              mode%np(i,j,k) = mode%np(i,j,k) - tend_n

              do s = 1, mode%nspecies
                mass = mode%q(i,j,k,s)
                tend_m  = f_scav_inc_m * max(0.0_field_r, mass)
                mode%qp(i,j,k,s) = mode%qp(i,j,k,s) - tend_m
                target_idx = inc_idx(mode%aerosol_type(s))
                qap_inc(i,j,k,target_idx) = qap_inc(i,j,k,target_idx) + tend_m
              end do
            end if
          end if
          end do
        end do
      end do
      ! Below-cloud
      !$acc parallel loop gang vector collapse(3) default(present) &
      !$acc private(mode_mean_mass, mode_mean_rho, mass, num, &
      !$acc         rainrate, mean_aerosol_radius, &
      !$acc         f_scav_blc_m, f_scav_blc_n, tend_m, tend_n, my_numb, &
      !$acc         target_idx) &
      !$acc async(2)
      do k = 1, kmax
        do j = 2, j1
          do i = 2, i1 
            if (qrmask(i,j,k) .and. sed_qr(i,j,k)*3600 > 0.01_field_r) then
              ! Mode mean properties
              mode_mean_mass = 0
              mode_mean_rho = 0

              do s = 1, mode%nspecies
                mass = mode%q(i,j,k,s)
                rho = rho_a(mode%aerosol_type(s))
                mode_mean_mass = mode_mean_mass + mass
                mode_mean_rho = mode_mean_rho + (mass / rho)
              end do

              mode_mean_rho = mode_mean_mass / (mode_mean_rho + eps0)

              if (mode_mean_mass > 0 .and. mode_mean_rho > 0) then
                ! Compute mean cloud droplet size and rain rate.
                ! Make sure both stay within the bounds of the lookup table
                mean_cloud_droplet_size = (3 * ql(i,j,k) * rhof(k) / &
                  (4 * pi * Nc(i,j,k) * rhow)) * 1E6
                mean_cloud_droplet_size = min(max(mean_cloud_droplet_size, &
                                                  0.001E-3_field_r), &
                                                  999.999_field_r)

                rainrate = min(max(sed_qr(i,j,k) * 3600, 0.01001_field_r), &
                               99.999_field_r)

                ! Compute mean aerosol radius in this mode
                num = mode%n(i,j,k)
                mean_aerosol_radius = 0.5 * (6 * mode_mean_mass / &
                  (pi * num * mode_mean_rho + eps0)) &
                  **(1.0_field_r / 3) &
                  * exp((-3 * mode%log_sigma_g * mode%log_sigma_g) / 2)

                mean_aerosol_radius = min(0.9999E3_field_r, &
                                          mean_aerosol_radius * 1E6_field_r)
                mean_aerosol_radius = max(mean_aerosol_radius, 1.001E-3_field_r)

                ! Compute how much aerosol is washed out (number and mass)
                f_scav_blc_m = LT2_get_col_inline(blc_tab_m, 1, &
                                           log(rainrate), &
                                           log(mean_aerosol_radius))

                f_scav_blc_n = LT2_get_col_inline(blc_tab_n, 1, &
                                           log(rainrate), &
                                           log(mean_aerosol_radius))

                f_scav_blc_m = merge(1 / delt, f_scav_blc_m, &
                  f_scav_blc_m * delt > 1 .or. f_scav_blc_n * delt > 1)
                f_scav_blc_n = merge(1 / delt, f_scav_blc_n, &
                  f_scav_blc_m * delt > 1 .or. f_scav_blc_n * delt > 1)

                ! Remove aerosol from the free modes
                tend_n = f_scav_blc_n * max(0.0_field_r, num)
                mode%np(i,j,k) = mode%np(i,j,k) - tend_n

                do s = 1, mode%nspecies
                  mass = mode%q(i,j,k,s)
                  tend_m  = f_scav_blc_m * max(0.0_field_r, mass)
                  target_idx = inc_idx(mode%aerosol_type(s))
                  mode%qp(i,j,k,s) = mode%qp(i,j,k,s) - tend_m
                  qap_inr(i,j,k,target_idx) = qap_inr(i,j,k,target_idx) + tend_m
                end do
              end if
            end if
          end do
        end do
      end do

    end do

    call timer_toc(routine)

  end subroutine scavenging

    elemental function erfinv(x) result(p)
      !$acc routine seq
      real(field_r), intent(in) :: x
      real(field_r) :: w, p
  ! TODO: in the following subroutines aero_xxx, make sure that the arguments are similarly named
  !> Moves aerosol mass from the in-cloud category to the in-rain category.
  !! \param qc Cloud water content.
  !! \param qcp Tendency of cloud water due to microphysical processes.
  subroutine aerosol_cloud_to_rain(qc, qcp)

    use modmicrodata, only: qa_c => qa_inc, qap_c => qap_inc, qap_r => qap_inr

    real(field_r), intent(in) :: qc(2:,2:,:)
    real(field_r), intent(in) :: qcp(2:,2:,:)

    character(len=*), parameter :: routine = modname//'aerosol_cloud_to_rain'

    integer       :: i, j, k, s
    real(field_r) :: dqadt

    call timer_tic(routine, 1)

    !$acc parallel loop collapse(4) default(present) private(dqadt)
    do s = 1, n_species_active
      do k = 1, kmax
        do j = 2, j1
          do i = 2, i1
            dqadt = qcp(i,j,k) / qc(i,j,k) * qa_c(i,j,k,s)
            dqadt = merge(dqadt, 0.0_field_r, qc(i,j,k) > qcmin)
            qap_c(i,j,k,s) = qap_c(i,j,k,s) - dqadt
            qap_r(i,j,k,s) = qap_r(i,j,k,s) + dqadt
          end do
        end do
      end do
    end do

    call timer_toc(routine)

  end subroutine aerosol_cloud_to_rain

      ! Safeguard for x close to 0 or 1
      if ( abs( x ) <= 1E-15 ) then
        p = huge(1.0_field_r)
        return
      else if ( abs( abs( x ) - 1.0_field_r )  <= 1E-15 ) then
        p = - huge(1.0_field_r)
        return
      end if

      w = -log( ( 1.0 - x ) * ( 1.0 + x ) )

      if ( w < 6.250000 ) then
          w = w - 3.125000;
          p = -3.6444120640178196996e-21
          p = -1.685059138182016589e-19 + p * w
          p = 1.2858480715256400167e-18 + p * w
          p = 1.115787767802518096e-17 + p * w
          p = -1.333171662854620906e-16 + p * w
          p = 2.0972767875968561637e-17 + p * w
          p = 6.6376381343583238325e-15 + p * w
          p = -4.0545662729752068639e-14 + p * w
          p = -8.1519341976054721522e-14 + p * w
          p = 2.6335093153082322977e-12 + p * w
          p = -1.2975133253453532498e-11 + p * w
          p = -5.4154120542946279317e-11 + p * w
          p = 1.051212273321532285e-09 + p * w
          p = -4.1126339803469836976e-09 + p * w
          p = -2.9070369957882005086e-08 + p * w
          p = 4.2347877827932403518e-07 + p * w
          p = -1.3654692000834678645e-06 + p * w
          p = -1.3882523362786468719e-05 + p * w
          p = 0.0001867342080340571352 + p * w
          p = -0.00074070253416626697512 + p * w
          p = -0.0060336708714301490533 + p * w
          p = 0.24015818242558961693 + p * w
          p = 1.6536545626831027356 + p * w
      else if ( w < 16.00000 ) then
          w = sqrt( w ) - 3.250000;
          p = 2.2137376921775787049e-09
          p = 9.0756561938885390979e-08 + p * w
          p = -2.7517406297064545428e-07 + p * w
          p = 1.8239629214389227755e-08 + p * w
          p = 1.5027403968909827627e-06 + p * w
          p = -4.013867526981545969e-06 + p * w
          p = 2.9234449089955446044e-06 + p * w
          p = 1.2475304481671778723e-05 + p * w
          p = -4.7318229009055733981e-05 + p * w
          p = 6.8284851459573175448e-05 + p * w
          p = 2.4031110387097893999e-05 + p * w
          p = -0.0003550375203628474796 + p * w
          p = 0.00095328937973738049703 + p * w
          p = -0.0016882755560235047313 + p * w
          p = 0.0024914420961078508066 + p * w
          p = -0.0037512085075692412107 + p * w
          p = 0.005370914553590063617 + p * w
          p = 1.0052589676941592334 + p * w
          p = 3.0838856104922207635 + p * w
      else
          w = sqrt( w ) - 5.000000
          p = -2.7109920616438573243e-11
          p = -2.5556418169965252055e-10 + p * w
          p = 1.5076572693500548083e-09 + p * w
          p = -3.7894654401267369937e-09 + p * w
          p = 7.6157012080783393804e-09 + p * w
          p = -1.4960026627149240478e-08 + p * w
          p = 2.9147953450901080826e-08 + p * w
          p = -6.7711997758452339498e-08 + p * w
          p = 2.2900482228026654717e-07 + p * w
          p = -9.9298272942317002539e-07 + p * w
          p = 4.5260625972231537039e-06 + p * w
          p = -1.9681778105531670567e-05 + p * w
          p = 7.5995277030017761139e-05 + p * w
          p = -0.00021503011930044477347 + p * w
          p = -0.00013871931833623122026 + p * w
          p = 1.0103004648645343977 + p * w
          p = 4.8499064014085844221 + p * w
        end if
      p = p * x;      
    end function erfinv

end module modaerosol
