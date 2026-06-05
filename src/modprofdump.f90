!> Dumps a vertical profile at a user-selected horizontal location.
module modprofdump

  use modfields,        only: u0, v0, w0, thl0, qt0, ql0, sv0
  use modglobal,        only: i1, j1, kmax, itot, jtot, imax, jmax, nsv, &
                              ifnamopt, checknamelisterror
  use modmpi,           only: myidx, myidy, myid, d_mpi_bcast, commwrld, mpierr
  use modprecision,     only: field_r
  use modnetcdf_file_t, only: profiles_file_t
  use modstat_nc_files, only: add_output_file, is_sampling_timestep
  use modtracers,       only: tracer_prop
  use fortran_support,  only: nnml_output
  use modlogging,       only: finish

  implicit none

  private

  character(len=*), parameter :: modname = 'modprofdump'

  public :: profdump_read_namelist
  public :: initprofdump
  public :: profdump

  type(profiles_file_t) :: ofile    !< Output file object.
  integer               :: ofile_id !< File ID in output file manager.

  logical       :: lprofdump = .false.            !< Enable profile dumping
  integer       :: i_glob = 1                     !< Global i-index of sampled column (1..itot)
  integer       :: j_glob = 1                     !< Global j-index of sampled column (1..jtot)
  real          :: dt_sample = 0.0                !< Sampling interval [s]
  character(80) :: profdump_fname = 'profdump.nc' !< Output filename

  logical :: is_owner = .false.
  integer :: i_loc = 0
  integer :: j_loc = 0

  character(len=16), allocatable :: sv_names(:)

contains

  !> Read profile dump namelist and broadcast settings.
  subroutine profdump_read_namelist(nml_filename)

    character(len=*), intent(in) :: nml_filename

    integer :: ierr

    namelist /profdump/ lprofdump, i_glob, j_glob, dt_sample

    if (myid == 0) then
      open(ifnamopt, file=nml_filename, status='old', iostat=ierr)
      read(ifnamopt, profdump, iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'profdump')
      write(nnml_output, profdump)
      close(ifnamopt)
    end if

    call d_mpi_bcast(lprofdump,       1, 0, commwrld, mpierr)
    call d_mpi_bcast(i_glob,      1, 0, commwrld, mpierr)
    call d_mpi_bcast(j_glob,      1, 0, commwrld, mpierr)
    call d_mpi_bcast(dt_sample,     1, 0, commwrld, mpierr)

  end subroutine profdump_read_namelist

  !> Initialize profile output and register file/variables.
  subroutine initprofdump()

    character(len=*), parameter :: routine = modname//'/initprofdump'

    integer :: owner_x, owner_y
    integer :: n

    if (.not. lprofdump) return

    if (dt_sample <= 0.0_field_r) then
      call finish(routine, 'dt_sample should be positive when lprofdump=.true.')
    end if

    if (i_glob < 1 .or. i_glob > itot) then
      call finish(routine, 'i_glob should be in [1,itot]')
    end if

    if (j_glob < 1 .or. j_glob > jtot) then
      call finish(routine, 'j_glob should be in [1,jtot]')
    end if

    owner_x = (i_glob - 1) / imax
    owner_y = (j_glob - 1) / jmax

    i_loc = mod(i_glob - 1, imax) + 2
    j_loc = mod(j_glob - 1, jmax) + 2

    is_owner = (myidx == owner_x) .and. (myidy == owner_y)

    if (.not. is_owner) return

    if (i_loc > i1 .or. j_loc > j1) then
      call finish(routine, 'local indices for profile point are out of bounds')
    end if

    ofile = profiles_file_t(trim(profdump_fname), nz=kmax)
    call add_output_file(ofile, dt_sample, ofile_id)

    call ofile%add_var('u', 'West-East velocity', 'm/s', 'tt')
    call ofile%add_var('v', 'South-North velocity', 'm/s', 'tt')
    call ofile%add_var('w', 'Vertical velocity', 'm/s', 'mt')
    call ofile%add_var('thl', 'Liquid water potential temperature', 'K', 'tt')
    call ofile%add_var('qt', 'Total water specific humidity', 'kg/kg', 'tt')
    call ofile%add_var('ql', 'Liquid water specific humidity', 'kg/kg', 'tt')

    if (nsv > 0) then
      allocate(sv_names(nsv))
      do n = 1, nsv
        sv_names(n) = trim(tracer_prop(n)%tracname)
        call ofile%add_var(trim(sv_names(n)), trim(tracer_prop(n)%traclong), &
                           trim(tracer_prop(n)%unit), 'tt')
      end do
    end if

  end subroutine initprofdump

  !> Sample the configured column and store it in the output buffer.
  subroutine profdump()

    integer :: n

    real(field_r), pointer :: profile(:)

    if (.not. is_owner) return

    if (.not. is_sampling_timestep(ofile_id)) return

    call ofile%get_pointer('u', profile)
    profile(:) = u0(i_loc, j_loc, 1:kmax)

    call ofile%get_pointer('v', profile)
    profile(:) = v0(i_loc, j_loc, 1:kmax)

    call ofile%get_pointer('w', profile)
    profile(:) = w0(i_loc, j_loc, 1:kmax)

    call ofile%get_pointer('thl', profile)
    profile(:) = thl0(i_loc, j_loc, 1:kmax)

    call ofile%get_pointer('qt', profile)
    profile(:) = qt0(i_loc, j_loc, 1:kmax)

    call ofile%get_pointer('ql', profile)
    profile(:) = ql0(i_loc, j_loc, 1:kmax)

    do n = 1, nsv
      call ofile%get_pointer(trim(sv_names(n)), profile)
      profile(:) = sv0(i_loc, j_loc, 1:kmax, n)
    end do

  end subroutine profdump

end module modprofdump
