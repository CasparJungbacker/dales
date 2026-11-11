!> FFT-based Poisson solver with various FFT backends.
module modpois_fft_solver

  use modglobal,        only: itot, jtot, imax, jmax, kmax, ih, jh, i1, j1, dxi, dyi, PI
  use modfft_backend_t, only: fft_backend_t, FORWARD_FFT, BACKWARD_FFT
  use modfftpack_backend_t, only: fftpack_backend_t
  use modmpi,           only: myidx, myidy, nprocx, nprocy
  use modprecision,     only: pois_r, longint
  use modtranspose,     only: t_transposer

  implicit none

  private

  character(len=*), parameter :: modname = 'modpois_fft_solver'

  public :: fft_solver_init
  public :: fft_solve
  public :: fft_solver_exit
  public :: print_slice

  ! FFT engines
  type(fftpack_backend_t) :: fft_fwd_x !< Forward FFT in the x-direction.
  type(fftpack_backend_t) :: fft_bwd_x !< Backward FFT in the x-direcction.
  type(fftpack_backend_t) :: fft_fwd_y !< Forward FFT in the y-direction.
  type(fftpack_backend_t) :: fft_bwd_y !< Backward FFT in the y-direction.

  type(t_transposer) :: transposer !< Transposer engine.

  real(pois_r), allocatable, target :: p_flat(:)
  real(pois_r), allocatable, target :: work(:)

  ! Helper pointers
  real(pois_r), pointer :: Fp(:,:,:)
  real(pois_r), pointer :: px(:,:,:)
  real(pois_r), pointer :: py(:,:,:)

  integer :: iony !< Number of cells in the x-direction per processor when p is y-aligned.
  integer :: jonx !< Number of cells in the y-direction per processor when p is y-aligned.
  integer :: konx !< Number of cells in the z-direction per processor when p is x or y-aligned.

  real(pois_r), allocatable :: xyrt(:,:) !< Eigenvalues of the Laplacian operator in the x-y plane.
  real(pois_r), allocatable :: d(:,:,:)  !< ?

  real(pois_r), allocatable :: a(:), b(:), c(:) !< Coefficients of the tridiagonal system in the z-direction.
  integer :: ps, pe, qs, qe

contains

  !> Initialize the FFT-based Poisson solver.
  !!
  !! Setup plans, allocate memory and generate Eigenvalues.
  !!
  !! @param[in] lperiodic Logical array indicating periodicity in each direction.
  !! @param[out] p Pointer to the pressure field.
  subroutine fft_solver_init(lperiodic, p)

    logical, intent(in) :: lperiodic(5)

    real(pois_r), pointer, intent(out) :: p(:,:,:)

    character(len=*), parameter :: routine = modname//'/fft_solver_init'

    integer(longint) :: sz !< Required size of the work space.

    ! Initialize the transposer engine.
    ! The transposer is wrapped in a derived type to easily support different backends (e.g. cuDecomp) in the future.
    transposer = t_transposer()

    ! TODO: why let this be handled by the transposer?
    iony = transposer%iony
    jonx = transposer%jonx
    konx = transposer%konx

    call fft_fwd_x%init(FORWARD_FFT, itot, jmax*konx, p_flat, p_flat)
    call fft_bwd_x%init(BACKWARD_FFT, itot, jmax*konx, p_flat, p_flat)
    call fft_fwd_y%init(FORWARD_FFT, jtot, iony*konx, p_flat, p_flat)
    call fft_bwd_y%init(BACKWARD_FFT, jtot, iony*konx, p_flat, p_flat)

    ! For a signal of length N, the FFT will result in floor(N/2) + 1 *complex* coefficients.
    ! Because we are dealing with real input data, the imaginary component of the first coefficient is always zero,
    ! and the imaginary component of the last coefficient is always zero if N is even.
    ! FFTW doesn't store these (at least when using the half-complex format), but other libraries might (like cuFFT).
    ! So, p_flat needs to be long enough to hold all complex coefficients of the longest signal.

    block

      integer :: n_coeff_x, n_coeff_y

      n_coeff_x = itot / 2 + 1 ! Number of *complex* coefficients.
      n_coeff_y = jtot / 2 + 1

      sz = max( &
        (imax + 2 * ih) * (jmax + 2 * jh) * kmax, & ! Before FFT; pressure data + halo's
        (n_coeff_x * 2) * jmax * konx, &            ! After FFT in x-direction.
        (n_coeff_y * 2) * (n_coeff_x * 2) * konx &  ! After FFT in y-direction.
      )

    end block

    allocate(p_flat(sz))

    p(2-ih:i1+ih,2-jh:j1+jh,1:kmax) => p_flat(1:(imax+2*ih)*(jmax+2*jh)*kmax)
    Fp(1:imax,1:jmax,1:konx) => p_flat(1:imax*jmax*konx)
    px(1:itot,1:jmax,1:konx) => p_flat(1:itot*jmax*konx)
    py(1:jtot,1:konx,1:iony) => p_flat(1:jtot*konx*iony)

    ! Work area can be used by the FFT engines and the transposer.
    sz = max( &
      fft_fwd_x%get_worksize(), &
      fft_bwd_x%get_worksize(), &
      fft_fwd_y%get_worksize(), &
      fft_bwd_y%get_worksize(), &
      transposer%get_buffer_size() &
    )

    allocate(work(sz))

    call fft_fwd_x%set_workspace(work)
    call fft_bwd_x%set_workspace(work)
    call fft_fwd_y%set_workspace(work)
    call fft_bwd_y%set_workspace(work)

    allocate( &
      xyrt(iony, jonx), &
      d(iony, jonx, kmax) &
    )

    ! Generate eigenvalues
    block

      integer :: i, j
      integer :: iv, jv

      real(pois_r) :: fac
      real(pois_r), allocatable :: xrt(:)
      real(pois_r), allocatable :: yrt(:)

      allocate(xrt(nprocy*iony), yrt(nprocx*jonx))

      xrt(:) = 0
      yrt(:) = 0

      if (lperiodic(1)) then 
        fac = 1.0_pois_r / (2 * itot)
      else
        fac = 1.0_pois_r / (4 * itot)
      end if

      ! The order of the eigenvalues depends on the FFT implementation.
      ! Here, we order the eigenvalues to match the output of FFTPACK.
      do i = 2, itot - 1, 2
        xrt(i) = -4 * dxi**2 * (sin(real(i, pois_r) * pi * fac))**2
        xrt(i+1) = xrt(i)
      end do

      xrt(1) = 0
      xrt(itot) = -4 * dxi**2

      if (lperiodic(3)) then 
        fac = 1.0_pois_r / (2 * jtot)
      else
        fac = 1.0_pois_r / (4 * jtot)
      end if

      do j = 2, jtot - 1, 2
        yrt(j) = -4 * dyi**2 * (sin(real(j, pois_r) * pi * fac))**2
        yrt(j+1) = yrt(j)
      end do

      yrt(1) = 0
      yrt(jtot) = -4 * dyi**2

      do j = 1, jonx
        do i = 1, iony
          iv = myidy * iony + i
          jv = myidx * jonx + j
          xyrt(i,j) = xrt(iv) + yrt(jv)
        end do
      end do

    end block

    allocate(a(kmax), b(kmax), c(kmax))
    ps = 1
    pe = iony
    qs = 1
    qe = jonx

  end subroutine fft_solver_init

  subroutine fft_solve(p)

    real(pois_r), pointer, intent(inout) :: p(:,:,:)

    character(len=*), parameter :: routine = modname//'/fft_solve'

    call transposer%z_to_x(p, px, work)
    call fft_fwd_x%exec(p_flat)
    call transposer%x_to_y(px, py, work)
    call fft_fwd_y%exec(p_flat)
    call transposer%y_to_z(py, Fp, work)

    call solmpj

    call transposer%z_to_y(Fp, py, work)
    call fft_bwd_y%exec(p_flat)
    call transposer%y_to_x(py, px, work)
    call fft_bwd_x%exec(p_flat)
    call transposer%x_to_z(px, p, work)

    ! Normalize the result of the inverse FFTs.
    p(2-ih:i1+ih,2-jh:j1+jh,1:kmax) = p(2-ih:i1+ih,2-jh:j1+jh,1:kmax) / real(itot*jtot, pois_r)

  end subroutine fft_solve

  subroutine generate_eigenvalues(fft_engine, lperiodic, xyrt)

    class(fft_backend_t), intent(in) :: fft_engine
    logical,              intent(in) :: lperiodic(5)

    real(pois_r), intent(out) :: xyrt(:,:)

    integer      :: i, j
    real(pois_r) :: fac

    real(pois_r) :: xrt(itot)
    real(pois_r) :: yrt(jtot)

    if (lperiodic(1)) then 
      fac = 1.0_pois_r / (2 * itot)
    else
      fac = 1.0_pois_r / (4 * itot)
    end if

    do i = 1, itot
      xrt(i) = -4 * dxi**2 * (sin(2 * real(i-1, pois_r) * pi * fac))**2
    end do

  end subroutine generate_eigenvalues



  subroutine fft_solver_exit()

    deallocate(p_flat, work)

  end subroutine fft_solver_exit

  subroutine solmpj
! version: working version, barrou's removed,
!          correct timing fft's
!          AAPC with MPI-provided routines
!          uses only 2 AAPC, using MPI-rovided routines,
!          to pre-distribute the arrays s.t. complete
!          2-D planes are present on each processor
!          uses ALLTOALL instead of ALLTOALLV
!          ONLY distribution in j-direction allowed

! NOTE: input array p is supposed to have the ip1ray distribution,
!       i.e. the entire range of the first index must be present on
!       each processor

!******************************************************************
!********************  FAST POISSON SOLVER ************************
!*****                                                        *****
!***               P_xx + P_yy + P_zz  =f(x,y,z)                ***
!****                                                         *****
!******************************************************************
!   FOURIER TRANSFORMS IN X AND Y DIRECTION   GIVE:
!   a^2 P + b^2 P + P_zz =  F(x,y,z) = FFT_i [ FTT_j (f(x,y,z))]

!   where a and b are the KNOWN eigenvalues, and P_zz is

!   P_zz =[ P_{i,j,k+1} - 2 P_{i,j,k} +P_{i,j,k-1} ] / (dz * dz)

!   a^2 P + b^2 P + P_zz =
!   [P_{i,j,k+1}-(2+ a^2 + b^2) P_{i,j,k} + P_{i,j,k-1}]/(dz*dz)=F(x,y,z)

!   The equation above results in a tridiagonal system in k which
!   can be solved with Gaussian elemination

    use modglobal, only : kmax,dzf,dzh
    use modfields, only : rhobf, rhobh
    implicit none

    real(pois_r) :: z,bbk
    integer :: i, j, k


  ! Generate tridiagonal matrix

    !$acc parallel loop default(present) async(1)
    do k=1,kmax
      ! SB fixed the coefficients
      a(k)=rhobh(k)  /(dzf(k)*dzh(k  ))
      c(k)=rhobh(k+1)/(dzf(k)*dzh(k+1))
      b(k)=-(a(k)+c(k))
    end do

    !$acc serial default(present) async(1)
    b(1   )=b(1)+a(1)        ! -c(1)
    a(1   )=0.
    b(kmax)=b(kmax)+c(kmax)  ! -a(kmax)
    c(kmax)=0.
    !$acc end serial

    ! SOLVE TRIDIAGONAL SYSTEMS WITH GAUSSIAN ELEMINATION
    ! a(i) x(i-1) + b(i) x(i) + c(i) x(i+1) = d(i)
    !
    ! a(i) => a(k) ~ rho/dz2
    ! b(i) => b(k)+rhobf(k)*xyrtij ~ -2 rho/dz2 -rho/dx2 -rho/dy2
    ! c(i) => c(k) ~ rho/dz2
    ! x(i) => Pij(k)
    ! d(i) => Pij(k) overwritten in-place

    ! Upward sweep i=1
    ! c'(1) = c(1) / b(1)
    ! d'(1) = d(1) / b(1)

    !$acc parallel loop collapse(2) default(present) private(z) async(1)
    do j=qs,qe
      do i=ps,pe
        z        = 1./(b(1)+rhobf(1)*xyrt(i,j))
        d(i,j,1) = c(1)*z
        Fp(i,j,1) = Fp(i,j,1)*z
      end do
    end do

    ! Upward sweep i=2..(n-1)
    ! c'(i) = c(i) / [ b(i) - c'(i-1) a(i) ]
    ! d'(i) = [ d(i) - d'(i-1) a(i) ] / [ b(i) - c'(i-1) a(i) ]
    !$acc parallel loop collapse(2) default(present) private(bbk, z) async(1)
    do  j=qs,qe
      do  i=ps,pe
        !$acc loop seq
        do k=2,kmax-1
          bbk      = b(k)+rhobf(k)*xyrt(i,j)
          z        = 1./(bbk-a(k)*d(i,j,k-1))
          d(i,j,k) = c(k)*z
          Fp(i,j,k) = (Fp(i,j,k)-a(k)*Fp(i,j,k-1))*z
        end do
      end do
    end do

    ! Upward sweep i=n and backsubstitution i=n
    ! x(n) = d'(n)

    !$acc parallel loop collapse(2) default(present) private(bbk, z) async(1)
    do j=qs,qe
      do i=ps,pe
        bbk = b(kmax) + rhobf(kmax)*xyrt(i,j)
        z        = bbk-a(kmax)*d(i,j,kmax-1)
        if(z/=0.) then
          Fp(i,j,kmax) = (Fp(i,j,kmax)-a(kmax)*Fp(i,j,kmax-1))/z
        else
          Fp(i,j,kmax) =0.
        end if
      end do
    end do

    ! Backsubstitution i=n-1..1
    ! x(i) = d'(i) - c'(i) x(i+1)

    !$acc parallel loop collapse(2) default(present) async(1)
    do j=qs,qe
      do i=ps,pe
        !$acc loop seq
        do k=kmax-1,1,-1
          Fp(i,j,k) = Fp(i,j,k)-d(i,j,k)*Fp(i,j,k+1)
        end do
      end do
    end do

    !$acc wait

  end subroutine solmpj

  !> Pretty-print the slice arr(:,:,k) to the screen.
  !!
  !! Rows correspond to the first index, columns to the second index.
  !!
  !! @param[in] arr 3D array.
  !! @param[in] k Index along the last dimension.
  !! @param[in] name Optional label printed in the header.
  subroutine print_slice(arr, k, name)

    real(pois_r), intent(in) :: arr(:,:,:)
    integer, intent(in) :: k
    character(len=*), intent(in), optional :: name

    integer :: i, j

    if (k < 1 .or. k > size(arr, 3)) then
      write(*,'(a,i0,a,i0,a)') 'print_slice: k = ', k, ' out of range [1, ', size(arr, 3), ']'
      return
    end if

    if (present(name)) then
      write(*,'(a,a,i0,a,i0,a,i0,a)') trim(name), '(1:', size(arr, 1), ', 1:', size(arr, 2), ', ', k, ')'
    else
      write(*,'(a,i0,a,i0,a,i0,a)') 'slice (1:', size(arr, 1), ', 1:', size(arr, 2), ', ', k, ')'
    end if

    write(*,'(a6,*(i12))') 'i\j', (j, j = 1, size(arr, 2))
    do i = 1, size(arr, 1)
      write(*,'(i6,*(es12.4))') i, (arr(i, j, k), j = 1, size(arr, 2))
    end do

  end subroutine print_slice

end module modpois_fft_solver