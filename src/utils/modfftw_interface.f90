module modfftw_interface

  use, intrinsic :: iso_c_binding

  public

  include 'fftw3.f03'

  interface fftw_plan_many_r2r_if
    procedure :: d_fftw_plan_many_r2r
    procedure :: d_fftwf_plan_many_r2r
  end interface
  interface fftw_plan_guru_r2r_if
    procedure :: d_fftw_plan_guru_r2r
    procedure :: d_fftwf_plan_guru_r2r
  end interface
  interface fftw_execute_r2r_if
    procedure :: fftw_execute_r2r
    procedure :: fftwf_execute_r2r
  end interface

contains

  subroutine D_fftw_execute_r2r(p, in, out)
      implicit none
      type(C_PTR) :: p
      real(C_DOUBLE), pointer, contiguous :: in(:)
      real(C_DOUBLE), pointer, contiguous :: out(:)
      call fftw_execute_r2r(p, in, out)
  end subroutine

  subroutine D_fftwf_execute_r2r(p, in, out)
      implicit none
      type(C_PTR) :: p
      real(C_FLOAT), pointer, contiguous :: in(:)
      real(C_FLOAT), pointer, contiguous :: out(:)
      call fftwf_execute_r2r(p, in, out)
  end subroutine

  type(c_ptr) function D_fftw_plan_many_r2r(rank,n,howmany,in,inembed,istride,idist,out,onembed,ostride,odist,kind,flags)
      implicit none
      integer(C_INT) :: rank
      integer(C_INT), intent(in) :: n(:)
      integer(C_INT) :: howmany
      real(C_DOUBLE), pointer :: in(:)
      integer(C_INT), intent(in) :: inembed(:)
      integer(C_INT) :: istride
      integer(C_INT) :: idist
      real(C_DOUBLE), pointer :: out(:)
      integer(C_INT), intent(in) :: onembed(:)
      integer(C_INT) :: ostride
      integer(C_INT) :: odist
      integer(C_FFTW_R2R_KIND), intent(in) :: kind(:)
      integer(C_INT) :: flags

      D_fftw_plan_many_r2r = fftw_plan_many_r2r(rank,n,howmany,in,inembed,istride,idist,out,onembed,ostride,odist,kind,flags)
  end function

  type(c_ptr) function D_fftwf_plan_many_r2r(rank,n,howmany,in,inembed,istride,idist,out,onembed,ostride,odist,kind,flags)
      implicit none
      integer(C_INT) :: rank
      integer(C_INT), intent(in) :: n(:)
      integer(C_INT) :: howmany
      real(C_FLOAT), pointer :: in(:)
      integer(C_INT), intent(in) :: inembed(:)
      integer(C_INT) :: istride
      integer(C_INT) :: idist
      real(C_FLOAT), pointer :: out(:)
      integer(C_INT), intent(in) :: onembed(:)
      integer(C_INT) :: ostride
      integer(C_INT) :: odist
      integer(C_FFTW_R2R_KIND), intent(in) :: kind(:)
      integer(C_INT) :: flags

      D_fftwf_plan_many_r2r = fftwf_plan_many_r2r(rank,n,howmany,in,inembed,istride,idist,out,onembed,ostride,odist,kind,flags)
  end function

    type(C_PTR) function d_fftw_plan_guru_r2r(rank,dims,howmany_rank,howmany_dims,in,out,kind,flags)
      integer(C_INT) :: rank
      type(fftw_iodim), intent(in) :: dims(:)
      integer(C_INT) :: howmany_rank
      type(fftw_iodim), intent(in) :: howmany_dims(:)
      real(C_DOUBLE) :: in(:)
      real(C_DOUBLE) :: out(:)
      integer(C_FFTW_R2R_KIND), intent(in) :: kind(:)
      integer(C_INT) :: flags
      d_fftw_plan_guru_r2r = fftw_plan_guru_r2r(rank,dims,howmany_rank,howmany_dims,in,out,kind,flags)
    end function

    type(C_PTR) function d_fftwf_plan_guru_r2r(rank,dims,howmany_rank,howmany_dims,in,out,kind,flags)
      integer(C_INT) :: rank
      type(fftw_iodim), intent(in) :: dims(:)
      integer(C_INT) :: howmany_rank
      type(fftw_iodim), intent(in) :: howmany_dims(:)
      real(C_FLOAT) :: in(:)
      real(C_FLOAT) :: out(:)
      integer(C_FFTW_R2R_KIND), intent(in) :: kind(:)
      integer(C_INT) :: flags

      !fftw_iodim and fftwf_iodim are exactly the same, fortran needs some convincing
      type(fftwf_iodim) :: dimsf(size(dims))
      type(fftwf_iodim) :: howmany_dimsf(size(howmany_dims))
      dimsf = transfer(dims,dimsf)
      howmany_dimsf = transfer(howmany_dims,dimsf)

      d_fftwf_plan_guru_r2r = fftwf_plan_guru_r2r(rank,dimsf,howmany_rank,howmany_dimsf,in,out,kind,flags)
    end function

end module modfftw_interface