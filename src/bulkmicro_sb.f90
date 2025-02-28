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
!  Copyright 1993-2024 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!
!> Kernels for Seifert-Beheng microphysics
module bulkmicro_sb
  use modaerosol,   only: mode_t, iINC, iINR, iACS, iCOS, n_species_active, &
                          rho_a, inc_idx, rho_a, aerosol_get_index_in_cloud, &
                          aerosol_get_index_in_mode, aerosol_get_type_in_cloud
  use modglobal,    only: ih, jh, i1, j1, k1, nsv, rlv, cp, eps1, pi, rv, &
                          mygamma21, mygamma251
  use modmath,      only: inv_sqrt_two
  use modmicrodata, only: Nc_0, pirhow, qrmin, iqr, inr, rhow, eps0, qcmin, &
                          ncmin, l_mur_cst, mur_cst
  use modprecision, only: field_r
  use modtimer,     only: timer_tic, timer_toc

  implicit none

  private

  character(*),  parameter :: modname = "bulkmicro_sb"

  ! Constants
  ! TODO, maybe read these from a namelist.
  real(field_r), parameter :: &
    a_tvsb = 9.65,    & !< Coefficient in terminal velocity param.
    avf = 0.78,       & !< Constant in ventilation factor.
    b_tvsb = 9.8,     & !< Coefficient in terminal velocity param.
    bvf = 0.308,      & !< Constant in ventilation factor.
    c_Nevap = 0.7,    & !< Coefficient for evaporation.
    c_tvsb = 600.,    & !< Coefficient in terminal velocity param.
    D_eq = 1.1e-3,    & !< Parameter for break-up.
    Dv = 2.4e-5,      & !< Diffusivity of water vapor [m^2/s].
    Dvcmax = 79.2-6,   & !< Max mean diameter of cw.
    !Dvcmax = 50E-6,   & !< Max mean diameter of cw.
    D_s = Dvcmax,     & !< Diameter separating the cloud and precipitation parts of the DSD.
    k_1 = 4.0e2,      & !< k_1 + k_2: coefficient for phi function in autoconversion rate SB2006.
    k_2 = 0.7,        & !< See k_1.
    kappa_r = 60.7,   & !< See eq. 11 in SB2006.
    k_br = 1000.,     & !< Parameter for break-up.
    k_c = 10.58e9,    & !< Long Kernel coefficient SB2006 (k'cc).
    k_cc = 4.44e9,    & !< Cloud selfcollection efficiency SB2006.
    k_l = 5.e-5,      & !< Coefficient for phi function in accretion rate.
    k_r = 5.25,       & !< Kernel SB2006.
    k_rr = 7.12,      & !< See eq. 11 in SB2006.
    Kt = 2.5e-2,      & !< Conductivity of heat [J/(sKm)].
    nu_a = 1.41e-5,   & !< Kinematic viscosity of air.
    rho0 = 1.225,     & !< Reference air density
    Sc_num = 0.71,    & !< Schmidt number.
    sig_gr = 1.5,     & !< GSD of rain drop DSD.
    wfallmax = 9.9,   & !< Terminal velocity (?)
    xcmin = 4.2e-15,  & !< Min mean mass of cw (D = 2.0e-6 m).
    xcmax = 2.6E-10,  & !< Max mean mass of cw.
    !xcmax = 6.5E-11,  & !< Max mean mass of cw.
    xrmin = xcmax,    & !< Min mean mass of pw.
    xrmax = 5.0e-6,   & !< Max mean maxx of pw.
    x_s = xcmax         !< Drop mass separating the cloud and precipitation parts of the DSD.

  ! Procedures
  public :: do_bulkmicro_sb

contains

  subroutine do_bulkmicro_sb
    use modaerosol,       only: modes, iINC, iINR, laerosol
    use modmicrodata,     only: qr, Nc, Nr, iqr, iNr, iNc, thlpmcr, qtpmcr, qcbase, &
                                qcroof, qrbase, qrroof, qcmask, qrmask, qrp, Nrp, &
                                Dvr, xr, lbdr, mur, delt, l_lognormal, l_mur_cst, &
                                mur_cst, precep, sed_qr, Ncp, qap_inc, qap_inr, &
                                qa_inc, qa_inr
    use modfields,        only: rhof, ql0, exnf, qvsl, tmp0, esl, svm, qt0, sv0
    use modglobal,        only: dzf
    use modbulkmicrostat, only: bulkmicrotend

    call calculate_rain_parameters(sv0(:,:,:,iNr), sv0(:,:,:,iqr), rhof, l_mur_cst, mur_cst, qrbase, &
                                   qrroof, qrmask, xr, Dvr, mur, lbdr)
    call bulkmicrotend
    call autoconversion(ql0, sv0(:,:,:,iqr), sv0(:,:,:,iNc), qa_inc, qa_inr, exnf, rhof, qcbase, qcroof, delt, &
                        thlpmcr, qtpmcr, qrp, Nrp, laerosol, Ncp, qap_inc, qap_inr)
    call bulkmicrotend
    call accretion(ql0, qr, Nc, Nr, qa_inc, qa_inr, exnf, rhof, &
                   qcbase, qcroof, qrbase, qrroof, laerosol, thlpmcr, qtpmcr, &
                   qrp, Ncp, Nrp, qap_inc, qap_inr)
    call bulkmicrotend
    call evaporation(ql0, qt0, qr, svm(:,:,:,iqr), svm(:,:,:,iNr), qvsl, tmp0, &
                     esl, qa_inr, exnf, rhof, Nr, qrbase, &
                     qrroof, laerosol, delt, qrp, Nrp, qtpmcr, thlpmcr, qap_inr, &
                     modes(iACS), modes(iCOS))
    call bulkmicrotend
#ifdef DALES_GPU
    call sedimentation_rain_gpu(sv0(:,:,:,iqr), sv0(:,:,:,iNr), rhof, dzf, qrbase, qrroof, qrmask, &
                                l_lognormal, l_mur_cst, mur_cst, delt, Dvr, &
                                lbdr, mur, xr, qrp, Nrp, precep, &
                                laerosol=laerosol, m_inr=modes(iINR), &
                                sed_qr_=sed_qr)
#else
    call sedimentation_rain(qr, Nr, rhof, dzf, qrbase, qrroof, qrmask, &
                            l_lognormal, delt, &
                            qrp, Nrp, precep, laerosol, qa_inr, qap_inr, sed_qr, n_species_active)
#endif
    call bulkmicrotend

  end subroutine do_bulkmicro_sb

  !function mur() result(mu)

  !end function mur

  !function lbdr() result(lambda)

  !end function lbdr

  !> Calculate rain DSD integral properties and parameters.
  !!
  !! \param nr rain drop number concentration.
  !! \param qr rain water mixing ratio.
  !! \param rhof Density at full levels.
  !! \param l_mur_cst Switch for selecting constant $\mu$.
  !! \param mur_cst Constant $\mu$ value.
  !! \param qrbase Lowest level with rain.
  !! \param qrroof Highest level with rain.
  !! \param qrmask Rain mask.
  !! \param xr Mean mass of rain drops.
  !! \param Dvr Rain water mean diameter.
  !! \param mur DSD $\mu$ parameter.
  !! \param lbdr DSD $\lambda$ parameter.
  subroutine calculate_rain_parameters(Nr, qr, rhof, l_mur_cst, mur_cst, qrbase, &
                                       qrroof, qrmask, xr, Dvr, mur, lbdr)
    real(field_r), intent(in)  :: Nr(2:i1,2:j1,1:k1)
    real(field_r), intent(in)  :: qr(2:i1,2:j1,1:k1)
    real(field_r), intent(in)  :: rhof(1:k1)

    logical,       intent(in)  :: l_mur_cst
    real(field_r), intent(in)  :: mur_cst

    integer,       intent(in)  :: qrbase, qrroof
    logical,       intent(in)  :: qrmask(2:i1,2:j1,1:k1)

    real(field_r), intent(out) :: xr(2:i1,2:j1,1:k1)
    real(field_r), intent(out) :: Dvr(2:i1,2:j1,1:k1)
    real(field_r), intent(out) :: mur(2:i1,2:j1,1:k1)
    real(field_r), intent(out) :: lbdr(2:i1,2:j1,1:k1)

    integer       :: i, j, k

    if (qrbase > qrroof) return

    call timer_tic('bulkmicro_sb/calculate_rain_parameters', 1)

    if (l_mur_cst) then
      !$acc parallel loop collapse(3) default(present)
      do k = qrbase, qrroof
        do j = 2, j1
          do i = 2, i1
            mur(i,j,k) = mur_cst
          end do
        end do
      end do
    else
      ! mur = f(Dv)
      !$acc parallel loop collapse(3) default(present)
      do k = qrbase, qrroof
        do j = 2, j1
          do i = 2, i1
            if (qrmask(i,j,k)) then
              mur(i,j,k) = min(30.0_field_r, &
                               -1 + 0.008_field_r / (qr(i,j,k) * rhof(k))**0.6_field_r)  ! G09b
            end if
          end do
        end do
      end do
    end if

    !$acc parallel loop collapse(3) default(present)
    do k = qrbase, qrroof
      do j = 2, j1
        do i = 2, i1
          if (qrmask(i,j,k)) then
            xr(i,j,k) = rhof(k) * qr(i,j,k) / Nr(i,j,k)

            ! to ensure xr is within bounds
            xr (i,j,k) = min(max(xr(i,j,k), xrmin), xrmax)
            Dvr(i,j,k) = (xr(i,j,k) / pirhow)**(1./3.)
            lbdr(i,j,k) = ((mur(i,j,k) + 3) * (mur(i,j,k) + 2) * (mur(i,j,k) + 1))**(1.0_field_r/3) / Dvr(i,j,k)
          end if
        end do
      end do
    end do

    call timer_toc('bulkmicro_sb/calculate_rain_parameters')

  end subroutine calculate_rain_parameters

  function calc_mur(qr, rho) result(mu)

    real(field_r), intent(in) :: qr, rho
    real(field_r)             :: mu
    !$acc routine seq

    if (.not. l_mur_cst) then
      mu = min(30.0_field_r, -1 + 0.008_field_r / (qr * rho**0.6_field_r))
    else
      mu = mur_cst
    end if

  end function calc_mur

  function calc_xr(rho, qr, nr) result(xr)

    real(field_r), intent(in) :: rho, qr, nr
    real(field_r)             :: xr
    !$acc routine seq

    xr = rho * qr / nr
    xr = min(max(xr, xrmin), xrmax)

  end function calc_xr

  !> Compute mean rain/cloud drop diameter given its mass.
  !! \param xr Mass of droplet.
  function calc_dvr(xr) result(d)

    real(field_r), intent(in) :: xr
    real(field_r)             :: d
    !$acc routine seq

    d = (xr / pirhow)**(1._field_r/3)

  end function calc_dvr

  function calc_lbdr(mur, dr) result(lambda)

    real(field_r), intent(in) :: mur
    real(field_r), intent(in) :: dr
    real(field_r)             :: lambda
    !$acc routine seq

    lambda = ((mur + 3) * (mur + 2) * (mur + 1))**(1.0_field_r/3) / dr

  end function calc_lbdr

  !> Calculate the autoconversion term.
  !!
  !! \param ql0 Liquid water mixing ratio.
  !! \param qr Rain water mixing ratio.
  !! \param Nc Cloud condensation nucleii.
  !! \param exnf Exner function at full levels.
  !! \param rhof Density at full levels.
  !! \param qcbase Lowest level with cloud.
  !! \param qcroof Highest level with cloud.
  !! \param qcmask Cloud mask.
  !! \param thlpmcr Tendency of $\theta_l$.
  !! \param qtpmcr Tendency of $\q_t$.
  !! \param qrp Tendency of rain water mixing ratio.
  !! \param Nrp Tendency of rain drop number concentration.
  !! \param laerosol Switch for applying autoconversion to aerosols.
  !! \param m_inc In-cloud aerosol mode.
  !! \param m_inr In-rain aerosol mode.
  subroutine autoconversion(ql, qr, Nc, qa_c, qa_r, exnf, rhof, qcbase, &
                            qcroof, delt, thlpmcr, qtpmcr, qrp, Nrp, &
                            laerosol, Ncp, qap_c, qap_r)
    real(field_r), intent(in)    :: ql(2:,2:,:)
    real(field_r), intent(in)    :: qr(2:,2:,:)
    real(field_r), intent(in)    :: Nc(2:,2:,:)
    real(field_r), intent(in)    :: qa_c(2:,2:,:,:)
    real(field_r), intent(in)    :: qa_r(2:,2:,:,:)
    real(field_r), intent(in)    :: exnf(:)
    real(field_r), intent(in)    :: rhof(:)
    integer,       intent(in)    :: qcbase, qcroof
    real(field_r), intent(in)    :: delt
    logical,       intent(in)    :: laerosol

    real(field_r), intent(inout) :: thlpmcr(2:,2:,:)
    real(field_r), intent(inout) :: qtpmcr(2-ih:,2-jh:,:)
    real(field_r), intent(inout) :: qrp(2:,2:,:)
    real(field_r), intent(inout) :: Nrp(2:,2:,:)
    real(field_r), intent(inout) :: Ncp(2:,2:,:)
    real(field_r), intent(inout) :: qap_c(2:,2:,:,:)
    real(field_r), intent(inout) :: qap_r(2:,2:,:,:)

    character(*), parameter :: routine = modname//"::autoconversion"

    integer       :: i, j, k, s, naer
    real(field_r) :: &
      au,   &
      tau,  & !< internal time scale
      phi,  & !< correction function (see SB2001)
      xc,   & !< mean mass of cloud water droplets
      nuc,  & !< width parameter of cloud DSD
      k_au, & !< Coefficient for autoconversion rate
      sc      !< Selfcollection rate
    real(field_r) :: ql_, Nc_
    if (qcbase > qcroof) return

    call timer_tic('bulkmicro_sb01/autoconversion', 1)

    k_au = k_c / (20 * x_s)

    !$acc parallel loop collapse(3) default(present)
    do k = qcbase, qcroof
      do j = 2, j1
        do i = 2, i1
          ql_ = ql(i,j,k)
          Nc_ = Nc(i,j,k)
          if (ql_ > qcmin .and. Nc_ > ncmin) then
              nuc = 1.58_field_r * (rhof(k) * ql_ * 1000.0_field_r) &
                    + 0.72_field_r - 1.0_field_r !G09a
              xc = rhof(k) * ql_ / (Nc_ + eps0)
              au = k_au * (nuc + 2) * (nuc + 4) / (nuc + 1)**2 &
                        * (ql_ * xc)**2 * rho0 ! *rho**2/rho/rho (= 1)

              tau = qr(i,j,k) / (ql_ + qr(i,j,k))
              phi = k_1 * tau**k_2 * (1 - tau**k_2)**3
              au = au * (1 + phi / (1 - tau)**2)

              ! Limit autoconversion to available cloud water
              au = min(ql_ / delt, au)

              qrp(i,j,k) = qrp(i,j,k) + au
              Nrp(i,j,k) = Nrp(i,j,k) + au / x_s * rhof(k)

              qtpmcr(i,j,k) = qtpmcr(i,j,k) - au
              thlpmcr(i,j,k) = thlpmcr(i,j,k) + (rlv / (cp * exnf(k))) * au
              
              if (laerosol) then
                ! When aerosols are enabled, we need to take selfcollection into account
                sc = -k_cc * ((nuc + 2) / (nuc + 1)) * rho0 / rhof(k) &
                     * (ql_ * rhof(k))**2

                Ncp(i,j,k) = Ncp(i,j,k) + sc - au / xc * rhof(k)

                ! Move aerosol mass
                do s = 1, n_species_active
                  qap_c(i,j,k,s) = qap_c(i,j,k,s) - au / ql_ * qa_c(i,j,k,s)
                  qap_r(i,j,k,s) = qap_r(i,j,k,s) + au / ql_ * qa_r(i,j,k,s)
                end do
              end if
           end if
        end do
      end do
    end do

    call timer_toc('bulkmicro_sb01/autoconversion')

  end subroutine autoconversion

  !> Calculate the accretion term.
  !!
  !! \param ql0 Liquid water mixing ratio.
  !! \param qr Rain water mixing ratio.
  !! \param Nr Rain drop number concentration.
  !! \param exnf Exner function at full levels.
  !! \param rhof Density at full levels.
  !! \param qcbase Lowest level with cloud.
  !! \param qcroof Highest level with cloud.
  !! \param qcmask Cloud mask.
  !! \param qrbase Lowest level with rain.
  !! \param qrroof Highest level with rain.
  !! \param qrmask Rain mask.
  !! \param Dvr Rain water mean diameter.
  !! \param lbdr DSD $\lambda$ parameter.
  !! \param thlpmcr Tendency of $\theta_l$.
  !! \param qtpmcr Tendency of total water mixing ratio.
  !! \param qrp Tendency of rain water mixing ratio.
  !! \param Nrp Tendency of rain drop number concentration.
  subroutine accretion(ql, qr, Nc, Nr, qa_c, qa_r, exnf, rhof, &
                       qcbase, qcroof, qrbase, qrroof, laerosol, thlpmcr, &
                       qtpmcr, qrp, Ncp, Nrp, qap_c, qap_r)

    real(field_r), intent(in)    :: ql(2:,2:,:)
    real(field_r), intent(in)    :: qr(2:,2:,:)
    real(field_r), intent(in)    :: Nc(2:,2:,:)
    real(field_r), intent(in)    :: Nr(2:,2:,:)
    real(field_r), intent(in)    :: qa_c(2:,2:,:,:)
    real(field_r), intent(in)    :: qa_r(2:,2:,:,:)
    real(field_r), intent(in)    :: exnf(1:k1)
    real(field_r), intent(in)    :: rhof(1:k1)
    integer,       intent(in)    :: qcbase, qcroof, qrbase, qrroof
    logical,       intent(in)    :: laerosol

    real(field_r), intent(inout) :: thlpmcr(2:,2:,:)
    real(field_r), intent(inout) :: qtpmcr(2-ih:,2-jh:,:)
    real(field_r), intent(inout) :: qrp(2:,2:,:)
    real(field_r), intent(inout) :: Ncp(2:,2:,:)
    real(field_r), intent(inout) :: Nrp(2:,2:,:)
    real(field_r), intent(inout) :: qap_c(2:,2:,:,:)
    real(field_r), intent(inout) :: qap_r(2:,2:,:,:)

    character(len=*), parameter :: routine = modname//"::accretion"

    integer :: i,j,k, s
    integer :: naer

    real(field_r) :: ac, sc, br
    real(field_r) :: xc
    real(field_r) :: phi     !  correction function (see SB2001)
    real(field_r) :: phi_br
    real(field_r) :: tau     !  internal time scale
    real(field_r) :: q_c, q_r
    real(field_r) :: xr, dvr, mur, lbdr

    if (max(qrbase, qcbase) > min(qrroof, qcroof)) return

    call timer_tic('bulkmicro_sb/accretion', 1)

    if (laerosol) naer = size(qa_c, dim=4)

    !$acc parallel loop gang vector collapse(3) default(present) &
    !$acc private(q_c, q_r, tau, phi, ac, xc)
    do k = max(qrbase,qcbase), min(qrroof, qcroof)
      do j = 2, j1
        do i = 2, i1
          q_c = ql(i,j,k)
          q_r = qr(i,j,k)
          if (q_r > qrmin .and. q_c > qcmin) then
            tau = q_r / (q_c + q_r)
            phi = (tau / (tau + k_l))**4
            ac = k_r * rhof(k) * q_c * q_r * phi &
                 * (1.225_field_r / rhof(k))**0.5_field_r

            qrp(i,j,k) = qrp(i,j,k) + ac
            qtpmcr(i,j,k) = qtpmcr(i,j,k) - ac
            thlpmcr(i,j,k) = thlpmcr(i,j,k) + (rlv / (cp * exnf(k))) * ac

            if (laerosol) then
              xc = rhof(k) * q_c / (Nc(i,j,k) + eps0)
              Ncp(i,j,k) = Ncp(i,j,k) - ac / xc

              do s = 1, n_species_active
                qap_c(i,j,k,s) = qap_c(i,j,k,s) - ac / q_c * qa_c(i,j,k,s)
                qap_r(i,j,k,s) = qap_r(i,j,k,s) + ac / q_c * qa_r(i,j,k,s)
              end do
            end if
          end if
        end do
      end do
    end do

    if (qrbase > qrroof) return

    ! TODO: CJ: maybe put this in another subroutine
    !$acc parallel loop gang vector collapse(3) default(present) &
    !$acc private(q_r, sc, phi_br, br)
    do k = qrbase, qrroof
      do j = 2, j1
        do i = 2, i1
          q_r = qr(i,j,k)
          if (q_r > qrmin) then
            xr = calc_xr(rhof(k), q_r, Nr(i,j,k))
            dvr = calc_dvr(xr)
            mur = calc_mur(q_r, rhof(k))
            lbdr = calc_lbdr(mur, dvr)

            sc = k_rr *rhof(k)* q_r * Nr(i,j,k)  &
                 * (1 + kappa_r/lbdr*pirhow**(1./3.))**(-9.)* (1.225/rhof(k))**0.5
            if (Dvr > 0.30E-3) then
              phi_br = k_br * (dvr-D_eq)
              br = (phi_br + 1.) * sc
            else
              br = 0.
            end if

            Nrp(i,j,k) = Nrp(i,j,k) - sc + br
          end if
        end do
      end do
    end do

    call timer_toc('bulkmicro_sb01/accretion')

  end subroutine accretion

  !> Calculate the evaporation term.
  !!
  !! \param ql0 Liquid water mixing ratio.
  !! \param qt0 Total water mixing ratio.
  !! \param qrm Rain water mixing ratio at previous time step.
  !! \param Nrm Rain drop number concentration at previous time step.
  !! \param qvsl Saturation humidity over liquid.
  !! \param tmp0 Temperature.
  !! \param esl Saturation vapor pressure over liquid.
  !! \param exnf Exner function at full levels.
  !! \param rhof Density at full levels.
  !! \param Nr Rain drop number concentration.
  !! \param qrbase Lowest level with rain.
  !! \param qrroof Highest level with rain.
  !! \param qrmask Rain mask.
  !! \param Dvr Rain water mean diameter.
  !! \param lbdr DSD $\lambda$ parameter.
  !! \param mur DSD $\mu$ parameter.
  !! \param xr Mean mass of rain drops.
  !! \param qrp Tendency of rain water mixing ratio.
  !! \param Nrp Tendency of rain drop number concentration.
  !! \param delt Time step size.
  !! \param qtpmcr Tendency of total water mixing ratio.
  !! \param thlpmcr Tendency of $\theta_l$.
  subroutine evaporation(ql, qt, qr, qrm, Nrm, qvsl, tmp0, esl, &
                         qa_r, exnf, rhof, Nr, qrbase, qrroof, laerosol, &
                         delt, qrp, Nrp, qtpmcr, thlpmcr, qap_r, m_acs, m_cos)

    real(field_r), intent(in)    :: ql(2-ih:i1+ih,2-jh:j1+jh,1:k1)
    real(field_r), intent(in)    :: qt(2-ih:i1+ih,2-jh:j1+jh,1:k1)
    real(field_r), intent(in)    :: qr(2:i1,2:j1,1:k1)
    real(field_r), intent(in)    :: qrm(2-ih:i1+ih,2-jh:j1+jh,1:k1)
    real(field_r), intent(in)    :: Nrm(2-ih:i1+ih,2-jh:j1+jh,1:k1)
    real(field_r), intent(in)    :: qvsl(2-ih:i1+ih,2-jh:j1+jh,1:k1)
    real(field_r), intent(in)    :: tmp0(2-ih:i1+ih,2-jh:j1+jh,1:k1)
    real(field_r), intent(in)    :: esl(2-ih:i1+ih,2-jh:j1+jh,1:k1)
    real(field_r), intent(in)    :: qa_r(2:,2:,:,:)
    real(field_r), intent(in)    :: exnf(1:k1)
    real(field_r), intent(in)    :: rhof(1:k1)
    real(field_r), intent(in)    :: Nr(2:i1,2:j1,1:k1)
    integer,       intent(in)    :: qrbase, qrroof
    logical,       intent(in)    :: laerosol
    real(field_r), intent(in)    :: delt

    real(field_r), intent(inout) :: qrp(2:i1,2:j1,1:k1)
    real(field_r), intent(inout) :: Nrp(2:i1,2:j1,1:k1)
    real(field_r), intent(inout) :: qtpmcr(2-ih:i1+ih,2-jh:j1+jh,1:k1)
    real(field_r), intent(inout) :: thlpmcr(2:i1,2:j1,1:k1)
    real(field_r), intent(inout) :: qap_r(2:,2:,:,:)
    type(mode_t),  intent(inout) :: m_acs
    type(mode_t),  intent(inout) :: m_cos

    integer       :: i,j,k,l
    integer       :: numel
    integer       :: idx
    real(field_r) :: F !< ventilation factor
    real(field_r) :: S !< super or undersaturation
    real(field_r) :: G !< cond/evap rate of a drop
    real(field_r) :: evap, Nevap
    real(field_r) :: xr, dvr, mur, lbdr
    real(field_r) :: dm, dm_fac, e, eps, evapt, f_evp, fm, fn, rho, v, dn
    integer :: itype, naer
    real(field_r) :: m_evp, v_evp, rho_evp

    character(*), parameter :: routine = modname//"::evaporation"

    real(field_r), parameter :: Dc = 1.0
    integer :: src_idx, target_idx

    if (qrbase > qrroof) return

    call timer_tic('bulkmicro_sb01/evaporation', 1)

    if(laerosol) naer = size(qa_r, dim=4)

    !$acc parallel loop collapse(3) default(present)
    do k = qrbase, qrroof
      do j = 2, j1
        do i = 2, i1
          if (qr(i,j,k) > qrmin) then
            xr = calc_xr(rhof(k), qr(i,j,k), Nr(i,j,k))
            dvr = calc_dvr(xr)
            mur = calc_mur(qr(i,j,k), rhof(k))
            lbdr = calc_lbdr(mur, dvr)

            numel = nint(mur * 100)
            F = avf * mygamma21(numel)*dvr +  &
               bvf*Sc_num**(1./3.)*(a_tvsb/nu_a)**0.5*mygamma251(numel)*dvr**(3./2.) * &
               (1. - (1.0_field_r/2)   * (b_tvsb / a_tvsb)    *(lbdr / (    c_tvsb + lbdr))**(mur + 2.5_field_r) &
                   - (1.0_field_r/8)   * (b_tvsb / a_tvsb)**2 *(lbdr / (2 * c_tvsb + lbdr))**(mur + 2.5_field_r) &
                   - (1.0_field_r/16)  * (b_tvsb / a_tvsb)**3 *(lbdr / (3 * c_tvsb + lbdr))**(mur + 2.5_field_r) &
                   - (5.0_field_r/128) * (b_tvsb / a_tvsb)**4 *(lbdr / (4 * c_tvsb + lbdr))**(mur + 2.5_field_r) )
            S = min(0.0_field_r, (qt(i,j,k) - ql(i,j,k)) / qvsl(i,j,k) - 1)
            G = (Rv * tmp0(i,j,k)) / (Dv * esl(i,j,k)) + rlv / (Kt * tmp0(i,j,k)) * (rlv / (Rv * tmp0(i,j,k)) - 1)
            G = 1/G

            evap = 2 * pi * Nr(i,j,k) * G * F * S / rhof(k)
            Nevap = c_Nevap * evap * rhof(k) / xr
            
            ! TODO: replace svm reference by qr and nr?
            !if (evap < -svm(i,j,k,iqr)/delt) then
            !  Nevap = - svm(i,j,k,inr)/delt
            !  evap  = - svm(i,j,k,iqr)/delt
            !end if
            if (evap < - qrm(i,j,k) / delt) then
              Nevap = - Nrm(i,j,k) / delt
              evap  = - qrm(i,j,k) / delt
            end if

            qrp(i,j,k) = qrp(i,j,k) + evap
            Nrp(i,j,k) = Nrp(i,j,k) + Nevap

            qtpmcr(i,j,k) = qtpmcr(i,j,k) - evap
            thlpmcr(i,j,k) = thlpmcr(i,j,k) + (rlv / (cp * exnf(k))) * evap

            if (laerosol) then
              ! Fraction of rain that evaporates
              f_evp = - evap / (qr(i,j,k) + eps0) * delt
              f_evp = max(min(f_evp, 1.0_field_r), 0.0_field_r)

              ! Correction factor from Gong et al. (2006)
              eps = (1 - exp(-2 * sqrt(f_evp)) * (1 + 2 * sqrt(f_evp) &
                    + 2 * f_evp + (4.0_field_r/3) * f_evp**(3.0_field_r/2))) &
                    * (1 - f_evp) + f_evp * f_evp

              ! Compute the mass and volume of the evaporated aerosol
              m_evp = 0
              v_evp = 0
              do l = 1, naer
                evapt = eps * f_evp * qa_r(i,j,k,l) / delt
                qap_r(i,j,k,l) = qap_r(i,j,k,l) - evapt
                
                ! Evaporation can't be a source
                evapt = max(0.0_field_r, evapt)

                m_evp = m_evp + evapt
                v_evp = v_evp + evapt / rho_a(l)
              end do

              rho_evp = m_evp / (v_evp + eps0)

              Nevap = max(0._field_r, -1 * Nevap)

              Dn = 1E6 * (6 * m_evp / (pi * Nevap * rho_evp * + eps0))**(1.0_field_r / 3) &
                   * exp(-(3.0_field_r / 2) * log(1.5_field_r)**2)
              Dm = Dn * exp(3 * log(1.5_field_r)**2)

              Fn = 0.5_field_r * erfc(-log(Dc/Dn) / log(1.5_field_r) * inv_sqrt_two)
              Fm = 0.5_field_r * erfc(-log(Dc/Dm) / log(1.5_field_r) * inv_sqrt_two)

              m_acs%np(i,j,k) = m_acs%np(i,j,k) + Fn * Nevap
              m_cos%np(i,j,k) = m_cos%np(i,j,k) + (1 - Fn) * Nevap

              do l = 1, naer
                evapt = eps * f_evp * qa_r(i,j,k,l) / delt
                itype = aerosol_get_type_in_cloud(l)
                target_idx = aerosol_get_index_in_mode(itype, m_acs)
                m_acs%qp(i,j,k,target_idx) = m_acs%qp(i,j,k,target_idx) + Fm * evapt
                target_idx = aerosol_get_index_in_mode(itype, m_cos)
                m_cos%qp(i,j,k,target_idx) = m_cos%qp(i,j,k,target_idx) + (1 - Fm) * evapt
              end do
            end if
          end if
        end do
      end do
    end do

    call timer_toc('bulkmicro_sb/evaporation')

  end subroutine evaporation

  !> Calculate the sedimentation term.
  !!
  !! \param qr Rain water mixing ratio.
  !! \param Nr Rain drop number concentration.
  !! \param rhof Density at full levels.
  !! \param dzf Thickness of vertical levels.
  !! \param qrbase Lowest level with rain.
  !! \param qrroof Highest level with rain.
  !! \param qrmask Rain mask.
  !! \param l_mur_cst Switch for selecting constant $\mu$.
  !! \param mur_cst Constant $\mu$ value.
  !! \param delt Time step size.
  !! \param Dvr Rain water mean diameter.
  !! \param lbdr DSD $\lambda$ parameter.
  !! \param mur DSD $\mu$ parameter.
  !! \param xr Mean mass of rain drops.
  !! \param qrp Tendency of rain water mixing ratio.
  !! \param Nrp Tendency of rain drop number concentration.
  !! \param precep Precipitation.
  subroutine sedimentation_rain(qr, Nr, rhof, dzf, qrbase, qrroof, qrmask, &
                                l_lognormal, delt, qrp, Nrp, precep, laerosol, qa_inr, qap_inr, &
                                sed_qr_, n_species_active)

    real(field_r), intent(in)    :: qr(2:i1,2:j1,1:k1)
    real(field_r), intent(in)    :: Nr(2:i1,2:j1,1:k1)
    real(field_r), intent(in)    :: rhof(1:k1)
    real(field_r), intent(in)    :: dzf(1:k1)

    integer,       intent(inout) :: qrbase
    integer,       intent(in)    :: qrroof
    logical,       intent(inout) :: qrmask(2:i1,2:j1,1:k1)

    logical,       intent(in)    :: l_lognormal
    real(field_r), intent(in)    :: delt

    real(field_r), intent(inout) :: qrp(2:i1,2:j1,1:k1)
    real(field_r), intent(inout) :: Nrp(2:i1,2:j1,1:k1)
    real(field_r), intent(out)   :: precep(2:i1,2:j1,1:k1)

    logical,       intent(in)    :: laerosol 
    real(field_r), intent(in)    :: qa_inr(2:,2:,:,:)
    real(field_r), intent(inout) :: qap_inr(2:,2:,:,:)
    real(field_r), intent(out)   :: sed_qr_(2:i1,2:j1,1:k1)
    integer,       intent(in)    :: n_species_active

    integer       :: i, j, k, jn, sedimbase, s
    integer       :: n_spl      !<  sedimentation time splitting loop
    real(field_r) :: pwcont
    real(field_r) :: delt_inv
    real(field_r) :: Dgr           !<  lognormal geometric diameter
    real(field_r) :: wfall_qr      !<  fall velocity for qr
    real(field_r) :: wfall_Nr      !<  fall velocity for Nr
    real(field_r) :: sed_qr
    real(field_r) :: sed_Nr
    real(field_r) :: xr, dvr, mur, lbdr

    real(field_r), allocatable :: qr_spl(:,:,:), Nr_spl(:,:,:), qa_spl(:,:,:,:)

    real(field_r) :: dt_spl

    call timer_tic('bulkmicro_sb01/sedimentation_rain', 1)

    precep(:,:,:) = 0 ! zero the precipitation flux field
                      ! the update below is not always performed

    if (qrbase > qrroof) return

    allocate(qr_spl(2:i1,2:j1,1:k1))
    allocate(Nr_spl(2:i1,2:j1,1:k1))

    if (laerosol) then
      allocate(qa_spl(2:i1,2:j1,1:k1,n_species_active))
    end if

    n_spl = ceiling(wfallmax * delt / minval(dzf))
    dt_spl = delt / real(n_spl, kind=field_r)

    do jn = 1, n_spl ! time splitting loop
      if (jn == 1) then
        qr_spl(:,:,:) = qr(:,:,:)
        Nr_spl(:,:,:) = Nr(:,:,:)
        if (laerosol) then
          do s = 1, n_species_active
            qa_spl(:,:,:,s) = qa_inr(:,:,:,s)
          end do
        end if
      else
        ! update parameters after the first iteration
        ! a new mask
        qrmask(:,:,:) = (qr_spl(:,:,:) > qrmin) .and. (Nr_spl(:,:,:) > 0)

        ! lower the rain base by one level to include the rain fall
        ! from the previous step
        qrbase = max(1, qrbase - 1)
      end if

      if (l_lognormal) then
        do k = qrbase,qrroof
          do j = 2, j1
            do i = 2, i1
              if (qrmask(i,j,k)) then
                xr = calc_xr(rhof(k), qr_spl(i,j,k), Nr_spl(i,j,k))
                dvr = calc_dvr(xr)

                ! correction for width of DSD
                Dgr = (exp(4.5_field_r * (log(sig_gr))**2))**(-1.0_field_r/3) * Dvr
                sed_qr = sed_flux(Nr_spl(i,j,k), Dgr, log(sig_gr)**2, D_s, 3)
                sed_Nr = 1 / pirhow * sed_flux(Nr_spl(i,j,k), Dgr, log(sig_gr)**2, D_s, 0)

                ! correction for the fact that pwcont.ne. qr_spl
                ! actually in this way for every grid box a fall velocity is determined
                pwcont = liq_cont(Nr_spl(i,j,k), Dgr, log(sig_gr)**2, D_s, 3)       ! note : kg m-3
                if (pwcont > eps1) then
                  sed_qr = (qr_spl(i,j,k) * rhof(k) / pwcont) * sed_qr
                  ! or:
                  ! qr_spl*(sed_qr/pwcont) = qr_spl*fallvel.
                end if

                qr_spl(i,j,k) = qr_spl(i,j,k) - sed_qr * dt_spl / (dzf(k) * rhof(k))
                Nr_spl(i,j,k) = Nr_spl(i,j,k) - sed_Nr * dt_spl / dzf(k)

                if (k > 1) then
                  qr_spl(i,j,k-1) = qr_spl(i,j,k-1) + sed_qr * dt_spl / (dzf(k-1) * rhof(k-1))
                  Nr_spl(i,j,k-1) = Nr_spl(i,j,k-1) + sed_Nr * dt_spl / dzf(k-1)
                end if
                if (jn == 1) then
                  precep(i,j,k) = sed_qr / rhof(k)   ! kg kg-1 m s-1
                end if
              end if ! qr_spl threshold statement
            end do
          end do
        end do
      else
        do k = qrbase, qrroof
          do j = 2, j1
            do i = 2, i1
              if (qrmask(i,j,k)) then
                xr = calc_xr(rhof(k), qr(i,j,k), Nr(i,j,k))
                dvr = calc_dvr(xr)
                mur = calc_mur(qr(i,j,k), rhof(k))
                lbdr = calc_lbdr(mur, dvr)

                wfall_qr = max(0._field_r, (a_tvsb - b_tvsb * (1 + c_tvsb / lbdr)**(-1 * (mur+4))))
                wfall_Nr = max(0._field_r, (a_tvsb - b_tvsb * (1 + c_tvsb / lbdr)**(-1 * (mur+1))))

                sed_qr  = wfall_qr * qr_spl(i,j,k) * rhof(k) ! m/s * kg/m3
                sed_Nr  = wfall_Nr * Nr_spl(i,j,k)

                qr_spl(i,j,k) = qr_spl(i,j,k) - sed_qr * dt_spl / (dzf(k) * rhof(k))
                Nr_spl(i,j,k) = Nr_spl(i,j,k) - sed_Nr * dt_spl / dzf(k)

                if (k .gt. 1) then
                  qr_spl(i,j,k-1) = qr_spl(i,j,k-1) + sed_qr * dt_spl / (dzf(k-1) * rhof(k-1))
                  Nr_spl(i,j,k-1) = Nr_spl(i,j,k-1) + sed_Nr * dt_spl / dzf(k-1)
                end if
                if (jn==1) then
                  precep(i,j,k) = sed_qr / rhof(k)   ! kg kg-1 m s-1
                end if

                if (laerosol) then
                  do s = 1, n_species_active
                    qa_spl(i,j,k,s) = qa_spl(i,j,k,s) - (sed_qr / qr_spl(i,j,k) * qa_spl(i,j,k,s)) * dt_spl / (dzf(k) * rhof(k))
                    if (k > 1) then
                      qa_spl(i,j,k-1,s) = qa_spl(i,j,k-1,s) + (sed_qr / qr_spl(i,j,k) * qa_spl(i,j,k,s)) * dt_spl / (dzf(k-1) * rhof(k-1))
                    end if
                  end do
                  ! Scavenging needs this, annoying to re-compute
                  sed_qr_(i,j,k) = sed_qr
                end if
              end if
            end do
          end do
        end do
      end if ! l_lognormal
    end do ! time splitting loop

    ! the last time splitting step lowered the base level
    ! and we still need to adjust for it
    qrbase = max(1, qrbase - 1)

    Nrp(:,:,qrbase:qrroof) = Nrp(:,:,qrbase:qrroof) + &
      (Nr_spl(:,:,qrbase:qrroof) - Nr(:,:,qrbase:qrroof))/delt

    qrp(:,:,qrbase:qrroof) = qrp(:,:,qrbase:qrroof) + &
      (qr_spl(:,:,qrbase:qrroof) - qr(:,:,qrbase:qrroof))/delt

    deallocate(qr_spl, Nr_spl)

    if (laerosol) then
      do s = 1, n_species_active
      qap_inr(:,:,qrbase:qrroof,s) = qap_inr(:,:,qrbase:qrroof,s) + &
        (qa_spl(:,:,qrbase:qrroof,s) - qa_inr(:,:,qrbase:qrroof,s)) / delt
      end do

      deallocate(qa_spl)
    end if

    call timer_toc('bulkmicro_sb01/sedimentation_rain')

  end subroutine sedimentation_rain

#ifdef DALES_GPU
  !> Calculate the sedimentation term.
  !!
  !! \param qr Rain water mixing ratio.
  !! \param Nr Rain drop number concentration.
  !! \param rhof Density at full levels.
  !! \param dzf Thickness of vertical levels.
  !! \param qrbase Lowest level with rain.
  !! \param qrroof Highest level with rain.
  !! \param qrmask Rain mask.
  !! \param l_mur_cst Switch for selecting constant $\mu$.
  !! \param mur_cst Constant $\mu$ value.
  !! \param delt Time step size.
  !! \param Dvr Rain water mean diameter.
  !! \param lbdr DSD $\lambda$ parameter.
  !! \param mur DSD $\mu$ parameter.
  !! \param xr Mean mass of rain drops.
  !! \param qrp Tendency of rain water mixing ratio.
  !! \param Nrp Tendency of rain drop number concentration.
  !! \param precep Precipitation.
  subroutine sedimentation_rain_gpu(qr, Nr, rhof, dzf, qrbase, qrroof, qrmask, &
                                    l_lognormal, l_mur_cst, mur_cst, delt, Dvr, lbdr, &
                                    mur, xr, qrp, Nrp, precep, laerosol, m_inr, sed_qr_)
    real(field_r), intent(in)    :: qr(2:i1,2:j1,1:k1)
    real(field_r), intent(in)    :: Nr(2:i1,2:j1,1:k1)
    real(field_r), intent(in)    :: rhof(1:k1)
    real(field_r), intent(in)    :: dzf(1:k1)

    integer,       intent(inout) :: qrbase
    integer,       intent(in)    :: qrroof
    logical,       intent(inout) :: qrmask(2:i1,2:j1,1:k1)

    logical,       intent(in)    :: l_lognormal, l_mur_cst
    real(field_r), intent(in)    :: mur_cst
    real(field_r), intent(in)    :: delt

    real(field_r), intent(inout) :: Dvr(2:i1,2:j1,1:k1)
    real(field_r), intent(inout) :: lbdr(2:i1,2:j1,1:k1)
    real(field_r), intent(inout) :: mur(2:i1,2:j1,1:k1)
    real(field_r), intent(inout) :: xr(2:i1,2:j1,1:k1)

    real(field_r), intent(inout) :: qrp(2:i1,2:j1,1:k1)
    real(field_r), intent(inout) :: Nrp(2:i1,2:j1,1:k1)
    real(field_r), intent(out)   :: precep(2:i1,2:j1,1:k1)

    logical,       intent(in),    optional :: laerosol
    type(mode_t),  intent(inout), optional :: m_inr
    real(field_r), intent(out),   optional :: sed_qr_(2:i1,2:j1,1:k1)

    integer       :: i, j, k, jn, sedimbase, s
    integer       :: n_spl      !<  sedimentation time splitting loop
    real(field_r) :: pwcont
    real(field_r) :: delt_inv
    real(field_r) :: Dgr           !<  lognormal geometric diameter
    real(field_r) :: wfall_qr      !<  fall velocity for qr
    real(field_r) :: wfall_Nr      !<  fall velocity for Nr
    real(field_r) :: sed_qr
    real(field_r) :: sed_Nr

    real(field_r), allocatable :: qr_spl(:,:,:), Nr_spl(:,:,:)
    real(field_r), allocatable :: qr_tmp(:,:,:), Nr_tmp(:,:,:)
    real(field_r), allocatable :: qa_spl(:,:,:,:), qa_tmp(:,:,:,:)

    real(field_r), save :: dt_spl

    logical :: laerosol_

    if (present(laerosol)) then
      laerosol_ = laerosol
    else
      laerosol_ = .false.
    end if

    !$acc parallel loop collapse(3) default(present)
    do k = 1, k1
      do j = 2, j1
        do i = 2, i1
          precep(i,j,k) = 0.0
        end do
      end do
    end do

    if (qrbase > qrroof) return

    call timer_tic('bulkmicro_sb/sedimentation_rain', 1)

    allocate(qr_spl(2:i1,2:j1,1:k1))
    allocate(Nr_spl(2:i1,2:j1,1:k1))
    allocate(qr_tmp(2:i1,2:j1,1:k1))
    allocate(Nr_tmp(2:i1,2:j1,1:k1))
    allocate(qa_spl(2:i1,2:j1,1:k1,1:m_inr%nspecies))
    allocate(qa_tmp(2:i1,2:j1,1:k1,1:m_inr%nspecies))

    !$acc enter data create(qr_spl, Nr_spl, qr_tmp, Nr_tmp, qa_spl, qa_tmp)

    n_spl = ceiling(wfallmax * delt / minval(dzf))
    dt_spl = delt / real(n_spl, kind=field_r)

    do jn = 1, n_spl ! time splitting loop
      if (jn == 1) then
        !$acc parallel loop collapse(3) default(present)
        do k = 1, k1
          do j = 2, j1
            do i = 2, i1
              qr_spl(i,j,k) = qr(i,j,k)
              Nr_spl(i,j,k) = Nr(i,j,k)
              qr_tmp(i,j,k) = qr(i,j,k)
              Nr_tmp(i,j,k) = Nr(i,j,k)
              if (laerosol_) then
                do s = 1, m_inr%nspecies
                  qa_spl(i,j,k,s) = qa_inr(i,j,k,s)
                  qa_tmp(i,j,k,s) = qa_inr(i,j,k,s)
                end do
              end if
            end do
          end do
        end do
      else
        !Copy from tmp into spl
        !$acc parallel loop collapse(3) default(present)
        do k = 1, k1
          do j = 2, j1
            do i = 2, i1
              qr_spl(i,j,k) = qr_tmp(i,j,k)
              Nr_spl(i,j,k) = Nr_tmp(i,j,k)
              if (laerosol_) then
                do s = 1, m_inr%nspecies
                  qa_spl(i,j,k,s) = qa_tmp(i,j,k,s)
                end do
              end if

              ! Update mask
              qrmask(i,j,k) = (qr_spl(i,j,k) > qrmin .and. Nr_spl(i,j,k) > 0.0)
            end do
          end do
        end do

        ! lower the rain base by one level to include the rain fall
        ! from the previous step
        qrbase = max(1, qrbase - 1)

        call calculate_rain_parameters(Nr_spl, qr_spl, rhof, l_mur_cst, mur_cst, &
                                      qrbase, qrroof, qrmask, xr, Dvr, mur, lbdr)
      end if

      ! Compute precep
      if (jn == 1) then
        if (l_lognormal) then
          !$acc parallel loop collapse(3) default(present) private(Dgr)
          do k = qrbase, qrroof
            do j = 2, j1
              do i = 2, i1
                if (qrmask(i,j,k)) then
                  Dgr = (exp(4.5_field_r * (log(sig_gr))**2))**(-1.0_field_r/3) * Dvr(i,j,k)
                  sed_qr = sed_flux(Nr_spl(i,j,k), Dgr, log(sig_gr)**2, D_s, 3)
                  pwcont = liq_cont(Nr_spl(i,j,k), Dgr, log(sig_gr)**2, D_s, 3)
                  if (pwcont > eps1) then
                    sed_qr = (qr_spl(i,j,k) * rhof(k) / pwcont) * sed_qr
                  end if
                  precep(i,j,k) = sed_qr / rhof(k)   ! kg kg-1 m s-1
                end if
              end do
            end do
          end do
        else ! l_lognormal
          !$acc parallel loop collapse(3) default(present)
          do k = qrbase, qrroof
            do j = 2, j1
              do i = 2, i1
                if (qrmask(i,j,k)) then
                  wfall_qr = max( &
                    0.0_field_r, &
                    a_tvsb - b_tvsb * (1 + c_tvsb / lbdr(i,j,k))**(-1 * (mur(i,j,k) + 4)) &
                  )
                  sed_qr  = wfall_qr * qr_spl(i,j,k) * rhof(k)
                  precep(i,j,k) = sed_qr / rhof(k)   ! kg kg-1 m s-1
                end if
              end do
            end do
          end do
        end if ! l_lognormal
      end if ! jn == 1

      sedimbase = qrbase

      ! k qrbase if == 1
      if (qrbase == 1) then
        sedimbase = sedimbase + 1
        k = 1
          if (l_lognormal) then
            !$acc parallel loop collapse(2) default(present) private(Dgr)
            do j = 2, j1
              do i = 2, i1
                if (qrmask(i,j,k)) then
                  ! correction for width of DSD
                  Dgr = (exp(4.5_field_r * (log(sig_gr))**2))**(-1.0_field_r/3) * Dvr(i,j,k)
                  sed_qr = sed_flux(Nr_spl(i,j,k), Dgr, log(sig_gr)**2, D_s, 3)
                  sed_Nr = 1.0_field_r / pirhow * sed_flux(Nr_spl(i,j,k), Dgr, log(sig_gr)**2, D_s, 0)

                  ! correction for the fact that pwcont .ne. qr_spl
                  ! actually in this way for every grid box a fall velocity is determined
                  pwcont = liq_cont(Nr_spl(i,j,k), Dgr, log(sig_gr)**2, D_s, 3)       ! note : kg m-3
                  if (pwcont > eps1) then
                    sed_qr = (qr_spl(i,j,k) * rhof(k) / pwcont) * sed_qr
                    ! or:
                    ! qr_spl*(sed_qr/pwcont) = qr_spl*fallvel.
                  end if

                  qr_tmp(i,j,k) = qr_tmp(i,j,k) - sed_qr*dt_spl / (dzf(k) * rhof(k))
                  Nr_tmp(i,j,k) = Nr_tmp(i,j,k) - sed_Nr*dt_spl / dzf(k)
                end if
              end do
            end do
          else ! l_lognormal
            !$acc parallel loop collapse(2) default(present)
            do j = 2, j1
              do i = 2, i1
                if (qrmask(i,j,k)) then
                  wfall_qr = max(0.0_field_r, (a_tvsb - b_tvsb * (1 + c_tvsb / lbdr(i,j,k))**(-1 * (mur(i,j,k) + 4))))
                  wfall_Nr = max(0.0_field_r, (a_tvsb - b_tvsb * (1 + c_tvsb / lbdr(i,j,k))**(-1 * (mur(i,j,k) + 1))))

                  sed_qr  = wfall_qr*qr_spl(i,j,k)*rhof(k)
                  sed_qr_(i,j,k) = sed_qr
                  sed_Nr  = wfall_Nr*Nr_spl(i,j,k)

                  qr_tmp(i,j,k) = qr_tmp(i,j,k) - sed_qr*dt_spl/(dzf(k)*rhof(k))
                  Nr_tmp(i,j,k) = Nr_tmp(i,j,k) - sed_Nr*dt_spl/dzf(k)
                  if (laerosol_) then
                    do s = 1, m_inr%nspecies
                      qa_tmp(i,j,k,s) = qa_tmp(i,j,k,s) - &
                        (sed_qr / qr_tmp(i,j,k) * qa_tmp(i,j,k,s)) * dt_spl &
                        / (dzf(k) * rhof(k))
                    end do
                  end if
                end if
              end do
            end do
          end if ! l_lognormal
      end if ! qrbase == 1

      if (l_lognormal) then
        !$acc parallel loop collapse(3) default(present) private(Dgr)
        do k = sedimbase, qrroof
          do j = 2, j1
            do i = 2, i1
              if (qrmask(i,j,k)) then
                ! correction for width of DSD
                Dgr = (exp(4.5_field_r * (log(sig_gr))**2))**(-1.0_field_r/3) * Dvr(i,j,k)
                sed_qr = sed_flux(Nr_spl(i,j,k),Dgr,log(sig_gr)**2,D_s,3)
                sed_Nr = 1.0_field_r / pirhow * sed_flux(Nr_spl(i,j,k), Dgr, log(sig_gr)**2, D_s, 0)

                ! correction for the fact that pwcont .ne. qr_spl
                ! actually in this way for every grid box a fall velocity is determined
                pwcont = liq_cont(Nr_spl(i,j,k), Dgr, log(sig_gr)**2, D_s, 3)       ! note : kg m-3
                if (pwcont > eps1) then
                  sed_qr = (qr_spl(i,j,k) * rhof(k) / pwcont) * sed_qr
                  ! or:
                  ! qr_spl*(sed_qr/pwcont) = qr_spl*fallvel.
                end if

                !$acc atomic update
                qr_tmp(i,j,k) = qr_tmp(i,j,k) - sed_qr * dt_spl / (dzf(k) * rhof(k))
                !$acc atomic update
                Nr_tmp(i,j,k) = Nr_tmp(i,j,k) - sed_Nr * dt_spl / dzf(k)

                !$acc atomic update
                qr_tmp(i,j,k-1) = qr_tmp(i,j,k-1) + sed_qr*dt_spl / (dzf(k-1) * rhof(k-1))
                !$acc atomic update
                Nr_tmp(i,j,k-1) = Nr_tmp(i,j,k-1) + sed_Nr*dt_spl / dzf(k-1)
              end if
            end do
          end do
        end do
      else
        !$acc parallel loop collapse(3) default(present)
        do k = sedimbase, qrroof
          do j = 2, j1
            do i = 2, i1
              if (qrmask(i,j,k)) then
                wfall_qr = max(0.0_field_r, (a_tvsb - b_tvsb * (1 + c_tvsb / lbdr(i,j,k))**(-1 * (mur(i,j,k) + 4))))
                wfall_Nr = max(0.0_field_r, (a_tvsb - b_tvsb * (1 + c_tvsb / lbdr(i,j,k))**(-1 * (mur(i,j,k) + 1))))

                sed_qr  = wfall_qr * qr_spl(i,j,k) * rhof(k)
                sed_qr_(i,j,k) = sed_qr
                sed_Nr  = wfall_Nr * Nr_spl(i,j,k)

                !$acc atomic update
                qr_tmp(i,j,k) = qr_tmp(i,j,k) - sed_qr * dt_spl / (dzf(k) * rhof(k))
                !$acc atomic update
                Nr_tmp(i,j,k) = Nr_tmp(i,j,k) - sed_Nr * dt_spl / dzf(k)

                !$acc atomic update
                qr_tmp(i,j,k-1) = qr_tmp(i,j,k-1) + sed_qr * dt_spl / (dzf(k-1) * rhof(k-1))
                !$acc atomic update
                Nr_tmp(i,j,k-1) = Nr_tmp(i,j,k-1) + sed_Nr * dt_spl / dzf(k-1)

                if (laerosol_) then
                  do s = 1, m_inr%nspecies
                    !$acc atomic update
                    qa_tmp(i,j,k,s) = qa_tmp(i,j,k,s) - &
                      (sed_qr / qr_tmp(i,j,k) * qa_tmp(i,j,k,s)) * dt_spl &
                      / (dzf(k) * rhof(k))
                    !$acc atomic update
                    qa_tmp(i,j,k-1,s) = qa_tmp(i,j,k-1,s) + &
                      (sed_qr / qr_tmp(i,j,k) * qa_tmp(i,j,k,s)) * dt_spl &
                      / (dzf(k-1) * rhof(k-1))
                  end do
                end if
              end if
            end do
          end do
        end do
      end if ! l_lognormal

    end do ! time splitting loop

    ! the last time splitting step lowered the base level
    ! and we still need to adjust for it
    qrbase = max(1,qrbase-1)

    delt_inv = 1.0 / delt

    !$acc parallel loop collapse(3) default(present)
    do k = qrbase, qrroof
      do j = 2, j1
        do i = 2, i1
          Nrp(i,j,k) = Nrp(i,j,k) + (Nr_tmp(i,j,k) - Nr(i,j,k)) * delt_inv
          qrp(i,j,k) = qrp(i,j,k) + (qr_tmp(i,j,k) - qr(i,j,k)) * delt_inv
          if (laerosol_) then
            do s = 1, m_inr%nspecies
              m_inr%tend(i,j,k,s+1) = m_inr%tend(i,j,k,s+1) + &
                (qa_spl(i,j,k,s) - m_inr%conc(i,j,k,s+1)) * delt_inv
            end do
          end if
        end do
      end do
    end do

    !$acc exit data delete(qr_spl, Nr_spl, qr_tmp, Nr_tmp, qa_spl, qa_tmp)

    deallocate(qr_spl, Nr_spl, qr_tmp, Nr_tmp, qa_spl, qa_tmp)

    call timer_toc('bulkmicro_sb01/sedimentation_rain')

  end subroutine sedimentation_rain_gpu
#endif

  real function sed_flux(Nin, Din, sig2, Ddiv, nnn)
  !*********************************************************************
  ! Function to calculate numerically the analytical solution of the
  ! sedimentation flux between Dmin and Dmax based on
  ! Feingold et al 1986 eq 17 -20.
  ! fall velocity is determined by alfa* D^beta with alfa+ beta taken as
  ! specified in Rogers and Yau 1989 Note here we work in D and in SI
  ! (in Roger+Yau in cm units + radius)
  ! flux is multiplied outside sed_flux with 1/rho_air to get proper
  ! kg/kg m/s units
  !
  ! M.C. van Zanten    August 2005
  !*********************************************************************

    real(field_r), intent(in) :: Nin, Ddiv
    real(field_r), intent(in) :: Din, sig2
    integer,       intent(in) :: nnn
    !para. def. lognormal DSD (sig2 = ln^2 sigma_g), D sep. droplets from drops
    !,power of of D in integral

    real(field_r), parameter ::   &
      C = rhow*pi/6., &
      D_intmin = 1e-6, &
      D_intmax = 4.3e-3

    real(field_r) :: &
      alfa,        & ! constant in fall velocity relation
      beta,        & ! power in fall vel. rel.
      D_min,       & ! min integration limit
      D_max,       & ! max integration limit
      flux           ![kg m^-2 s^-1]

    flux = 0.0_field_r

    if (Din < Ddiv) then
      alfa = 3.e5*100  ![1/ms]
      beta = 2
      D_min = D_intmin
      D_max = Ddiv
      flux = C*Nin*alfa*erfint(beta,Din,D_min,D_max,sig2,nnn)
    else
      ! fall speed ~ D^2
      alfa = 3.e5*100 ![1/m 1/s]
      beta = 2
      D_min = Ddiv
      D_max = 133e-6
      flux = flux + C*Nin*alfa*erfint(beta,Din,D_min,D_max,sig2,nnn)

      ! fall speed ~ D
      alfa = 4e3     ![1/s]
      beta = 1
      D_min = 133e-6
      D_max = 1.25e-3
      flux = flux + C*Nin*alfa*erfint(beta,Din,D_min,D_max,sig2,nnn)

      ! fall speed ~ sqrt(D)
      alfa = 1.4e3 *0.1  ![m^.5 1/s]
      beta = .5
      D_min = 1.25e-3
      D_max = D_intmax
      flux = flux + C*Nin*alfa*erfint(beta,Din,D_min,D_max,sig2,nnn)
    end if
    sed_flux = flux
  end function sed_flux

  real function liq_cont(Nin,Din,sig2,Ddiv,nnn)
  !$acc routine seq
  !*********************************************************************
  ! Function to calculate numerically the analytical solution of the
  ! liq. water content between Dmin and Dmax based on
  ! Feingold et al 1986 eq 17 -20.
  !
  ! M.C. van Zanten    September 2005
  !*********************************************************************
    use modglobal, only : pi,rhow
    implicit none

    real(field_r), intent(in) :: Nin, Ddiv
    real(field_r), intent(in) :: Din, sig2
    integer, intent(in) :: nnn
    !para. def. lognormal DSD (sig2 = ln^2 sigma_g), D sep. droplets from drops
    !,power of of D in integral

    real(field_r), parameter :: beta = 0           &
                      ,C = pi/6.*rhow     &
                      ,D_intmin = 80e-6    &   ! value of start of rain D
                      ,D_intmax = 3e-3         !4.3e-3    !  value is now max value for sqrt fall speed rel.

    real(field_r) ::  D_min        & ! min integration limit
            ,D_max        & ! max integration limit
            ,sn

    sn = sign(0.5_field_r, Din - Ddiv)
    D_min = (0.5 - sn) * D_intmin + (0.5 + sn) * Ddiv
    D_max = (0.5 - sn) * Ddiv     + (0.5 + sn) * D_intmax

    liq_cont = C*Nin*erfint(beta,Din,D_min,D_max,sig2,nnn)
  end function liq_cont

  real function erfint(beta, D, D_min, D_max, sig2,nnn )
  !$acc routine seq
  !*********************************************************************
  ! Function to calculate erf(x) approximated by a polynomial as
  ! specified in 7.1.27 in Abramowitz and Stegun
  ! NB phi(x) = 0.5(erf(0.707107*x)+1) but 1 disappears by substraction
  !
  !*********************************************************************
    implicit none
    real(field_r), intent(in) :: beta, D, D_min, D_max, sig2
    integer, intent(in) :: nnn

    real(field_r), parameter :: eps = 1e-10       
    !                  ,a1 = 0.278393    & !a1 till a4 constants in polynomial fit to the error
    !                  ,a2 = 0.230389    & !function 7.1.27 in Abramowitz and Stegun
    !                  ,a3 = 0.000972    &
    !                  ,a4 = 0.078108
    real(field_r) :: nn, ymin, ymax, erfymin, erfymax, D_inv

    D_inv = 1./(eps + D)
    nn = beta + nnn

    ymin = 0.707107*(log(D_min*D_inv) - nn*sig2)/(sqrt(sig2))
    ymax = 0.707107*(log(D_max*D_inv) - nn*sig2)/(sqrt(sig2))

    !erfymin = 1.-1./((1.+a1*abs(ymin) + a2*abs(ymin)**2 + a3*abs(ymin)**3 +a4*abs(ymin)**4)**4)
    !erfymax = 1.-1./((1.+a1*abs(ymax) + a2*abs(ymax)**2 + a3*abs(ymax)**3 +a4*abs(ymax)**4)**4)
    erfymin = erf(abs(ymin))
    erfymax = erf(abs(ymax))

    erfymin = sign(erfymin, ymin)
    erfymax = sign(erfymax, ymax)

    erfint = max(0., D**nn*exp(0.5*nn**2*sig2)*0.5*(erfymax-erfymin))
  end function erfint

end module bulkmicro_sb
