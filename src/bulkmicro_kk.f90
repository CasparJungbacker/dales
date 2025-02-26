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
!> Kernels for Khairoutdinov-Kogan microphysics.
module bulkmicro_kk
  use modglobal,    only: i1, ih, j1, jh, k1, rlv, cp, pi, rv, nsv
  use modmicrodata, only: pirhow, qrmin, Nc_0, eps0, qcmin, ncmin
  use modprecision, only: field_r
  use modtimer,     only: timer_tic, timer_toc
  use modaerosol,   only: mode_t, idx_tab, iACS, iCOS, rho_a, &
                          aerosol_get_type_in_cloud, aerosol_get_index_in_cloud, &
                          aerosol_get_index_in_mode
  use modmath,      only: inv_sqrt_two

  implicit none

  private

  character(*), parameter :: modname = "bulkmicro_kk"

  real(field_r), parameter :: &
    c_evap = 0.87,  & !< Coefficient for evaporation.
    D0 = 50e-6,     & !< Diameter separating cloud and precipitation parts of the DSD.
    Dv = 2.4e-5,    & !< Diffusivity of water vapor [m2/s].
    Kt = 2.5e-2,    & !< Conductivity of heat [J/(sKm)].
    wfallmax = 9.9, & !< Terminal fall velocity.
    xrmax = 5.2e-7    !< Max mean mass of pw.
  
  public :: do_bulkmicro_kk

contains

  !> Calculate microphysical source terms according to Khairoutdinov and Kogan (2000).
  subroutine do_bulkmicro_kk
    use modmicrodata,     only: qr, Nr, iqr, iNr, thlpmcr, qtpmcr, qcbase, &
                                qcroof, qrbase, qrroof, qcmask, qrmask, qrp, Nrp, &
                                Dvr, xr, delt, precep, Nc, sed_qr, Ncp, qa_c => qa_inc, qa_r => qa_inr, &
                                qap_r => qap_inr, qap_c => qap_inc
    use modfields,        only: rhof, ql0, exnf, qvsl, tmp0, esl, svm, qt0
    use modglobal,        only: dzf
    use modbulkmicrostat, only: bulkmicrotend
    use modaerosol,       only: modes, iACS, iCOS, iINR, iINC, laerosol

    call calculate_rain_parameters(Nr, qr, rhof, qrbase, qrroof, qrmask, Dvr, &
                                   xr)
    call bulkmicrotend
    call autoconversion(ql0, Nc, qa_c, qa_r, rhof, exnf, qcbase, qcroof, laerosol, &
                        delt, thlpmcr, qtpmcr, Ncp, qrp, Nrp, qap_c, qap_r)
    call bulkmicrotend
    call accretion(ql0, qr, Nc, rhof, exnf, qa_c, qa_r, qcbase, qcroof, qrbase, &
                   qrroof, laerosol, thlpmcr, qtpmcr, qrp, Ncp, qap_c, qap_r)
    call bulkmicrotend
    call evaporation(ql0, qt0, qr, Nr, qvsl, esl, tmp0, svm(:,:,:,iqr), &
                     svm(:,:,:,iNr), qa_r, rhof, exnf, qrbase, qrroof, Dvr, &
                     xr, delt, thlpmcr, qtpmcr, qrp, Nrp, qap_r, laerosol, &
                     modes(iACS), modes(iCOS))
    call bulkmicrotend
#ifdef DALES_GPU
    call sedimentation_rain_gpu(qr, Nr, rhof, dzf, qrbase, qrroof, qrmask, delt, &
                                Dvr, xr, qrp, Nrp, precep, laerosol, modes(iINR), &
                                sed_qr)
#else
    call sedimentation_rain(qr, Nr, qa_r, rhof, dzf, qrroof, delt, laerosol, qrmask, &
                            Dvr, xr, qrp, Nrp, qap_r, qrbase, precep)
#endif
    call bulkmicrotend

  end subroutine do_bulkmicro_kk

  !> Calculate rain DSD integral properties and parameters.
  !!
  !! \param nr rain drop number concentration.
  !! \param qr rain water mixing ratio.
  !! \param rhof Density at full levels.
  !! \param qrbase Lowest level with rain.
  !! \param qrroof Highest level with rain.
  !! \param qrmask Rain mask.
  !! \param Dvr Rain water mean diameter.
  !! \param xr Mean mass of rain drops.
  subroutine calculate_rain_parameters(Nr, qr, rhof, qrbase, qrroof, qrmask, &
                                       Dvr, xr)
    real(field_r), intent(in)  :: Nr(2:i1,2:j1,1:k1)
    real(field_r), intent(in)  :: qr(2:i1,2:j1,1:k1)
    real(field_r), intent(in)  :: rhof(1:k1)

    integer,       intent(in)  :: qrbase, qrroof
    logical,       intent(in)  :: qrmask(2:i1,2:j1,1:k1)
    
    real(field_r), intent(out) :: xr(2:i1,2:j1,1:k1)
    real(field_r), intent(out) :: Dvr(2:i1,2:j1,1:k1)

    integer :: i, j, k

    if (qrbase > qrroof) return

    call timer_tic('bulkmicro_kk/calculate_rain_parameters', 1)

    !$acc parallel loop collapse(3) default(present)
    do k = qrbase, qrroof
      do j = 2, j1
        do i = 2, i1
          if (qrmask(i,j,k)) then
            xr(i,j,k) = rhof(k) * qr(i,j,k) / Nr(i,j,k)

            ! to ensure x_pw is within bounds
            xr(i,j,k) = min(xr(i,j,k), xrmax)
            Dvr(i,j,k) = (xr(i,j,k) / pirhow)**(1.0_field_r/3)
          endif
        enddo
      enddo
    enddo

    call timer_toc('bulkmicro_kk/calculate_rain_parameters')

  end subroutine calculate_rain_parameters

  !> Calculate the autoconversion term.
  !!
  !! \param ql0 Liquid water mixing ratio.
  !! \param rhof Density at full levels.
  !! \param exnf Exner function at full levels.
  !! \param qcbase Lowest level with cloud.
  !! \param qcroof Highest level with cloud.
  !! \param qcmask Cloud mask.
  !! \param thlpmcr Tendency of $\theta_l$.
  !! \param qtpmcr Tendency of $\q_t$.
  !! \param qrp Tendency of rain water mixing ratio.
  !! \param Nrp Tendency of rain drop number concentration.
  subroutine autoconversion(ql, nc, qa_c, qa_r, rhof, exnf, qcbase, qcroof, &
                            laerosol, delt, thlpmcr, qtpmcr, ncp, qrp, nrp, &
                            qap_c, qap_r)

    real(field_r), intent(in)    :: ql(2:,2:,:)
    real(field_r), intent(in)    :: nc(2:,2:,:)
    real(field_r), intent(in)    :: qa_c(2:,2:,:,:)
    real(field_r), intent(in)    :: qa_r(2:,2:,:,:)
    real(field_r), intent(in)    :: rhof(:)
    real(field_r), intent(in)    :: exnf(:)
    integer,       intent(in)    :: qcbase, qcroof
    logical,       intent(in)    :: laerosol
    real(field_r), intent(in)    :: delt

    real(field_r), intent(inout) :: thlpmcr(2:,2:,:)
    real(field_r), intent(inout) :: qtpmcr(2-ih:,2-jh:,1:)
    real(field_r), intent(inout) :: ncp(2:,2:,:)
    real(field_r), intent(inout) :: qrp(2:,2:,:)
    real(field_r), intent(inout) :: nrp(2:,2:,:)
    real(field_r), intent(inout) :: qap_c(2:,2:,:,:)
    real(field_r), intent(inout) :: qap_r(2:,2:,:,:)

    character(len=*), parameter :: routine = modname//"autoconversion"

    integer       :: i, j, k, s, naer
    real(field_r) :: q_c, n_c
    real(field_r) :: au, xc

    if (qcbase > qcroof) return

    call timer_tic(routine, 1)

    naer = size(qa_c, dim=4)

    !$acc parallel loop collapse(3) default(present) private(au, xc)
    do k = qcbase, qcroof
      do j = 2, j1
        do i = 2, i1
          q_c = ql(i,j,k) 
          n_c = nc(i,j,k)
          if (q_c > qcmin .and. n_c > ncmin) then
            au = 1350 * q_c**(2.47_field_r) * (n_c / 1E6)**(-1.79_field_r)
            au = min(q_c / delt, au)
            qrp(i,j,k) = qrp(i,j,k) + au
            qtpmcr(i,j,k) = qtpmcr(i,j,k) - au
            thlpmcr(i,j,k) = thlpmcr(i,j,k) + (rlv / (cp * exnf(k))) * au
            Nrp(i,j,k) = Nrp(i,j,k) + au * rhof(k) / (pirhow * D0**3)

            if (laerosol) then
              xc = rhof(k) * q_c / (n_c + eps0)
              ncp(i,j,k) = ncp(i,j,k) - au / xc * rhof(k)
              do s = 1, naer
                qap_c(i,j,k,s) = qap_c(i,j,k,s) - au / q_c * qa_c(i,j,k,s)
                qap_r(i,j,k,s) = qap_r(i,j,k,s) + au / q_c * qa_r(i,j,k,s)
              end do
            end if
          end if
        end do
      end do
    end do

    call timer_toc(routine)

  end subroutine autoconversion

  !> Calculate the accretion term.
  !!
  !! \param ql0 Liquid water mixing ratio.
  !! \param qr Rain water mixing ratio.
  !! \param exnf Exner function at full levels.
  !! \param qcbase Lowest level with cloud.
  !! \param qcroof Highest level with cloud.
  !! \param qcmask Cloud mask.
  !! \param qrbase Lowest level with rain.
  !! \param qrroof Highest level with rain.
  !! \param qrmask Rain mask.
  !! \param thlpmcr Tendency of $\theta_l$.
  !! \param qtpmcr Tendency of total water mixing ratio.
  !! \param qrp Tendency of rain water mixing ratio.
  subroutine accretion(ql, qr, nc, rhof, exnf, qa_c, qa_r, qcbase, qcroof, qrbase, &
                       qrroof, laerosol, thlpmcr, qtpmcr, qrp, ncp, qap_c, qap_r)

    real(field_r), intent(in)    :: ql(2-ih:,2-jh:,:)
    real(field_r), intent(in)    :: qr(2:,2:,:)
    real(field_r), intent(in)    :: nc(2:,2:,:)
    real(field_r), intent(in)    :: rhof(:)
    real(field_r), intent(in)    :: exnf(:)
    real(field_r), intent(in)    :: qa_c(2:,2:,:,:)
    real(field_r), intent(in)    :: qa_r(2:,2:,:,:)
    integer,       intent(in)    :: qcbase, qcroof, qrbase, qrroof
    logical,       intent(in)    :: laerosol

    real(field_r), intent(inout) :: thlpmcr(2:,2:,:)
    real(field_r), intent(inout) :: qtpmcr(2-ih:,2-jh:,:)
    real(field_r), intent(inout) :: qrp(2:,2:,:)
    real(field_r), intent(inout) :: ncp(2:,2:,:)
    real(field_r), intent(inout) :: qap_c(2:,2:,:,:)
    real(field_r), intent(inout) :: qap_r(2:,2:,:,:)

    character(len=*), parameter :: routine = modname//"accretion"
    
    integer       :: i, j, k, s, naer
    real(field_r) :: ac, xc
    real(field_r) :: q_c, q_r

    if (max(qrbase, qcbase) > min(qcroof, qcroof)) return

    call timer_tic(routine, 1)

    naer = size(qa_c, dim=4)

    !$acc parallel loop collapse(3) default(present) private(ac, xc, q_c, q_r)
    do k = max(qrbase, qcbase), min(qcroof, qrroof)
      do j = 2, j1
        do i = 2, i1
          q_c = ql(i,j,k)
          q_r = qr(i,j,k)
          ! CJ: removed the check for Nr to prevent passing yet another array
          ! -> is this bad?
          if (q_c > qcmin .and. q_r > qrmin) then
            ac = 67 * (q_c * q_r)**1.15_field_r
            qrp(i,j,k) = qrp(i,j,k) + ac
            qtpmcr(i,j,k) = qtpmcr(i,j,k) - ac
            thlpmcr(i,j,k) = thlpmcr(i,j,k) + (rlv / (cp * exnf(k))) * ac

            if (laerosol) then
              xc = rhof(k) * q_c / nc(i,j,k)
              ncp(i,j,k) = ncp(i,j,k) - ac / xc

              do s = 1, naer
                qap_c(i,j,k,s) = qap_c(i,j,k,s) - ac / q_c * qa_c(i,j,k,s)
                qap_r(i,j,k,s) = qap_r(i,j,k,s) + ac / q_c * qa_r(i,j,k,s)
              end do
            end if
          end if
        end do
      end do
    end do

    call timer_toc(routine)

  end subroutine accretion

  !> Calculate the evaporation term.
  !!
  !! \param ql0 Liquid water mixing ratio.
  !! \param qt0 Total water mixing ratio.
  !! \param qvsl 
  !! \param esl 
  !! \param tmp0 Temperature.
  !! \param qrm Rain water mixing ratio at previous time step.
  !! \param Nrm Rain drop number concentration at previous time step.
  !! \param Nr Rain drop number concentration.
  !! \param rhof Density at full levels.
  !! \param exnf Exner function at full levels.
  !! \param qrbase Lowest level with rain.
  !! \param qrroof Highest level with rain.
  !! \param qrmask Rain mask.
  !! \param Dvr Rain water mean diameter.
  !! \param xr Mean mass of rain drops.
  !! \param delt Time step size.
  !! \param thlpmcr Tendency of $\theta_l$.
  !! \param qtpmcr Tendency of total water mixing ratio.
  !! \param qrp Tendency of rain water mixing ratio.
  !! \param Nrp Tendency of rain drop number concentration.
  subroutine evaporation(ql, qt, qr, Nr, qvsl, esl, T, qrm, Nrm, qa_r, rhof, exnf, &
                         qrbase, qrroof, Dvr, xr, delt, thlpmcr, &
                         qtpmcr, qrp, Nrp, qap_r, laerosol, m_acs, m_cos)

    real(field_r), intent(in)    :: ql(2-ih:,2-jh:,:)
    real(field_r), intent(in)    :: qt(2-ih:,2-jh:,:)
    real(field_r), intent(in)    :: qr(2:,2:,:)
    real(field_r), intent(in)    :: Nr(2:,2:,:)
    real(field_r), intent(in)    :: qvsl(2-ih:,2-jh:,:)
    real(field_r), intent(in)    :: esl(2-ih:,2-jh:,:)
    real(field_r), intent(in)    :: T(2-ih:,2-jh:,:)
    real(field_r), intent(in)    :: qrm(2-ih:,2-jh:,:)
    real(field_r), intent(in)    :: Nrm(2-ih:,2-jh:,:)
    real(field_r), intent(in)    :: qa_r(2:,2:,:,:)
    real(field_r), intent(in)    :: rhof(:)
    real(field_r), intent(in)    :: exnf(:)
    integer,       intent(in)    :: qrbase, qrroof
    real(field_r), intent(in)    :: Dvr(2:,2:,:)
    real(field_r), intent(in)    :: xr(2:,2:,:)
    logical,       intent(in)    :: laerosol
    real(field_r), intent(in)    :: delt

    real(field_r), intent(inout) :: thlpmcr(2:,2:,:)
    real(field_r), intent(inout) :: qtpmcr(2-ih:,2-jh:,:)
    real(field_r), intent(inout) :: qrp(2:,2:,:)
    real(field_r), intent(inout) :: Nrp(2:,2:,:)
    real(field_r), intent(inout) :: qap_r(2:,2:,:,:)
    type(mode_t),  intent(inout) :: m_acs
    type(mode_t),  intent(inout) :: m_cos

    character(*), parameter :: routine = modname//"evaporation"

    integer       :: i, j, k
    real(field_r) :: S, G
    real(field_r) :: evap, Nevap
    real(field_r) :: E, V, evapt, rho, eps, f_evp
    real(field_r) :: m_evp, v_evp, rho_evp
    integer       :: idx, l, src_idx, target_idx, itype, naer
    real(field_r) :: dm, dm_fac, fm, fn, dn
    real(field_r), parameter :: Dc = 1.0

    if (qrbase > qrroof) return

    call timer_tic(routine, 1)

    if (laerosol) naer = size(qa_r, dim=4)

    !$acc parallel loop collapse(3) default(present) private(S, G, evap, Nevap)
    do k = qrbase, qrroof
      do j = 2, j1
        do i = 2, i1
          if (qr(i,j,k) > qrmin) then
            ! Consider computing qvsl, esl and T in-place
            S = min(0.0_field_r, (qt(i,j,k) - ql(i,j,k)) / qvsl(i,j,k) - 1)
            G = (Rv * T(i,j,k)) / (Dv * esl(i,j,k)) + rlv / &
                (Kt * T(i,j,k)) * (rlv / (Rv * T(i,j,k)) - 1)
            G = 1 / G

            evap = c_evap * 2 * pi * Dvr(i,j,k) * G * S * Nr(i,j,k) / rhof(k)
            Nevap = evap * rhof(k) / xr(i,j,k)

            ! Limit evaporation to available rain water
            if (evap < - qrm(i,j,k) / delt) then
              Nevap = - Nrm(i,j,k) / delt
              evap  = - qrm(i,j,k) / delt
            endif

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

    call timer_toc(routine)

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
  !! \param delt Time step size.
  !! \param Dvr Rain water mean diameter.
  !! \param xr Mean mass of rain drops.
  !! \param qrp Tendency of rain water mixing ratio.
  !! \param Nrp Tendency of rain drop number concentration.
  !! \param precep Precipitation.
  subroutine sedimentation_rain(qr, Nr, qa_inr, rhof, dzf, qrroof, delt, laerosol, &
                                qrmask, Dvr, xr, qrp, Nrp, qap_inr, qrbase, precep)

    real(field_r), intent(in)    :: qr(2:,2:,:)
    real(field_r), intent(in)    :: Nr(2:,2:,:)
    real(field_r), intent(in)    :: qa_inr(2:,2:,:,:)
    real(field_r), intent(in)    :: rhof(:)
    real(field_r), intent(in)    :: dzf(:)
    integer,       intent(in)    :: qrroof
    real(field_r), intent(in)    :: delt
    logical,       intent(in)    :: laerosol 

    logical,       intent(inout) :: qrmask(2:i1,2:j1,1:k1)
    real(field_r), intent(inout) :: Dvr(2:i1,2:j1,1:k1)
    real(field_r), intent(inout) :: xr(2:i1,2:j1,1:k1)
    real(field_r), intent(inout) :: qrp(2:i1,2:j1,1:k1)
    real(field_r), intent(inout) :: Nrp(2:i1,2:j1,1:k1)
    real(field_r), intent(inout) :: qap_inr(2:,2:,:,:)
    integer,       intent(inout) :: qrbase

    real(field_r), intent(out)   :: precep(2:,2:,:)

    character(*), parameter :: routine = modname//"sedimentation_rain"

    integer       :: i, j, k, jn, s, naer
    integer       :: n_spl      !<  sedimentation time splitting loop
    real(field_r) :: dt_spl
    real(field_r) :: sed_qr
    real(field_r) :: sed_Nr

    real(field_r), allocatable :: qr_spl(:,:,:), Nr_spl(:,:,:), qa_spl(:,:,:,:)

    precep(:,:,:) = 0 ! zero the precipitation flux field
                      ! the update below is not always performed

    if (qrbase > qrroof) return

    call timer_tic(modname, 1)

    if (laerosol) naer = size(qa_inr, dim=4)

    allocate(qr_spl(2:i1,2:j1,1:k1))
    allocate(Nr_spl(2:i1,2:j1,1:k1))
    if (laerosol) allocate(qa_spl(2:i1,2:j1,1:k1,naer))


    n_spl = ceiling(wfallmax * delt / minval(dzf))
    dt_spl = delt / real(n_spl, kind=field_r)

    do jn = 1, n_spl ! time splitting loop
      if (jn == 1) then
        qr_spl(:,:,:) = qr(:,:,:)
        Nr_spl(:,:,:) = Nr(:,:,:)
        if (laerosol) qa_spl(:,:,:,:) = qa_inr(:,:,:,:)
      else
        ! update parameters after the first iteration
        ! a new mask
        qrmask(:,:,:) = (qr_spl(:,:,:) > qrmin) .and. (Nr_spl(:,:,:) > 0)

        ! lower the rain base by one level to include the rain fall
        ! from the previous step
        qrbase = max(1, qrbase - 1)

        call calculate_rain_parameters(Nr_spl, qr_spl, rhof, qrbase, qrroof, qrmask, Dvr, xr)
      end if

      do k = qrbase, qrroof
        do j = 2, j1
          do i = 2, i1
            if (qr_spl(i,j,k) > qrmin .and. Nr_spl(i,j,k) > 0) then
              sed_qr = max(0.0_field_r, 0.006_field_r * 1E6 * Dvr(i,j,k) - 0.2_field_r) * qr_spl(i,j,k) * rhof(k)
              sed_Nr = max(0.0_field_r, 0.0035_field_r * 1E6 * Dvr(i,j,k) - 0.1_field_r) * Nr_spl(i,j,k)

              qr_spl(i,j,k) = qr_spl(i,j,k) - sed_qr * dt_spl / (dzf(k) * rhof(k))
              Nr_spl(i,j,k) = Nr_spl(i,j,k) - sed_Nr * dt_spl / dzf(k)

              if (k > 1) then
                qr_spl(i,j,k-1) = qr_spl(i,j,k-1) + sed_qr * dt_spl / (dzf(k-1) * rhof(k-1))
                Nr_spl(i,j,k-1) = Nr_spl(i,j,k-1) + sed_Nr * dt_spl / dzf(k-1)
              endif
              if (jn==1) then
                precep(i,j,k) = sed_qr / rhof(k)   ! kg kg-1 m s-1
              endif
              ! This can probably be moved outside the time-splitting loop.
              ! I.e., just take the total qr tendency and use that to diagnose the
              ! aerosol tendency
              if (laerosol) then
                do s = 1, naer
                  qa_spl(i,j,k,s) = qa_spl(i,j,k,s) - (sed_qr / qr_spl(i,j,k) * qa_spl(i,j,k,s)) * dt_spl / (dzf(k) * rhof(k))
                  if (k > 1) then
                    qa_spl(i,j,k-1,s) = qa_spl(i,j,k-1,s) + (sed_qr / qr_spl(i,j,k) * qa_spl(i,j,k,s)) * dt_spl / (dzf(k-1) * rhof(k-1))
                  end if
                end do
              end if
            end if
          end do
        end do
      end do    
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
      do s = 1, naer
        qap_inr(:,:,qrbase:qrroof,s) = qap_inr(:,:,qrbase:qrroof,s) + &
          (qa_spl(:,:,qrbase:qrroof,s) - qa_inr(:,:,qrbase:qrroof,s)) / delt
      end do

      deallocate(qa_spl)
    end if

    call timer_toc(routine)

  end subroutine sedimentation_rain

#ifdef DALES_GPU
  !> Calculate the sedimentation term. Optimized for GPU's.
  !!
  !! \param qr Rain water mixing ratio.
  !! \param Nr Rain drop number concentration.
  !! \param rhof Density at full levels.
  !! \param dzf Thickness of vertical levels.
  !! \param qrbase Lowest level with rain.
  !! \param qrroof Highest level with rain.
  !! \param qrmask Rain mask.
  !! \param delt Time step size.
  !! \param Dvr Rain water mean diameter.
  !! \param xr Mean mass of rain drops.
  !! \param qrp Tendency of rain water mixing ratio.
  !! \param Nrp Tendency of rain drop number concentration.
  !! \param precep Precipitation.
  subroutine sedimentation_rain_gpu(qr, Nr, rhof, dzf, qrbase, qrroof, qrmask, &
                                    delt, Dvr, xr, qrp, Nrp, precep, laerosol, &
                                    m_inr, sed_qr_)
    real(field_r), intent(in)    :: qr(2:i1,2:j1,1:k1)
    real(field_r), intent(in)    :: Nr(2:i1,2:j1,1:k1)
    real(field_r), intent(in)    :: rhof(1:k1)
    real(field_r), intent(in)    :: dzf(1:k1)

    integer,       intent(inout) :: qrbase
    integer,       intent(in)    :: qrroof
    logical,       intent(inout) :: qrmask(2:i1,2:j1,1:k1)

    real(field_r), intent(in)    :: delt

    real(field_r), intent(inout) :: Dvr(2:i1,2:j1,1:k1)
    real(field_r), intent(inout) :: xr(2:i1,2:j1,1:k1)

    real(field_r), intent(inout) :: qrp(2:i1,2:j1,1:k1)
    real(field_r), intent(inout) :: Nrp(2:i1,2:j1,1:k1)
    real(field_r), intent(out)   :: precep(2:i1,2:j1,1:k1)

    logical,       intent(in)    :: laerosol
    type(mode_t),  intent(inout) :: m_inr
    real(field_r), intent(out)   :: sed_qr_(2:i1,2:j1,1:k1)

    character(*), parameter :: routine = modname//"sedimentation_rain"

    integer       :: i, j, k, jn, sedimbase, s
    integer       :: n_spl      !<  sedimentation time splitting loop
    real(field_r) :: sed_qr
    real(field_r) :: sed_Nr
    real(field_r) :: dt_spl
    real(field_r) :: delt_inv

    real(field_r), allocatable :: qr_spl(:,:,:), Nr_spl(:,:,:)
    real(field_r), allocatable :: qr_tmp(:,:,:), Nr_tmp(:,:,:)
    real(field_r), allocatable :: qa_spl(:,:,:,:), qa_tmp(:,:,:,:)

    !$acc parallel loop collapse(3) default(present)
    do k = 1, k1
      do j = 2, j1
        do i = 2, i1
          precep(i,j,k) = 0.0
        end do
      end do
    end do

    if (qrbase > qrroof) return

    call timer_tic(routine, 1)

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
              if (laerosol) then
                do s = 1, m_inr%nspecies
                  qa_spl(i,j,k,s) = m_inr%conc(i,j,k,s+1)
                  qa_tmp(i,j,k,s) = m_inr%conc(i,j,k,s+1)
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
              if (laerosol) then
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

        call calculate_rain_parameters(Nr_spl, qr_spl, rhof, qrbase, qrroof, qrmask, Dvr, xr)
      end if

      ! Compute precep
      if (jn == 1) then
        !$acc parallel loop collapse(3) default(present)
        do k = qrbase, qrroof
          do j = 2, j1
            do i = 2, i1
              if (qrmask(i,j,k)) then
                precep(i,j,k) = max(0.0_field_r, 0.006_field_r * 1E6 * Dvr(i,j,k) - 0.2_field_r) * qr_spl(i,j,k)
              endif
            enddo
          enddo
        enddo
      end if ! jn == 1

      sedimbase = qrbase

      ! k qrbase if == 1
      if (qrbase == 1) then
        sedimbase = sedimbase + 1
        k = 1

        !$acc parallel loop collapse(2) default(present) private(sed_qr, sed_Nr)
        do j = 2, j1
          do i = 2, i1
            if (qrmask(i,j,k)) then
              sed_qr = max(0.0_field_r, 0.006_field_r *1E6 * Dvr(i,j,k) - 0.2_field_r) * qr_spl(i,j,k) * rhof(k)
              sed_qr_(i,j,k) = sed_qr
              sed_Nr = max(0.0_field_r, 0.0035_field_r *1E6 * Dvr(i,j,k) - 0.1_field_r) * Nr_spl(i,j,k)

              qr_tmp(i,j,k) = qr_tmp(i,j,k) - sed_qr * dt_spl / (dzf(k) * rhof(k))
              Nr_tmp(i,j,k) = Nr_tmp(i,j,k) - sed_Nr * dt_spl / dzf(k)
              if (laerosol) then
                do s = 1, m_inr%nspecies
                  qa_tmp(i,j,k,s) = qa_tmp(i,j,k,s) - &
                    (sed_qr / qr_tmp(i,j,k) * qa_tmp(i,j,k,s)) * dt_spl &
                    / (dzf(k) * rhof(k))
                end do
              end if
            endif
          enddo
        enddo
      end if ! qrbase == 1

      !$acc parallel loop collapse(3) default(present) private(sed_qr, sed_Nr)
      do k = sedimbase, qrroof
        do j = 2, j1
          do i = 2, i1
            if (qrmask(i,j,k)) then
              sed_qr = max(0.0_field_r, 0.006_field_r *1E6 * Dvr(i,j,k) - 0.2_field_r) * qr_spl(i,j,k) * rhof(k)
              sed_qr_(i,j,k) = sed_qr
              sed_Nr = max(0.0_field_r, 0.0035_field_r *1E6 * Dvr(i,j,k) - 0.1_field_r) * Nr_spl(i,j,k)

              !$acc atomic update
              qr_tmp(i,j,k) = qr_tmp(i,j,k) - sed_qr*dt_spl/(dzf(k)*rhof(k))
              !$acc atomic update
              Nr_tmp(i,j,k) = Nr_tmp(i,j,k) - sed_Nr*dt_spl/dzf(k)

              !$acc atomic update
              qr_tmp(i,j,k-1) = qr_tmp(i,j,k-1) + sed_qr*dt_spl/(dzf(k-1)*rhof(k-1))
              !$acc atomic update
              Nr_tmp(i,j,k-1) = Nr_tmp(i,j,k-1) + sed_Nr*dt_spl/dzf(k-1)
              if (laerosol) then
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
            endif
          enddo
        enddo
      enddo
    end do ! time splitting loop

    ! the last time splitting step lowered the base level
    ! and we still need to adjust for it
    qrbase = max(1, qrbase - 1)

    delt_inv = 1 / delt

    !$acc parallel loop collapse(3) default(present)
    do k = qrbase, qrroof
      do j = 2, j1
        do i = 2, i1
          Nrp(i,j,k) = Nrp(i,j,k) + (Nr_tmp(i,j,k) - Nr(i,j,k)) * delt_inv
          qrp(i,j,k) = qrp(i,j,k) + (qr_tmp(i,j,k) - qr(i,j,k)) * delt_inv
          if (laerosol) then
            do s = 1, m_inr%nspecies
              m_inr%tend(i,j,k,s+1) = m_inr%tend(i,j,k,s+1) + &
                (qa_spl(i,j,k,s) - m_inr%conc(i,j,k,s+1)) * delt_inv
            end do
          end if
        end do
      end do
    end do

    !$acc exit data delete(qr_spl, Nr_spl, qr_tmp, Nr_tmp)

    deallocate(qr_spl, Nr_spl, qr_tmp, Nr_tmp)

    call timer_toc(routine)

  end subroutine sedimentation_rain_gpu
#endif

end module bulkmicro_kk
