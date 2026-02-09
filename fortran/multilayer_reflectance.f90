! multilayer_reflectance.f90
! Fortran translation of multilayer reflectance calculations

module multilayer_reflectance
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    complex(dp), parameter :: zi = (0.0_dp, 1.0_dp)
    
contains

    ! TE Fresnel reflection coefficient
    function r_te(qi, qj) result(r)
        complex(dp), intent(in) :: qi, qj
        complex(dp) :: r
        r = (qi - qj) / (qi + qj)
    end function r_te

    ! TM Fresnel reflection coefficient
    function r_tm(qi, qj, ei, ej) result(r)
        complex(dp), intent(in) :: qi, qj
        real(dp), intent(in) :: ei, ej
        complex(dp) :: r
        r = (qi/ei - qj/ej) / (qi/ei + qj/ej)
    end function r_tm

    ! RF: reflection looking downward from target_layer
    function RF_multilayer(n_list, d_list, N, polarization, k0, kp, target_layer) result(r_eff)
        integer, intent(in) :: N, target_layer
        real(dp), intent(in) :: n_list(N), d_list(N-1), k0, kp
        character(len=*), intent(in) :: polarization
        complex(dp) :: r_eff
        
        real(dp) :: e_list(N), eps_up, eps_down, kp2
        complex(dp) :: q_up, q_down, xL, wn, wm, wp, wq, al, bl
        integer :: L
        
        kp2 = kp**2
        e_list = n_list**2
        r_eff = cmplx(0.0_dp, 0.0_dp, kind=dp)
        
        ! From bottom to target_layer
        do L = N, target_layer+1, -1
            eps_up = e_list(L-1)
            eps_down = e_list(L)
            
            q_up = sqrt(cmplx(kp2 - eps_up * k0**2, 0.0_dp, kind=dp)) / k0
            q_down = sqrt(cmplx(kp2 - eps_down * k0**2, 0.0_dp, kind=dp)) / k0
            xL = exp(-q_up * k0 * d_list(L-1))
            
            if (polarization == "TE" .or. polarization == "te") then
                wn = q_down / q_up
            else  ! TM
                wn = (q_down / q_up) * (eps_up / eps_down)
            end if
            
            wm = cmplx(1.0_dp, 0.0_dp, kind=dp)
            wp = wm + wn
            wq = wm - wn
            
            al = (wp + r_eff * wq) / 2.0_dp
            bl = (wq + r_eff * wp) / 2.0_dp
            r_eff = xL * (bl / al) * xL
        end do
    end function RF_multilayer

    ! RB: reflection looking upward from target_layer
    function RB_multilayer(n_list, d_list, N, polarization, k0, kp, target_layer) result(r_eff)
        integer, intent(in) :: N, target_layer
        real(dp), intent(in) :: n_list(N), d_list(N-1), k0, kp
        character(len=*), intent(in) :: polarization
        complex(dp) :: r_eff
        
        real(dp) :: e_list(N), eps_up, eps_down, kp2
        complex(dp) :: q_up, q_down, xL, wn, wm, wp, wq, al, bl
        integer :: L
        
        kp2 = kp**2
        e_list = n_list**2
        r_eff = cmplx(0.0_dp, 0.0_dp, kind=dp)
        
        ! From top to target_layer
        do L = 1, target_layer
            eps_up = e_list(L)
            eps_down = e_list(L+1)
            
            q_up = sqrt(cmplx(kp2 - eps_up * k0**2, 0.0_dp, kind=dp)) / k0
            q_down = sqrt(cmplx(kp2 - eps_down * k0**2, 0.0_dp, kind=dp)) / k0
            xL = exp(-q_up * k0 * d_list(L))
            
            if (polarization == "TE" .or. polarization == "te") then
                wn = q_up / q_down
            else  ! TM
                wn = (q_up / q_down) * (eps_down / eps_up)
            end if
            
            wm = cmplx(1.0_dp, 0.0_dp, kind=dp)
            wp = wm + wn
            wq = wm - wn
            
            bl = wq + wp * xL * r_eff * xL
            al = wp + wq * xL * r_eff * xL
            r_eff = bl / al
        end do
    end function RB_multilayer

end module multilayer_reflectance
