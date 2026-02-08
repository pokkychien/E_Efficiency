! te_greens.f90
! TE Green's function engine for multilayer structures

module te_greens
    use multilayer_reflectance
    implicit none
    
contains

    ! Level 1: compute q list for all layers
    subroutine compute_q_list(n_list, N, k0, kp, q_list)
        integer, intent(in) :: N
        real(dp), intent(in) :: n_list(N), k0, kp
        complex(dp), intent(out) :: q_list(N)
        real(dp) :: eps_list(N)
        integer :: i
        
        eps_list = n_list**2
        do i = 1, N
            q_list(i) = sqrt(cmplx(kp**2 - eps_list(i) * k0**2, 0.0_dp, kind=dp)) / k0
        end do
    end subroutine compute_q_list

    ! Compute RF for all layers
    subroutine compute_RF_all(n_list, d_list, N, k0, kp, polarization, RF_all)
        integer, intent(in) :: N
        real(dp), intent(in) :: n_list(N), d_list(N-1), k0, kp
        character(len=*), intent(in) :: polarization
        complex(dp), intent(out) :: RF_all(N)
        
        real(dp) :: eps(N), kp2
        complex(dp) :: q_up, q_down, wn, wp, wq, x_val, r_eff, num, den
        integer :: L
        
        eps = n_list**2
        kp2 = kp**2
        r_eff = cmplx(0.0_dp, 0.0_dp, kind=dp)
        
        ! Recursion from bottom to top
        do L = N-1, 1, -1
            q_up = sqrt(cmplx(kp2 - eps(L) * k0**2, 0.0_dp, kind=dp)) / k0
            q_down = sqrt(cmplx(kp2 - eps(L+1) * k0**2, 0.0_dp, kind=dp)) / k0
            x_val = exp(-q_up * k0 * d_list(L))
            
            if (polarization == "TE" .or. polarization == "te") then
                wn = q_down / q_up
            else  ! TM
                wn = (q_down / q_up) * (eps(L) / eps(L+1))
            end if
            
            wp = 1.0_dp + wn
            wq = 1.0_dp - wn
            
            num = wq + wp * r_eff
            den = wp + wq * r_eff
            r_eff = x_val * (num / den) * x_val
            RF_all(L) = r_eff
        end do
        
        RF_all(N) = cmplx(0.0_dp, 0.0_dp, kind=dp)
    end subroutine compute_RF_all

    ! Compute RB for all layers
    subroutine compute_RB_all(n_list, d_list, N, k0, kp, polarization, RB_all)
        integer, intent(in) :: N
        real(dp), intent(in) :: n_list(N), d_list(N-1), k0, kp
        character(len=*), intent(in) :: polarization
        complex(dp), intent(out) :: RB_all(N)
        
        real(dp) :: eps(N), kp2
        complex(dp) :: q_up, q_down, wn, wp, wq, x_val, r_eff, xr, num, den
        integer :: L
        
        eps = n_list**2
        kp2 = kp**2
        r_eff = cmplx(0.0_dp, 0.0_dp, kind=dp)
        
        RB_all(1) = cmplx(0.0_dp, 0.0_dp, kind=dp)
        
        ! Recursion from top to bottom
        do L = 1, N-1
            q_up = sqrt(cmplx(kp2 - eps(L) * k0**2, 0.0_dp, kind=dp)) / k0
            q_down = sqrt(cmplx(kp2 - eps(L+1) * k0**2, 0.0_dp, kind=dp)) / k0
            x_val = exp(-q_up * k0 * d_list(L))
            
            if (polarization == "TE" .or. polarization == "te") then
                wn = q_up / q_down
            else  ! TM
                wn = (q_up / q_down) * (eps(L+1) / eps(L))
            end if
            
            wp = 1.0_dp + wn
            wq = 1.0_dp - wn
            
            xr = x_val * r_eff * x_val
            num = wq + wp * xr
            den = wp + wq * xr
            r_eff = num / den
            RB_all(L+1) = r_eff
        end do
    end subroutine compute_RB_all

    ! Level 3: same-layer source amplitudes
    subroutine TE_f1yf2y_same_layer(q_list, R_down, R_up, N, layer_src, z_src, k0, f1y, f2y)
        integer, intent(in) :: N, layer_src
        complex(dp), intent(in) :: q_list(N), R_down(N), R_up(N)
        real(dp), intent(in) :: z_src, k0
        complex(dp), intent(out) :: f1y, f2y
        
        complex(dp) :: q, RF, RB, ez, emz, den
        
        q = q_list(layer_src)
        RF = R_down(layer_src)
        RB = R_up(layer_src)
        
        ez = exp(q * k0 * z_src)
        emz = exp(-q * k0 * z_src)
        den = 2.0_dp * q * k0 * (1.0_dp - RF * RB)
        
        f1y = (ez + emz * RB) / den
        f2y = (emz + ez * RF) / den
    end subroutine TE_f1yf2y_same_layer

    ! Level 4: propagate downward
    function propagate_down_TE(f1_src, q_list, R_down, N, layer_src, layer_obs, k0, d_list) result(f)
        integer, intent(in) :: N, layer_src, layer_obs
        complex(dp), intent(in) :: f1_src, q_list(N), R_down(N)
        real(dp), intent(in) :: k0, d_list(N-1)
        complex(dp) :: f
        
        complex(dp) :: ql, qlp, wn, wp, wq, Rnext, T, x_val
        real(dp) :: d_l
        integer :: l
        
        f = f1_src
        
        do l = layer_src, layer_obs-1
            d_l = d_list(l)
            ql = q_list(l)
            qlp = q_list(l+1)
            
            wn = qlp / ql
            wp = 1.0_dp + wn
            wq = 1.0_dp - wn
            Rnext = R_down(l+1)
            
            T = 1.0_dp / (wp + wq * Rnext)
            x_val = exp(-ql * k0 * d_l)
            f = T * x_val * f
        end do
    end function propagate_down_TE

    ! Level 4: propagate upward
    function propagate_up_TE(f2_src, q_list, R_up, N, layer_src, layer_obs, k0, d_list) result(f)
        integer, intent(in) :: N, layer_src, layer_obs
        complex(dp), intent(in) :: f2_src, q_list(N), R_up(N)
        real(dp), intent(in) :: k0, d_list(N-1)
        complex(dp) :: f
        
        complex(dp) :: ql, qlp1, wn, mp, mm, r_l, T, x_val
        real(dp) :: d_l
        integer :: l
        
        f = f2_src
        
        do l = layer_src-1, layer_obs, -1
            d_l = d_list(l)
            ql = q_list(l)
            qlp1 = q_list(l+1)
            
            x_val = exp(-ql * k0 * d_l)
            wn = ql / qlp1
            mp = 1.0_dp + wn
            mm = 1.0_dp - wn
            r_l = R_up(l)
            
            T = 1.0_dp / (mp + mm * x_val * r_l * x_val)
            f = x_val * T * f
        end do
    end function propagate_up_TE

    ! Level 5: assemble G_yy
    function gyy_TE(n_list, d_list, N, layer_src, z_src, layer_obs, z_obs, k0, kp) result(Gyy)
        integer, intent(in) :: N, layer_src, layer_obs
        real(dp), intent(in) :: n_list(N), d_list(N-1), z_src, z_obs, k0, kp
        complex(dp) :: Gyy
        
        complex(dp) :: q_list(N), R_down(N), R_up(N)
        complex(dp) :: f1y_src, f2y_src, q_obs, RF_obs, RB_obs
        complex(dp) :: q, RF, RB, direct, cavity, f1y_at_obs, f2y_at_obs
        
        call compute_q_list(n_list, N, k0, kp, q_list)
        call compute_RF_all(n_list, d_list, N, k0, kp, "TE", R_down)
        call compute_RB_all(n_list, d_list, N, k0, kp, "TE", R_up)
        
        call TE_f1yf2y_same_layer(q_list, R_down, R_up, N, layer_src, z_src, k0, f1y_src, f2y_src)
        
        q_obs = q_list(layer_obs)
        RF_obs = R_down(layer_obs)
        RB_obs = R_up(layer_obs)
        
        if (layer_obs == layer_src) then
            q = q_list(layer_src)
            RF = R_down(layer_src)
            RB = R_up(layer_src)
            
            direct = (1.0_dp / (2.0_dp * q * k0)) * exp(-q * k0 * abs(z_obs - z_src))
            cavity = exp(q * k0 * z_obs) * RF * f1y_src + exp(-q * k0 * z_obs) * RB * f2y_src
            Gyy = direct + cavity
            
        else if (layer_obs > layer_src) then
            f1y_at_obs = propagate_down_TE(f1y_src, q_list, R_down, N, layer_src, layer_obs, k0, d_list)
            Gyy = (exp(-q_obs * k0 * z_obs) + exp(q_obs * k0 * z_obs) * RF_obs) * f1y_at_obs
            
        else
            f2y_at_obs = propagate_up_TE(f2y_src, q_list, R_up, N, layer_src, layer_obs, k0, d_list)
            Gyy = (exp(q_obs * k0 * z_obs) + exp(-q_obs * k0 * z_obs) * RB_obs) * f2y_at_obs
        end if
    end function gyy_TE

end module te_greens

! Simple test program
program test_te_greens
    use te_greens
    implicit none
    
    integer, parameter :: N = 4
    real(dp) :: n_list(N), d_list(N-1)
    real(dp) :: z_src, z_obs, k0, kp, wl
    complex(dp) :: Gyy
    
    ! Setup
    n_list = [1.0_dp, 1.5_dp, 1.3_dp, 1.0_dp]
    d_list = [0.0_dp, 2000.0_dp, 1500.0_dp]
    
    z_src = 100.0_dp
    z_obs = 200.0_dp
    wl = 650.0_dp
    k0 = 2.0_dp * 3.14159265359_dp / wl
    kp = 0.3_dp * k0
    
    Gyy = gyy_TE(n_list, d_list, N, 1, z_src, 1, z_obs, k0, kp)
    
    print *, 'Gyy = ', real(Gyy), aimag(Gyy)
    
end program test_te_greens
