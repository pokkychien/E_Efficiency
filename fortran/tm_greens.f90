! tm_greens.f90
! TM Green's function engine for multilayer structures

module tm_greens
    use te_greens
    implicit none
    
contains

    ! Level 3: TM source amplitudes for g_xx
    subroutine TM_f1xf2x_same_layer(q_list, R_down, R_up, N, layer_src, z_src, k0, kp, f1x, f2x)
        integer, intent(in) :: N, layer_src
        complex(dp), intent(in) :: q_list(N), R_down(N), R_up(N)
        real(dp), intent(in) :: z_src, k0, kp
        complex(dp), intent(out) :: f1x, f2x
        
        complex(dp) :: q, RF, RB, ez, emz, den
        real(dp) :: e2
        
        q = q_list(layer_src)
        RF = R_down(layer_src)
        RB = R_up(layer_src)
        
        e2 = (kp**2 - (abs(q) * k0)**2) / (k0**2)
        ez = exp(q * k0 * z_src)
        emz = exp(-q * k0 * z_src)
        
        den = 2.0_dp * e2 * (k0 / q) * (1.0_dp - RF * RB)
        f1x = -(ez + emz * RB) / den
        f2x = -(emz + ez * RF) / den
    end subroutine TM_f1xf2x_same_layer

    ! TM derivative source amplitudes
    subroutine TM_df1xdf2x_same_layer(q_list, R_down, R_up, N, layer_src, z_src, k0, kp, df1x, df2x)
        integer, intent(in) :: N, layer_src
        complex(dp), intent(in) :: q_list(N), R_down(N), R_up(N)
        real(dp), intent(in) :: z_src, k0, kp
        complex(dp), intent(out) :: df1x, df2x
        
        complex(dp) :: q, RF, RB, ez, emz, den
        real(dp) :: e2
        
        q = q_list(layer_src)
        RF = R_down(layer_src)
        RB = R_up(layer_src)
        
        e2 = (kp**2 - (abs(q) * k0)**2) / (k0**2)
        ez = exp(q * k0 * z_src)
        emz = exp(-q * k0 * z_src)
        
        den = 2.0_dp * e2 * (k0 / q) * (1.0_dp - RF * RB)
        df1x = -(q * k0) * (ez - RB * emz) / den
        df2x = (q * k0) * (emz - RF * ez) / den
    end subroutine TM_df1xdf2x_same_layer

    ! TM z-component source amplitudes
    subroutine TM_f1zf2z_same_layer(q_list, R_down, R_up, N, layer_src, z_src, k0, kp, f1z, f2z)
        integer, intent(in) :: N, layer_src
        complex(dp), intent(in) :: q_list(N), R_down(N), R_up(N)
        real(dp), intent(in) :: z_src, k0, kp
        complex(dp), intent(out) :: f1z, f2z
        
        complex(dp) :: q, RF, RB, ez, emz, den
        real(dp) :: e2
        
        q = q_list(layer_src)
        RF = R_down(layer_src)
        RB = R_up(layer_src)
        
        e2 = (kp**2 - (abs(q) * k0)**2) / (k0**2)
        ez = exp(q * k0 * z_src)
        emz = exp(-q * k0 * z_src)
        
        den = (2.0_dp * k0**2 * e2 / (zi * kp)) * (1.0_dp - RF * RB)
        f1z = (ez - emz * RB) / den
        f2z = (-emz + ez * RF) / den
    end subroutine TM_f1zf2z_same_layer

    ! Propagate downward for TM
    function propagate_down_TM(f1_src, q_list, R_down, eps_list, N, layer_src, layer_obs, k0, d_list) result(f)
        integer, intent(in) :: N, layer_src, layer_obs
        complex(dp), intent(in) :: f1_src, q_list(N), R_down(N)
        real(dp), intent(in) :: eps_list(N), k0, d_list(N-1)
        complex(dp) :: f
        
        complex(dp) :: ql, qlp, wn, wp, wq, Rnext, T, x_val
        real(dp) :: d_l, eps_l, eps_lp
        integer :: l
        
        f = f1_src
        
        do l = layer_src, layer_obs-1
            d_l = d_list(l)
            ql = q_list(l)
            qlp = q_list(l+1)
            eps_l = eps_list(l)
            eps_lp = eps_list(l+1)
            
            wn = (qlp / ql) * (eps_l / eps_lp)
            wp = 1.0_dp + wn
            wq = 1.0_dp - wn
            Rnext = R_down(l+1)
            
            T = 1.0_dp / (wp + wq * Rnext)
            x_val = exp(-ql * k0 * d_l)
            f = T * x_val * f
        end do
    end function propagate_down_TM

    ! Propagate upward for TM
    function propagate_up_TM(f2_src, q_list, R_up, eps_list, N, layer_src, layer_obs, k0, d_list) result(f)
        integer, intent(in) :: N, layer_src, layer_obs
        complex(dp), intent(in) :: f2_src, q_list(N), R_up(N)
        real(dp), intent(in) :: eps_list(N), k0, d_list(N-1)
        complex(dp) :: f
        
        complex(dp) :: ql, qlp1, wn, mp, mm, r_l, T, x_val
        real(dp) :: d_l, eps_l, eps_lp
        integer :: l
        
        f = f2_src
        
        do l = layer_src-1, layer_obs, -1
            d_l = d_list(l)
            ql = q_list(l)
            qlp1 = q_list(l+1)
            eps_l = eps_list(l)
            eps_lp = eps_list(l+1)
            
            wn = (ql / qlp1) * (eps_lp / eps_l)
            mp = 1.0_dp + wn
            mm = 1.0_dp - wn
            r_l = R_up(l)
            x_val = exp(-ql * k0 * d_l)
            
            T = 1.0_dp / (mp + mm * x_val * r_l * x_val)
            f = x_val * T * f
        end do
    end function propagate_up_TM

    ! Level 5: g_xx component
    function gxx_TM(n_list, d_list, N, layer_src, z_src, layer_obs, z_obs, k0, kp) result(Gxx)
        integer, intent(in) :: N, layer_src, layer_obs
        real(dp), intent(in) :: n_list(N), d_list(N-1), z_src, z_obs, k0, kp
        complex(dp) :: Gxx
        
        real(dp) :: eps_list(N), eps_obs, eps
        complex(dp) :: q_list(N), R_down(N), R_up(N)
        complex(dp) :: f1x_src, f2x_src, q_obs, RF_obs, RB_obs
        complex(dp) :: q, RF, RB, direct, refl, dress, f1x_at_obs, f2x_at_obs
        
        eps_list = n_list**2
        
        call compute_q_list(n_list, N, k0, kp, q_list)
        call compute_RF_all(n_list, d_list, N, k0, kp, "TM", R_down)
        call compute_RB_all(n_list, d_list, N, k0, kp, "TM", R_up)
        
        call TM_f1xf2x_same_layer(q_list, R_down, R_up, N, layer_src, z_src, k0, kp, f1x_src, f2x_src)
        
        q_obs = q_list(layer_obs)
        eps_obs = eps_list(layer_obs)
        RF_obs = R_down(layer_obs)
        RB_obs = R_up(layer_obs)
        
        if (layer_obs == layer_src) then
            q = q_list(layer_src)
            eps = eps_list(layer_src)
            RF = R_down(layer_src)
            RB = R_up(layer_src)
            
            direct = (-q / (2.0_dp * k0 * eps)) * exp(-q * k0 * abs(z_obs - z_src))
            refl = exp(q * k0 * z_obs) * RF * f1x_src + exp(-q * k0 * z_obs) * RB * f2x_src
            Gxx = direct + refl
            
        else if (layer_obs > layer_src) then
            f1x_at_obs = propagate_down_TM(f1x_src, q_list, R_down, eps_list, N, layer_src, layer_obs, k0, d_list)
            dress = exp(-q_obs * k0 * z_obs) + exp(q_obs * k0 * z_obs) * RF_obs
            Gxx = dress * f1x_at_obs
            
        else
            f2x_at_obs = propagate_up_TM(f2x_src, q_list, R_up, eps_list, N, layer_src, layer_obs, k0, d_list)
            dress = exp(q_obs * k0 * z_obs) + exp(-q_obs * k0 * z_obs) * RB_obs
            Gxx = dress * f2x_at_obs
        end if
    end function gxx_TM

    ! g_zx component
    function gzx_TM(n_list, d_list, N, layer_src, z_src, layer_obs, z_obs, k0, kp) result(gzx)
        integer, intent(in) :: N, layer_src, layer_obs
        real(dp), intent(in) :: n_list(N), d_list(N-1), z_src, z_obs, k0, kp
        complex(dp) :: gzx
        
        real(dp) :: eps_list(N), sgn
        complex(dp) :: q_list(N), R_down(N), R_up(N)
        complex(dp) :: f1x_src, f2x_src, q_obs, RF_obs, RB_obs
        complex(dp) :: q, RF, RB, f1x_at_obs, f2x_at_obs
        
        eps_list = n_list**2
        
        call compute_q_list(n_list, N, k0, kp, q_list)
        call compute_RF_all(n_list, d_list, N, k0, kp, "TM", R_down)
        call compute_RB_all(n_list, d_list, N, k0, kp, "TM", R_up)
        
        call TM_f1xf2x_same_layer(q_list, R_down, R_up, N, layer_src, z_src, k0, kp, f1x_src, f2x_src)
        
        q_obs = q_list(layer_obs)
        RF_obs = R_down(layer_obs)
        RB_obs = R_up(layer_obs)
        
        if (layer_obs == layer_src) then
            q = q_list(layer_src)
            RF = R_down(layer_src)
            RB = R_up(layer_src)
            
            sgn = sign(1.0_dp, z_obs - z_src)
            
            gzx = -zi * kp / (2.0_dp * k0**2) * sgn * exp(-q * k0 * abs(z_obs - z_src)) &
                  + (zi * kp / q) * (-exp(q * k0 * z_obs) * RF * f1x_src &
                                     + exp(-q * k0 * z_obs) * RB * f2x_src)
            
        else if (layer_obs > layer_src) then
            f1x_at_obs = propagate_down_TM(f1x_src, q_list, R_down, eps_list, N, layer_src, layer_obs, k0, d_list)
            gzx = (-zi * kp / (q_obs * k0)) * (-exp(-q_obs * k0 * z_obs) * f1x_at_obs &
                                                + exp(q_obs * k0 * z_obs) * RF_obs * f1x_at_obs)
            
        else
            f2x_at_obs = propagate_up_TM(f2x_src, q_list, R_up, eps_list, N, layer_src, layer_obs, k0, d_list)
            gzx = (-zi * kp / (q_obs * k0)) * (exp(q_obs * k0 * z_obs) * f2x_at_obs &
                                                - exp(-q_obs * k0 * z_obs) * RB_obs * f2x_at_obs)
        end if
    end function gzx_TM

    ! d/dz_obs of G_xx for TM polarization
    function dgxx_dzobs_TM(n_list, d_list, N, layer_src, z_src, layer_obs, z_obs, k0, kp) result(dGxx)
        integer, intent(in) :: N, layer_src, layer_obs
        real(dp), intent(in) :: n_list(N), d_list(N-1), z_src, z_obs, k0, kp
        complex(dp) :: dGxx
        
        real(dp) :: eps_list(N), dz, s
        complex(dp) :: q_list(N), R_down(N), R_up(N)
        complex(dp) :: f1x_src, f2x_src, q_obs, RF_obs, RB_obs
        complex(dp) :: q, RF, RB, a, ddirect, drefl, ddress
        complex(dp) :: f1x_at_obs, f2x_at_obs, a_obs, eps
        
        eps_list = n_list**2
        
        call compute_q_list(n_list, N, k0, kp, q_list)
        call compute_RF_all(n_list, d_list, N, k0, kp, "TM", R_down)
        call compute_RB_all(n_list, d_list, N, k0, kp, "TM", R_up)
        
        call TM_f1xf2x_same_layer(q_list, R_down, R_up, N, layer_src, z_src, k0, kp, f1x_src, f2x_src)
        
        q_obs = q_list(layer_obs)
        RF_obs = R_down(layer_obs)
        RB_obs = R_up(layer_obs)
        
        ! Case A: same layer
        if (layer_obs == layer_src) then
            q = q_list(layer_src)
            eps = cmplx(eps_list(layer_src), 0.0_dp, kind=dp)
            RF = R_down(layer_src)
            RB = R_up(layer_src)
            
            a = q * k0
            dz = z_obs - z_src
            s = sign(1.0_dp, dz)
            
            ! ddirect = (q^2/(2*eps)) * exp(-a*|dz|) * sign(dz)
            ddirect = (q*q/(2.0_dp*eps)) * exp(-a * abs(dz)) * s
            
            ! drefl = a * (exp(+a*z_obs)*RF*f1x - exp(-a*z_obs)*RB*f2x)
            drefl = a * (exp(a * z_obs) * RF * f1x_src - exp(-a * z_obs) * RB * f2x_src)
            
            dGxx = ddirect + drefl
            
        ! Case B: obs below src
        else if (layer_obs > layer_src) then
            f1x_at_obs = propagate_down_TM(f1x_src, q_list, R_down, eps_list, N, layer_src, layer_obs, k0, d_list)
            
            a_obs = q_obs * k0
            ! ddress = -a_obs*exp(-a_obs*z_obs) + a_obs*exp(+a_obs*z_obs)*RF_obs
            ddress = (-a_obs) * exp(-a_obs * z_obs) + a_obs * exp(a_obs * z_obs) * RF_obs
            
            dGxx = ddress * f1x_at_obs
            
        ! Case C: obs above src
        else
            f2x_at_obs = propagate_up_TM(f2x_src, q_list, R_up, eps_list, N, layer_src, layer_obs, k0, d_list)
            
            a_obs = q_obs * k0
            ! ddress = a_obs*exp(+a_obs*z_obs) - a_obs*exp(-a_obs*z_obs)*RB_obs
            ddress = a_obs * exp(a_obs * z_obs) + (-a_obs) * exp(-a_obs * z_obs) * RB_obs
            
            dGxx = ddress * f2x_at_obs
        end if
    end function dgxx_dzobs_TM

end module tm_greens
