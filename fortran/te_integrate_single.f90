!======================================================================
! Single-file Fortran implementation of TE Green's function integration
! No modules, no contains - only external subroutines and functions
! Compile: gfortran -O2 -o te_single te_integrate_single.f90
!======================================================================

program te_integrate_single
    implicit none
    
    ! Double precision kind parameter
    integer, parameter :: dp = selected_real_kind(15, 307)
    
    ! Physical constants
    real(dp), parameter :: pi = 3.14159265358979323846_dp
    complex(dp), parameter :: zi = (0.0_dp, 1.0_dp)
    
    ! Layer parameters (4 layers: 0=substrate, 1,2=middle, 3=superstrate)
    integer, parameter :: num_layers = 4
    complex(dp) :: n_list(num_layers)
    real(dp) :: d_list(num_layers-1)  ! thicknesses of layers 1 to num_layers-2
    
    ! Source and observation
    integer :: layer_src, layer_obs
    real(dp) :: z_src, z_obs
    real(dp) :: wl, k0
    
    ! Integration parameters
    real(dp) :: k_parallel_max
    integer :: num_k
    
    ! Rho grid
    integer, parameter :: num_rho = 500
    real(dp) :: rho_min, rho_max
    real(dp) :: rho_list(num_rho)
    complex(dp) :: G_rho(num_rho)
    
    ! Loop variables
    integer :: i
    
    !------------------------------------------------------------------
    ! Setup parameters (matching Python)
    !------------------------------------------------------------------
    n_list(1) = (1.0_dp, 0.0_dp)   ! layer 0 (substrate)
    n_list(2) = (1.0_dp, 0.0_dp)   ! layer 1
    n_list(3) = (1.0_dp, 0.0_dp)   ! layer 2
    n_list(4) = (1.0_dp, 0.0_dp)   ! layer 3 (superstrate)
    
    d_list(1) = 0.0_dp             ! d[0] not used
    d_list(2) = 2.0_dp             ! d[1] = thickness of layer 1
    d_list(3) = 1.5_dp             ! d[2] = thickness of layer 2
    
    ! Source in layer 1, observation in layer 1
    layer_src = 2   ! Fortran index (Python layer 1)
    layer_obs = 2   ! Fortran index (Python layer 1)
    z_src = -0.35_dp
    z_obs = -0.20_dp
    
    wl = 0.65_dp
    k0 = 2.0_dp * pi / wl
    
    ! Integration range
    k_parallel_max = 3.5_dp * k0
    num_k = 2001
    
    ! Rho range
    rho_min = 0.01_dp
    rho_max = 5.0_dp
    
    !------------------------------------------------------------------
    ! Print setup
    !------------------------------------------------------------------
    print '(A)', '=== Single-file Fortran TE Integration ==='
    print '(A, 4F8.3)', 'n_list: ', real(n_list)
    print '(A, 3F8.3)', 'd_list: ', d_list
    print '(A,I3,A,F8.4)', 'layer_src=', layer_src, ', z_src=', z_src
    print '(A,I3,A,F8.4)', 'layer_obs=', layer_obs, ', z_obs=', z_obs
    print '(A,F10.4)', 'wavelength: ', wl
    print '(A,ES12.4E2)', 'k0: ', k0
    print '(A,ES12.4E2)', 'k_parallel_max: ', k_parallel_max
    print '(A,I6)', 'num_k: ', num_k
    print *
    
    !------------------------------------------------------------------
    ! Generate rho grid
    !------------------------------------------------------------------
    do i = 1, num_rho
        rho_list(i) = rho_min + (rho_max - rho_min) * real(i-1, dp) / real(num_rho-1, dp)
    end do
    
    !------------------------------------------------------------------
    ! Compute G(rho) for each rho
    !------------------------------------------------------------------
    print '(A)', 'Computing G(rho)...'
    
    do i = 1, num_rho
        call gyy_TE_rho_full(n_list, d_list, num_layers, &
                             layer_src, z_src, layer_obs, z_obs, &
                             k0, rho_list(i), k_parallel_max, num_k, &
                             zi, G_rho(i))
        
        if (mod(i, 50) == 0) then
            print '(A,I5,A,I5)', 'Progress: ', i, ' / ', num_rho
        end if
    end do
    
    !------------------------------------------------------------------
    ! Save results
    !------------------------------------------------------------------
    call save_results_ext(rho_list, G_rho, num_rho)
    
    print '(A)', 'Results saved to: gyy_rho_fortran.dat'
    print '(A)', 'Done!'

end program te_integrate_single


!======================================================================
! Fresnel reflection coefficient for TE polarization
!======================================================================
function r_te_ext(n1, n2, q1, q2) result(r)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    
    complex(dp), intent(in) :: n1, n2, q1, q2
    complex(dp) :: r
    complex(dp) :: num, den
    
    num = q1 - q2
    den = q1 + q2
    
    if (abs(den) < 1.0e-30_dp) then
        r = (0.0_dp, 0.0_dp)
    else
        r = num / den
    end if
end function r_te_ext


!======================================================================
! Compute q = sqrt(n^2 * k0^2 - k_parallel^2) for all layers
!======================================================================
subroutine compute_q_list(n_list_in, num_layers_in, k0_in, k_parallel_in, q_list_out)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    
    integer, intent(in) :: num_layers_in
    complex(dp), intent(in) :: n_list_in(num_layers_in)
    real(dp), intent(in) :: k0_in, k_parallel_in
    complex(dp), intent(out) :: q_list_out(num_layers_in)
    
    integer :: j
    complex(dp) :: arg
    
    do j = 1, num_layers_in
        arg = n_list_in(j)**2 * k0_in**2 - k_parallel_in**2
        q_list_out(j) = sqrt(arg)
        ! Ensure Im(q) >= 0 for evanescent waves
        if (aimag(q_list_out(j)) < 0.0_dp) then
            q_list_out(j) = -q_list_out(j)
        end if
    end do
end subroutine compute_q_list


!======================================================================
! RF_multilayer: Forward reflection from layer j looking down
!======================================================================
subroutine RF_multilayer(n_list_in, d_list_in, q_list_in, num_layers_in, j_in, RF_out)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    complex(dp), parameter :: zi_local = (0.0_dp, 1.0_dp)
    
    integer, intent(in) :: num_layers_in, j_in
    complex(dp), intent(in) :: n_list_in(num_layers_in), q_list_in(num_layers_in)
    real(dp), intent(in) :: d_list_in(num_layers_in-1)
    complex(dp), intent(out) :: RF_out
    
    integer :: m
    complex(dp) :: r_jm1_j, phase, num, den
    complex(dp) :: r_te_ext
    external :: r_te_ext
    
    ! Start from bottom (layer 0, index 1)
    RF_out = (0.0_dp, 0.0_dp)
    
    do m = 2, j_in
        r_jm1_j = r_te_ext(n_list_in(m-1), n_list_in(m), q_list_in(m-1), q_list_in(m))
        phase = exp(2.0_dp * zi_local * q_list_in(m-1) * d_list_in(m-1))
        num = r_jm1_j + RF_out * phase
        den = 1.0_dp + r_jm1_j * RF_out * phase
        if (abs(den) < 1.0e-30_dp) then
            RF_out = (0.0_dp, 0.0_dp)
        else
            RF_out = num / den
        end if
    end do
end subroutine RF_multilayer


!======================================================================
! RB_multilayer: Backward reflection from layer j looking up
!======================================================================
subroutine RB_multilayer(n_list_in, d_list_in, q_list_in, num_layers_in, j_in, RB_out)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    complex(dp), parameter :: zi_local = (0.0_dp, 1.0_dp)
    
    integer, intent(in) :: num_layers_in, j_in
    complex(dp), intent(in) :: n_list_in(num_layers_in), q_list_in(num_layers_in)
    real(dp), intent(in) :: d_list_in(num_layers_in-1)
    complex(dp), intent(out) :: RB_out
    
    integer :: m
    complex(dp) :: r_jp1_j, phase, num, den
    complex(dp) :: r_te_ext
    external :: r_te_ext
    
    ! Start from top (layer N-1, index num_layers)
    RB_out = (0.0_dp, 0.0_dp)
    
    do m = num_layers_in-1, j_in, -1
        r_jp1_j = r_te_ext(n_list_in(m+1), n_list_in(m), q_list_in(m+1), q_list_in(m))
        phase = exp(2.0_dp * zi_local * q_list_in(m+1) * d_list_in(m))
        num = r_jp1_j + RB_out * phase
        den = 1.0_dp + r_jp1_j * RB_out * phase
        if (abs(den) < 1.0e-30_dp) then
            RB_out = (0.0_dp, 0.0_dp)
        else
            RB_out = num / den
        end if
    end do
end subroutine RB_multilayer


!======================================================================
! Compute RF for all layers
!======================================================================
subroutine compute_RF_all(n_list_in, d_list_in, q_list_in, num_layers_in, RF_all_out)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    
    integer, intent(in) :: num_layers_in
    complex(dp), intent(in) :: n_list_in(num_layers_in), q_list_in(num_layers_in)
    real(dp), intent(in) :: d_list_in(num_layers_in-1)
    complex(dp), intent(out) :: RF_all_out(num_layers_in)
    
    integer :: j
    
    RF_all_out(1) = (0.0_dp, 0.0_dp)
    do j = 2, num_layers_in
        call RF_multilayer(n_list_in, d_list_in, q_list_in, num_layers_in, j, RF_all_out(j))
    end do
end subroutine compute_RF_all


!======================================================================
! Compute RB for all layers
!======================================================================
subroutine compute_RB_all(n_list_in, d_list_in, q_list_in, num_layers_in, RB_all_out)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    
    integer, intent(in) :: num_layers_in
    complex(dp), intent(in) :: n_list_in(num_layers_in), q_list_in(num_layers_in)
    real(dp), intent(in) :: d_list_in(num_layers_in-1)
    complex(dp), intent(out) :: RB_all_out(num_layers_in)
    
    integer :: j
    
    RB_all_out(num_layers_in) = (0.0_dp, 0.0_dp)
    do j = 1, num_layers_in-1
        call RB_multilayer(n_list_in, d_list_in, q_list_in, num_layers_in, j, RB_all_out(j))
    end do
end subroutine compute_RB_all


!======================================================================
! TE source amplitudes f1y, f2y in source layer
!======================================================================
subroutine TE_f1yf2y_same_layer_ext(q_list_in, RF_all_in, RB_all_in, num_layers_in, &
                                     layer_src_in, z_src_in, k0_in, f1y_out, f2y_out)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    
    integer, intent(in) :: num_layers_in, layer_src_in
    complex(dp), intent(in) :: q_list_in(num_layers_in)
    complex(dp), intent(in) :: RF_all_in(num_layers_in), RB_all_in(num_layers_in)
    real(dp), intent(in) :: z_src_in, k0_in
    complex(dp), intent(out) :: f1y_out, f2y_out
    
    complex(dp) :: q, RF, RB, ez, emz, den
    
    q  = q_list_in(layer_src_in)
    RF = RF_all_in(layer_src_in)
    RB = RB_all_in(layer_src_in)
    
    ez  = exp(+q * k0_in * z_src_in)
    emz = exp(-q * k0_in * z_src_in)
    den = 2.0_dp * q * k0_in * (1.0_dp - RF * RB)
    
    f1y_out = (ez + emz * RB) / den
    f2y_out = (emz + ez * RF) / den
end subroutine TE_f1yf2y_same_layer_ext


!======================================================================
! Propagate downward (layer_obs > layer_src)
!======================================================================
function propagate_down_TE_ext(f1_src_in, q_list_in, RF_all_in, num_layers_in, &
                                layer_src_in, layer_obs_in, k0_in, d_list_in) result(f_out)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    
    integer, intent(in) :: num_layers_in, layer_src_in, layer_obs_in
    complex(dp), intent(in) :: f1_src_in, q_list_in(num_layers_in), RF_all_in(num_layers_in)
    real(dp), intent(in) :: k0_in, d_list_in(num_layers_in-1)
    complex(dp) :: f_out
    
    complex(dp) :: ql, qlp, wn, wp, wq, Rnext, T, x_val
    real(dp) :: d_l
    integer :: l
    
    f_out = f1_src_in
    
    do l = layer_src_in, layer_obs_in - 1
        d_l = d_list_in(l)
        ql  = q_list_in(l)
        qlp = q_list_in(l + 1)
        
        wn = qlp / ql
        wp = 1.0_dp + wn
        wq = 1.0_dp - wn
        Rnext = RF_all_in(l + 1)
        
        T = 1.0_dp / (wp + wq * Rnext)
        x_val = exp(-ql * k0_in * d_l)
        f_out = T * x_val * f_out
    end do
end function propagate_down_TE_ext


!======================================================================
! Propagate upward (layer_obs < layer_src)
!======================================================================
function propagate_up_TE_ext(f2_src_in, q_list_in, RB_all_in, num_layers_in, &
                              layer_src_in, layer_obs_in, k0_in, d_list_in) result(f_out)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    
    integer, intent(in) :: num_layers_in, layer_src_in, layer_obs_in
    complex(dp), intent(in) :: f2_src_in, q_list_in(num_layers_in), RB_all_in(num_layers_in)
    real(dp), intent(in) :: k0_in, d_list_in(num_layers_in-1)
    complex(dp) :: f_out
    
    complex(dp) :: ql, qlp1, wn, mp, mm, r_l, T, x_val
    real(dp) :: d_l
    integer :: l
    
    f_out = f2_src_in
    
    do l = layer_src_in - 1, layer_obs_in, -1
        d_l = d_list_in(l)
        ql   = q_list_in(l)
        qlp1 = q_list_in(l + 1)
        
        x_val = exp(-ql * k0_in * d_l)
        wn = ql / qlp1
        mp = 1.0_dp + wn
        mm = 1.0_dp - wn
        r_l = RB_all_in(l)
        
        T = 1.0_dp / (mp + mm * x_val * r_l * x_val)
        f_out = x_val * T * f_out
    end do
end function propagate_up_TE_ext


!======================================================================
! gyy_TE: TE Green's function component (with cross-layer support)
!======================================================================
function gyy_TE_ext(q_list_in, RF_all_in, RB_all_in, d_list_in, num_layers_in, &
                layer_src_in, z_src_in, layer_obs_in, z_obs_in, k0_in) result(G)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    
    integer, intent(in) :: num_layers_in, layer_src_in, layer_obs_in
    complex(dp), intent(in) :: q_list_in(num_layers_in)
    complex(dp), intent(in) :: RF_all_in(num_layers_in), RB_all_in(num_layers_in)
    real(dp), intent(in) :: d_list_in(num_layers_in-1)
    real(dp), intent(in) :: z_src_in, z_obs_in, k0_in
    complex(dp) :: G
    
    ! External functions
    complex(dp) :: propagate_down_TE_ext, propagate_up_TE_ext
    external :: propagate_down_TE_ext, propagate_up_TE_ext
    
    complex(dp) :: f1y_src, f2y_src
    complex(dp) :: q, RF, RB, q_obs, RF_obs, RB_obs
    complex(dp) :: direct, cavity
    complex(dp) :: f1y_at_obs, f2y_at_obs
    
    ! Get source amplitudes
    call TE_f1yf2y_same_layer_ext(q_list_in, RF_all_in, RB_all_in, num_layers_in, &
                                   layer_src_in, z_src_in, k0_in, f1y_src, f2y_src)
    
    q_obs  = q_list_in(layer_obs_in)
    RF_obs = RF_all_in(layer_obs_in)
    RB_obs = RB_all_in(layer_obs_in)
    
    ! Case A: same layer
    if (layer_obs_in == layer_src_in) then
        q  = q_list_in(layer_src_in)
        RF = RF_all_in(layer_src_in)
        RB = RB_all_in(layer_src_in)
        
        direct = (1.0_dp / (2.0_dp * q * k0_in)) * exp(-q * k0_in * abs(z_obs_in - z_src_in))
        cavity = exp(+q * k0_in * z_obs_in) * RF * f1y_src + &
                 exp(-q * k0_in * z_obs_in) * RB * f2y_src
        G = direct + cavity
        
    ! Case B: obs below src (layer_obs > layer_src)
    else if (layer_obs_in > layer_src_in) then
        f1y_at_obs = propagate_down_TE_ext(f1y_src, q_list_in, RF_all_in, num_layers_in, &
                                            layer_src_in, layer_obs_in, k0_in, d_list_in)
        G = (exp(-q_obs * k0_in * z_obs_in) + exp(+q_obs * k0_in * z_obs_in) * RF_obs) * f1y_at_obs
        
    ! Case C: obs above src (layer_obs < layer_src)
    else
        f2y_at_obs = propagate_up_TE_ext(f2y_src, q_list_in, RB_all_in, num_layers_in, &
                                          layer_src_in, layer_obs_in, k0_in, d_list_in)
        G = (exp(+q_obs * k0_in * z_obs_in) + exp(-q_obs * k0_in * z_obs_in) * RB_obs) * f2y_at_obs
    end if
end function gyy_TE_ext


!======================================================================
! dgyy_dzobs_TE: Derivative of gyy_TE with respect to z_obs (with cross-layer)
!======================================================================
function dgyy_dzobs_TE_ext(q_list_in, RF_all_in, RB_all_in, d_list_in, num_layers_in, &
                       layer_src_in, z_src_in, layer_obs_in, z_obs_in, k0_in) result(dG)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    
    integer, intent(in) :: num_layers_in, layer_src_in, layer_obs_in
    complex(dp), intent(in) :: q_list_in(num_layers_in)
    complex(dp), intent(in) :: RF_all_in(num_layers_in), RB_all_in(num_layers_in)
    real(dp), intent(in) :: d_list_in(num_layers_in-1)
    real(dp), intent(in) :: z_src_in, z_obs_in, k0_in
    complex(dp) :: dG
    
    ! External functions
    complex(dp) :: propagate_down_TE_ext, propagate_up_TE_ext
    external :: propagate_down_TE_ext, propagate_up_TE_ext
    
    complex(dp) :: f1y_src, f2y_src
    complex(dp) :: q, RF, RB, a, q_obs, RF_obs, RB_obs, a_obs
    complex(dp) :: ddirect, dcavity, ddress
    complex(dp) :: f1y_at_obs, f2y_at_obs
    real(dp) :: dz, s
    
    ! Get source amplitudes
    call TE_f1yf2y_same_layer_ext(q_list_in, RF_all_in, RB_all_in, num_layers_in, &
                                   layer_src_in, z_src_in, k0_in, f1y_src, f2y_src)
    
    q_obs  = q_list_in(layer_obs_in)
    RF_obs = RF_all_in(layer_obs_in)
    RB_obs = RB_all_in(layer_obs_in)
    
    ! Case A: same layer
    if (layer_obs_in == layer_src_in) then
        q  = q_list_in(layer_src_in)
        RF = RF_all_in(layer_src_in)
        RB = RB_all_in(layer_src_in)
        
        a  = q * k0_in
        dz = z_obs_in - z_src_in
        
        if (dz >= 0.0_dp) then
            s = 1.0_dp
        else
            s = -1.0_dp
        end if
        
        ! ddirect = -0.5 * exp(-a*|dz|) * sign(dz)
        ddirect = (-0.5_dp) * exp(-a * abs(dz)) * s
        
        ! dcavity = a * (exp(+a*z_obs)*RF*f1y - exp(-a*z_obs)*RB*f2y)
        dcavity = a * (exp(+a * z_obs_in) * RF * f1y_src - exp(-a * z_obs_in) * RB * f2y_src)
        
        dG = ddirect + dcavity
        
    ! Case B: obs below src (layer_obs > layer_src)
    else if (layer_obs_in > layer_src_in) then
        f1y_at_obs = propagate_down_TE_ext(f1y_src, q_list_in, RF_all_in, num_layers_in, &
                                            layer_src_in, layer_obs_in, k0_in, d_list_in)
        a_obs = q_obs * k0_in
        ! ddress = -a_obs*exp(-a_obs*z) + a_obs*exp(+a_obs*z)*RF_obs
        ddress = (-a_obs) * exp(-a_obs * z_obs_in) + a_obs * exp(+a_obs * z_obs_in) * RF_obs
        dG = ddress * f1y_at_obs
        
    ! Case C: obs above src (layer_obs < layer_src)
    else
        f2y_at_obs = propagate_up_TE_ext(f2y_src, q_list_in, RB_all_in, num_layers_in, &
                                          layer_src_in, layer_obs_in, k0_in, d_list_in)
        a_obs = q_obs * k0_in
        ! ddress = a_obs*exp(+a_obs*z) - a_obs*exp(-a_obs*z)*RB_obs
        ddress = a_obs * exp(+a_obs * z_obs_in) + (-a_obs) * exp(-a_obs * z_obs_in) * RB_obs
        dG = ddress * f2y_at_obs
    end if
end function dgyy_dzobs_TE_ext


!======================================================================
! Trapezoidal integration for complex arrays
!======================================================================
function trapz_complex_ext(y_arr, dx_in, n_in) result(integral)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    
    integer, intent(in) :: n_in
    complex(dp), intent(in) :: y_arr(n_in)
    real(dp), intent(in) :: dx_in
    complex(dp) :: integral
    
    integer :: i
    
    integral = (0.0_dp, 0.0_dp)
    do i = 1, n_in-1
        integral = integral + (y_arr(i) + y_arr(i+1))
    end do
    integral = integral * dx_in * 0.5_dp
end function trapz_complex_ext


!======================================================================
! Main integration routine: gyy_TE_rho_full
!======================================================================
subroutine gyy_TE_rho_full(n_list_in, d_list_in, num_layers_in, &
                            layer_src_in, z_src_in, layer_obs_in, z_obs_in, &
                            k0_in, rho_in, k_parallel_max_in, num_k_in, &
                            zi_in, result_out)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    
    integer, intent(in) :: num_layers_in, layer_src_in, layer_obs_in, num_k_in
    complex(dp), intent(in) :: n_list_in(num_layers_in)
    real(dp), intent(in) :: d_list_in(num_layers_in-1)
    real(dp), intent(in) :: z_src_in, z_obs_in, k0_in, rho_in, k_parallel_max_in
    complex(dp), intent(in) :: zi_in
    complex(dp), intent(out) :: result_out
    
    ! External function declarations
    complex(dp) :: gyy_TE_ext, dgyy_dzobs_TE_ext, trapz_complex_ext
    external :: gyy_TE_ext, dgyy_dzobs_TE_ext, trapz_complex_ext
    
    ! Local arrays
    real(dp), allocatable :: k_parallel(:)
    complex(dp), allocatable :: q_list_local(:), RF_all_local(:), RB_all_local(:)
    complex(dp), allocatable :: integrand(:,:)
    
    ! DBESJ variables
    real(8) :: bessel_x, bessel_y(3)
    integer :: bessel_nz
    
    real(dp) :: dk, kp, J0_val, J2_val
    complex(dp) :: gyy_val, dgyy_val
    complex(dp) :: I1, I2, I3, I4
    integer :: i
    
    ! Allocate arrays
    allocate(k_parallel(num_k_in))
    allocate(q_list_local(num_layers_in))
    allocate(RF_all_local(num_layers_in))
    allocate(RB_all_local(num_layers_in))
    allocate(integrand(num_k_in, 4))
    
    ! Generate k_parallel grid (avoid k=0)
    dk = k_parallel_max_in / real(num_k_in - 1, dp)
    do i = 1, num_k_in
        k_parallel(i) = real(i - 1, dp) * dk
        if (i == 1) k_parallel(i) = 1.0e-10_dp
    end do
    
    ! Compute integrands
    do i = 1, num_k_in
        kp = k_parallel(i)
        
        ! Compute q_list, RF_all, RB_all
        call compute_q_list(n_list_in, num_layers_in, k0_in, kp, q_list_local)
        call compute_RF_all(n_list_in, d_list_in, q_list_local, num_layers_in, RF_all_local)
        call compute_RB_all(n_list_in, d_list_in, q_list_local, num_layers_in, RB_all_local)
        
        ! Get Green's functions
        gyy_val = gyy_TE_ext(q_list_local, RF_all_local, RB_all_local, d_list_in, num_layers_in, &
                        layer_src_in, z_src_in, layer_obs_in, z_obs_in, k0_in)
        dgyy_val = dgyy_dzobs_TE_ext(q_list_local, RF_all_local, RB_all_local, d_list_in, num_layers_in, &
                                 layer_src_in, z_src_in, layer_obs_in, z_obs_in, k0_in)
        
        ! Bessel functions using DBESJ (compute J0, J1, J2 in one call)
        bessel_x = kp * rho_in
        if (bessel_x < 1.0d-30) then
            J0_val = 1.0_dp
            J2_val = 0.0_dp
        else
            call DBESJ(bessel_x, 0.0d0, 3, bessel_y, bessel_nz)
            J0_val = bessel_y(1)  ! J0
            J2_val = bessel_y(3)  ! J2
        end if
        
        ! 4 integrands (TE mode only, matching Python I1-I4)
        ! I5-I10 require TM mode (gxx, dgxx, gzx) - not yet implemented
        integrand(i, 1) = kp * J0_val * gyy_val                           ! I1: kp * J0 * Gyy
        integrand(i, 2) = kp * J0_val * conjg(dgyy_val)                   ! I2: kp * J0 * conj(DGyy)
        integrand(i, 3) = kp * J2_val * gyy_val                           ! I3: kp * J2 * Gyy
        integrand(i, 4) = kp * J2_val * conjg(dgyy_val)                   ! I4: kp * J2 * conj(DGyy)
    end do
    
    ! Integrate each component (TE mode: I1-I4 only)
    I1 = trapz_complex_ext(integrand(:,1), dk, num_k_in)
    I2 = trapz_complex_ext(integrand(:,2), dk, num_k_in)
    I3 = trapz_complex_ext(integrand(:,3), dk, num_k_in)
    I4 = trapz_complex_ext(integrand(:,4), dk, num_k_in)
    
    ! Combine TE contribution (I5-I10 TM mode not implemented)
    result_out = I1 + I2 + I3 + I4
    
    ! Deallocate
    deallocate(k_parallel, q_list_local, RF_all_local, RB_all_local, integrand)
    
end subroutine gyy_TE_rho_full


!======================================================================
! Save results to file
!======================================================================
subroutine save_results_ext(rho_list_in, G_rho_in, n_in)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    
    integer, intent(in) :: n_in
    real(dp), intent(in) :: rho_list_in(n_in)
    complex(dp), intent(in) :: G_rho_in(n_in)
    
    integer :: i, unit_num
    
    unit_num = 20
    open(unit=unit_num, file='gyy_rho_single.dat', status='replace')
    
    write(unit_num, '(A)') '# rho  Re(G)  Im(G)  |G|'
    do i = 1, n_in
        write(unit_num, '(4ES18.10E2)') rho_list_in(i), real(G_rho_in(i)), &
                                        aimag(G_rho_in(i)), abs(G_rho_in(i))
    end do
    
    close(unit_num)
end subroutine save_results_ext
