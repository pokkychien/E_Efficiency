!======================================================================
! Single-file Fortran implementation of TE Green's function integration
! No modules - only subroutines and functions
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
    call save_results(rho_list, G_rho, num_rho)
    
    print '(A)', 'Results saved to: gyy_rho_fortran.dat'
    print '(A)', 'Done!'

contains

    !==================================================================
    ! Fresnel reflection coefficient for TE polarization
    !==================================================================
    function r_te(n1, n2, q1, q2) result(r)
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
    end function r_te
    
    !==================================================================
    ! Compute q = sqrt(n^2 * k0^2 - k_parallel^2) for all layers
    !==================================================================
    subroutine compute_q_list(n_list, num_layers, k0, k_parallel, q_list)
        integer, intent(in) :: num_layers
        complex(dp), intent(in) :: n_list(num_layers)
        real(dp), intent(in) :: k0, k_parallel
        complex(dp), intent(out) :: q_list(num_layers)
        
        integer :: j
        complex(dp) :: arg
        
        do j = 1, num_layers
            arg = n_list(j)**2 * k0**2 - k_parallel**2
            q_list(j) = sqrt(arg)
            ! Ensure Im(q) >= 0 for evanescent waves
            if (aimag(q_list(j)) < 0.0_dp) then
                q_list(j) = -q_list(j)
            end if
        end do
    end subroutine compute_q_list
    
    !==================================================================
    ! RF_multilayer: Forward reflection from layer j looking down
    !==================================================================
    subroutine RF_multilayer(n_list, d_list, q_list, num_layers, j, RF)
        integer, intent(in) :: num_layers, j
        complex(dp), intent(in) :: n_list(num_layers), q_list(num_layers)
        real(dp), intent(in) :: d_list(num_layers-1)
        complex(dp), intent(out) :: RF
        
        integer :: m
        complex(dp) :: r_jm1_j, phase, num, den
        
        ! Start from bottom (layer 0, index 1)
        RF = (0.0_dp, 0.0_dp)
        
        do m = 2, j
            r_jm1_j = r_te(n_list(m-1), n_list(m), q_list(m-1), q_list(m))
            phase = exp(2.0_dp * zi * q_list(m-1) * d_list(m-1))
            num = r_jm1_j + RF * phase
            den = 1.0_dp + r_jm1_j * RF * phase
            if (abs(den) < 1.0e-30_dp) then
                RF = (0.0_dp, 0.0_dp)
            else
                RF = num / den
            end if
        end do
    end subroutine RF_multilayer
    
    !==================================================================
    ! RB_multilayer: Backward reflection from layer j looking up
    !==================================================================
    subroutine RB_multilayer(n_list, d_list, q_list, num_layers, j, RB)
        integer, intent(in) :: num_layers, j
        complex(dp), intent(in) :: n_list(num_layers), q_list(num_layers)
        real(dp), intent(in) :: d_list(num_layers-1)
        complex(dp), intent(out) :: RB
        
        integer :: m
        complex(dp) :: r_jp1_j, phase, num, den
        
        ! Start from top (layer N-1, index num_layers)
        RB = (0.0_dp, 0.0_dp)
        
        do m = num_layers-1, j, -1
            r_jp1_j = r_te(n_list(m+1), n_list(m), q_list(m+1), q_list(m))
            phase = exp(2.0_dp * zi * q_list(m+1) * d_list(m))
            num = r_jp1_j + RB * phase
            den = 1.0_dp + r_jp1_j * RB * phase
            if (abs(den) < 1.0e-30_dp) then
                RB = (0.0_dp, 0.0_dp)
            else
                RB = num / den
            end if
        end do
    end subroutine RB_multilayer
    
    !==================================================================
    ! Compute RF for all layers
    !==================================================================
    subroutine compute_RF_all(n_list, d_list, q_list, num_layers, RF_all)
        integer, intent(in) :: num_layers
        complex(dp), intent(in) :: n_list(num_layers), q_list(num_layers)
        real(dp), intent(in) :: d_list(num_layers-1)
        complex(dp), intent(out) :: RF_all(num_layers)
        
        integer :: j
        
        RF_all(1) = (0.0_dp, 0.0_dp)
        do j = 2, num_layers
            call RF_multilayer(n_list, d_list, q_list, num_layers, j, RF_all(j))
        end do
    end subroutine compute_RF_all
    
    !==================================================================
    ! Compute RB for all layers
    !==================================================================
    subroutine compute_RB_all(n_list, d_list, q_list, num_layers, RB_all)
        integer, intent(in) :: num_layers
        complex(dp), intent(in) :: n_list(num_layers), q_list(num_layers)
        real(dp), intent(in) :: d_list(num_layers-1)
        complex(dp), intent(out) :: RB_all(num_layers)
        
        integer :: j
        
        RB_all(num_layers) = (0.0_dp, 0.0_dp)
        do j = 1, num_layers-1
            call RB_multilayer(n_list, d_list, q_list, num_layers, j, RB_all(j))
        end do
    end subroutine compute_RB_all
    
    !==================================================================
    ! gyy_TE: TE Green's function component
    !==================================================================
    function gyy_TE(q_list, RF_all, RB_all, d_list, num_layers, &
                    layer_src, z_src, layer_obs, z_obs, k0) result(G)
        integer, intent(in) :: num_layers, layer_src, layer_obs
        complex(dp), intent(in) :: q_list(num_layers)
        complex(dp), intent(in) :: RF_all(num_layers), RB_all(num_layers)
        real(dp), intent(in) :: d_list(num_layers-1)
        real(dp), intent(in) :: z_src, z_obs, k0
        complex(dp) :: G
        
        complex(dp) :: q_j, RF_j, RB_j, d_j
        complex(dp) :: denom, prefactor
        complex(dp) :: term1, term2, term3, term4
        complex(dp) :: exp_sum, exp_diff
        
        integer :: j
        
        j = layer_src  ! Assuming layer_src == layer_obs
        
        q_j = q_list(j)
        RF_j = RF_all(j)
        RB_j = RB_all(j)
        d_j = d_list(j)
        
        denom = 1.0_dp - RF_j * RB_j * exp(2.0_dp * zi * q_j * d_j)
        prefactor = zi / (2.0_dp * q_j) / denom
        
        exp_sum = exp(zi * q_j * (z_obs + z_src))
        exp_diff = exp(zi * q_j * abs(z_obs - z_src))
        
        ! Four terms
        term1 = exp_diff
        term2 = RF_j * exp(zi * q_j * 2.0_dp * d_j) * exp(zi * q_j * (-z_obs - z_src))
        term3 = RB_j * exp_sum
        term4 = RF_j * RB_j * exp(zi * q_j * 2.0_dp * d_j) * &
                exp(zi * q_j * (z_obs - z_src) * sign(1.0_dp, z_src - z_obs))
        
        G = prefactor * (term1 + term2 + term3 + term4)
    end function gyy_TE
    
    !==================================================================
    ! dgyy_dzobs_TE: Derivative of gyy_TE with respect to z_obs
    !==================================================================
    function dgyy_dzobs_TE(q_list, RF_all, RB_all, d_list, num_layers, &
                           layer_src, z_src, layer_obs, z_obs, k0) result(dG)
        integer, intent(in) :: num_layers, layer_src, layer_obs
        complex(dp), intent(in) :: q_list(num_layers)
        complex(dp), intent(in) :: RF_all(num_layers), RB_all(num_layers)
        real(dp), intent(in) :: d_list(num_layers-1)
        real(dp), intent(in) :: z_src, z_obs, k0
        complex(dp) :: dG
        
        complex(dp) :: q_j, RF_j, RB_j, d_j
        complex(dp) :: denom, prefactor
        complex(dp) :: term1, term2, term3, term4
        real(dp) :: sgn
        
        integer :: j
        
        j = layer_src
        
        q_j = q_list(j)
        RF_j = RF_all(j)
        RB_j = RB_all(j)
        d_j = d_list(j)
        
        denom = 1.0_dp - RF_j * RB_j * exp(2.0_dp * zi * q_j * d_j)
        prefactor = zi / (2.0_dp * q_j) / denom
        
        if (z_obs >= z_src) then
            sgn = 1.0_dp
        else
            sgn = -1.0_dp
        end if
        
        term1 = zi * q_j * sgn * exp(zi * q_j * abs(z_obs - z_src))
        term2 = -zi * q_j * RF_j * exp(zi * q_j * 2.0_dp * d_j) * &
                exp(zi * q_j * (-z_obs - z_src))
        term3 = zi * q_j * RB_j * exp(zi * q_j * (z_obs + z_src))
        term4 = zi * q_j * sgn * RF_j * RB_j * exp(zi * q_j * 2.0_dp * d_j) * &
                exp(zi * q_j * (z_obs - z_src) * (-sgn))
        
        dG = prefactor * (term1 + term2 + term3 + term4)
    end function dgyy_dzobs_TE
    
    !==================================================================
    ! Bessel J0 series expansion
    !==================================================================
    function J0_series(x) result(j0)
        real(dp), intent(in) :: x
        real(dp) :: j0
        
        real(dp) :: term, x2
        integer :: k
        
        j0 = 1.0_dp
        term = 1.0_dp
        x2 = (x / 2.0_dp)**2
        
        do k = 1, 50
            term = -term * x2 / real(k, dp)**2
            j0 = j0 + term
            if (abs(term) < 1.0e-15_dp * abs(j0)) exit
        end do
    end function J0_series
    
    !==================================================================
    ! Bessel J2 series expansion
    !==================================================================
    function J2_series(x) result(j2)
        real(dp), intent(in) :: x
        real(dp) :: j2
        
        real(dp) :: term, x2
        integer :: k
        
        if (abs(x) < 1.0e-30_dp) then
            j2 = 0.0_dp
            return
        end if
        
        j2 = (x / 2.0_dp)**2 / 2.0_dp
        term = j2
        x2 = (x / 2.0_dp)**2
        
        do k = 1, 50
            term = -term * x2 / (real(k, dp) * real(k + 2, dp))
            j2 = j2 + term
            if (abs(term) < 1.0e-15_dp * abs(j2)) exit
        end do
    end function J2_series
    
    !==================================================================
    ! Trapezoidal integration for complex arrays
    !==================================================================
    function trapz_complex(y, dx, n) result(integral)
        integer, intent(in) :: n
        complex(dp), intent(in) :: y(n)
        real(dp), intent(in) :: dx
        complex(dp) :: integral
        
        integer :: i
        
        integral = (0.0_dp, 0.0_dp)
        do i = 1, n-1
            integral = integral + (y(i) + y(i+1))
        end do
        integral = integral * dx * 0.5_dp
    end function trapz_complex
    
    !==================================================================
    ! Main integration routine: gyy_TE_rho_full
    !==================================================================
    subroutine gyy_TE_rho_full(n_list, d_list, num_layers, &
                                layer_src, z_src, layer_obs, z_obs, &
                                k0, rho, k_parallel_max, num_k, &
                                zi_in, result)
        integer, intent(in) :: num_layers, layer_src, layer_obs, num_k
        complex(dp), intent(in) :: n_list(num_layers)
        real(dp), intent(in) :: d_list(num_layers-1)
        real(dp), intent(in) :: z_src, z_obs, k0, rho, k_parallel_max
        complex(dp), intent(in) :: zi_in
        complex(dp), intent(out) :: result
        
        ! Local arrays
        real(dp), allocatable :: k_parallel(:)
        complex(dp), allocatable :: q_list(:), RF_all(:), RB_all(:)
        complex(dp), allocatable :: integrand(:,:)
        
        real(dp) :: dk, kp, J0_val, J2_val
        complex(dp) :: gyy_val, dgyy_val, q_j
        complex(dp) :: I1, I2, I3, I4, I5, I6, I7, I8, I9, I10
        integer :: i
        
        ! Allocate arrays
        allocate(k_parallel(num_k))
        allocate(q_list(num_layers))
        allocate(RF_all(num_layers))
        allocate(RB_all(num_layers))
        allocate(integrand(num_k, 10))
        
        ! Generate k_parallel grid (avoid k=0)
        dk = k_parallel_max / real(num_k - 1, dp)
        do i = 1, num_k
            k_parallel(i) = real(i - 1, dp) * dk
            if (i == 1) k_parallel(i) = 1.0e-10_dp
        end do
        
        ! Compute integrands
        do i = 1, num_k
            kp = k_parallel(i)
            
            ! Compute q_list, RF_all, RB_all
            call compute_q_list(n_list, num_layers, k0, kp, q_list)
            call compute_RF_all(n_list, d_list, q_list, num_layers, RF_all)
            call compute_RB_all(n_list, d_list, q_list, num_layers, RB_all)
            
            ! Get Green's functions
            gyy_val = gyy_TE(q_list, RF_all, RB_all, d_list, num_layers, &
                            layer_src, z_src, layer_obs, z_obs, k0)
            dgyy_val = dgyy_dzobs_TE(q_list, RF_all, RB_all, d_list, num_layers, &
                                     layer_src, z_src, layer_obs, z_obs, k0)
            
            q_j = q_list(layer_src)
            
            ! Bessel functions
            J0_val = J0_series(kp * rho)
            J2_val = J2_series(kp * rho)
            
            ! 10 integrands (matching Python)
            integrand(i, 1) = kp * J0_val * gyy_val                           ! Integral 1
            integrand(i, 2) = kp * J2_val * gyy_val                           ! Integral 2
            integrand(i, 3) = kp * J0_val * dgyy_val                          ! Integral 3
            integrand(i, 4) = kp * J2_val * dgyy_val                          ! Integral 4
            integrand(i, 5) = kp**3 / q_j**2 * J0_val * gyy_val               ! Integral 5
            integrand(i, 6) = kp**3 / q_j**2 * J2_val * gyy_val               ! Integral 6
            integrand(i, 7) = kp**3 / q_j**2 * J0_val * dgyy_val              ! Integral 7
            integrand(i, 8) = kp**3 / q_j**2 * J2_val * dgyy_val              ! Integral 8
            integrand(i, 9) = kp**2 / q_j * J0_val * gyy_val                  ! Integral 9
            integrand(i, 10) = kp**2 / q_j * J2_val * gyy_val                 ! Integral 10
        end do
        
        ! Integrate each component
        I1 = trapz_complex(integrand(:,1), dk, num_k)
        I2 = trapz_complex(integrand(:,2), dk, num_k)
        I3 = trapz_complex(integrand(:,3), dk, num_k)
        I4 = trapz_complex(integrand(:,4), dk, num_k)
        I5 = trapz_complex(integrand(:,5), dk, num_k)
        I6 = trapz_complex(integrand(:,6), dk, num_k)
        I7 = trapz_complex(integrand(:,7), dk, num_k)
        I8 = trapz_complex(integrand(:,8), dk, num_k)
        I9 = trapz_complex(integrand(:,9), dk, num_k)
        I10 = trapz_complex(integrand(:,10), dk, num_k)
        
        ! Combine (simplified - just sum for demonstration)
        result = I1 + I2 + I3 + I4 + I5 + I6 + I7 + I8 + I9 + I10
        
        ! Deallocate
        deallocate(k_parallel, q_list, RF_all, RB_all, integrand)
        
    end subroutine gyy_TE_rho_full
    
    !==================================================================
    ! Save results to file
    !==================================================================
    subroutine save_results(rho_list, G_rho, n)
        integer, intent(in) :: n
        real(dp), intent(in) :: rho_list(n)
        complex(dp), intent(in) :: G_rho(n)
        
        integer :: i, unit_num
        
        unit_num = 20
        open(unit=unit_num, file='gyy_rho_single.dat', status='replace')
        
        write(unit_num, '(A)') '# rho  Re(G)  Im(G)  |G|'
        do i = 1, n
            write(unit_num, '(4ES18.10E2)') rho_list(i), real(G_rho(i)), &
                                            aimag(G_rho(i)), abs(G_rho(i))
        end do
        
        close(unit_num)
    end subroutine save_results

end program te_integrate_single
