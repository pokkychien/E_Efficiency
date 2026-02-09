! te_integrate.f90
! k_parallel integration for multilayer Green's functions
! Fortran translation of te_integrate_plot.py :: gyy_TE_rho()

module te_integrate
    use multilayer_reflectance, only: dp, zi
    use te_greens
    use tm_greens
    use bessel_functions
    implicit none
    
contains

    ! Main integration function: compute the combined Green's function integral
    ! This is the Fortran equivalent of gyy_TE_rho() in te_integrate_plot.py
    function gyy_TE_rho_full(n_list, d_list, N, layer_src, z_src, layer_obs, z_obs, &
                        k0, rho, k_parallel_max, num_k) result(G_result)
        integer, intent(in) :: N, layer_src, layer_obs, num_k
        real(dp), intent(in) :: n_list(N), d_list(N-1), z_src, z_obs
        real(dp), intent(in) :: k0, rho, k_parallel_max
        complex(dp) :: G_result
        
        real(dp) :: dk, kp_val, pref
        real(dp), allocatable :: kps(:)
        complex(dp), allocatable :: Gyykp(:), DGyykp(:), Gxxkp(:), DGxxkp(:), Gzxkp(:)
        complex(dp), allocatable :: J0_arr(:), J2_arr(:)
        complex(dp), allocatable :: integrand_1(:), integrand_2(:), integrand_3(:), integrand_4(:)
        complex(dp), allocatable :: integrand_5(:), integrand_6(:), integrand_7(:), integrand_8(:)
        complex(dp), allocatable :: integrand_9(:), integrand_10(:)
        complex(dp) :: I1, I2, I3, I4, I5, I6, I7, I8, I9, I10
        real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
        integer :: i
        
        ! Allocate arrays
        allocate(kps(num_k))
        allocate(Gyykp(num_k), DGyykp(num_k), Gxxkp(num_k), DGxxkp(num_k), Gzxkp(num_k))
        allocate(J0_arr(num_k), J2_arr(num_k))
        allocate(integrand_1(num_k), integrand_2(num_k), integrand_3(num_k), integrand_4(num_k))
        allocate(integrand_5(num_k), integrand_6(num_k), integrand_7(num_k), integrand_8(num_k))
        allocate(integrand_9(num_k), integrand_10(num_k))
        
        ! --- midpoint grid (avoid hitting branch point exactly) ---
        dk = k_parallel_max / real(num_k, dp)
        do i = 1, num_k
            kps(i) = (real(i, dp) - 0.5_dp) * dk
        end do
        
        ! --- Compute Green's functions at each kp ---
        do i = 1, num_k
            kp_val = kps(i)
            
            ! TE: Gyy and dGyy/dz_obs
            Gyykp(i) = gyy_TE(n_list, d_list, N, layer_src, z_src, layer_obs, z_obs, k0, kp_val)
            DGyykp(i) = dgyy_dzobs_TE(n_list, d_list, N, layer_src, z_src, layer_obs, z_obs, k0, kp_val)
            
            ! TM: Gxx, dGxx/dz_obs, Gzx
            Gxxkp(i) = gxx_TM(n_list, d_list, N, layer_src, z_src, layer_obs, z_obs, k0, kp_val)
            DGxxkp(i) = dgxx_dzobs_TM(n_list, d_list, N, layer_src, z_src, layer_obs, z_obs, k0, kp_val)
            Gzxkp(i) = gzx_TM(n_list, d_list, N, layer_src, z_src, layer_obs, z_obs, k0, kp_val)
        end do
        
        ! --- Bessel functions ---
        do i = 1, num_k
            J0_arr(i) = J0_series(cmplx(kps(i) * rho, 0.0_dp, kind=dp))
            J2_arr(i) = J2_series(cmplx(kps(i) * rho, 0.0_dp, kind=dp))
        end do
        
        ! --- Form integrands (matching Python exactly) ---
        pref = 1.0_dp / pi
        
        do i = 1, num_k
            integrand_1(i) = kps(i) * J0_arr(i) * Gyykp(i)
            integrand_2(i) = kps(i) * J0_arr(i) * conjg(DGyykp(i))
            integrand_3(i) = kps(i) * J2_arr(i) * Gyykp(i)
            integrand_4(i) = kps(i) * J2_arr(i) * conjg(DGyykp(i))
            integrand_5(i) = kps(i) * J0_arr(i) * Gxxkp(i)
            integrand_6(i) = kps(i) * J0_arr(i) * conjg(DGxxkp(i))
            integrand_7(i) = kps(i) * J2_arr(i) * Gxxkp(i)
            integrand_8(i) = kps(i) * J2_arr(i) * conjg(DGxxkp(i))
            integrand_9(i) = kps(i)**2 * J0_arr(i) * conjg(Gzxkp(i))
            integrand_10(i) = kps(i)**2 * J2_arr(i) * conjg(Gzxkp(i))
        end do
        
        ! --- Integrate using trapezoidal rule ---
        I1 = pref * trapz_complex(integrand_1, kps, num_k)
        I2 = pref * trapz_complex(integrand_2, kps, num_k)
        I3 = pref * trapz_complex(integrand_3, kps, num_k)
        I4 = pref * trapz_complex(integrand_4, kps, num_k)
        I5 = pref * trapz_complex(integrand_5, kps, num_k)
        I6 = pref * trapz_complex(integrand_6, kps, num_k)
        I7 = pref * trapz_complex(integrand_7, kps, num_k)
        I8 = pref * trapz_complex(integrand_8, kps, num_k)
        I9 = pref * trapz_complex(integrand_9, kps, num_k)
        I10 = pref * trapz_complex(integrand_10, kps, num_k)
        
        ! --- Final combination (matching Python exactly) ---
        G_result = I1*I2 + I3*I4 + I5*I6 + I7*I8 &
                 + zi*I5*I9 + zi*I5*I10 &
                 + I5*I2 + I7*I4 + I1*I6 + I3*I8 &
                 + zi*I1*I9 + zi*I1*I10
        
        ! Deallocate
        deallocate(kps, Gyykp, DGyykp, Gxxkp, DGxxkp, Gzxkp)
        deallocate(J0_arr, J2_arr)
        deallocate(integrand_1, integrand_2, integrand_3, integrand_4)
        deallocate(integrand_5, integrand_6, integrand_7, integrand_8)
        deallocate(integrand_9, integrand_10)
        
    end function gyy_TE_rho_full

    ! Trapezoidal integration (like numpy.trapz)
    function trapz_complex(y, x, n) result(integral)
        integer, intent(in) :: n
        complex(dp), intent(in) :: y(n)
        real(dp), intent(in) :: x(n)
        complex(dp) :: integral
        integer :: i
        
        integral = cmplx(0.0_dp, 0.0_dp, kind=dp)
        do i = 1, n-1
            integral = integral + 0.5_dp * (y(i) + y(i+1)) * (x(i+1) - x(i))
        end do
    end function trapz_complex

    ! Compute G_result for an array of rho values and write to file
    subroutine compute_and_save(n_list, d_list, N, layer_src, z_src, layer_obs, z_obs, &
                                k0, rho_min, rho_max, n_rho, k_parallel_max, num_k, filename)
        integer, intent(in) :: N, layer_src, layer_obs, n_rho, num_k
        real(dp), intent(in) :: n_list(N), d_list(N-1), z_src, z_obs
        real(dp), intent(in) :: k0, rho_min, rho_max, k_parallel_max
        character(len=*), intent(in) :: filename
        
        real(dp) :: rho, drho
        complex(dp) :: G_val
        integer :: i, unit_num
        
        if (n_rho > 1) then
            drho = (rho_max - rho_min) / real(n_rho - 1, dp)
        else
            drho = 0.0_dp
        end if
        unit_num = 20
        
        open(unit=unit_num, file=filename, status='replace', action='write')
        write(unit_num, '(A)') '# rho, Re(G), Im(G)'
        
        do i = 1, n_rho
            rho = rho_min + real(i-1, dp) * drho
            G_val = gyy_TE_rho_full(n_list, d_list, N, layer_src, z_src, layer_obs, z_obs, &
                               k0, rho, k_parallel_max, num_k)
            write(unit_num, '(3ES20.12)') rho, real(G_val), aimag(G_val)
            
            ! Progress output
            if (mod(i, 50) == 0) then
                print '(A,I5,A,I5)', 'Progress: ', i, ' / ', n_rho
            end if
        end do
        
        close(unit_num)
        print '(A,A)', 'Results saved to: ', trim(filename)
        
    end subroutine compute_and_save

end module te_integrate

! Main test program
program main_te_integrate
    use multilayer_reflectance, only: dp
    use te_integrate
    implicit none
    
    integer, parameter :: N = 4
    real(dp) :: n_list(N), d_list(N-1)
    real(dp) :: z_src, z_obs, k0, wl, k_parallel_max
    real(dp) :: rho_min, rho_max
    integer :: n_rho, num_k, layer_src, layer_obs
    real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
    
    ! --- Setup geometry (matching Python demo_plot_TE) ---
    n_list = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
    d_list = [0.0_dp, 2.0_dp, 1.5_dp]
    
    ! Source and observation positions
    ! NOTE: Fortran uses 1-based indexing, Python uses 0-based
    ! Python layer_src=0 corresponds to Fortran layer_src=1
    layer_src = 1
    layer_obs = 1
    z_src = -0.35_dp
    z_obs = -0.20_dp
    
    ! Wavelength
    wl = 0.650_dp
    k0 = 2.0_dp * pi / wl
    
    ! Integration parameters
    k_parallel_max = 3.5_dp * k0
    num_k = 2001
    
    ! Output range
    rho_min = 0.0_dp
    rho_max = 1.30_dp
    n_rho = 500
    
    print '(A)', '=== Fortran TE Integration Demo ==='
    print '(A,4F8.3)', 'n_list: ', n_list
    print '(A,3F8.3)', 'd_list: ', d_list
    print '(A,I3,A,F8.4)', 'layer_src=', layer_src, ', z_src=', z_src
    print '(A,I3,A,F8.4)', 'layer_obs=', layer_obs, ', z_obs=', z_obs
    print '(A,F8.4)', 'wavelength: ', wl
    print '(A,ES12.4)', 'k0: ', k0
    print '(A,ES12.4)', 'k_parallel_max: ', k_parallel_max
    print '(A,I6)', 'num_k: ', num_k
    print '(A)', ''
    print '(A)', 'Computing G(rho)...'
    
    call compute_and_save(n_list, d_list, N, layer_src, z_src, layer_obs, z_obs, &
                          k0, rho_min, rho_max, n_rho, k_parallel_max, num_k, &
                          'gyy_rho_fortran.dat')
    
    print '(A)', ''
    print '(A)', 'Done! Use plot_fortran_results.py to visualize.'
    
end program main_te_integrate
