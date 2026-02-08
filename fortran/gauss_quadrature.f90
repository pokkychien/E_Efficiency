! gauss_quadrature.f90
! Gauss-Legendre quadrature nodes and weights

module gauss_quadrature
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    
contains

    ! Get Gauss-Legendre nodes and weights for n-point quadrature on [-1, 1]
    ! Uses recurrence relations to compute nodes and weights
    subroutine gauss_legendre(n, xi, wi)
        integer, intent(in) :: n
        real(dp), intent(out) :: xi(n), wi(n)
        
        real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
        real(dp), parameter :: tol = 1.0e-14_dp
        integer :: i, j, m
        real(dp) :: z, z1, p1, p2, p3, pp, dz
        
        m = (n + 1) / 2
        
        do i = 1, m
            ! Initial guess for ith root
            z = cos(pi * (real(i, dp) - 0.25_dp) / (real(n, dp) + 0.5_dp))
            
            ! Newton-Raphson iteration
            do j = 1, 100
                p1 = 1.0_dp
                p2 = 0.0_dp
                
                ! Evaluate Legendre polynomial using recurrence
                do k = 1, n
                    p3 = p2
                    p2 = p1
                    p1 = ((2.0_dp * real(k, dp) - 1.0_dp) * z * p2 - &
                          (real(k, dp) - 1.0_dp) * p3) / real(k, dp)
                end do
                
                ! Derivative of Legendre polynomial
                pp = real(n, dp) * (z * p1 - p2) / (z * z - 1.0_dp)
                
                ! Newton update
                z1 = z
                z = z1 - p1 / pp
                
                if (abs(z - z1) < tol) exit
            end do
            
            ! Store nodes and weights
            xi(i) = -z
            xi(n + 1 - i) = z
            wi(i) = 2.0_dp / ((1.0_dp - z * z) * pp * pp)
            wi(n + 1 - i) = wi(i)
        end do
    end subroutine gauss_legendre

    ! Transform nodes from [-1, 1] to [a, b] and scale weights
    subroutine transform_interval(xi, wi, n, a, b, x_ab, w_ab)
        integer, intent(in) :: n
        real(dp), intent(in) :: xi(n), wi(n), a, b
        real(dp), intent(out) :: x_ab(n), w_ab(n)
        
        real(dp) :: scale
        integer :: i
        
        scale = (b - a) / 2.0_dp
        
        do i = 1, n
            x_ab(i) = a + (xi(i) + 1.0_dp) * scale
            w_ab(i) = wi(i) * scale
        end do
    end subroutine transform_interval

    ! Integrate function f on [a, b] using n-point Gauss-Legendre quadrature
    ! f_vals should be function values at the transformed nodes
    function integrate_gauss(f_vals, wi, n) result(integral)
        integer, intent(in) :: n
        real(dp), intent(in) :: f_vals(n), wi(n)
        real(dp) :: integral
        integer :: i
        
        integral = 0.0_dp
        do i = 1, n
            integral = integral + wi(i) * f_vals(i)
        end do
    end function integrate_gauss

    ! Complex version
    function integrate_gauss_complex(f_vals, wi, n) result(integral)
        integer, intent(in) :: n
        complex(dp), intent(in) :: f_vals(n)
        real(dp), intent(in) :: wi(n)
        complex(dp) :: integral
        integer :: i
        
        integral = cmplx(0.0_dp, 0.0_dp, kind=dp)
        do i = 1, n
            integral = integral + wi(i) * f_vals(i)
        end do
    end function integrate_gauss_complex

end module gauss_quadrature
