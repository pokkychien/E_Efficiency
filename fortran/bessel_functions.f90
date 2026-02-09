! bessel_functions.f90
! Minimal Bessel function implementations for J0 and J2

module bessel_functions
    use multilayer_reflectance, only: dp
    implicit none
    
contains

    ! J0 Bessel function using series expansion
    function J0_series(x) result(J0_val)
        complex(dp), intent(in) :: x
        complex(dp) :: J0_val
        
        complex(dp) :: x2_over_4, term
        real(dp) :: max_abs
        integer :: m, M_max
        
        x2_over_4 = (x * x) / 4.0_dp
        
        ! Adaptive truncation based on magnitude
        max_abs = abs(x)
        if (max_abs < 5.0_dp) then
            M_max = 40
        else if (max_abs < 20.0_dp) then
            M_max = 80
        else
            M_max = 140
        end if
        
        J0_val = cmplx(1.0_dp, 0.0_dp, kind=dp)
        term = cmplx(1.0_dp, 0.0_dp, kind=dp)
        
        do m = 1, M_max
            term = term * (-x2_over_4) / real(m * m, dp)
            J0_val = J0_val + term
        end do
    end function J0_series

    ! J2 Bessel function using series expansion
    function J2_series(x) result(J2_val)
        complex(dp), intent(in) :: x
        complex(dp) :: J2_val
        
        complex(dp) :: x2_over_4, term
        real(dp) :: max_abs
        integer :: m, M_max
        
        x2_over_4 = (x * x) / 4.0_dp
        
        ! Adaptive truncation
        max_abs = abs(x)
        if (max_abs < 5.0_dp) then
            M_max = 50
        else if (max_abs < 20.0_dp) then
            M_max = 120
        else
            M_max = 220
        end if
        
        ! First term (m=0): (x/2)^2 / (0! * 2!) = x^2 / 8
        term = (x * x) / 8.0_dp
        J2_val = term
        
        ! Recurrence: term_{m+1} = term_m * [-(x^2/4)] / [(m+1)(m+3)]
        do m = 0, M_max - 1
            term = term * (-x2_over_4) / real((m + 1) * (m + 3), dp)
            J2_val = J2_val + term
        end do
    end function J2_series

    ! Array version of J0
    subroutine J0_array(x_arr, n, J0_arr)
        integer, intent(in) :: n
        complex(dp), intent(in) :: x_arr(n)
        complex(dp), intent(out) :: J0_arr(n)
        integer :: i
        
        do i = 1, n
            J0_arr(i) = J0_series(x_arr(i))
        end do
    end subroutine J0_array

    ! Array version of J2
    subroutine J2_array(x_arr, n, J2_arr)
        integer, intent(in) :: n
        complex(dp), intent(in) :: x_arr(n)
        complex(dp), intent(out) :: J2_arr(n)
        integer :: i
        
        do i = 1, n
            J2_arr(i) = J2_series(x_arr(i))
        end do
    end subroutine J2_array

end module bessel_functions
