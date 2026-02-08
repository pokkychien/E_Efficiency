! split_integration.f90
! Split integration around poles for k_parallel integrals

module split_integration
    use gauss_quadrature
    implicit none
    
contains

    ! Find pole positions in k_parallel space
    ! Poles occur at kp = k0 * sqrt(eps) for each layer
    subroutine find_poles(n_list, N, k0, poles, n_poles)
        integer, intent(in) :: N
        real(dp), intent(in) :: n_list(N), k0
        real(dp), intent(out) :: poles(N)
        integer, intent(out) :: n_poles
        
        integer :: i
        real(dp) :: eps
        
        n_poles = 0
        do i = 1, N
            eps = n_list(i)**2
            if (eps > 0.0_dp) then
                n_poles = n_poles + 1
                poles(n_poles) = k0 * sqrt(eps)
            end if
        end do
    end subroutine find_poles

    ! Create integration intervals avoiding poles
    ! Given [a, b] and pole positions, create sub-intervals with safety margin
    subroutine create_intervals(a, b, poles, n_poles, safety_factor, &
                                 intervals, n_intervals)
        real(dp), intent(in) :: a, b, safety_factor
        integer, intent(in) :: n_poles
        real(dp), intent(in) :: poles(n_poles)
        integer, intent(out) :: n_intervals
        real(dp), intent(out) :: intervals(2, n_poles+1)  ! (start, end) for each interval
        
        real(dp) :: sorted_poles(n_poles), delta
        integer :: i, j
        real(dp) :: temp, current_start
        
        ! Sort poles
        sorted_poles = poles(1:n_poles)
        do i = 1, n_poles - 1
            do j = i + 1, n_poles
                if (sorted_poles(j) < sorted_poles(i)) then
                    temp = sorted_poles(i)
                    sorted_poles(i) = sorted_poles(j)
                    sorted_poles(j) = temp
                end if
            end do
        end do
        
        ! Create intervals
        n_intervals = 0
        current_start = a
        
        do i = 1, n_poles
            delta = safety_factor * sorted_poles(i)
            
            ! Check if pole is within [a, b]
            if (sorted_poles(i) > a .and. sorted_poles(i) < b) then
                ! Add interval before pole
                if (current_start < sorted_poles(i) - delta) then
                    n_intervals = n_intervals + 1
                    intervals(1, n_intervals) = current_start
                    intervals(2, n_intervals) = sorted_poles(i) - delta
                end if
                
                ! Skip pole region
                current_start = sorted_poles(i) + delta
            end if
        end do
        
        ! Add final interval
        if (current_start < b) then
            n_intervals = n_intervals + 1
            intervals(1, n_intervals) = current_start
            intervals(2, n_intervals) = b
        end if
        
        ! Handle case where no poles or no valid intervals
        if (n_intervals == 0) then
            n_intervals = 1
            intervals(1, 1) = a
            intervals(2, 1) = b
        end if
    end subroutine create_intervals

    ! Integrate over multiple intervals using Gauss quadrature
    ! This is the main interface function
    subroutine integrate_split_gauss(intervals, n_intervals, n_gauss, &
                                     kp_all, w_all, n_total)
        integer, intent(in) :: n_intervals, n_gauss
        real(dp), intent(in) :: intervals(2, n_intervals)
        integer, intent(out) :: n_total
        real(dp), intent(out) :: kp_all(n_intervals * n_gauss)
        real(dp), intent(out) :: w_all(n_intervals * n_gauss)
        
        real(dp) :: xi(n_gauss), wi(n_gauss)
        real(dp) :: x_local(n_gauss), w_local(n_gauss)
        integer :: i, j, idx
        
        ! Get base Gauss-Legendre nodes and weights on [-1, 1]
        call gauss_legendre(n_gauss, xi, wi)
        
        ! Transform to each interval
        idx = 0
        do i = 1, n_intervals
            call transform_interval(xi, wi, n_gauss, &
                                   intervals(1, i), intervals(2, i), &
                                   x_local, w_local)
            
            do j = 1, n_gauss
                idx = idx + 1
                kp_all(idx) = x_local(j)
                w_all(idx) = w_local(j)
            end do
        end do
        
        n_total = idx
    end subroutine integrate_split_gauss

end module split_integration
