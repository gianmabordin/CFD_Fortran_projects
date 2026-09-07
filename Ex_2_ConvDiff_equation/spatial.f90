module mod_spatial
    use mod_param, only: rp
    implicit none

contains

    ! =====================================================================
    ! Scheme I: First-order Forward (FW1) for 1st derivative
    ! =====================================================================
    function deriv_FW1(phi, dx) result(dphi)
        real(rp), dimension(:), intent(in) :: phi
        real(rp), intent(in)               :: dx
        real(rp), dimension(size(phi))     :: dphi
        integer                            :: i, N

        N = size(phi)
        
        ! Internal nodes
        do i = 1, N - 1
            dphi(i) = (phi(i+1) - phi(i)) / dx
        end do
        
        ! Periodic boundary for the last node (i+1 becomes 1)
        dphi(N) = (phi(1) - phi(N)) / dx
    end function deriv_FW1

    ! =====================================================================
    ! Scheme II: First-order Backward (BW1) for 1st derivative
    ! =====================================================================
    function deriv_BW1(phi, dx) result(dphi)
        real(rp), dimension(:), intent(in) :: phi
        real(rp), intent(in)               :: dx
        real(rp), dimension(size(phi))     :: dphi
        integer                            :: i, N

        N = size(phi)
        
        ! Periodic boundary for the first node (i-1 becomes N)
        dphi(1) = (phi(1) - phi(N)) / dx
        
        ! Internal nodes
        do i = 2, N
            dphi(i) = (phi(i) - phi(i-1)) / dx
        end do
    end function deriv_BW1

    ! =====================================================================
    ! Scheme III: Second-order Centered (CS2) for 1st derivative
    ! =====================================================================
    function deriv_CS2(phi, dx) result(dphi)
        real(rp), dimension(:), intent(in) :: phi
        real(rp), intent(in)               :: dx
        real(rp), dimension(size(phi))     :: dphi
        integer                            :: i, N

        N = size(phi)
        
        ! Periodic boundary for the first node (left neighbor is N)
        dphi(1) = (phi(2) - phi(N)) / (2.0_rp * dx)
        
        ! Internal nodes
        do i = 2, N - 1
            dphi(i) = (phi(i+1) - phi(i-1)) / (2.0_rp * dx)
        end do
        
        ! Periodic boundary for the last node (right neighbor is 1)
        dphi(N) = (phi(1) - phi(N-1)) / (2.0_rp * dx)
    end function deriv_CS2

    ! =====================================================================
    ! Scheme A: Fourth-order Centered (CS4) for 2nd derivative (Diffusion)
    ! Formula: (-phi_{i-2} + 16*phi_{i-1} - 30*phi_i + 16*phi_{i+1} - phi_{i+2}) / (12*dx^2)
    ! =====================================================================
    function deriv2_CS4(phi, dx) result(d2phi)
        real(rp), dimension(:), intent(in) :: phi
        real(rp), intent(in)               :: dx
        real(rp), dimension(size(phi))     :: d2phi
        integer                            :: i, N
        real(rp)                           :: dx212

        N = size(phi)
        dx212 = 12.0_rp * (dx**2)  ! Precompute the denominator for efficiency
        
        ! 1. Left boundary (i = 1) -> needs i-1=N, i-2=N-1
        d2phi(1) = (-phi(N-1) + 16.0_rp*phi(N) - 30.0_rp*phi(1) + 16.0_rp*phi(2) - phi(3)) / dx212
        
        ! 2. Second node (i = 2) -> needs i-2=N
        d2phi(2) = (-phi(N) + 16.0_rp*phi(1) - 30.0_rp*phi(2) + 16.0_rp*phi(3) - phi(4)) / dx212

        ! 3. Internal nodes
        do i = 3, N - 2
            d2phi(i) = (-phi(i-2) + 16.0_rp*phi(i-1) - 30.0_rp*phi(i) + 16.0_rp*phi(i+1) - phi(i+2)) / dx212
        end do

        ! 4. Second-to-last node (i = N-1) -> needs i+2=1
        d2phi(N-1) = (-phi(N-3) + 16.0_rp*phi(N-2) - 30.0_rp*phi(N-1) + 16.0_rp*phi(N) - phi(1)) / dx212
        
        ! 5. Right boundary (i = N) -> needs i+1=1, i+2=2
        d2phi(N) = (-phi(N-2) + 16.0_rp*phi(N-1) - 30.0_rp*phi(N) + 16.0_rp*phi(1) - phi(2)) / dx212
        
    end function deriv2_CS4

end module mod_spatial