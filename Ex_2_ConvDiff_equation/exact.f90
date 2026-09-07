module mod_exact
    use mod_param, only: rp, L, c, nu, x0, sigma
    implicit none

contains

    ! =====================================================================
    ! Exact solution for the Gaussian initial condition
    ! =====================================================================
    function exact_gaussian(x, t) result(phi_val)
        real(rp), intent(in) :: x
        real(rp), intent(in) :: t
        real(rp)             :: phi_val
        real(rp)             :: dist, denominator
        
        ! Calculate the distance from the center, moving with velocity c.
        dist = x - x0 - c * t
        
        ! Periodic wrapping (modulo arithmetic) to keep the pulse inside [0, L].
        ! This correctly handles the periodic boundary conditions.
        dist = dist - L * anint(dist / L)
        
        ! Calculate the denominator of the exact solution formula
        denominator = sigma**2 + 2.0_rp * nu * t
        
        phi_val = (sigma / sqrt(denominator)) * &
                  exp( -(dist**2) / (2.0_rp * denominator) )
    end function exact_gaussian

    ! =====================================================================
    ! Exact solution for the Harmonic initial condition
    ! =====================================================================
    function exact_harmonic(x, t) result(phi_val)
        real(rp), intent(in) :: x
        real(rp), intent(in) :: t
        real(rp)             :: phi_val
        real(rp), parameter  :: pi = acos(-1.0_rp)
        real(rp)             :: k
        
        ! Wave number
        k = 2.0_rp * pi / L
        
        ! Formula from the assignment
        phi_val = exp(-nu * (k**2) * t) * sin(k * (x - c * t))
    end function exact_harmonic

end module mod_exact