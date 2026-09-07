module mod_rk3
    use mod_param, only: rp, c, nu, dx, dt
    use mod_spatial, only: deriv_FW1, deriv_BW1, deriv_CS2, deriv2_CS4
    implicit none

contains

    ! =====================================================================
    ! Right-Hand Side (RHS) calculator: L(phi) = -c * dphi/dx + nu * d2phi/dx2
    ! =====================================================================
    function calc_rhs(phi, scheme_conv) result(rhs)
        real(rp), dimension(:), intent(in) :: phi
        character(len=3), intent(in)       :: scheme_conv
        real(rp), dimension(size(phi))     :: rhs
        real(rp), dimension(size(phi))     :: dphi_dx, d2phi_dx2

        ! 1. Calculate Convection term (First derivative)
        if (abs(c) > 0.0_rp) then
            select case (scheme_conv)
                case ('FW1')
                    dphi_dx = deriv_FW1(phi, dx)
                case ('BW1')
                    dphi_dx = deriv_BW1(phi, dx)
                case ('CS2')
                    dphi_dx = deriv_CS2(phi, dx)
                case default
                    print *, "ERROR: Unknown convection scheme!"
                    stop
            end select
        else 
            dphi_dx = 0.0_rp
        end if

        ! 2. Calculate Diffusion term (Second derivative)
        ! If nu is 0, we can skip it to save calculation time
        if (nu > 0.0_rp) then
            d2phi_dx2 = deriv2_CS4(phi, dx)
        else
            d2phi_dx2 = 0.0_rp
        end if

        ! 3. Assemble the full physical equation
        rhs = -c * dphi_dx + nu * d2phi_dx2

    end function calc_rhs

    ! =====================================================================
    ! 3rd Order Runge-Kutta (Coherent with Professor's Table - Heun's RK3)
    ! =====================================================================
    function rk3_step(phi_n, scheme_conv) result(phi_next)
        real(rp), dimension(:), intent(in) :: phi_n
        character(len=3), intent(in)       :: scheme_conv
        real(rp), dimension(size(phi_n))   :: phi_next
        real(rp), dimension(size(phi_n))   :: phi_1, phi_2
        real(rp), dimension(size(phi_n))   :: k1, k2, k3

        ! Stage 1 (uses alpha_1 = 0)
        k1 = calc_rhs(phi_n, scheme_conv)
        phi_1 = phi_n + (1.0_rp / 3.0_rp) * dt * k1   ! alpha_2 = 1/3

        ! Stage 2 
        k2 = calc_rhs(phi_1, scheme_conv)
        phi_2 = phi_n + (2.0_rp / 3.0_rp) * dt * k2   ! alpha_3 = 2/3

        ! Stage 3
        k3 = calc_rhs(phi_2, scheme_conv)

        ! Final Assembly using beta coefficients (beta_1 = 1/4, beta_2 = 0, beta_3 = 3/4)
        phi_next = phi_n + dt * (0.25_rp * k1 + 0.75_rp * k3)

    end function rk3_step

end module mod_rk3