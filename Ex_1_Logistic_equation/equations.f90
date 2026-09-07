module mod_equations
    use mod_param, only: rp
    implicit none

    ! New variables in the module
    real(rp) :: r_val
    real(rp) :: K_val


contains

    subroutine pass_var(r_in, K_in)
        
        real(rp), intent(in) :: r_in, K_in
        
        r_val = r_in
        K_val = K_in

    end subroutine pass_var

   
    function f(t, phi) result(val)
        real(rp), intent(in) :: t, phi
        real(rp)             :: val
        
        val = r_val * phi * (1.0_rp - (phi / K_val))

    end function f

    
    function dfdy(t, phi) result(val)
        real(rp), intent(in) :: t, phi
        real(rp)             :: val
        
        val = r_val * (1.0_rp - ((2.0_rp * phi) / K_val))

    end function dfdy


    function exact_sol(t, phi_0) result(val)
        real(rp), intent(in) :: t, phi_0
        real(rp)             :: val
        
        val = K_val / (1.0_rp + ((K_val - phi_0) / phi_0) * exp(-r_val * t))
        
    end function exact_sol

end module mod_equations