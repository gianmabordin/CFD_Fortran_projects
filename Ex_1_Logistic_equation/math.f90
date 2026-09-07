module mod_math
    use mod_param, only: rp, R_func_template
    implicit none

contains

    subroutine newton_raphson(R_func, dR_func, x_guess, x_root)
        
        procedure(R_func_template) :: R_func, dR_func

        ! Declarations
        real(rp), intent(in)  :: x_guess
        real(rp), intent(out) :: x_root
        
        integer  :: iter
        integer, parameter  :: max_iter = 50
        real(rp), parameter :: tolerance = epsilon(1.0_rp) * 100.0_rp
        real(rp) :: x_curr, x_next, error_val, deriv
        
        ! Alghoritm
        x_curr = x_guess  
        
        do iter = 1, max_iter
            deriv = dR_func(x_curr)
            if (abs(deriv) < tiny(1.0_rp)) then
                print *, "Newton-Raphson: zero derivative at iteration", iter
                stop
            end if
            x_next = x_curr - ( R_func(x_curr) / deriv )
            error_val = abs(x_next - x_curr)
            if (error_val < tolerance * (abs(x_next) + 1.0_rp)) then
                x_root = x_next
                return  
            end if
            x_curr = x_next
        end do
        
        print *, "CRITICAL ERROR: Newton-Raphson doesn't converge!"
        stop
        
    end subroutine newton_raphson


end module mod_math