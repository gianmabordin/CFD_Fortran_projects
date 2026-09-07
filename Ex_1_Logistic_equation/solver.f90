module mod_solver
    use mod_param, only: rp, func_template
    use mod_math, only: newton_raphson
    implicit none

contains

    ! IMPLICIT EULER SOLVER
    subroutine compute_imp_euler_nonlin(f, dfdy, exact_sol, phi_0, dt, t_0, N, phi_n, filename, write_output)
        
        procedure(func_template) :: f, dfdy, exact_sol

        character(len=*), intent(in) :: filename
        integer :: iunit

        logical, intent(in), optional :: write_output
        logical :: do_write

        real(rp), intent(in)  :: phi_0, dt, t_0
        integer,  intent(in)  :: N
        real(rp), intent(out) :: phi_n
        
        real(rp) :: t_current, phi_old, phi_new, phi_exact_new
        integer  :: i

        if (present(write_output)) then
            do_write = write_output
        else
            do_write = .true. 
        end if
        
        phi_old   = phi_0
        t_current = t_0

        if (do_write) open(newunit=iunit, file=filename, status='replace')
        if (do_write) write(iunit, *) t_0, phi_0, exact_sol(t_0, phi_0)

        ! Time cycle
        do i = 1, N
            t_current = t_current + dt
            
            call newton_raphson(Residual_Imp_Euler, Jacobian_Imp_Euler, phi_old, phi_new)
            phi_old = phi_new

            phi_exact_new = exact_sol(t_current, phi_0)

            if (do_write) write(iunit, *) t_current, phi_new, phi_exact_new
        end do

        if (do_write) close(iunit)
        
        phi_n = phi_old
        
    ! Residual and derivative
    contains 
    
       
        function Residual_Imp_Euler(phi_guess) result(val)
            real(rp), intent(in) :: phi_guess
            real(rp)             :: val
            val = phi_guess - phi_old - dt * f(t_current, phi_guess)
        end function Residual_Imp_Euler

    
        function Jacobian_Imp_Euler(phi_guess) result(val)
            real(rp), intent(in) :: phi_guess
            real(rp)             :: val
            val = 1.0_rp - dt * dfdy(t_current, phi_guess)
        end function Jacobian_Imp_Euler
    

    end subroutine compute_imp_euler_nonlin


    ! CRANK-NICOLSON SOLVER
     subroutine compute_crank_nic_nonlin(f, dfdy, exact_sol, phi_0, dt, t_0, N, phi_n, filename, write_output)
        
        procedure(func_template) :: f, dfdy, exact_sol

        character(len=*), intent(in) :: filename
        integer :: iunit

        logical, intent(in), optional :: write_output
        logical :: do_write

        real(rp), intent(in)  :: phi_0, dt, t_0
        integer,  intent(in)  :: N
        real(rp), intent(out) :: phi_n
        
        real(rp) :: t_current, phi_old, phi_new, phi_exact_new
        integer  :: i

        if (present(write_output)) then
            do_write = write_output
        else
            do_write = .true. 
        end if
        
        phi_old   = phi_0
        t_current = t_0

        if (do_write) open(newunit=iunit, file=filename, status='replace')
        if (do_write) write(iunit, *) t_0, phi_0, exact_sol(t_0, phi_0)

        ! Time cycle
        do i = 1, N
            t_current = t_current + dt
            
            call newton_raphson(Residual_Crank_Nic, Jacobian_Crank_Nic, phi_old, phi_new)
            phi_old = phi_new

            phi_exact_new = exact_sol(t_current, phi_0)

            if (do_write) write(iunit, *) t_current, phi_new, phi_exact_new
        end do

        if (do_write) close(iunit)
        
        phi_n = phi_old

    ! Residual and derivative
    contains 
    
       
        function Residual_Crank_Nic(phi_guess) result(val)
            real(rp), intent(in) :: phi_guess
            real(rp)             :: val
            val = phi_guess - phi_old - 0.5_rp * dt * ( f(t_current, phi_guess) + f(t_current - dt, phi_old) )
        end function Residual_Crank_Nic

    
        function Jacobian_Crank_Nic(phi_guess) result(val)
            real(rp), intent(in) :: phi_guess
            real(rp)             :: val
            val = 1.0_rp - 0.5_rp * dt * dfdy(t_current, phi_guess)
        end function Jacobian_Crank_Nic
    

    end subroutine compute_crank_nic_nonlin


end module mod_solver