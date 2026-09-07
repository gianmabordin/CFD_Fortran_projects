module mod_simulator
    use mod_param,   only: rp, Nx, c, nu, dx, dt, T_final, ic_type
    use mod_exact,   only: exact_gaussian, exact_harmonic
    use mod_rk3,     only: rk3_step
    implicit none


contains

    ! =========================================================
    ! INTERNAL SUBROUTINE: Runs the time loop and saves data
    ! =========================================================
    subroutine run_simulation(filename, scheme_conv, num_saves)
        character(len=*), intent(in) :: filename
        character(len=3), intent(in) :: scheme_conv
        integer, intent(in)          :: num_saves
        
        real(rp), dimension(Nx) :: x, phi
        real(rp)                :: t, phi_ex
        real(rp)                :: cfl_conv, cfl_diff
        integer                 :: i, n, n_steps, save_interval
        integer                 :: file_unit

        ! 1. Initialize spatial grid (cell centers)
        do i = 1, Nx
            x(i) = (real(i, rp) - 0.5_rp) * dx
        end do

        ! 2. Initialize phi at t=0
        t = 0.0_rp
        do i = 1, Nx
            if (ic_type == 1) then
                phi(i) = exact_gaussian(x(i), t)
            else
                phi(i) = exact_harmonic(x(i), t)
            end if
        end do

        ! 3. Setup time loop and saving intervals
        n_steps = nint(T_final / dt)
        ! Prevent division by zero if num_saves is 1 and saves the code if num saves>steps
        save_interval = max(1, n_steps / max(1, num_saves - 1))

        ! 3.5 Check Stability Conditions (CFL)
        if(abs(c) > 0.0_rp) then
            cfl_conv = abs(c * dt / dx)
            print *, "   |-- Courant (Convection) : ", cfl_conv
        else
            cfl_conv = 0.0_rp
        end if
        
        if (nu > 0.0_rp) then
            cfl_diff = abs(4.0_rp * nu * dt / (dx**2))
            print *, "   |-- Fourier (Diffusion)  : ", cfl_diff
        else
            cfl_diff = 0.0_rp
        end if
        
        ! Security warnings
        if (cfl_conv > 1.0_rp) then
            print *, "   [WARNING] Courant > 1! Possible instability"
        end if
        if (cfl_diff > 1.0_rp) then
            print *, "   [WARNING] Fourier > 1! Possible instability"
        end if

        ! 4. Open file for writing
        open(newunit=file_unit, file=filename, status='replace')
        write(file_unit, *) "Variables: x, t, phi_numerical, phi_exact"

        ! Save Initial Condition (t=0)
        do i = 1, Nx
            write(file_unit, '(4(E15.7, 2X))') x(i), t, phi(i), phi(i)
        end do
        write(file_unit, *) "" ! Empty line to separate blocks for plotting tools

        ! 5. Time Integration Loop
        do n = 1, n_steps
            ! Advance one time step using Runge-Kutta
            phi = rk3_step(phi, scheme_conv)
            t = t + dt

            !Safety switch (to stop calculation in case of NaN)
            if (maxval(abs(phi)) > 1.0e6_rp) then
                print *, "[WARNING] Simulation exploded at time t =", t
                print *, "Forced interruption to avoid corrupted data"
                exit 
            end if

            ! Save data at specific intervals or at the very last step
            if (mod(n, save_interval) == 0 .or. n == n_steps) then
                do i = 1, Nx
                    ! Calculate exact solution for comparison
                    if (ic_type == 1) then
                        phi_ex = exact_gaussian(x(i), t)
                    else
                        phi_ex = exact_harmonic(x(i), t)
                    end if
                    
                    ! Write to file
                    write(file_unit, '(4(E15.7, 2X))') x(i), t, phi(i), phi_ex
                end do
                write(file_unit, *) "" ! Empty line separator
            end if
        end do

        close(file_unit)
    end subroutine run_simulation

end module mod_simulator