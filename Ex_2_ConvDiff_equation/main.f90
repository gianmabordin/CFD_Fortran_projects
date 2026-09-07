program main
    use mod_param
    use mod_simulator, only: run_simulation
    implicit none

    print *, "========================================"
    print *, " EXERCISE 2: Convection-Diffusion Solver"
    print *, "========================================"

    ! =========================================================
    ! TASK 1A: Pure Convection (Gaussian)
    ! =========================================================
    print *, "-> Running Task 1A (CS2, Gaussian)..."
    Nx = 100
    dx = L / real(Nx, rp)
    dt = 0.01_rp
    T_final = 1.0_rp
    c = 0.1_rp
    nu = 0.0_rp
    ic_type = 1 ! 1 = Gaussian

    ! Call the solver
    call run_simulation("task1a.txt", "CS2", 6)


    ! =========================================================
    ! TASK 1B: Pure Convection (Harmonic)
    ! =========================================================
    print *, "-> Running Task 1B (Harmonic)..."
    Nx = 50
    dx = L / real(Nx, rp)
    dt = 0.25_rp
    T_final = 100.0_rp
    c = 0.1_rp
    nu = 0.0_rp
    ic_type = 2 ! 2 = Harmonic

    print *, "   - FW1 Scheme..."
    call run_simulation("task1b_fw1.txt", "FW1", 101)
    
    print *, "   - BW1 Scheme..."
    call run_simulation("task1b_bw1.txt", "BW1", 6)
    
    print *, "   - CS2 Scheme..."
    call run_simulation("task1b_cs2.txt", "CS2", 6)

    ! =========================================================================
    ! TASK 2: Pure Diffusion
    ! ν = 0.01, c = 0.0, Gaussian IC, CS4 Scheme
    ! =========================================================================
    print *, ""
    print *, "-> Running Task 2 (Pure Diffusion, CS4)..."
    
    Nx      = 50
    dx      = L / real(Nx, rp)
    dt      = 1.0e-3_rp       
    T_final = 3.0_rp 
    c       = 0.0_rp         
    nu      = 0.01_rp         
    ic_type = 1 ! 1 = Gaussian
    
    call run_simulation("task2.txt", "CS4", 6)

    ! =========================================================================
    ! TASK 3: Convection-Diffusion
    ! ν = 0.01, c = 0.2, Gaussian IC
    ! =========================================================================
    print *, ""
    print *, "-> Running Task 3 (Convection-Diffusion)..."
    
    Nx      = 100
    dx      = L / real(Nx, rp)
    dt      = 0.001_rp        ! Scelto per soddisfare CFL e Fourier
    T_final = 2.5_rp          ! L'onda farà mezzo giro del dominio (c*T = 0.5)
    c       = 0.2_rp
    nu      = 0.01_rp
    ic_type = 1               ! 1 = Gaussian Initial Condition
    
    ! Chiamiamo il CS2 per la convezione (la diffusione userà il CS4)
    call run_simulation("task3.txt", "CS2", 5)


    print *, "========================================"
    print *, " ALL TASKS COMPLETED SUCCESSFULLY!"
    print *, "========================================"

end program main