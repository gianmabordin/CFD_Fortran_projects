program main
    use mod_param, only: rp
    use mod_solver, only: compute_imp_euler_nonlin, compute_crank_nic_nonlin       
    use mod_equations, only : f, dfdy, exact_sol, pass_var  
    implicit none

    ! Variables
    integer, parameter :: len_vec = 5
    real(rp) :: phi_0, t_0, T
    real(rp), dimension(len_vec) :: dt
    real(rp) :: r_val, K_val
    integer  :: N, i_dt, iunit
    
    ! Results
    real(rp) :: phi_n_imp, phi_n_cn, phi_exact

    ! Data
    t_0    = 0.0_rp
    T      = 10.0_rp       
    dt     = [1.0_rp, 1.0e-1_rp, 1.0e-2_rp, 1.0e-3_rp, 1.0e-4_rp]   
    phi_0  = 0.5_rp 
    r_val  = 1.0_rp
    K_val  = 10.0_rp
    
    call pass_var(r_val, K_val)

    ! Task 1
    
    i_dt = 1
    N = int((T-t_0) / dt(i_dt))
    
    call compute_imp_euler_nonlin(f, dfdy, exact_sol, phi_0, dt(i_dt), t_0, N, phi_n_imp, 'Task_1eul.txt')
    call compute_crank_nic_nonlin(f, dfdy, exact_sol, phi_0, dt(i_dt), t_0, N, phi_n_cn, 'Task_1cn.txt')

    ! Task 2

    open(newunit=iunit, file='task_2_errors.txt', status='replace')

    do i_dt = 1, len_vec

        N = int((T-t_0) / dt(i_dt))
       
        call compute_imp_euler_nonlin(f, dfdy, exact_sol, phi_0, dt(i_dt), t_0, N, phi_n_imp, 'temp_eul.txt', write_output=.false.)
        call compute_crank_nic_nonlin(f, dfdy, exact_sol, phi_0, dt(i_dt), t_0, N, phi_n_cn, 'temp_cn.txt', write_output=.false.)

        write(iunit, *) dt(i_dt), abs(exact_sol(T, phi_0) - phi_n_imp), abs(exact_sol(T, phi_0) - phi_n_cn)

    end do 

    close(iunit)




end program main