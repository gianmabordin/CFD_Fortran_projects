program main_adaptive_limit
    use iso_fortran_env, only: real64
    implicit none
    
    integer, parameter :: rp = real64
    
    real(rp) :: a, b, x_k1, f_k, f_k1, h, sum_integral, pi_calc, pi_exact
    real(rp) :: error_pi, error_old, tol
    integer  :: N, k

    a = 0.0_rp
    b = 1.0_rp
    pi_exact = acos(-1.0_rp)
    
    ! Tolleranza target: 10 volte la machine precision
    tol = 10.0_rp * epsilon(1.0_rp)
    
    ! Inizializzazione
    N = 10
    error_pi = 1.0_rp 
    error_old = 100.0_rp ! Valore iniziale fittizio molto alto

    print *, "Target tolerance: ", tol
    print *, "---------------------------------------------------"

    ! Ciclo principale
    do while (error_pi > tol .and. N <= 1000000000)
        
        h = (b - a) / real(N, rp)
        sum_integral = 0.0_rp
        f_k = sqrt(max(0.0_rp,1.0_rp - a**2))
        
        do k = 0, N - 1
            x_k1 = a + (k + 1) * h
            f_k1 = sqrt(1.0_rp - x_k1**2)
            sum_integral = sum_integral + 0.5_rp * (f_k + f_k1)
            f_k = f_k1
        end do

        pi_calc = 4.0_rp * h * sum_integral 
        error_pi = abs(pi_calc - pi_exact)

        print *, "N: ", N, " | Error: ", error_pi
        
        ! --- MECCANISMO DI CONTROLLO ROUND-OFF ---
        if (error_pi >= error_old) then
            print *, "---------------------------------------------------"
            print *, ">>> ATTENZIONE: Limite di arrotondamento raggiunto!"
            print *, ">>> L'errore ha smesso di decrescere. Interruzione."
            exit
        end if
        
        ! Aggiorniamo la memoria per il prossimo giro
        error_old = error_pi
        
        ! Incremento più dolce (raddoppio i nodi anziché x10)
        N = N * 2
        
    end do

    print *, "---------------------------------------------------"
    print *, "Final Nodes N:          ", N / 2  ! Stampo il nodo dell'iterazione precedente (la migliore)
    print *, "Calculated value of pi: ", pi_calc
    print *, "Absolute error:         ", error_old ! Stampo l'errore minimo raggiunto

end program main_adaptive_limit