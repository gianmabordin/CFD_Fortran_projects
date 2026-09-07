program main_mpi
    use mpi
    use iso_fortran_env, only: real64, int64
    implicit none
    
    integer, parameter :: rp = real64
    integer, parameter :: i8 = int64  ! Interi a 64 bit per griglie enormi
    
    integer :: ierr, rank, size
    integer(i8) :: N, N_local, k
    character(len=32) :: arg
    
    real(rp) :: a, b, h, global_sum, pi_calc, pi_exact, err_pi
    real(rp) :: a_local, x_k1, f_k, f_k1, local_sum
    
    ! Variabili per il cronometro
    real(rp) :: t_start, t_end, t_local, t_global

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)

    ! --- LETTURA DA RIGA DI COMANDO ---
    ! Se specifichi un numero lanciando il programma, usa quello. Altrimenti usa 10 milioni.
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg)
        read(arg, *) N
    else
        N = 10000000_i8
    end if

    a = 0.0_rp
    b = 1.0_rp
    h = (b - a) / real(N, rp)
    pi_exact = acos(-1.0_rp)
    
    N_local = N / size
    a_local = a + rank * N_local * h

    ! ==========================================
    ! SINCRONIZZAZIONE E PARTENZA CRONOMETRO
    ! ==========================================
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    t_start = MPI_Wtime()

    local_sum = 0.0_rp
    f_k = sqrt(max(0.0_rp, 1.0_rp - a_local**2))
    
    do k = 0, N_local - 1
        x_k1 = a_local + (k + 1) * h
        f_k1 = sqrt(max(0.0_rp, 1.0_rp - x_k1**2))
        local_sum = local_sum + 0.5_rp * (f_k + f_k1)
        f_k = f_k1
    end do

    call MPI_Reduce(local_sum, global_sum, 1, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    ! ==========================================
    ! FINE CRONOMETRO E CALCOLO TEMPO MASSIMO
    ! ==========================================
    t_end = MPI_Wtime()
    t_local = t_end - t_start
    
    ! Troviamo chi ci ha messo di più
    call MPI_Reduce(t_local, t_global, 1, MPI_DOUBLE_PRECISION, &
                    MPI_MAX, 0, MPI_COMM_WORLD, ierr)

    if (rank == 0) then
        pi_calc = 4.0_rp * h * global_sum 
        err_pi = abs(pi_calc - pi_exact)
        ! Stampiamo un output formattato per essere letto facilmente
        print "(A, I2, A, I12, A, F10.6, A, A, F16.14, A, ES12.4)", &
      "Processi: ", size, " | N Totale: ", N, " | Tempo: ", t_global, " sec", &
      " | Pi: ", pi_calc, " | Err: ", err_pi
    end if

    call MPI_Finalize(ierr)

end program main_mpi