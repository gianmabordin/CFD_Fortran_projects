program main
    use iso_fortran_env, only: real32, real64
    implicit none
    
    integer, parameter :: rp = real64  ! Double precision
    
    ! variables definition
    real(rp) :: a, b, x_k1, f_k, f_k1, h, sum_integral, pi_calc, pi_exact, error_pi
    integer  :: N, k

    ! variables initialization
    a = 0.0_rp
    b = 1.0_rp
    N = 10000000

    ! increment value
    h = (b - a) / real(N, rp)
    sum_integral = 0.0_rp

    ! cycle for calculation
    f_k = sqrt(max(0.0_rp, 1.0_rp - a**2))
    
    do k = 0, N - 1
        x_k1 = a + (k + 1) * h
        f_k1 = sqrt(max(0.0_rp, 1.0_rp - x_k1**2))

        sum_integral = sum_integral + 0.5_rp * (f_k + f_k1)

        f_k = f_k1
    end do

    pi_calc = 4.0_rp * h * sum_integral 
    pi_exact = acos(-1.0_rp)
    error_pi = abs(pi_calc - pi_exact)

    print *, "Nodes N:                ", N
    print *, "Calculated value of pi: ", pi_calc
    print *, "Exact value of pi:      ", pi_exact
    print *, "Absolute error:         ", error_pi

end program main