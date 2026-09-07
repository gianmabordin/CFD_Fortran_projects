module mod_param
    use iso_fortran_env, only: real32, real64
    implicit none

    ! IMPORTANT: use comment in order to select the precision desired
    integer,parameter :: rp=real64  ! Double precision
    !integer,parameter :: rp=real32   ! Single precision

    

    abstract interface
        
        function func_template(t, phi) result(val)
            import :: rp
            real(rp), intent(in) :: t, phi
            real(rp)             :: val
        end function func_template

        function R_func_template(x) result(val)
            import :: rp
            real(rp), intent(in) :: x
            real(rp)             :: val
        end function R_func_template

    end interface

    
end module mod_param