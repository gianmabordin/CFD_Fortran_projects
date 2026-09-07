module mod_types
  !
  use iso_fortran_env, only: real64
  !
  implicit none
  !
  public
  !
  ! ============================================================
  ! Define the real precision used throughout the code.
  !
  ! The parameter rp is used as a shorthand for the floating-point
  ! kind, ensuring consistency and portability across all modules.
  ! In this case, double precision (real64) is selected.
  ! ============================================================
  integer, parameter :: rp = real64
  !
end module mod_types