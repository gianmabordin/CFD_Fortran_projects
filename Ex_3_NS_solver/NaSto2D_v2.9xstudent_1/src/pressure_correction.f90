module mod_pressure_correction
  !
  use mod_types,  only: rp
  use mod_common, only: nx, ny
  use mod_common, only: dx, dy
  !
  implicit none
  !
  private
  !
  public :: PC_ComputeDivergence
  public :: PC_BuildRHS
  public :: PC_CorrectVelocity
  !
contains
  !
  ! ============================================================
  ! Compute the cell-centered divergence of a staggered velocity field.
  ! The divergence is evaluated over each control volume using the
  ! face-normal velocity components stored on the staggered grid.
  ! ============================================================
  subroutine PC_ComputeDivergence(div, velx, vely)
    !
    real(rp), intent(out) :: div  (0:nx+1,0:ny+1)
    real(rp), intent(in)  :: velx(0:nx+1,0:ny+1)
    real(rp), intent(in)  :: vely(0:nx+1,0:ny+1)
    !
    integer :: i, j
    !
    div = 0.0_rp
    !
    do j = 1, ny
       do i = 1, nx
          div(i,j) = (velx(i,j) - velx(i-1,j)) / dx + &
                     (vely(i,j) - vely(i,j-1)) / dy
                     
       end do
    end do
    !
    return
  end subroutine PC_ComputeDivergence
  !
  ! ============================================================
  ! Build the right-hand side of the Poisson equation used in the
  ! pressure-correction step. The RHS is proportional to the
  ! divergence of the intermediate velocity field.
  ! ============================================================
  subroutine PC_BuildRHS(rhs_p, ustar, vstar, dt_rk3)
    !
    real(rp), intent(out) :: rhs_p (0:nx+1,0:ny+1)
    real(rp), intent(in)  :: ustar (0:nx+1,0:ny+1)
    real(rp), intent(in)  :: vstar (0:nx+1,0:ny+1)
    real(rp), intent(in)  :: dt_rk3
    !
    call PC_ComputeDivergence(rhs_p, ustar, vstar)
    !
    rhs_p(1:nx,1:ny) = rhs_p(1:nx,1:ny) / dt_rk3
    !
    return
  end subroutine PC_BuildRHS
  !
  ! ============================================================
  ! Correct the intermediate velocity field using the pressure
  ! solution of the current RK3 stage. The pressure gradient is
  ! applied on the staggered locations of u and v.
  ! ============================================================
  subroutine PC_CorrectVelocity(u, v, ustar, vstar, p_star, dt_rk3)
    !
    real(rp), intent(out) :: u     (0:nx+1,0:ny+1)
    real(rp), intent(out) :: v     (0:nx+1,0:ny+1)
    real(rp), intent(in)  :: ustar (0:nx+1,0:ny+1)
    real(rp), intent(in)  :: vstar (0:nx+1,0:ny+1)
    real(rp), intent(in)  :: p_star(0:nx+1,0:ny+1)
    real(rp), intent(in)  :: dt_rk3
    !
    integer :: i, j
    !
    u = ustar
    v = vstar
    !
    do j = 1, ny
       do i = 1, nx
          u(i,j) = ustar(i,j) - dt_rk3 * (p_star(i+1,j) - p_star(i,j)) / dx
       end do
    end do
    !
    do j = 1, ny
       do i = 1, nx
          v(i,j) = vstar(i,j) - dt_rk3 * (p_star(i,j+1) - p_star(i,j)) / dy
       end do
    end do
    !
    return
  end subroutine PC_CorrectVelocity
  !
end module mod_pressure_correction
