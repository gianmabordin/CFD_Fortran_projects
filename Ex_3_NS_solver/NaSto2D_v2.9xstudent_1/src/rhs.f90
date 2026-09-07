module mod_rhs
  !
  use mod_types,  only: rp
  !
  use mod_common, only: nx, ny
  use mod_common, only: dx, dy
  use mod_common, only: nu, rho
  use mod_common, only: p
  !
  implicit none
  !
  private
  !
  public :: RHS_ComputeMomentumPressure
  public :: RHS_ComputeMomentumConvDiff
  !
contains
  !
  ! ============================================================
  ! Compute the pressure contribution to the momentum right-hand side.
  ! Pressure gradients are evaluated at the staggered locations of
  ! the velocity components and returned in array form.
  !
  ! u(i,j) is stored on vertical cell faces.
  ! v(i,j) is stored on horizontal cell faces.
  !
  ! The output arrays are initialized everywhere, but only the
  ! physical staggered locations are used by the solver.
  ! ============================================================
  subroutine RHS_ComputeMomentumPressure(rhs_u_out, rhs_v_out)
    !
    real(rp), intent(out) :: rhs_u_out(0:nx+1,0:ny+1)
    real(rp), intent(out) :: rhs_v_out(0:nx+1,0:ny+1)
    !
    integer :: i, j
    !
    rhs_u_out = 0.0_rp
    rhs_v_out = 0.0_rp
    !
    do j = 1, ny
       do i = 1, nx
          rhs_u_out(i,j) = -(p(i+1,j) - p(i,j)) / (dx * rho)
       end do
    end do
    !
    do j = 1, ny
       do i = 1, nx
          rhs_v_out(i,j) = -(p(i,j+1) - p(i,j)) / (dy * rho)
       end do
    end do
    !
    return
  end subroutine RHS_ComputeMomentumPressure
  !
  ! ============================================================
  ! Compute the convective and diffusive contributions to the
  ! momentum right-hand side on the staggered grid.
  !
  ! The returned terms are:
  !   rhs_u_out = -conv_u + diff_u
  !   rhs_v_out = -conv_v + diff_v
  !
  ! Convective fluxes are evaluated in conservative form, while
  ! viscous terms are approximated with second-order central
  ! differences.
  ! ============================================================
  subroutine RHS_ComputeMomentumConvDiff(rhs_u_out, rhs_v_out, u, v)
    !
    real(rp), intent(out) :: rhs_u_out(0:nx+1,0:ny+1)
    real(rp), intent(out) :: rhs_v_out(0:nx+1,0:ny+1)
    real(rp), intent(in)  :: u        (0:nx+1,0:ny+1)
    real(rp), intent(in)  :: v        (0:nx+1,0:ny+1)
    !
    integer :: i, j
    !
    real(rp) :: duu_dx, duv_dy
    real(rp) :: dvu_dx, dvv_dy
    real(rp) :: d2u_dx2, d2u_dy2
    real(rp) :: d2v_dx2, d2v_dy2
    !
    real(rp) :: conv_x, conv_y
    real(rp) :: diff_x, diff_y
    !
    rhs_u_out = 0.0_rp
    rhs_v_out = 0.0_rp
    !
    ! ----------------------------------------------------------
    ! Compute the u-momentum contribution on u-staggered faces.
    ! The convective term includes the x-flux of u^2 and the
    ! y-flux of uv, while the diffusive term is the Laplacian of u.
    ! ----------------------------------------------------------
    do j = 1, ny
       do i = 1, nx
          !
          duu_dx = ( ((u(i+1,j) + u(i,j)) * 0.5_rp)**2 - &
                     ((u(i,j) + u(i-1,j)) * 0.5_rp)**2 ) / dx
                  
          !
          duv_dy = ( ((u(i,j)   + u(i,j+1)) * 0.5_rp) * ((v(i,j)   + v(i+1,j))   * 0.5_rp) - &
                     ((u(i,j)   + u(i,j-1)) * 0.5_rp) * ((v(i,j-1) + v(i+1,j-1)) * 0.5_rp) ) / dy
          !
          d2u_dx2 = (u(i+1,j) - 2.0_rp * u(i,j) + u(i-1,j)) / (dx * dx)
          d2u_dy2 = (u(i,j+1) - 2.0_rp * u(i,j) + u(i,j-1)) / (dy * dy)
          !
          conv_x = duu_dx + duv_dy
          diff_x = nu * (d2u_dx2 + d2u_dy2)
          !
          rhs_u_out(i,j) = -conv_x + diff_x
          !
       end do
    end do
    !
    ! ----------------------------------------------------------
    ! Compute the v-momentum contribution on v-staggered faces.
    ! The convective term includes the x-flux of vu and the
    ! y-flux of v^2, while the diffusive term is the Laplacian of v.
    ! ----------------------------------------------------------
    do j = 1, ny
       do i = 1, nx
          !
          dvu_dx = ( ((u(i,j)   + u(i,j+1))   * 0.5_rp) * ((v(i,j)   + v(i+1,j)) * 0.5_rp) - &
                     ((u(i-1,j) + u(i-1,j+1)) * 0.5_rp) * ((v(i-1,j) + v(i,j))   * 0.5_rp) ) / dx
                   
          !
          dvv_dy = ( ((v(i,j+1) + v(i,j)) * 0.5_rp)**2 - &
                     ((v(i,j) + v(i,j-1)) * 0.5_rp)**2 ) / dy
                   
          !
          d2v_dx2 = (v(i+1,j) - 2.0_rp * v(i,j) + v(i-1,j)) / (dx * dx)
          d2v_dy2 = (v(i,j+1) - 2.0_rp * v(i,j) + v(i,j-1)) / (dy * dy)
          !
          conv_y = dvu_dx + dvv_dy
          diff_y = nu * (d2v_dx2 + d2v_dy2)
          !
          rhs_v_out(i,j) = -conv_y + diff_y
          !
       end do
    end do
    !
    return
  end subroutine RHS_ComputeMomentumConvDiff
  !
end module mod_rhs
