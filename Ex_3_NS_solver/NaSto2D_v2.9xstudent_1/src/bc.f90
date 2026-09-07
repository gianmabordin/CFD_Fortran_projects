module mod_bc
  !
  use mod_types,  only: rp
  use mod_common, only: nx, ny
  use mod_common, only: dx, dy
  !
  use mod_common, only: BC_DIRICHLET
  use mod_common, only: BC_NEUMANN
  use mod_common, only: BC_PERIODIC
  !
  use mod_common, only: bc_u_west,  bc_u_east,  bc_u_south,  bc_u_north
  use mod_common, only: bc_v_west,  bc_v_east,  bc_v_south,  bc_v_north
  use mod_common, only: bc_p_west,  bc_p_east,  bc_p_south,  bc_p_north
  !
  use mod_common, only: bc_u_west_val,  bc_u_east_val
  use mod_common, only: bc_u_south_val, bc_u_north_val
  use mod_common, only: bc_v_west_val,  bc_v_east_val
  use mod_common, only: bc_v_south_val, bc_v_north_val
  use mod_common, only: bc_p_west_val,  bc_p_east_val
  use mod_common, only: bc_p_south_val, bc_p_north_val
  !
  implicit none
  !
  private
  !
  public :: BC_ApplyVelocity
  public :: BC_ApplyPressure
  public :: BC_ApplyVelocityAfterProjection
  !
contains
  !
  ! ============================================================
  ! Apply boundary conditions to both velocity components.
  ! This is the standard velocity update used before projection,
  ! when all velocity boundary values and ghost cells are refreshed.
  ! ============================================================
  subroutine BC_ApplyVelocity(u, v)
    !
    real(rp), intent(inout) :: u(0:nx+1,0:ny+1)
    real(rp), intent(inout) :: v(0:nx+1,0:ny+1)
    !
    call BC_ApplyU(u)
    call BC_ApplyV(v)
    !
  end subroutine BC_ApplyVelocity
  !
  ! ============================================================
  ! Apply velocity boundary conditions after the pressure projection.
  ! Dirichlet conditions are re-imposed on normal velocity components,
  ! while Neumann normal boundaries are left unchanged so that the
  ! projected field is not modified unnecessarily.
  ! ============================================================
  subroutine BC_ApplyVelocityAfterProjection(u, v)
    !
    real(rp), intent(inout) :: u(0:nx+1,0:ny+1)
    real(rp), intent(inout) :: v(0:nx+1,0:ny+1)
    !
    integer :: i, j
    !
    ! ----------------------------------------------------------
    ! u on west/east boundaries: normal velocity component.
    ! Re-impose Dirichlet values, keep Neumann faces unchanged,
    ! and update periodic boundaries when needed.
    ! ----------------------------------------------------------
    do j = 1, ny
       select case (bc_u_west)
       case (BC_DIRICHLET)
          u(0,j) = bc_u_west_val(j)
       case (BC_PERIODIC)
          u(0,j) = u(nx,j)
       case default
          ! BC_NEUMANN: leave u(0,j) unchanged
       end select
       !
       select case (bc_u_east)
       case (BC_DIRICHLET)
          u(nx,j) = bc_u_east_val(j)
       case (BC_PERIODIC)
          u(nx,j) = u(0,j)
       case default
          ! BC_NEUMANN: leave u(nx,j) unchanged
       end select
    end do
    !
    ! ----------------------------------------------------------
    ! Update the outer ghost layer on the east side for u.
    ! This keeps storage consistent with the selected boundary type.
    ! ----------------------------------------------------------
    do j = 1, ny
       select case (bc_u_east)
       case (BC_DIRICHLET)
          u(nx+1,j) = 2.0_rp * u(nx,j) - u(nx-1,j)
       case (BC_NEUMANN)
          u(nx+1,j) = u(nx,j) + dx * bc_u_east_val(j)
       case (BC_PERIODIC)
          u(nx+1,j) = u(1,j)
       end select
    end do
    !
    ! ----------------------------------------------------------
    ! u on south/north boundaries: tangential velocity component.
    ! Here ghost values are updated in the standard way.
    ! ----------------------------------------------------------
    do i = 1, nx-1
       select case (bc_u_south)
       case (BC_DIRICHLET)
          u(i,0) = 2.0_rp * bc_u_south_val(i) - u(i,1)
       case (BC_NEUMANN)
          u(i,0) = u(i,1) - dy * bc_u_south_val(i)
       case (BC_PERIODIC)
          u(i,0) = u(i,ny)
       end select
       !
       select case (bc_u_north)
       case (BC_DIRICHLET)
          u(i,ny+1) = 2.0_rp * bc_u_north_val(i) - u(i,ny)
       case (BC_NEUMANN)
          u(i,ny+1) = u(i,ny) + dy * bc_u_north_val(i)
       case (BC_PERIODIC)
          u(i,ny+1) = u(i,1)
       end select
    end do
    !
    ! Complete corner and outer-corner ghost values by copying
    ! from the closest available interior or boundary-adjacent value.
    u(0,0)       = u(0,1)
    u(0,ny+1)    = u(0,ny)
    u(nx,0)      = u(nx,1)
    u(nx,ny+1)   = u(nx,ny)
    u(nx+1,0)    = u(nx+1,1)
    u(nx+1,ny+1) = u(nx+1,ny)
    !
    ! ----------------------------------------------------------
    ! v on west/east boundaries: tangential velocity component.
    ! Ghost values are updated in the standard way.
    ! ----------------------------------------------------------
    do j = 1, ny-1
       select case (bc_v_west)
       case (BC_DIRICHLET)
          v(0,j) = 2.0_rp * bc_v_west_val(j) - v(1,j)
       case (BC_NEUMANN)
          v(0,j) = v(1,j) - dx * bc_v_west_val(j)
       case (BC_PERIODIC)
          v(0,j) = v(nx,j)
       end select
       !
       select case (bc_v_east)
       case (BC_DIRICHLET)
          v(nx+1,j) = 2.0_rp * bc_v_east_val(j) - v(nx,j)
       case (BC_NEUMANN)
          v(nx+1,j) = v(nx,j) + dx * bc_v_east_val(j)
       case (BC_PERIODIC)
          v(nx+1,j) = v(1,j)
       end select
    end do
    !
    ! ----------------------------------------------------------
    ! v on south/north boundaries: normal velocity component.
    ! Re-impose Dirichlet values, keep Neumann faces unchanged,
    ! and update periodic boundaries when needed.
    ! ----------------------------------------------------------
    do i = 1, nx
       select case (bc_v_south)
       case (BC_DIRICHLET)
          v(i,0) = bc_v_south_val(i)
       case (BC_PERIODIC)
          v(i,0) = v(i,ny)
       case default
          ! BC_NEUMANN: leave v(i,0) unchanged
       end select
       !
       select case (bc_v_north)
       case (BC_DIRICHLET)
          v(i,ny) = bc_v_north_val(i)
       case (BC_PERIODIC)
          v(i,ny) = v(i,0)
       case default
          ! BC_NEUMANN: leave v(i,ny) unchanged
       end select
    end do
    !
    ! ----------------------------------------------------------
    ! Update the outer ghost layer on the north side for v.
    ! This keeps storage consistent with the selected boundary type.
    ! ----------------------------------------------------------
    do i = 1, nx
       select case (bc_v_north)
       case (BC_DIRICHLET)
          v(i,ny+1) = 2.0_rp * v(i,ny) - v(i,ny-1)
       case (BC_NEUMANN)
          v(i,ny+1) = v(i,ny) + dy * bc_v_north_val(i)
       case (BC_PERIODIC)
          v(i,ny+1) = v(i,1)
       end select
    end do
    !
    ! Complete corner and outer-corner ghost values by copying
    ! from the closest available interior or boundary-adjacent value.
    v(0,0)       = v(1,0)
    v(nx+1,0)    = v(nx,0)
    v(0,ny)      = v(1,ny)
    v(nx+1,ny)   = v(nx,ny)
    v(0,ny+1)    = v(1,ny+1)
    v(nx+1,ny+1) = v(nx,ny+1)
    !
    return
  end subroutine BC_ApplyVelocityAfterProjection
  !
  ! ============================================================
  ! Apply boundary conditions to the horizontal velocity field u.
  ! This routine updates physical boundary values and ghost cells
  ! according to the selected Dirichlet, Neumann, or periodic type.
  ! ============================================================
  subroutine BC_ApplyU(u)
    !
    real(rp), intent(inout) :: u(0:nx+1,0:ny+1)
    !
    integer :: i, j
    !
    do j = 1, ny
       select case (bc_u_west)
       case (BC_DIRICHLET)
          u(0,j) = bc_u_west_val(j)
       case (BC_NEUMANN)
          u(0,j) = u(1,j) - dx * bc_u_west_val(j)
       case (BC_PERIODIC)
          u(0,j) = u(nx,j)
       end select
       !
       select case (bc_u_east)
       case (BC_DIRICHLET)
          u(nx,j) = bc_u_east_val(j)
       case (BC_NEUMANN)
          u(nx,j) = u(nx-1,j) + dx * bc_u_east_val(j)
       case (BC_PERIODIC)
          u(nx,j) = u(0,j)
       end select
    end do
    !
    do j = 1, ny
       select case (bc_u_east)
       case (BC_DIRICHLET)
          u(nx+1,j) = 2.0_rp * u(nx,j) - u(nx-1,j)
       case (BC_NEUMANN)
          u(nx+1,j) = u(nx,j) + dx * bc_u_east_val(j)
       case (BC_PERIODIC)
          u(nx+1,j) = u(1,j)
       end select
    end do
    !
    do i = 1, nx-1
       select case (bc_u_south)
       case (BC_DIRICHLET)
          u(i,0) = 2.0_rp * bc_u_south_val(i) - u(i,1)
       case (BC_NEUMANN)
          u(i,0) = u(i,1) - dy * bc_u_south_val(i)
       case (BC_PERIODIC)
          u(i,0) = u(i,ny)
       end select
       !
       select case (bc_u_north)
       case (BC_DIRICHLET)
          u(i,ny+1) = 2.0_rp * bc_u_north_val(i) - u(i,ny)
       case (BC_NEUMANN)
          u(i,ny+1) = u(i,ny) + dy * bc_u_north_val(i)
       case (BC_PERIODIC)
          u(i,ny+1) = u(i,1)
       end select
    end do
    !
    ! Complete corner and outer-corner ghost values.
    u(0,0)       = u(0,1)
    u(0,ny+1)    = u(0,ny)
    u(nx,0)      = u(nx,1)
    u(nx,ny+1)   = u(nx,ny)
    u(nx+1,0)    = u(nx+1,1)
    u(nx+1,ny+1) = u(nx+1,ny)
    !
    return
  end subroutine BC_ApplyU
  !
  ! ============================================================
  ! Apply boundary conditions to the vertical velocity field v.
  ! This routine updates physical boundary values and ghost cells
  ! according to the selected Dirichlet, Neumann, or periodic type.
  ! ============================================================
  subroutine BC_ApplyV(v)
    !
    real(rp), intent(inout) :: v(0:nx+1,0:ny+1)
    !
    integer :: i, j
    !
    do j = 1, ny-1
       select case (bc_v_west)
       case (BC_DIRICHLET)
          v(0,j) = 2.0_rp * bc_v_west_val(j) - v(1,j)
       case (BC_NEUMANN)
          v(0,j) = v(1,j) - dx * bc_v_west_val(j)
       case (BC_PERIODIC)
          v(0,j) = v(nx,j)
       end select
       !
       select case (bc_v_east)
       case (BC_DIRICHLET)
          v(nx+1,j) = 2.0_rp * bc_v_east_val(j) - v(nx,j)
       case (BC_NEUMANN)
          v(nx+1,j) = v(nx,j) + dx * bc_v_east_val(j)
       case (BC_PERIODIC)
          v(nx+1,j) = v(1,j)
       end select
    end do
    !
    do i = 1, nx
       select case (bc_v_south)
       case (BC_DIRICHLET)
          v(i,0) = bc_v_south_val(i)
       case (BC_NEUMANN)
          v(i,0) = v(i,1) - dy * bc_v_south_val(i)
       case (BC_PERIODIC)
          v(i,0) = v(i,ny)
       end select
       !
       select case (bc_v_north)
       case (BC_DIRICHLET)
          v(i,ny) = bc_v_north_val(i)
       case (BC_NEUMANN)
          v(i,ny) = v(i,ny-1) + dy * bc_v_north_val(i)
       case (BC_PERIODIC)
          v(i,ny) = v(i,0)
       end select
    end do
    !
    do i = 1, nx
       select case (bc_v_north)
       case (BC_DIRICHLET)
          v(i,ny+1) = 2.0_rp * v(i,ny) - v(i,ny-1)
       case (BC_NEUMANN)
          v(i,ny+1) = v(i,ny) + dy * bc_v_north_val(i)
       case (BC_PERIODIC)
          v(i,ny+1) = v(i,1)
       end select
    end do
    !
    ! Complete corner and outer-corner ghost values.
    v(0,0)       = v(1,0)
    v(nx+1,0)    = v(nx,0)
    v(0,ny)      = v(1,ny)
    v(nx+1,ny)   = v(nx,ny)
    v(0,ny+1)    = v(1,ny+1)
    v(nx+1,ny+1) = v(nx,ny+1)
    !
    return
  end subroutine BC_ApplyV
  !
  ! ============================================================
  ! Apply boundary conditions to the pressure field.
  ! Pressure is stored at cell centers, so boundary conditions are
  ! enforced through ghost-cell updates around the domain.
  ! ============================================================
  subroutine BC_ApplyPressure(p)
    !
    real(rp), intent(inout) :: p(0:nx+1,0:ny+1)
    !
    integer :: i, j
    !
    do j = 1, ny
       select case (bc_p_west)
       case (BC_DIRICHLET)
          p(0,j) = 2.0_rp * bc_p_west_val(j) - p(1,j)
       case (BC_NEUMANN)
          p(0,j) = p(1,j) - dx * bc_p_west_val(j)
       case (BC_PERIODIC)
          p(0,j) = p(nx,j)
       end select
       !
       select case (bc_p_east)
       case (BC_DIRICHLET)
          p(nx+1,j) = 2.0_rp * bc_p_east_val(j) - p(nx,j)
       case (BC_NEUMANN)
          p(nx+1,j) = p(nx,j) + dx * bc_p_east_val(j)
       case (BC_PERIODIC)
          p(nx+1,j) = p(1,j)
       end select
    end do
    !
    do i = 1, nx
       select case (bc_p_south)
       case (BC_DIRICHLET)
          p(i,0) = 2.0_rp * bc_p_south_val(i) - p(i,1)
       case (BC_NEUMANN)
          p(i,0) = p(i,1) - dy * bc_p_south_val(i)
       case (BC_PERIODIC)
          p(i,0) = p(i,ny)
       end select
       !
       select case (bc_p_north)
       case (BC_DIRICHLET)
          p(i,ny+1) = 2.0_rp * bc_p_north_val(i) - p(i,ny)
       case (BC_NEUMANN)
          p(i,ny+1) = p(i,ny) + dy * bc_p_north_val(i)
       case (BC_PERIODIC)
          p(i,ny+1) = p(i,1)
       end select
    end do
    !
    ! Complete corner ghost values.
    p(0,0)         = p(0,1)
    p(0,ny+1)      = p(0,ny)
    p(nx+1,0)      = p(nx+1,1)
    p(nx+1,ny+1)   = p(nx+1,ny)
    !
    return
  end subroutine BC_ApplyPressure
  !
end module mod_bc
