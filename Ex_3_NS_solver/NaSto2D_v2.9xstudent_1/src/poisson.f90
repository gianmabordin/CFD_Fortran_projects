module mod_poisson2d
  !
  use mod_types,  only: rp
  !
  use mod_common, only: nx, ny
  use mod_common, only: dx, dy
  !
  use mod_common, only: BC_DIRICHLET
  use mod_common, only: BC_NEUMANN
  use mod_common, only: BC_PERIODIC
  !
  use mod_common, only: poisson_maxit
  use mod_common, only: poisson_tol
  !
  use mod_common, only: poisson_iter
  use mod_common, only: poisson_flag
  use mod_common, only: poisson_res0
  use mod_common, only: poisson_res
  use mod_common, only: poisson_converged
  !
  use mod_common, only: bc_p_west
  use mod_common, only: bc_p_east
  use mod_common, only: bc_p_south
  use mod_common, only: bc_p_north
  !
  use mod_common, only: bc_p_west_val
  use mod_common, only: bc_p_east_val
  use mod_common, only: bc_p_south_val
  use mod_common, only: bc_p_north_val
  !
  use mod_common, only: poisson_active
  use mod_common, only: poisson_diagA
  use mod_common, only: poisson_r
  use mod_common, only: poisson_rhat
  use mod_common, only: poisson_p
  use mod_common, only: poisson_v
  use mod_common, only: poisson_s
  use mod_common, only: poisson_t
  use mod_common, only: poisson_phat
  use mod_common, only: poisson_shat
  use mod_common, only: poisson_Ax
  use mod_common, only: poisson_xt
  !
  implicit none
  !
  private
  !
  public :: POISSON2D_SolveBiCGStab
  public :: POISSON2D_CheckSetup
  !
contains
  !
  logical function IsPureNeumann()
    !
    logical :: west_free, east_free, south_free, north_free
    !
    west_free  = (bc_p_west  == BC_NEUMANN) .or. (bc_p_west  == BC_PERIODIC)
    east_free  = (bc_p_east  == BC_NEUMANN) .or. (bc_p_east  == BC_PERIODIC)
    south_free = (bc_p_south == BC_NEUMANN) .or. (bc_p_south == BC_PERIODIC)
    north_free = (bc_p_north == BC_NEUMANN) .or. (bc_p_north == BC_PERIODIC)
    !
    IsPureNeumann = west_free .and. east_free .and. south_free .and. north_free
    !
  end function
  !
  subroutine POISSON2D_CheckSetup()
    !
    if (nx < 1 .or. ny < 1) then
       stop 'STOP ERROR: invalid nx, ny in POISSON2D_CheckSetup'
    end if
    !
    if (.not. allocated(bc_p_west_val))  stop 'STOP ERROR: bc_p_west_val not allocated'
    if (.not. allocated(bc_p_east_val))  stop 'STOP ERROR: bc_p_east_val not allocated'
    if (.not. allocated(bc_p_south_val)) stop 'STOP ERROR: bc_p_south_val not allocated'
    if (.not. allocated(bc_p_north_val)) stop 'STOP ERROR: bc_p_north_val not allocated'
    !
    if (size(bc_p_west_val)  /= ny) stop 'STOP ERROR: wrong size for bc_p_west_val'
    if (size(bc_p_east_val)  /= ny) stop 'STOP ERROR: wrong size for bc_p_east_val'
    if (size(bc_p_south_val) /= nx) stop 'STOP ERROR: wrong size for bc_p_south_val'
    if (size(bc_p_north_val) /= nx) stop 'STOP ERROR: wrong size for bc_p_north_val'
    !
    if (bc_p_west  /= BC_DIRICHLET .and. bc_p_west  /= BC_NEUMANN .and. &
        bc_p_west  /= BC_PERIODIC) stop 'STOP ERROR: invalid bc_p_west'
    if (bc_p_east  /= BC_DIRICHLET .and. bc_p_east  /= BC_NEUMANN .and. &
        bc_p_east  /= BC_PERIODIC) stop 'STOP ERROR: invalid bc_p_east'
    if (bc_p_south /= BC_DIRICHLET .and. bc_p_south /= BC_NEUMANN .and. &
        bc_p_south /= BC_PERIODIC) stop 'STOP ERROR: invalid bc_p_south'
    if (bc_p_north /= BC_DIRICHLET .and. bc_p_north /= BC_NEUMANN .and. &
        bc_p_north /= BC_PERIODIC) stop 'STOP ERROR: invalid bc_p_north'
    !
    if (.not. allocated(poisson_active)) stop 'STOP ERROR: poisson_active not allocated'
    if (.not. allocated(poisson_diagA))  stop 'STOP ERROR: poisson_diagA not allocated'
    if (.not. allocated(poisson_r))      stop 'STOP ERROR: poisson_r not allocated'
    if (.not. allocated(poisson_rhat))   stop 'STOP ERROR: poisson_rhat not allocated'
    if (.not. allocated(poisson_p))      stop 'STOP ERROR: poisson_p not allocated'
    if (.not. allocated(poisson_v))      stop 'STOP ERROR: poisson_v not allocated'
    if (.not. allocated(poisson_s))      stop 'STOP ERROR: poisson_s not allocated'
    if (.not. allocated(poisson_t))      stop 'STOP ERROR: poisson_t not allocated'
    if (.not. allocated(poisson_phat))   stop 'STOP ERROR: poisson_phat not allocated'
    if (.not. allocated(poisson_shat))   stop 'STOP ERROR: poisson_shat not allocated'
    if (.not. allocated(poisson_Ax))     stop 'STOP ERROR: poisson_Ax not allocated'
    if (.not. allocated(poisson_xt))     stop 'STOP ERROR: poisson_xt not allocated'
    !
  end subroutine POISSON2D_CheckSetup
  !
  subroutine POISSON2D_SolveBiCGStab(p_star, rhs)
    !
    real(rp), intent(inout) :: p_star(0:nx+1,0:ny+1)
    real(rp), intent(in)    :: rhs   (0:nx+1,0:ny+1)
    !
    integer  :: k
    real(rp) :: rho_old, rho_new
    real(rp) :: alpha, omega, beta
    real(rp) :: den
    logical  :: pure_neumann
    !
    call POISSON2D_CheckSetup()
    !
    pure_neumann = IsPureNeumann()
    !
    call BuildActiveMask(poisson_active)
    !
    if (count(poisson_active(1:nx,1:ny)) <= 0) then
       stop 'STOP ERROR: no active dofs in POISSON2D_SolveBiCGStab'
    end if
    !
    call ApplyBC(p_star)
    call ApplyGaugeFixPointIfNeeded(p_star)
    !
    call BuildDiagPrecond(poisson_diagA, poisson_active)
    !
    call ApplyOperator(poisson_Ax, p_star, poisson_active)
    poisson_r = rhs - poisson_Ax
    call ZeroInactive(poisson_r, poisson_active)
    call ApplyGaugeFixPointToVectorIfNeeded(poisson_r)
    !
    poisson_rhat = poisson_r
    poisson_p    = 0.0_rp
    poisson_v    = 0.0_rp
    !
    rho_old = 1.0_rp
    alpha   = 1.0_rp
    omega   = 1.0_rp
    !
    poisson_iter      = 0
    poisson_flag      = 1
    poisson_converged = .false.
    poisson_res0      = sqrt(max(DotMasked(poisson_r, poisson_r, poisson_active), 0.0_rp))
    poisson_res       = poisson_res0
    !
    if (poisson_res0 <= poisson_tol) then
       poisson_flag      = 0
       poisson_converged = .true.
       return
    end if
    !
    do k = 1, poisson_maxit
       !
       rho_new = DotMasked(poisson_rhat, poisson_r, poisson_active)
       !
       if (abs(rho_new) <= tiny(1.0_rp)) then
          poisson_iter = k - 1
          poisson_flag = 2
          exit
       end if
       !
       if (k == 1) then
          poisson_p = poisson_r
       else
          if (abs(omega) <= tiny(1.0_rp)) then
             poisson_iter = k - 1
             poisson_flag = 5
             exit
          end if
          beta = (rho_new / rho_old) * (alpha / omega)
          poisson_p = poisson_r + beta * (poisson_p - omega * poisson_v)
       end if
       !
       call ZeroInactive(poisson_p, poisson_active)
       call ApplyGaugeFixPointToVectorIfNeeded(poisson_p)
       !
       call ApplyJacobiPrecond(poisson_phat, poisson_p, poisson_diagA, poisson_active)
       !
       call ApplyGaugeFixPointToVectorIfNeeded(poisson_phat)
       !
       call ApplyOperator(poisson_v, poisson_phat, poisson_active)
       den = DotMasked(poisson_rhat, poisson_v, poisson_active)
       !
       if (abs(den) <= tiny(1.0_rp)) then
          poisson_iter = k - 1
          poisson_flag = 3
          exit
       end if
       !
       alpha = rho_new / den
       !
       poisson_s = poisson_r - alpha * poisson_v
       call ZeroInactive(poisson_s, poisson_active)
       call ApplyGaugeFixPointToVectorIfNeeded(poisson_s)
       !
       poisson_res = sqrt(max(DotMasked(poisson_s, poisson_s, poisson_active), 0.0_rp))
       !
       if (poisson_res <= poisson_tol) then
          p_star = p_star + alpha * poisson_phat
          call ApplyBC(p_star)
          call ApplyGaugeFixPointIfNeeded(p_star)
          poisson_iter      = k
          poisson_flag      = 0
          poisson_converged = .true.
          exit
       end if
       !
       call ApplyJacobiPrecond(poisson_shat, poisson_s, poisson_diagA, poisson_active)
       !
       call ApplyGaugeFixPointToVectorIfNeeded(poisson_shat)
       !
       call ApplyOperator(poisson_t, poisson_shat, poisson_active)
       den = DotMasked(poisson_t, poisson_t, poisson_active)
       !
       if (abs(den) <= tiny(1.0_rp)) then
          poisson_iter = k - 1
          poisson_flag = 4
          exit
       end if
       !
       omega = DotMasked(poisson_t, poisson_s, poisson_active) / den
       !
       if (abs(omega) <= tiny(1.0_rp)) then
          poisson_iter = k - 1
          poisson_flag = 5
          exit
       end if
       !
       p_star = p_star + alpha * poisson_phat + omega * poisson_shat
       call ApplyBC(p_star)
       call ApplyGaugeFixPointIfNeeded(p_star)
       !
       poisson_r = poisson_s - omega * poisson_t
       call ZeroInactive(poisson_r, poisson_active)
       call ApplyGaugeFixPointToVectorIfNeeded(poisson_r)
       !
       poisson_res  = sqrt(max(DotMasked(poisson_r, poisson_r, poisson_active), 0.0_rp))
       poisson_iter = k
       !
       if (poisson_res <= poisson_tol) then
          poisson_flag      = 0
          poisson_converged = .true.
          exit
       end if
       !
       rho_old = rho_new
       !
    end do
    !
    call ApplyBC(p_star)
    call ApplyGaugeFixPointIfNeeded(p_star)
    !
  end subroutine POISSON2D_SolveBiCGStab
  !
  subroutine ApplyOperator(y, x, is_active)
    !
    real(rp), intent(out) :: y(0:nx+1,0:ny+1)
    real(rp), intent(in)  :: x(0:nx+1,0:ny+1)
    logical, intent(in)   :: is_active(0:nx+1,0:ny+1)
    !
    integer :: i, j
    real(rp) :: idx2, idy2
    !
    poisson_xt = x
    !
    call ApplyBC(poisson_xt)
    call ApplyGaugeFixPointIfNeeded(poisson_xt)
    !
    idx2 = 1.0_rp / (dx*dx)
    idy2 = 1.0_rp / (dy*dy)
    !
    y = 0.0_rp
    !
    do j = 1, ny
       do i = 1, nx
          if (is_active(i,j)) then
             y(i,j) = ( &
                  (poisson_xt(i+1,j) - 2.0_rp*poisson_xt(i,j) + poisson_xt(i-1,j)) * idx2 + &
                  (poisson_xt(i,j+1) - 2.0_rp*poisson_xt(i,j) + poisson_xt(i,j-1)) * idy2 )
          end if
       end do
    end do
    !
    call ZeroInactive(y, is_active)
    call ApplyGaugeFixPointToVectorIfNeeded(y)
    !
  end subroutine ApplyOperator
  !
  subroutine ApplyBC(p_star)
    !
    real(rp), intent(inout) :: p_star(0:nx+1,0:ny+1)
    !
    integer :: i, j
    !
    do j = 1, ny
       select case (bc_p_west)
       case (BC_DIRICHLET)
          p_star(0,j) = 2.0_rp * bc_p_west_val(j) - p_star(1,j)
       case (BC_NEUMANN)
          p_star(0,j) = p_star(1,j) - dx * bc_p_west_val(j)
       case (BC_PERIODIC)
          p_star(0,j) = p_star(nx,j)
       end select
       !
       select case (bc_p_east)
       case (BC_DIRICHLET)
          p_star(nx+1,j) = 2.0_rp * bc_p_east_val(j) - p_star(nx,j)
       case (BC_NEUMANN)
          p_star(nx+1,j) = p_star(nx,j) + dx * bc_p_east_val(j)
       case (BC_PERIODIC)
          p_star(nx+1,j) = p_star(1,j)
       end select
    end do
    !
    do i = 1, nx
       select case (bc_p_south)
       case (BC_DIRICHLET)
          p_star(i,0) = 2.0_rp * bc_p_south_val(i) - p_star(i,1)
       case (BC_NEUMANN)
          p_star(i,0) = p_star(i,1) - dy * bc_p_south_val(i)
       case (BC_PERIODIC)
          p_star(i,0) = p_star(i,ny)
       end select
       !
       select case (bc_p_north)
       case (BC_DIRICHLET)
          p_star(i,ny+1) = 2.0_rp * bc_p_north_val(i) - p_star(i,ny)
       case (BC_NEUMANN)
          p_star(i,ny+1) = p_star(i,ny) + dy * bc_p_north_val(i)
       case (BC_PERIODIC)
          p_star(i,ny+1) = p_star(i,1)
       end select
    end do
    !
    p_star(0,0)       = 0.5_rp * (p_star(1,0)      + p_star(0,1))
    p_star(0,ny+1)    = 0.5_rp * (p_star(1,ny+1)   + p_star(0,ny))
    p_star(nx+1,0)    = 0.5_rp * (p_star(nx,0)     + p_star(nx+1,1))
    p_star(nx+1,ny+1) = 0.5_rp * (p_star(nx,ny+1)  + p_star(nx+1,ny))
    !
  end subroutine ApplyBC
  !
  subroutine BuildActiveMask(is_active)
    !
    logical, intent(out) :: is_active(0:nx+1,0:ny+1)
    !
    integer :: i, j
    !
    is_active = .false.
    !
    do j = 1, ny
       do i = 1, nx
          is_active(i,j) = .true.
       end do
    end do
    !
    if (IsPureNeumann()) then
       is_active(1,1) = .false.
    end if
    !
  end subroutine BuildActiveMask
  !
  subroutine ApplyGaugeFixPointIfNeeded(a)
    !
    real(rp), intent(inout) :: a(0:nx+1,0:ny+1)
    !
    if (IsPureNeumann()) then
       a(1,1) = 0.0_rp
    end if
    !
  end subroutine ApplyGaugeFixPointIfNeeded
  !
  subroutine ApplyGaugeFixPointToVectorIfNeeded(a)
    !
    real(rp), intent(inout) :: a(0:nx+1,0:ny+1)
    !
    if (IsPureNeumann()) then
       a(1,1) = 0.0_rp
    end if
    !
  end subroutine ApplyGaugeFixPointToVectorIfNeeded
  !
  subroutine BuildDiagPrecond(diagA, is_active)
    !
    real(rp), intent(out) :: diagA(0:nx+1,0:ny+1)
    logical, intent(in)   :: is_active(0:nx+1,0:ny+1)
    !
    integer :: i, j
    real(rp) :: idx2, idy2
    real(rp) :: val
    !
    idx2 = 1.0_rp / (dx*dx)
    idy2 = 1.0_rp / (dy*dy)
    val  = + ( 2.0_rp * idx2 + 2.0_rp * idy2 )
    !
    diagA = 1.0_rp
    !
    do j = 1, ny
       do i = 1, nx
          if (is_active(i,j)) then
             diagA(i,j) = val
          end if
       end do
    end do
    !
  end subroutine BuildDiagPrecond
  !
  subroutine ApplyJacobiPrecond(z, r, diagA, is_active)
    !
    real(rp), intent(out) :: z(0:nx+1,0:ny+1)
    real(rp), intent(in)  :: r(0:nx+1,0:ny+1)
    real(rp), intent(in)  :: diagA(0:nx+1,0:ny+1)
    logical, intent(in)   :: is_active(0:nx+1,0:ny+1)
    !
    integer :: i, j
    !
    z = 0.0_rp
    !
    do j = 1, ny
       do i = 1, nx
          if (is_active(i,j)) then
             z(i,j) = r(i,j) / diagA(i,j)
          end if
       end do
    end do
    !
  end subroutine ApplyJacobiPrecond
  !
  subroutine ZeroInactive(a, is_active)
    !
    real(rp), intent(inout) :: a(0:nx+1,0:ny+1)
    logical, intent(in)     :: is_active(0:nx+1,0:ny+1)
    !
    where (.not. is_active) a = 0.0_rp
    !
  end subroutine ZeroInactive
  !
  real(rp) function DotMasked(a, b, is_active)
    !
    real(rp), intent(in) :: a(0:nx+1,0:ny+1)
    real(rp), intent(in) :: b(0:nx+1,0:ny+1)
    logical, intent(in)  :: is_active(0:nx+1,0:ny+1)
    !
    DotMasked = sum(a*b, mask=is_active)
    !
  end function DotMasked
  !
end module mod_poisson2d
