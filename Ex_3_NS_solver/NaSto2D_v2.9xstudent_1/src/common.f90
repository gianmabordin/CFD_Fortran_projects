module mod_common
  !
  use mod_types, only: rp
  !
  implicit none
  !
  private
  !
  ! ============================================================
  ! Integer flags used to identify the boundary-condition type
  ! applied to a field on each side of the domain.
  ! ============================================================
  integer, parameter, public :: BC_DIRICHLET = 1
  integer, parameter, public :: BC_NEUMANN   = 2
  integer, parameter, public :: BC_PERIODIC  = 3
  !
  ! ============================================================
  ! Global grid and domain parameters.
  ! nx and ny are the number of physical cells, while lx and ly
  ! are the domain lengths in the x and y directions.
  ! dx and dy are the corresponding uniform grid spacings.
  ! ============================================================
  integer, public :: nx = 0
  integer, public :: ny = 0
  !
  real(rp), public :: lx = 0.0_rp
  real(rp), public :: ly = 0.0_rp
  real(rp), public :: dx = 0.0_rp
  real(rp), public :: dy = 0.0_rp
  !
  ! Cell-center coordinates for physical cells only.
  real(rp), allocatable, public :: xc(:)
  real(rp), allocatable, public :: yc(:)
  !
  ! Optional staggered coordinates for physical lines and faces.
  ! These are useful for post-processing and interpretation of the
  ! staggered-grid arrangement used by the solver.
  real(rp), allocatable, public :: xu(:)
  real(rp), allocatable, public :: yu(:)
  real(rp), allocatable, public :: xv(:)
  real(rp), allocatable, public :: yv(:)
  !
  ! ============================================================
  ! Main flow fields including one layer of ghost cells.
  !
  ! Physical cells are stored in:
  !   i = 1:nx, j = 1:ny
  !
  ! Ghost layers are stored at:
  !   i = 0 and nx+1, j = 0 and ny+1
  ! ============================================================
  real(rp), allocatable, public :: u(:,:)
  real(rp), allocatable, public :: v(:,:)
  real(rp), allocatable, public :: p(:,:)
  !
  ! Pressure correction of the current projection stage.
  real(rp), allocatable, public :: p_star(:,:)
  !
  ! Predictor fields and Poisson right-hand side used during the
  ! fractional-step time advancement.
  real(rp), allocatable, public :: ustar(:,:)
  real(rp), allocatable, public :: vstar(:,:)
  real(rp), allocatable, public :: rhs_p(:,:)
  !
  ! RK3 stage storage fields.
  real(rp), allocatable, public :: urk1(:,:)
  real(rp), allocatable, public :: vrk1(:,:)
  real(rp), allocatable, public :: urk2(:,:)
  real(rp), allocatable, public :: vrk2(:,:)
  !
  ! ============================================================
  ! Work arrays for the momentum right-hand side.
  ! These arrays store convective, diffusive, and combined RHS
  ! contributions, as well as the previous-stage history needed
  ! by the RK3 scheme.
  ! ============================================================
  real(rp), allocatable, public :: conv_u(:,:)
  real(rp), allocatable, public :: conv_v(:,:)
  real(rp), allocatable, public :: diff_u(:,:)
  real(rp), allocatable, public :: diff_v(:,:)
  real(rp), allocatable, public :: rhs_u(:,:)
  real(rp), allocatable, public :: rhs_v(:,:)
  real(rp), allocatable, public :: rhss_u(:,:)
  real(rp), allocatable, public :: rhss_v(:,:)
  !
  ! ============================================================
  ! Work arrays, controls, and diagnostics for the Poisson solver.
  ! These variables are used internally by the BiCGStab algorithm
  ! employed in the pressure-correction step.
  ! ============================================================
  integer, public :: poisson_maxit = 5000
  real(rp), public :: poisson_tol  = 1.0e-10_rp
  !
  integer, public :: poisson_iter = 0
  integer, public :: poisson_flag = 0
  real(rp), public :: poisson_res0 = 0.0_rp
  real(rp), public :: poisson_res  = 0.0_rp
  logical, public :: poisson_converged = .false.
  !
  logical,  allocatable, public :: poisson_active(:,:)
  real(rp), allocatable, public :: poisson_diagA(:,:)
  real(rp), allocatable, public :: poisson_r(:,:)
  real(rp), allocatable, public :: poisson_rhat(:,:)
  real(rp), allocatable, public :: poisson_p(:,:)
  real(rp), allocatable, public :: poisson_v(:,:)
  real(rp), allocatable, public :: poisson_s(:,:)
  real(rp), allocatable, public :: poisson_t(:,:)
  real(rp), allocatable, public :: poisson_phat(:,:)
  real(rp), allocatable, public :: poisson_shat(:,:)
  real(rp), allocatable, public :: poisson_Ax(:,:)
  real(rp), allocatable, public :: poisson_xt(:,:)
  !
  ! ============================================================
  ! Physical and numerical parameters of the simulation.
  ! rho is the density, nu is the kinematic viscosity, dt is the
  ! time step, and dt_rk3 is the effective time increment of the
  ! current RK3 stage.
  ! ============================================================
  real(rp), public :: rho = 1.0_rp
  real(rp), public :: nu  = 1.0e-3_rp
  real(rp), public :: dt  = 1.0e-3_rp
  real(rp), public :: dt_rk3 = 0.0_rp
  !
  real(rp), public :: time     = 0.0_rp
  real(rp), public :: time_end = 0.0_rp
  !
  integer, public :: istep  = 0
  integer, public :: nsteps = 0
  integer, public :: it_rk3 = 0
  !
  ! ============================================================
  ! Boundary-condition type assigned to each side of the domain
  ! for the velocity and pressure fields.
  ! ============================================================
  integer, public :: bc_u_west  = BC_DIRICHLET
  integer, public :: bc_u_east  = BC_DIRICHLET
  integer, public :: bc_u_south = BC_DIRICHLET
  integer, public :: bc_u_north = BC_DIRICHLET
  !
  integer, public :: bc_v_west  = BC_DIRICHLET
  integer, public :: bc_v_east  = BC_DIRICHLET
  integer, public :: bc_v_south = BC_DIRICHLET
  integer, public :: bc_v_north = BC_DIRICHLET
  !
  integer, public :: bc_p_west  = BC_NEUMANN
  integer, public :: bc_p_east  = BC_NEUMANN
  integer, public :: bc_p_south = BC_NEUMANN
  integer, public :: bc_p_north = BC_NEUMANN
  !
  ! ============================================================
  ! Boundary-condition values associated with the NS fields.
  ! These values are stored only on the physical boundary lines:
  !
  !   west/east   -> indexed with j = 1:ny
  !   south/north -> indexed with i = 1:nx
  ! ============================================================
  real(rp), allocatable, public :: bc_u_west_val(:)
  real(rp), allocatable, public :: bc_u_east_val(:)
  real(rp), allocatable, public :: bc_u_south_val(:)
  real(rp), allocatable, public :: bc_u_north_val(:)
  !
  real(rp), allocatable, public :: bc_v_west_val(:)
  real(rp), allocatable, public :: bc_v_east_val(:)
  real(rp), allocatable, public :: bc_v_south_val(:)
  real(rp), allocatable, public :: bc_v_north_val(:)
  !
  real(rp), allocatable, public :: bc_p_west_val(:)
  real(rp), allocatable, public :: bc_p_east_val(:)
  real(rp), allocatable, public :: bc_p_south_val(:)
  real(rp), allocatable, public :: bc_p_north_val(:)
  !
  public :: CMM_InitGrid
  public :: CMM_AllocFields
  public :: CMM_FreeFields
  !
contains
  !
  ! ============================================================
  ! Initialize the grid and the coordinate arrays from the user
  ! input. This routine sets the domain size, computes dx and dy,
  ! and builds the cell-centered and staggered coordinates.
  ! ============================================================
  subroutine CMM_InitGrid(nx_in, ny_in, lx_in, ly_in)
    !
    integer, intent(in)  :: nx_in, ny_in
    real(rp), intent(in) :: lx_in, ly_in
    !
    integer :: i
    !
    if (nx_in <= 0 .or. ny_in <= 0) then
       stop 'STOP ERROR: nx and ny must be > 0 in CMM_InitGrid'
    end if
    !
    nx = nx_in
    ny = ny_in
    lx = lx_in
    ly = ly_in
    !
    dx = lx / real(nx, rp)
    dy = ly / real(ny, rp)
    !
    if (allocated(xc)) deallocate(xc)
    if (allocated(yc)) deallocate(yc)
    if (allocated(xu)) deallocate(xu)
    if (allocated(yu)) deallocate(yu)
    if (allocated(xv)) deallocate(xv)
    if (allocated(yv)) deallocate(yv)
    !
    allocate(xc(1:nx))
    allocate(yc(1:ny))
    !
    allocate(xu(0:nx))
    allocate(yu(1:ny))
    !
    allocate(xv(1:nx))
    allocate(yv(0:ny))
    !
    do i = 1, nx
       xc(i) = (real(i, rp) - 0.5_rp) * dx
    end do
    !
    do i = 1, ny
       yc(i) = (real(i, rp) - 0.5_rp) * dy
    end do
    !
    do i = 0, nx
       xu(i) = real(i, rp) * dx
    end do
    !
    do i = 1, ny
       yu(i) = (real(i, rp) - 0.5_rp) * dy
    end do
    !
    do i = 1, nx
       xv(i) = (real(i, rp) - 0.5_rp) * dx
    end do
    !
    do i = 0, ny
       yv(i) = real(i, rp) * dy
    end do
    !
  end subroutine CMM_InitGrid
  !
  ! ============================================================
  ! Allocate all main solver fields and work arrays.
  ! This routine assumes that the grid has already been initialized.
  ! All arrays are reset to zero and solver counters are reinitialized.
  ! ============================================================
  subroutine CMM_AllocFields()
    !
    if (nx <= 0 .or. ny <= 0) then
       stop 'STOP ERROR: call CMM_InitGrid before CMM_AllocFields'
    end if
    !
    call CMM_FreeFields()
    !
    ! ------------------------------------------------------------
    ! Allocate main flow fields including ghost layers.
    ! ------------------------------------------------------------
    allocate(u     (0:nx+1,0:ny+1))
    allocate(v     (0:nx+1,0:ny+1))
    allocate(p     (0:nx+1,0:ny+1))
    allocate(p_star(0:nx+1,0:ny+1))
    !
    allocate(ustar (0:nx+1,0:ny+1))
    allocate(vstar (0:nx+1,0:ny+1))
    allocate(rhs_p (0:nx+1,0:ny+1))
    !
    allocate(urk1  (0:nx+1,0:ny+1))
    allocate(vrk1  (0:nx+1,0:ny+1))
    allocate(urk2  (0:nx+1,0:ny+1))
    allocate(vrk2  (0:nx+1,0:ny+1))
    !
    ! ------------------------------------------------------------
    ! Allocate momentum RHS and workspace arrays.
    ! ------------------------------------------------------------
    allocate(conv_u(0:nx+1,0:ny+1))
    allocate(conv_v(0:nx+1,0:ny+1))
    allocate(diff_u(0:nx+1,0:ny+1))
    allocate(diff_v(0:nx+1,0:ny+1))
    allocate(rhs_u (0:nx+1,0:ny+1))
    allocate(rhs_v (0:nx+1,0:ny+1))
    allocate(rhss_u(0:nx+1,0:ny+1))
    allocate(rhss_v(0:nx+1,0:ny+1))
    !
    ! ------------------------------------------------------------
    ! Allocate Poisson solver workspaces.
    ! ------------------------------------------------------------
    allocate(poisson_active(0:nx+1,0:ny+1))
    allocate(poisson_diagA (0:nx+1,0:ny+1))
    allocate(poisson_r     (0:nx+1,0:ny+1))
    allocate(poisson_rhat  (0:nx+1,0:ny+1))
    allocate(poisson_p     (0:nx+1,0:ny+1))
    allocate(poisson_v     (0:nx+1,0:ny+1))
    allocate(poisson_s     (0:nx+1,0:ny+1))
    allocate(poisson_t     (0:nx+1,0:ny+1))
    allocate(poisson_phat  (0:nx+1,0:ny+1))
    allocate(poisson_shat  (0:nx+1,0:ny+1))
    allocate(poisson_Ax    (0:nx+1,0:ny+1))
    allocate(poisson_xt    (0:nx+1,0:ny+1))
    !
    ! ------------------------------------------------------------
    ! Allocate boundary-condition value arrays.
    ! ------------------------------------------------------------
    allocate(bc_u_west_val (1:ny))
    allocate(bc_u_east_val (1:ny))
    allocate(bc_u_south_val(1:nx))
    allocate(bc_u_north_val(1:nx))
    !
    allocate(bc_v_west_val (1:ny))
    allocate(bc_v_east_val (1:ny))
    allocate(bc_v_south_val(1:nx))
    allocate(bc_v_north_val(1:nx))
    !
    allocate(bc_p_west_val (1:ny))
    allocate(bc_p_east_val (1:ny))
    allocate(bc_p_south_val(1:nx))
    allocate(bc_p_north_val(1:nx))
    !
    u      = 0.0_rp
    v      = 0.0_rp
    p      = 0.0_rp
    p_star = 0.0_rp
    ustar  = 0.0_rp
    vstar  = 0.0_rp
    rhs_p  = 0.0_rp
    !
    urk1   = 0.0_rp
    vrk1   = 0.0_rp
    urk2   = 0.0_rp
    vrk2   = 0.0_rp
    !
    conv_u = 0.0_rp
    conv_v = 0.0_rp
    diff_u = 0.0_rp
    diff_v = 0.0_rp
    rhs_u  = 0.0_rp
    rhs_v  = 0.0_rp
    rhss_u = 0.0_rp
    rhss_v = 0.0_rp
    !
    poisson_active = .false.
    poisson_diagA  = 0.0_rp
    poisson_r      = 0.0_rp
    poisson_rhat   = 0.0_rp
    poisson_p      = 0.0_rp
    poisson_v      = 0.0_rp
    poisson_s      = 0.0_rp
    poisson_t      = 0.0_rp
    poisson_phat   = 0.0_rp
    poisson_shat   = 0.0_rp
    poisson_Ax     = 0.0_rp
    poisson_xt     = 0.0_rp
    !
    poisson_iter      = 0
    poisson_flag      = 0
    poisson_res0      = 0.0_rp
    poisson_res       = 0.0_rp
    poisson_converged = .false.
    !
    it_rk3 = 0
    dt_rk3 = 0.0_rp
    !
    bc_u_west_val  = 0.0_rp
    bc_u_east_val  = 0.0_rp
    bc_u_south_val = 0.0_rp
    bc_u_north_val = 0.0_rp
    !
    bc_v_west_val  = 0.0_rp
    bc_v_east_val  = 0.0_rp
    bc_v_south_val = 0.0_rp
    bc_v_north_val = 0.0_rp
    !
    bc_p_west_val  = 0.0_rp
    bc_p_east_val  = 0.0_rp
    bc_p_south_val = 0.0_rp
    bc_p_north_val = 0.0_rp
    !
  end subroutine CMM_AllocFields
  !
  ! ============================================================
  ! Deallocate all allocatable solver arrays if they are present.
  ! This routine can be called safely before a new allocation step
  ! or at the end of the simulation.
  ! ============================================================
  subroutine CMM_FreeFields()
    !
    if (allocated(u))      deallocate(u)
    if (allocated(v))      deallocate(v)
    if (allocated(p))      deallocate(p)
    if (allocated(p_star)) deallocate(p_star)
    if (allocated(ustar))  deallocate(ustar)
    if (allocated(vstar))  deallocate(vstar)
    if (allocated(rhs_p))  deallocate(rhs_p)
    !
    if (allocated(urk1))   deallocate(urk1)
    if (allocated(vrk1))   deallocate(vrk1)
    if (allocated(urk2))   deallocate(urk2)
    if (allocated(vrk2))   deallocate(vrk2)
    !
    if (allocated(conv_u)) deallocate(conv_u)
    if (allocated(conv_v)) deallocate(conv_v)
    if (allocated(diff_u)) deallocate(diff_u)
    if (allocated(diff_v)) deallocate(diff_v)
    if (allocated(rhs_u))  deallocate(rhs_u)
    if (allocated(rhs_v))  deallocate(rhs_v)
    if (allocated(rhss_u)) deallocate(rhss_u)
    if (allocated(rhss_v)) deallocate(rhss_v)
    !
    if (allocated(poisson_active)) deallocate(poisson_active)
    if (allocated(poisson_diagA))  deallocate(poisson_diagA)
    if (allocated(poisson_r))      deallocate(poisson_r)
    if (allocated(poisson_rhat))   deallocate(poisson_rhat)
    if (allocated(poisson_p))      deallocate(poisson_p)
    if (allocated(poisson_v))      deallocate(poisson_v)
    if (allocated(poisson_s))      deallocate(poisson_s)
    if (allocated(poisson_t))      deallocate(poisson_t)
    if (allocated(poisson_phat))   deallocate(poisson_phat)
    if (allocated(poisson_shat))   deallocate(poisson_shat)
    if (allocated(poisson_Ax))     deallocate(poisson_Ax)
    if (allocated(poisson_xt))     deallocate(poisson_xt)
    !
    if (allocated(bc_u_west_val))  deallocate(bc_u_west_val)
    if (allocated(bc_u_east_val))  deallocate(bc_u_east_val)
    if (allocated(bc_u_south_val)) deallocate(bc_u_south_val)
    if (allocated(bc_u_north_val)) deallocate(bc_u_north_val)
    !
    if (allocated(bc_v_west_val))  deallocate(bc_v_west_val)
    if (allocated(bc_v_east_val))  deallocate(bc_v_east_val)
    if (allocated(bc_v_south_val)) deallocate(bc_v_south_val)
    if (allocated(bc_v_north_val)) deallocate(bc_v_north_val)
    !
    if (allocated(bc_p_west_val))  deallocate(bc_p_west_val)
    if (allocated(bc_p_east_val))  deallocate(bc_p_east_val)
    if (allocated(bc_p_south_val)) deallocate(bc_p_south_val)
    if (allocated(bc_p_north_val)) deallocate(bc_p_north_val)
    !
  end subroutine CMM_FreeFields
  !
end module mod_common
