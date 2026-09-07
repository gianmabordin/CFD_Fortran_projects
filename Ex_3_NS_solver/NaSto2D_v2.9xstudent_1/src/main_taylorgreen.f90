program main
  !
  use mod_types,               only: rp
  !
  use mod_common,              only: nx, ny
  use mod_common,              only: dx, dy
  use mod_common,              only: dt, dt_rk3, nu, rho
  use mod_common,              only: time, time_end
  use mod_common,              only: istep, nsteps, it_rk3
  !
  use mod_common,              only: u, v, p, p_star
  use mod_common,              only: ustar, vstar
  use mod_common,              only: rhs_p
  !
  use mod_common,              only: BC_DIRICHLET, BC_NEUMANN, BC_PERIODIC
  !
  use mod_common,              only: bc_u_west, bc_u_east, bc_u_south, bc_u_north
  use mod_common,              only: bc_v_west, bc_v_east, bc_v_south, bc_v_north
  use mod_common,              only: bc_p_west, bc_p_east, bc_p_south, bc_p_north
  !
  use mod_common,              only: bc_u_west_val, bc_u_east_val, bc_u_south_val, bc_u_north_val
  use mod_common,              only: bc_v_west_val, bc_v_east_val, bc_v_south_val, bc_v_north_val
  use mod_common,              only: bc_p_west_val, bc_p_east_val, bc_p_south_val, bc_p_north_val
  !
  use mod_common,              only: poisson_maxit, poisson_tol
  use mod_common,              only: poisson_res, poisson_converged, poisson_iter
  !
  use mod_common,              only: CMM_InitGrid, CMM_AllocFields, CMM_FreeFields
  !
  use mod_common,              only: xc, yc, xu, yu, xv, yv
  !
  use mod_bc,                  only: BC_ApplyVelocity
  use mod_bc,                  only: BC_ApplyPressure
  !
  use mod_pressure_correction, only: PC_ComputeDivergence
  !
  use mod_io,                  only: IO_WriteVTK
  use mod_io,                  only: IO_SaveRestart
  use mod_io,                  only: IO_LoadRestart
  use mod_io,                  only: IO_WriteFieldsASCII
  !
  use mod_timeint,             only: TIMEINT_AdvanceRK3FractionalStep
  !
  implicit none
  !
  integer :: log_every, log_header_every
  integer :: restart_every, vtk_every, ascii_every
  integer :: it_restart
  !
  real(rp) :: div_star_max, div_corr_max
  real(rp) :: cfl_adv, cfl_diff
  real(rp) :: umax, vmax
  !
  integer :: i,j
  real(rp) :: U0
  !
  logical :: need_final_restart, need_final_vtk, need_final_ascii
  logical :: use_restart
  !
  character(len=256) :: fname
  character(len=256) :: out_dir
  character(len=256) :: restart_file
  !
  integer :: istep_start
  !
  ! ------------------------------------------------------------
  ! Initialize grid, allocate arrays, and define the main
  ! physical, numerical, and output-control parameters.
  ! ------------------------------------------------------------
  call CMM_InitGrid(64, 64, 2.0_rp*acos(-1.0_rp), 2.0_rp*acos(-1.0_rp))
  call CMM_AllocFields()
  !
  rho = 1.0_rp
  nu  = 1.0e-2_rp
  dt  = 1.0e-2_rp
  !
  time      = 0.0_rp
  time_end  = 75.0_rp
  nsteps    = int(time_end / dt)
  istep     = 0
  it_rk3    = 0
  dt_rk3    = 0.0_rp
  !
  poisson_maxit = 5000
  poisson_tol   = 1.0e-13_rp
  !
  log_every        = 1
  log_header_every = 50
  restart_every    = 100
  vtk_every        = 100
  ascii_every      = 100
  !
  out_dir = './data'
  !
  it_restart  = 0
  !
  U0         = 1.0_rp
  !
  ! ------------------------------------------------------------
  ! Define the boundary-condition type on each side of the domain
  ! and assign the corresponding boundary values.
  ! ------------------------------------------------------------
  bc_u_west  = BC_PERIODIC
  bc_v_west  = BC_PERIODIC
  bc_u_east  = BC_PERIODIC
  bc_v_east  = BC_PERIODIC
  !
  bc_u_south = BC_PERIODIC
  bc_u_north = BC_PERIODIC
  bc_v_south = BC_PERIODIC
  bc_v_north = BC_PERIODIC
  !
  bc_p_west  = BC_PERIODIC
  bc_p_east  = BC_PERIODIC
  bc_p_south = BC_PERIODIC
  bc_p_north = BC_PERIODIC
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
  ! ------------------------------------------------------------
  ! Initialize the solution either from scratch or from a restart
  ! file, then enforce boundary conditions on the loaded state.
  ! ------------------------------------------------------------
  !
  use_restart = .false.
  !
  if(it_restart>0) use_restart = .true.
  !
  write(restart_file,'(A,"/field_",I6.6,".bin")') trim(out_dir), it_restart
  !
  u      = 0.0_rp
  v      = 0.0_rp
  p      = 0.0_rp
  p_star = 0.0_rp
  ustar  = 0.0_rp
  vstar  = 0.0_rp
  rhs_p  = 0.0_rp
  !
  if (use_restart) then
     call IO_LoadRestart(trim(restart_file))
     call BC_ApplyVelocity(u, v)
     call BC_ApplyPressure(p)
     istep_start = istep + 1
  else
     do j = 1, ny
        do i = 0, nx
           u(i,j) = U0 * sin(xu(i)) * cos(yu(j))
        end do
     end do
     !
     do j = 0, ny
        do i = 1, nx
           v(i,j) = -U0 * cos(xv(i)) * sin(yv(j))
        end do
     end do
     !
     do j = 1, ny
        do i = 1, nx
           p(i,j) = 0.25_rp * U0 * U0 * (cos(2.0_rp*xc(i)) + cos(2.0_rp*yc(j)))
        end do
     end do
     !
     call BC_ApplyVelocity(u, v)
     call BC_ApplyPressure(p)
     istep_start = 1
  end if
  !
  ! ------------------------------------------------------------
  ! Print a summary of the simulation settings and start the log.
  ! ------------------------------------------------------------
  write(*,'(A)') '===================================================================================='
  write(*,'(A)') 'NaSto'
  write(*,'(A)') 'Fractional-step RK3 with projection'
  write(*,'(A)') '===================================================================================='
  write(*,'(A,I8)')     'nx              = ', nx
  write(*,'(A,I8)')     'ny              = ', ny
  write(*,'(A,ES12.4)') 'dx              = ', dx
  write(*,'(A,ES12.4)') 'dy              = ', dy
  write(*,'(A,ES12.4)') 'dt              = ', dt
  write(*,'(A,ES12.4)') 'nu              = ', nu
  write(*,'(A,ES12.4)') 'rho             = ', rho
  write(*,'(A,ES12.4)') 'time            = ', time
  write(*,'(A,ES12.4)') 'time_end        = ', time_end
  write(*,'(A,I8)')     'istep           = ', istep
  write(*,'(A,I8)')     'nsteps          = ', nsteps
  write(*,'(A,I8)')     'log_every       = ', log_every
  write(*,'(A,I8)')     'restart_every   = ', restart_every
  write(*,'(A,I8)')     'vtk_every       = ', vtk_every
  write(*,'(A,I8)')     'ascii_every     = ', ascii_every
  write(*,'(A,I8)')     'poisson_maxit   = ', poisson_maxit
  write(*,'(A,ES12.4)') 'poisson_tol     = ', poisson_tol
  write(*,'(A,L1)')     'use_restart     = ', use_restart
  write(*,'(A,A)')      'restart_file    = ', trim(restart_file)
  write(*,'(A,A)')      'out_dir         = ', trim(out_dir)
  write(*,'(A)') '===================================================================================='
  call PrintLogHeader()
  !
  if (.not. use_restart) then
     call WriteOutputSet(istep, .true., .true., .true.)
  end if
  !
  ! ------------------------------------------------------------
  ! Main time-integration loop.
  ! At each step the code applies boundary conditions, advances
  ! the solution, computes diagnostics, prints logs, and writes
  ! output files at the selected frequency.
  ! ------------------------------------------------------------
  do istep = istep_start, nsteps
     !
     call BC_ApplyVelocity(u, v)
     call BC_ApplyPressure(p)
     !
     call TIMEINT_AdvanceRK3FractionalStep(u, v, p)
     !
     call PC_ComputeDivergence(rhs_p, ustar, vstar)
     div_star_max = maxval(abs(rhs_p(2:nx-1,2:ny-1)))
     !
     call PC_ComputeDivergence(rhs_p, u, v)
     div_corr_max = maxval(abs(rhs_p(2:nx-1,2:ny-1)))
     !
     umax = maxval(abs(u))
     vmax = maxval(abs(v))
     !
     cfl_adv  = dt * (umax / dx + vmax / dy)
     cfl_diff = nu * dt * (1.0_rp/(dx*dx) + 1.0_rp/(dy*dy))
     !
     time = time + dt
     !
     if (mod(istep, log_header_every) == 0) call PrintLogHeader()
     !
     if (mod(istep, log_every) == 0) then
        write(*,'(I6,2X,ES10.3,2X,ES10.3,2X,ES10.3,2X,ES9.2,2X,ES9.2,2X,L6,2X,I6,2X,ES10.3)') &
             istep, time, div_star_max, div_corr_max, cfl_adv, cfl_diff, &
             poisson_converged, poisson_iter, poisson_res
     end if
     !
     if (mod(istep, restart_every) == 0 .or. &
         mod(istep, vtk_every)     == 0 .or. &
         mod(istep, ascii_every)   == 0) then
        call WriteOutputSet(istep,                           &
                            mod(istep, restart_every) == 0,  &
                            mod(istep, vtk_every)     == 0,  &
                            mod(istep, ascii_every)   == 0)
     end if
     !
  end do
  !
  ! ------------------------------------------------------------
  ! Write a final set of output files and release memory.
  ! ------------------------------------------------------------
  need_final_restart = .true.
  need_final_vtk     = .true.
  need_final_ascii   = .true.
  !
  call WriteOutputSet(istep, need_final_restart, need_final_vtk, need_final_ascii)
  !
  call CMM_FreeFields()
  !
contains
  !
  ! ============================================================
  ! Print the header of the step-by-step log written to screen.
  ! ============================================================
  subroutine PrintLogHeader()
    !
    write(*,'(A)') '------------------------------------------------------------------------------------------------'
    write(*,'(A)') '   It        time        div*     div_corr    CFL_adv   CFL_diff   p_conv   p_iter     p_res'
    write(*,'(A)') '------------------------------------------------------------------------------------------------'
    !
  end subroutine PrintLogHeader
  !
  ! ============================================================
  ! Write the selected output files for a given step.
  ! Depending on the logical flags, the routine can generate a
  ! restart file, a VTK file, and an ASCII field file.
  ! ============================================================
  subroutine WriteOutputSet(step_id, do_restart, do_vtk, do_ascii)
    !
    integer, intent(in) :: step_id
    logical, intent(in) :: do_restart
    logical, intent(in) :: do_vtk
    logical, intent(in) :: do_ascii
    !
    if (do_restart) then
       write(fname,'(A,"/field_",I6.6,".bin")') trim(out_dir), step_id
       call IO_SaveRestart(trim(fname))
    end if
    !
    if (do_vtk) then
       write(fname,'(A,"/field_",I6.6,".vtk")') trim(out_dir), step_id
       call IO_WriteVTK(trim(fname))
    end if
    !
    if (do_ascii) then
       write(fname,'(A,"/field_",I6.6,".dat")') trim(out_dir), step_id
       call IO_WriteFieldsASCII(trim(fname))
    end if
    !
  end subroutine WriteOutputSet
  !
end program main
