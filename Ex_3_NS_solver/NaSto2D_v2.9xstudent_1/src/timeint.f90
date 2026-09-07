module mod_timeint
  !
  use mod_types,  only: rp
  !
  use mod_common, only: nx, ny
  use mod_common, only: dt, dt_rk3
  use mod_common, only: it_rk3
  !
  use mod_common, only: u, v, p
  use mod_common, only: ustar, vstar, p_star
  !
  use mod_common, only: rhs_u, rhs_v
  use mod_common, only: rhss_u, rhss_v
  use mod_common, only: rhs_p
  !
  use mod_rhs,                 only: RHS_ComputeMomentumPressure
  use mod_rhs,                 only: RHS_ComputeMomentumConvDiff
  !
  use mod_bc,                  only: BC_ApplyVelocity
  use mod_bc,                  only: BC_ApplyPressure
  use mod_bc,                  only: BC_ApplyVelocityAfterProjection
  !
  use mod_pressure_correction, only: PC_BuildRHS
  use mod_pressure_correction, only: PC_CorrectVelocity
  !
  use mod_poisson2d,           only: POISSON2D_SolveBiCGStab
  !
  implicit none
  !
  private
  !
  public :: TIMEINT_AdvanceRK3FractionalStep
  !
  real(rp), parameter :: cff_rk3(2,3) = reshape( (/ &
       32.0_rp/60.0_rp,   0.0_rp,         &
       25.0_rp/60.0_rp, -17.0_rp/60.0_rp, &
       45.0_rp/60.0_rp, -25.0_rp/60.0_rp  /), (/2,3/) )
  !
contains
  !
  ! ============================================================
  ! Advance one full time step with a three-stage RK3 fractional-step
  ! scheme. Each stage includes pressure contribution, convective and
  ! diffusive terms, projection of the predictor field, and update of
  ! the accumulated pressure.
  !
  ! ustar and vstar are the stage predictor velocities.
  ! p_star is the stage pressure correction.
  ! ============================================================
  subroutine TIMEINT_AdvanceRK3FractionalStep(u_inout, v_inout, p_inout)
    !
    real(rp), intent(inout) :: u_inout(0:nx+1,0:ny+1)
    real(rp), intent(inout) :: v_inout(0:nx+1,0:ny+1)
    real(rp), intent(inout) :: p_inout(0:nx+1,0:ny+1)
    !
    integer  :: i, j
    real(rp) :: fct_1, fct_2, fct_3
    !
    do it_rk3 = 1, 3
       !
       fct_1  = cff_rk3(1,it_rk3) * dt
       fct_2  = cff_rk3(2,it_rk3) * dt
       fct_3  = fct_1 + fct_2
       dt_rk3 = fct_3
       !
       ! -------------------------------------------------------
       ! Apply boundary conditions to the current stage fields
       ! before computing the momentum right-hand side.
       ! -------------------------------------------------------
       call BC_ApplyVelocity(u_inout, v_inout)
       call BC_ApplyPressure(p_inout)
       !
       ! -------------------------------------------------------
       ! Add the pressure-gradient contribution to the stage
       ! predictor velocity.
       ! -------------------------------------------------------
       call RHS_ComputeMomentumPressure(rhs_u, rhs_v)
       !
       do j = 1, ny
          do i = 1, nx
             ustar(i,j) = u_inout(i,j) + fct_3 * rhs_u(i,j)
             vstar(i,j) = v_inout(i,j) + fct_3 * rhs_v(i,j)
          end do
       end do
       !
       ! -------------------------------------------------------
       ! Add the convective and diffusive contribution using the
       ! RK3 coefficients and the stored RHS history.
       ! -------------------------------------------------------
       call RHS_ComputeMomentumConvDiff(rhs_u, rhs_v, u_inout, v_inout)
       !
       do j = 1, ny
          do i = 1, nx
             ustar(i,j) = ustar(i,j) + fct_1 * rhs_u(i,j) + fct_2 * rhss_u(i,j)
             vstar(i,j) = vstar(i,j) + fct_1 * rhs_v(i,j) + fct_2 * rhss_v(i,j)
          end do
       end do
       !
       ! Store the current convective-diffusive RHS so it can be
       ! reused in the following RK stage.
       do j = 1, ny
          do i = 1, nx
             rhss_u(i,j) = rhs_u(i,j)
             rhss_v(i,j) = rhs_v(i,j)
          end do
       end do
       !
       ! -------------------------------------------------------
       ! Apply boundary conditions to the predictor field before
       ! solving the pressure-correction problem.
       ! -------------------------------------------------------
       call BC_ApplyVelocity(ustar, vstar)
       !
       ! -------------------------------------------------------
       ! Project the predictor velocity to obtain a discrete
       ! divergence-free stage velocity.
       ! -------------------------------------------------------
       call ProjectStageVelocity(u_inout, v_inout, ustar, vstar, p_star, dt_rk3)
       !
       call BC_ApplyVelocityAfterProjection(u_inout, v_inout)
       !
       ! Accumulate the stage pressure correction into the total
       ! pressure field and update its boundary conditions.
       p_inout = p_inout + p_star
       call BC_ApplyPressure(p_inout)
       !
    end do
    !
  end subroutine TIMEINT_AdvanceRK3FractionalStep
  !
  ! ============================================================
  ! Project the predictor velocity of one RK stage.
  ! The routine builds the Poisson right-hand side from the
  ! predictor divergence, solves for the stage pressure correction,
  ! and then corrects the velocity field.
  !
  ! Input:
  !   uin, vin   : stage predictor velocity
  !
  ! Output:
  !   uout, vout : projected stage velocity
  !   pstage     : stage pressure correction
  ! ============================================================
  subroutine ProjectStageVelocity(uout, vout, uin, vin, pstage, dt_stage)
    !
    real(rp), intent(out)   :: uout  (0:nx+1,0:ny+1)
    real(rp), intent(out)   :: vout  (0:nx+1,0:ny+1)
    real(rp), intent(in)    :: uin   (0:nx+1,0:ny+1)
    real(rp), intent(in)    :: vin   (0:nx+1,0:ny+1)
    real(rp), intent(inout) :: pstage(0:nx+1,0:ny+1)
    real(rp), intent(in)    :: dt_stage
    !
    pstage = 0.0_rp
    !
    call PC_BuildRHS(rhs_p, uin, vin, dt_stage)
    !
    call POISSON2D_SolveBiCGStab(pstage, rhs_p)
    !
    call BC_ApplyPressure(pstage)
    !
    call PC_CorrectVelocity(uout, vout, uin, vin, pstage, dt_stage)
    !
    return
  end subroutine ProjectStageVelocity
  !
end module mod_timeint
