module mod_param
    use iso_fortran_env, only: real32, real64
    implicit none
    
    ! The 'save' attribute ensures these variables remain in memory 
    ! when the module is accessed from different parts of the code.
    save

    ! =====================================================================
    ! 1. PRECISION
    ! =====================================================================
    ! IMPORTANT: use comment in order to select the precision desired
    integer, parameter :: rp = real64  ! Double precision
    !integer, parameter :: rp = real32   ! Single precision

    ! =====================================================================
    ! 2. DOMAIN PARAMETERS (Constants)
    ! =====================================================================
    real(rp), parameter :: L = 1.0_rp
    

    ! =====================================================================
    ! 3. WORKING VARIABLES (Dynamically modified by the main program)
    ! =====================================================================
    integer  :: Nx             ! Current number of grid cells (e.g., 100 or 50)
    real(rp) :: dx             ! Current spatial step size
    real(rp) :: dt             ! Current time step size
    real(rp) :: T_final        ! Final simulation time
    
    ! Physical parameters of the PDE
    real(rp) :: c              ! Convection velocity
    real(rp) :: nu             ! Diffusion coefficient
    
    ! Parameters for the Gaussian initial condition
    real(rp) :: x0 = 0.5_rp    ! Center of the Gaussian (L/2)
    real(rp) :: sigma = 0.1_rp ! Standard deviation (width)
    
    ! Initial condition selector:
    ! 1 = Gaussian, 2 = Harmonic
    integer  :: ic_type        

end module mod_param