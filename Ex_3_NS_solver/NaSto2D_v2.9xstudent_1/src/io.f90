module mod_io
  !
  use mod_types,  only: rp
  use mod_common, only: nx, ny
  use mod_common, only: lx, ly, dx, dy
  use mod_common, only: time, istep
  use mod_common, only: u, v, p
  !
  implicit none
  !
  private
  !
  public :: IO_SaveRestart
  public :: IO_LoadRestart
  public :: IO_WriteVTK
  public :: IO_WriteFieldsASCII
  !
contains
  !
  ! ============================================================
  ! Save a restart file containing the full simulation state.
  ! The file stores grid information, current time and step,
  ! and the full field arrays including ghost cells.
  ! ============================================================
  subroutine IO_SaveRestart(filename)
    !
    character(len=*), intent(in) :: filename
    !
    integer :: iu
    integer :: ios
    !
    open(newunit=iu, file=filename, form='unformatted', status='replace', &
         action='write', iostat=ios)
    !
    if (ios /= 0) then
       stop 'STOP ERROR: cannot open restart file for writing in IO_SaveRestart'
    end if
    !
    write(iu) nx, ny
    write(iu) lx, ly, dx, dy
    write(iu) time, istep
    write(iu) u
    write(iu) v
    write(iu) p
    !
    close(iu)
    !
    return
  end subroutine IO_SaveRestart
  !
  ! ============================================================
  ! Load a restart file and restore the simulation state.
  ! This routine assumes that the arrays are already allocated
  ! consistently with the current grid size nx x ny.
  ! ============================================================
  subroutine IO_LoadRestart(filename)
    !
    character(len=*), intent(in) :: filename
    !
    integer :: iu
    integer :: ios
    integer :: nx_file, ny_file
    real(rp) :: lx_file, ly_file, dx_file, dy_file
    !
    open(newunit=iu, file=filename, form='unformatted', status='old', &
         action='read', iostat=ios)
    !
    if (ios /= 0) then
       stop 'STOP ERROR: cannot open restart file for reading in IO_LoadRestart'
    end if
    !
    read(iu) nx_file, ny_file
    read(iu) lx_file, ly_file, dx_file, dy_file
    !
    if (nx_file /= nx .or. ny_file /= ny) then
       stop 'STOP ERROR: restart grid size mismatch in IO_LoadRestart'
    end if
    !
    read(iu) time, istep
    read(iu) u
    read(iu) v
    read(iu) p
    !
    close(iu)
    !
    return
  end subroutine IO_LoadRestart
  !
  ! ============================================================
  ! Write a legacy VTK ASCII file for visualization in ParaView.
  ! The output contains both cell-centered fields and point-wise
  ! interpolated fields so that scalar and vector data can be
  ! displayed more smoothly.
  ! ============================================================
  subroutine IO_WriteVTK(filename)
    !
    character(len=*), intent(in) :: filename
    !
    integer :: iu
    integer :: ios
    integer :: i, j
    integer :: ic, jc
    integer :: ncnt
    real(rp) :: xnode, ynode
    real(rp) :: uc, vc
    real(rp), allocatable :: pc(:,:), uc_cell(:,:), vc_cell(:,:)
    real(rp), allocatable :: pp(:,:), up(:,:), vp(:,:)
    !
    open(newunit=iu, file=filename, form='formatted', status='replace', &
         action='write', iostat=ios)
    !
    if (ios /= 0) then
       stop 'STOP ERROR: cannot open VTK file for writing in IO_WriteVTK'
    end if
    !
    allocate(pc(1:nx,1:ny))
    allocate(uc_cell(1:nx,1:ny))
    allocate(vc_cell(1:nx,1:ny))
    !
    allocate(pp(0:nx,0:ny))
    allocate(up(0:nx,0:ny))
    allocate(vp(0:nx,0:ny))
    !
    ! ----------------------------------------------------------
    ! Build cell-centered fields from the staggered variables.
    ! Pressure is already cell-centered, while velocity is
    ! averaged from face locations to cell centers.
    ! ----------------------------------------------------------
    do j = 1, ny
       do i = 1, nx
          pc(i,j)      = p(i,j)
          uc_cell(i,j) = 0.5_rp * (u(i,j) + u(i-1,j))
          vc_cell(i,j) = 0.5_rp * (v(i,j) + v(i,j-1))
       end do
    end do
    !
    ! ----------------------------------------------------------
    ! Build point-wise fields by averaging neighboring cells.
    ! This is only for visualization, so ParaView can display
    ! smoother nodal-looking scalar and vector fields.
    ! ----------------------------------------------------------
    do j = 0, ny
       do i = 0, nx
          pp(i,j) = 0.0_rp
          up(i,j) = 0.0_rp
          vp(i,j) = 0.0_rp
          ncnt    = 0
          !
          do jc = j, j+1
             do ic = i, i+1
                if (ic >= 1 .and. ic <= nx .and. jc >= 1 .and. jc <= ny) then
                   pp(i,j) = pp(i,j) + pc(ic,jc)
                   up(i,j) = up(i,j) + uc_cell(ic,jc)
                   vp(i,j) = vp(i,j) + vc_cell(ic,jc)
                   ncnt    = ncnt + 1
                end if
             end do
          end do
          !
          if (ncnt > 0) then
             pp(i,j) = pp(i,j) / real(ncnt, rp)
             up(i,j) = up(i,j) / real(ncnt, rp)
             vp(i,j) = vp(i,j) / real(ncnt, rp)
          end if
       end do
    end do
    !
    ! ----------------------------------------------------------
    ! Write the VTK file using the legacy STRUCTURED_GRID format.
    ! Grid points are placed at the Cartesian mesh nodes.
    ! ----------------------------------------------------------
    write(iu,'(A)') '# vtk DataFile Version 3.0'
    write(iu,'(A)') '2D incompressible NS output'
    write(iu,'(A)') 'ASCII'
    write(iu,'(A)') 'DATASET STRUCTURED_GRID'
    write(iu,'(A,3I8)') 'DIMENSIONS ', nx+1, ny+1, 1
    !
    write(iu,'(A,I12,A)') 'POINTS ', (nx+1)*(ny+1), ' double'
    do j = 0, ny
       ynode = real(j, rp) * dy
       do i = 0, nx
          xnode = real(i, rp) * dx
          write(iu,'(3ES24.16)') xnode, ynode, 0.0_rp
       end do
    end do
    !
    ! ----------------------------------------------------------
    ! Write CELL_DATA using the physical cell-centered fields.
    ! ----------------------------------------------------------
    write(iu,'(A,I12)') 'CELL_DATA ', nx*ny
    !
    write(iu,'(A)') 'SCALARS pressure double 1'
    write(iu,'(A)') 'LOOKUP_TABLE default'
    do j = 1, ny
       do i = 1, nx
          write(iu,'(ES24.16)') pc(i,j)
       end do
    end do
    !
    write(iu,'(A)') 'VECTORS velocity double'
    do j = 1, ny
       do i = 1, nx
          uc = uc_cell(i,j)
          vc = vc_cell(i,j)
          write(iu,'(3ES24.16)') uc, vc, 0.0_rp
       end do
    end do
    !
    write(iu,'(A)') 'SCALARS u_center double 1'
    write(iu,'(A)') 'LOOKUP_TABLE default'
    do j = 1, ny
       do i = 1, nx
          write(iu,'(ES24.16)') uc_cell(i,j)
       end do
    end do
    !
    write(iu,'(A)') 'SCALARS v_center double 1'
    write(iu,'(A)') 'LOOKUP_TABLE default'
    do j = 1, ny
       do i = 1, nx
          write(iu,'(ES24.16)') vc_cell(i,j)
       end do
    end do
    !
    ! ----------------------------------------------------------
    ! Write POINT_DATA using interpolated nodal values.
    ! These fields are mainly useful for smoother visualization.
    ! ----------------------------------------------------------
    write(iu,'(A,I12)') 'POINT_DATA ', (nx+1)*(ny+1)
    !
    write(iu,'(A)') 'SCALARS pressure_point double 1'
    write(iu,'(A)') 'LOOKUP_TABLE default'
    do j = 0, ny
       do i = 0, nx
          write(iu,'(ES24.16)') pp(i,j)
       end do
    end do
    !
    write(iu,'(A)') 'VECTORS velocity_point double'
    do j = 0, ny
       do i = 0, nx
          write(iu,'(3ES24.16)') up(i,j), vp(i,j), 0.0_rp
       end do
    end do
    !
    close(iu)
    !
    deallocate(pc, uc_cell, vc_cell)
    deallocate(pp, up, vp)
    !
    return
  end subroutine IO_WriteVTK
  !
  ! ============================================================
  ! Write the main flow fields to a simple ASCII table.
  ! This output is mainly intended for students and for quick
  ! post-processing in Python, Matlab, Julia, or spreadsheets.
  ! ============================================================
  subroutine IO_WriteFieldsASCII(filename)
    !
    character(len=*), intent(in) :: filename
    !
    integer :: iu
    integer :: ios
    integer :: i, j
    real(rp) :: xc, yc
    real(rp) :: uc, vc
    !
    open(newunit=iu, file=filename, form='formatted', status='replace', &
         action='write', iostat=ios)
    !
    if (ios /= 0) then
       stop 'STOP ERROR: cannot open ASCII field file for writing in IO_WriteFieldsASCII'
    end if
    !
    write(iu,'(A)') '# ============================================================'
    write(iu,'(A)') '# 2D incompressible NS fields'
    write(iu,'(A)') '#'
    write(iu,'(A,ES24.16)') '# time = ', time
    write(iu,'(A,I0)')       '# istep = ', istep
    write(iu,'(A,I0)')       '# nx = ', nx
    write(iu,'(A,I0)')       '# ny = ', ny
    write(iu,'(A,ES24.16)')  '# dx = ', dx
    write(iu,'(A,ES24.16)')  '# dy = ', dy
    write(iu,'(A)') '#'
    write(iu,'(A)') '# Columns:'
    write(iu,'(A)') '# i  j  x  y  p  uc  vc'
    write(iu,'(A)') '# ============================================================'
    !
    do j = 1, ny
       yc = (real(j, rp) - 0.5_rp) * dy
       do i = 1, nx
          xc = (real(i, rp) - 0.5_rp) * dx
          uc = 0.5_rp * (u(i,j) + u(i-1,j))
          vc = 0.5_rp * (v(i,j) + v(i,j-1))
          !
          write(iu,'(2(I8,1X),5(ES24.16,1X))') i, j, xc, yc, p(i,j), uc, vc
       end do
    end do
    !
    close(iu)
    !
    return
  end subroutine IO_WriteFieldsASCII
  !
end module mod_io
