program hexagonal_periodic_grid

   use cell_indexing
   use write_netcdf

   implicit none

   real (kind=8), parameter :: pi = 3.141592653589793
   real (kind=8), parameter :: ONE = 1.0_8
   real (kind=8), parameter :: TWO = 2.0_8
   real (kind=8), parameter :: THREE = 3.0_8
   real (kind=8), parameter :: FOUR = 4.0_8
   real (kind=8), parameter :: SIX = 6.0_8

   integer, allocatable, dimension(:) :: indexToCellID, indexToEdgeID, indexToVertexID
   integer, allocatable, dimension(:) :: nEdgesOnCell, nEdgesOnEdge
   integer, allocatable, dimension(:,:) :: cellsOnCell, edgesOnCell, verticesOnCell
   integer, allocatable, dimension(:,:) :: cellsOnEdge, edgesOnEdge, verticesOnEdge
   integer, allocatable, dimension(:,:) :: edgesOnVertex, cellsOnVertex, uBC
   real (kind=8), allocatable, dimension(:) :: areaTriangle, areaCell, angleEdge
   real (kind=8), allocatable, dimension(:) :: dcEdge, dvEdge
   real (kind=8), allocatable, dimension(:) :: latCell, lonCell, xCell, yCell, zCell
   real (kind=8), allocatable, dimension(:) :: latEdge, lonEdge, xEdge, yEdge, zEdge
   real (kind=8), allocatable, dimension(:) :: latVertex, lonVertex, xVertex, yVertex, zVertex
   real (kind=8), allocatable, dimension(:,:) :: weightsOnEdge, kiteAreasOnVertex
   real (kind=8), allocatable, dimension(:) :: fEdge, fVertex, h_s
   real (kind=8), allocatable, dimension(:,:,:) :: u, v, h, vh, circulation, vorticity, ke
   real (kind=8), allocatable, dimension(:,:,:,:) :: tracers

   integer :: i, j, np, iCell
   integer :: nCells, nEdges, nVertices
   integer :: iRow, iCol, ii, jj
   integer :: nprocx, nprocy
   real (kind=8) :: r
   character (len=32) :: decomp_fname

   call cell_indexing_read_nl()

   nCells = nx*ny
   nEdges = 3*nCells
   nVertices = 2*nCells

   allocate(indexToCellID(nCells))
   allocate(indexToEdgeID(nEdges))
   allocate(indexToVertexID(nVertices))

   allocate(nEdgesOnCell(nCells))
   allocate(cellsOnCell(maxEdges, nCells))
   allocate(edgesOnCell(maxEdges, nCells))
   allocate(verticesOnCell(maxEdges, nCells))

   allocate(nEdgesOnEdge(nEdges))
   allocate(cellsOnEdge(2,nEdges))
   allocate(verticesOnEdge(2,nEdges))
   allocate(edgesOnEdge(2*maxEdges,nEdges))
   allocate(weightsOnEdge(2*maxEdges,nEdges))

   allocate(edgesOnVertex(3,nVertices))
   allocate(cellsOnVertex(3,nVertices))
   allocate(kiteAreasOnVertex(3,nVertices))

   allocate(areaTriangle(nVertices))
   allocate(areaCell(nCells))

   allocate(dcEdge(nEdges))
   allocate(dvEdge(nEdges))
   allocate(angleEdge(nEdges))

   allocate(latCell(nCells))
   allocate(lonCell(nCells))
   allocate(xCell(nCells))
   allocate(yCell(nCells))
   allocate(zCell(nCells))
   allocate(latEdge(nEdges))
   allocate(lonEdge(nEdges))
   allocate(xEdge(nEdges))
   allocate(yEdge(nEdges))
   allocate(zEdge(nEdges))
   allocate(latVertex(nVertices))
   allocate(lonVertex(nVertices))
   allocate(xVertex(nVertices))
   allocate(yVertex(nVertices))
   allocate(zVertex(nVertices))

   allocate(fEdge(nEdges))
   allocate(fVertex(nVertices))
   allocate(h_s(nCells))
   allocate(uBC(nVertLevels, nEdges))

   allocate(u(nVertLevels,nEdges,1))
   allocate(v(nVertLevels,nEdges,1))
   allocate(vh(nVertLevels,nEdges,1))
   allocate(h(nVertLevels,nCells,1))
   allocate(circulation(nVertLevels,nVertices,1))
   allocate(vorticity(nVertLevels,nVertices,1))
   allocate(ke(nVertLevels,nCells,1))
   allocate(tracers(nTracers,nVertLevels,nCells,1))

   do iRow = 1, ny
   do iCol = 1, nx
      iCell = cellIdx(iCol,iRow)
      nEdgesOnCell(iCell) = 6
      do j=1,maxEdges
         cellsOnCell(j,iCell) = cellOnCell(iCol,iRow,j)
         edgesOnCell(j,iCell) = edgeOnCell(iCell,j)
         verticesOnCell(j,iCell) = vertexOnCell(iCell,j)
      end do
      do j=1,3
         cellsOnEdge(2,edgesOnCell(j,iCell)) = iCell     ! Edges owned by this cell
      end do
      do j=4,6
         cellsOnEdge(1,edgesOnCell(j,iCell)) = iCell 
      end do
      verticesOnEdge(1,edgesOnCell(1,iCell)) = verticesOnCell(2,iCell)  ! For edges owned by this cell
      verticesOnEdge(2,edgesOnCell(1,iCell)) = verticesOnCell(1,iCell)
      verticesOnEdge(1,edgesOnCell(2,iCell)) = verticesOnCell(3,iCell)
      verticesOnEdge(2,edgesOnCell(2,iCell)) = verticesOnCell(2,iCell)
      verticesOnEdge(1,edgesOnCell(3,iCell)) = verticesOnCell(4,iCell)
      verticesOnEdge(2,edgesOnCell(3,iCell)) = verticesOnCell(3,iCell)

      edgesOnEdge(1,edgesOnCell(4,iCell)) = edgesOnCell(5,iCell)
      edgesOnEdge(2,edgesOnCell(4,iCell)) = edgesOnCell(6,iCell)
      edgesOnEdge(3,edgesOnCell(4,iCell)) = edgesOnCell(1,iCell)
      edgesOnEdge(4,edgesOnCell(4,iCell)) = edgesOnCell(2,iCell)
      edgesOnEdge(5,edgesOnCell(4,iCell)) = edgesOnCell(3,iCell)

      edgesOnEdge(1,edgesOnCell(5,iCell)) = edgesOnCell(6,iCell)
      edgesOnEdge(2,edgesOnCell(5,iCell)) = edgesOnCell(1,iCell)
      edgesOnEdge(3,edgesOnCell(5,iCell)) = edgesOnCell(2,iCell)
      edgesOnEdge(4,edgesOnCell(5,iCell)) = edgesOnCell(3,iCell)
      edgesOnEdge(5,edgesOnCell(5,iCell)) = edgesOnCell(4,iCell)

      edgesOnEdge(1,edgesOnCell(6,iCell)) = edgesOnCell(1,iCell)
      edgesOnEdge(2,edgesOnCell(6,iCell)) = edgesOnCell(2,iCell)
      edgesOnEdge(3,edgesOnCell(6,iCell)) = edgesOnCell(3,iCell)
      edgesOnEdge(4,edgesOnCell(6,iCell)) = edgesOnCell(4,iCell)
      edgesOnEdge(5,edgesOnCell(6,iCell)) = edgesOnCell(5,iCell)

      edgesOnEdge(6,edgesOnCell(1,iCell)) = edgesOnCell(2,iCell)
      edgesOnEdge(7,edgesOnCell(1,iCell)) = edgesOnCell(3,iCell)
      edgesOnEdge(8,edgesOnCell(1,iCell)) = edgesOnCell(4,iCell)
      edgesOnEdge(9,edgesOnCell(1,iCell)) = edgesOnCell(5,iCell)
      edgesOnEdge(10,edgesOnCell(1,iCell)) = edgesOnCell(6,iCell)

      edgesOnEdge(6,edgesOnCell(2,iCell)) = edgesOnCell(3,iCell)
      edgesOnEdge(7,edgesOnCell(2,iCell)) = edgesOnCell(4,iCell)
      edgesOnEdge(8,edgesOnCell(2,iCell)) = edgesOnCell(5,iCell)
      edgesOnEdge(9,edgesOnCell(2,iCell)) = edgesOnCell(6,iCell)
      edgesOnEdge(10,edgesOnCell(2,iCell)) = edgesOnCell(1,iCell)

      edgesOnEdge(6,edgesOnCell(3,iCell)) = edgesOnCell(4,iCell)
      edgesOnEdge(7,edgesOnCell(3,iCell)) = edgesOnCell(5,iCell)
      edgesOnEdge(8,edgesOnCell(3,iCell)) = edgesOnCell(6,iCell)
      edgesOnEdge(9,edgesOnCell(3,iCell)) = edgesOnCell(1,iCell)
      edgesOnEdge(10,edgesOnCell(3,iCell)) = edgesOnCell(2,iCell)

      weightsOnEdge(1,edgesOnCell(4,iCell)) = ONE / THREE
      weightsOnEdge(2,edgesOnCell(4,iCell)) = ONE / SIX
      weightsOnEdge(3,edgesOnCell(4,iCell)) = 0.0
      weightsOnEdge(4,edgesOnCell(4,iCell)) = ONE / SIX
      weightsOnEdge(5,edgesOnCell(4,iCell)) = ONE / THREE

      weightsOnEdge(1,edgesOnCell(5,iCell)) = ONE / THREE
      weightsOnEdge(2,edgesOnCell(5,iCell)) = -ONE / SIX
      weightsOnEdge(3,edgesOnCell(5,iCell)) = 0.0
      weightsOnEdge(4,edgesOnCell(5,iCell)) = ONE / SIX
      weightsOnEdge(5,edgesOnCell(5,iCell)) = -ONE / THREE

      weightsOnEdge(1,edgesOnCell(6,iCell)) = -ONE / THREE
      weightsOnEdge(2,edgesOnCell(6,iCell)) = -ONE / SIX
      weightsOnEdge(3,edgesOnCell(6,iCell)) = 0.0
      weightsOnEdge(4,edgesOnCell(6,iCell)) = -ONE / SIX
      weightsOnEdge(5,edgesOnCell(6,iCell)) = -ONE / THREE

      weightsOnEdge(6,edgesOnCell(1,iCell)) = ONE / THREE
      weightsOnEdge(7,edgesOnCell(1,iCell)) = ONE / SIX
      weightsOnEdge(8,edgesOnCell(1,iCell)) = 0.0
      weightsOnEdge(9,edgesOnCell(1,iCell)) = ONE / SIX
      weightsOnEdge(10,edgesOnCell(1,iCell)) = ONE / THREE

      weightsOnEdge(6,edgesOnCell(2,iCell)) = ONE / THREE
      weightsOnEdge(7,edgesOnCell(2,iCell)) = -ONE / SIX
      weightsOnEdge(8,edgesOnCell(2,iCell)) = 0.0
      weightsOnEdge(9,edgesOnCell(2,iCell)) = ONE / SIX
      weightsOnEdge(10,edgesOnCell(2,iCell)) = -ONE / THREE

      weightsOnEdge(6,edgesOnCell(3,iCell)) = -ONE / THREE
      weightsOnEdge(7,edgesOnCell(3,iCell)) = -ONE / SIX
      weightsOnEdge(8,edgesOnCell(3,iCell)) = 0.0
      weightsOnEdge(9,edgesOnCell(3,iCell)) = -ONE / SIX
      weightsOnEdge(10,edgesOnCell(3,iCell)) = -ONE / THREE

      cellsOnVertex(3,verticesOnCell(2,iCell)) = iCell
      cellsOnVertex(1,verticesOnCell(4,iCell)) = iCell
      cellsOnVertex(2,verticesOnCell(6,iCell)) = iCell
      cellsOnVertex(1,verticesOnCell(1,iCell)) = iCell
      cellsOnVertex(2,verticesOnCell(3,iCell)) = iCell
      cellsOnVertex(3,verticesOnCell(5,iCell)) = iCell

      edgesOnVertex(1,verticesOnCell(1,iCell)) = edgesOnCell(1,iCell)
      edgesOnVertex(1,verticesOnCell(2,iCell)) = edgesOnCell(1,iCell)
      edgesOnVertex(3,verticesOnCell(3,iCell)) = edgesOnCell(2,iCell)
      edgesOnVertex(3,verticesOnCell(2,iCell)) = edgesOnCell(2,iCell)
      edgesOnVertex(2,verticesOnCell(3,iCell)) = edgesOnCell(3,iCell)
      edgesOnVertex(2,verticesOnCell(4,iCell)) = edgesOnCell(3,iCell)
   end do
   end do

   weightsOnEdge(:,:) = weightsOnEdge(:,:) * ONE/sqrt(THREE)

   do iRow = 1, ny
   do iCol = 1, nx
      iCell = cellIdx(iCol, iRow)
      indexToCellID(iCell) = iCell
      areaCell = dc*dc*sqrt(THREE) / TWO
      latCell(iCell) = 0.0
      lonCell(iCell) = 0.0

      if (mod(iRow,2) == 1) then
         xCell(iCell) = dc*real(iCol) - 0.5*dc
         yCell(iCell) = dc*real(iRow)*sqrt(THREE) / TWO
         zCell(iCell) = 0.0
      else
         xCell(iCell) = dc*real(iCol)
         yCell(iCell) = dc*real(iRow)*sqrt(THREE) / TWO
         zCell(iCell) = 0.0
      end if

      latEdge(iCell) = 0.0
      lonEdge(iCell) = 0.0
      xEdge(edgesOnCell(1,iCell)) = xCell(iCell) - 0.5*dc
      yEdge(edgesOnCell(1,iCell)) = yCell(iCell)
      zEdge(edgesOnCell(1,iCell)) = 0.0

      xEdge(edgesOnCell(2,iCell)) = xCell(iCell) - 0.5 * dc * cos(pi/THREE)
      yEdge(edgesOnCell(2,iCell)) = yCell(iCell) - 0.5 * dc * sin(pi/THREE)
      zEdge(edgesOnCell(2,iCell)) = 0.0

      xEdge(edgesOnCell(3,iCell)) = xCell(iCell) + 0.5 * dc * cos(pi/THREE)
      yEdge(edgesOnCell(3,iCell)) = yCell(iCell) - 0.5 * dc * sin(pi/THREE)
      zEdge(edgesOnCell(3,iCell)) = 0.0

      latVertex(iCell) = 0.0
      lonVertex(iCell) = 0.0
      xVertex(verticesOnCell(1,iCell)) = xCell(iCell) - 0.5*dc
      yVertex(verticesOnCell(1,iCell)) = yCell(iCell) + dc * sqrt(THREE) / SIX
      zVertex(verticesOnCell(1,iCell)) = 0.0

      xVertex(verticesOnCell(2,iCell)) = xCell(iCell) - 0.5*dc
      yVertex(verticesOnCell(2,iCell)) = yCell(iCell) - dc * sqrt(THREE) / SIX
      zVertex(verticesOnCell(2,iCell)) = 0.0

      angleEdge(edgesOnCell(1,iCell)) = 0.0
      angleEdge(edgesOnCell(2,iCell)) = pi / THREE
      angleEdge(edgesOnCell(3,iCell)) = TWO*pi/THREE
   end do
   end do

   do i=1,nEdges
      indexToEdgeID(i) = i
      nEdgesOnEdge(i) = 10 
      dcEdge(i) = dc
      dvEdge(i) = dcEdge(i) * sqrt(THREE)/THREE
   end do

   do i=1,nVertices
      indexToVertexID(i) = i
      areaTriangle(i) = dc*dc*sqrt(THREE)/FOUR
      do j=1,3
         kiteAreasOnVertex(j,i) = dc*dc*sqrt(THREE)/(TWO*SIX)
      end do
   end do


   !
   ! fill in initial conditions below
   ! NOTE: these initial conditions will likely be removed
   !   from the grid.nc files at some point (soon).
   ! Initialize fields in grid
   !

   fEdge(:) = 1.0e-4
   fVertex(:) = 1.0e-4

   h_s(:) = 0.0
   u(:,:,:) = 0.0
   v(:,:,:) = 0.0
   vh(:,:,:) = 0.0
   circulation(:,:,:) = 0.0
   vorticity(:,:,:) = 0.0
   ke(:,:,:) = 0.0
   tracers(:,:,:,:) = 0.0
   h(:,:,:) = 1.0

   do i=1,nCells
      r = sqrt((xCell(i) - (nx/2)*(10.0*dc))**2.0 + (yCell(i) - (ny/2)*(10.0*dc))**2.0)
      if (r < 10.0*10.0*dc) then
         tracers(1,1,i,1) = (20.0 / 2.0) * (1.0 + cos(pi*r/(10.0*10.0*dc))) + 0.0
         h(1,i,1) = 1.0 + 0.1*cos(pi*r/(20.0*10.0*dc))
      else
         tracers(1,1,i,1) = 0.0
         h(1,i,1) = 1.0
      end if
   end do

   !
   ! Write grid to grid.nc file
   !
   call write_netcdf_init( nCells, nEdges, nVertices, maxEdges, nVertLevels, nTracers, vertexDegree, dc, nx, ny )
 
   call write_netcdf_fields( 1, &
                             latCell, lonCell, xCell, yCell, zCell, indexToCellID, &
                             latEdge, lonEdge, xEdge, yEdge, zEdge, indexToEdgeID, &
                             latVertex, lonVertex, xVertex, yVertex, zVertex, indexToVertexID, &
                             cellsOnEdge, &
                             nEdgesOnCell, &
                             nEdgesOnEdge, &
                             edgesOnCell, &
                             edgesOnEdge, &
                             weightsOnEdge, &
                             dvEdge, &
                             dcEdge, &
                             angleEdge, &
                             areaCell, &
                             areaTriangle, &
                             cellsOnCell, &
                             verticesOnCell, &
                             verticesOnEdge, &
                             edgesOnVertex, &
                             cellsOnVertex, &
                             kiteAreasOnVertex, &
                             fEdge, &
                             fVertex, &
                             h_s, &
                             uBC, &
                             u, &
                             v, &
                             h, &
                             vh, &
                             circulation, &
                             vorticity, &
                             ke, &
                             tracers &
                            )

   call write_netcdf_finalize()

   !
   ! Write a graph.info file to be partitioned by kmetis
   !
   np = 1
   do while (nproc(np) > 0)
      call decompose_nproc(nproc(np), nprocx, nprocy)
      if (nproc(np) < 10) then
         write(decomp_fname,'(a,i1)') 'graph.info.part.',nproc(np)
      else if (nproc(np) < 100) then
         write(decomp_fname,'(a,i2)') 'graph.info.part.',nproc(np)
      else if (nproc(np) < 1000) then
         write(decomp_fname,'(a,i3)') 'graph.info.part.',nproc(np)
      else if (nproc(np) < 10000) then
         write(decomp_fname,'(a,i4)') 'graph.info.part.',nproc(np)
      end if
      indexToCellID(:) = -1
      do iRow = 1, ny
      do iCol = 1, nx
         iCell = cellIdx(iCol, iRow)
         ii = nprocx*real(iCol-1)/real(nx)
         jj = nprocy*real(iRow-1)/real(ny)
         indexToCellID(iCell) = jj*nprocx+ii
      end do
      end do
      open(21,file=trim(decomp_fname),status='unknown')
      do i=1,nCells
         write(21,*) indexToCellID(i)
      end do
      close(21)
      np = np + 1
   end do

end program hexagonal_periodic_grid


subroutine enforce_uBC(u, uBC, xCell, yCell, zCell, nCells, nEdges, nVertLevels, dc)
! this suboutine provides a hook into uBC. the uBC field is read into the ocean
! model and used to enforce boundary conditions on the velocity field.
! uBC is written to the grid.nc file, even if the forward model does not use it.

real (kind=8), intent(in) :: dc
real (kind=8), intent(inout), dimension(nVertLevels, nEdges, 1) :: u
real (kind=8), intent(in), dimension(nCells) :: xCell, yCell, zCell
integer, intent(inout), dimension(nVertLevels, nEdges) :: uBC

uBC = -10

end subroutine enforce_uBC


subroutine decompose_nproc(nproc, nprocx, nprocy)

   implicit none

   integer, intent(in) :: nproc
   integer, intent(out) :: nprocx, nprocy

   do nprocx=int(sqrt(real(nproc))),1,-1
      nprocy = nproc / nprocx
      if (nprocy == ceiling(real(nproc)/real(nprocx))) return 
   end do

end subroutine decompose_nproc
