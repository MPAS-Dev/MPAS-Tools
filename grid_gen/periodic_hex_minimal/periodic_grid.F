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

   integer, allocatable, dimension(:,:) :: verticesOnCell, cellsOnVertex
   real (kind=8), allocatable, dimension(:) :: xCell, yCell, zCell
   real (kind=8), allocatable, dimension(:) :: xVertex, yVertex, zVertex

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

   call write_netcdf_init( nCells, nEdges, nVertices, maxEdges, vertexDegree, dc, nx, ny )

   allocate(verticesOnCell(maxEdges, nCells))

   do iRow = 1, ny
   do iCol = 1, nx
      iCell = cellIdx(iCol,iRow)
      do j=1,maxEdges
         verticesOnCell(j,iCell) = vertexOnCell(iCell,j)
      end do
    end do
    end do

   allocate(cellsOnVertex(3,nVertices))

   do iRow = 1, ny
   do iCol = 1, nx
      iCell = cellIdx(iCol,iRow)
      cellsOnVertex(3,verticesOnCell(2,iCell)) = iCell
      cellsOnVertex(1,verticesOnCell(4,iCell)) = iCell
      cellsOnVertex(2,verticesOnCell(6,iCell)) = iCell
      cellsOnVertex(1,verticesOnCell(1,iCell)) = iCell
      cellsOnVertex(2,verticesOnCell(3,iCell)) = iCell
      cellsOnVertex(3,verticesOnCell(5,iCell)) = iCell
   end do
   end do

   call write_netcdf_fields( cellsOnVertex=cellsOnVertex )

   allocate(xCell(nCells))
   allocate(yCell(nCells))
   allocate(zCell(nCells))
   allocate(xVertex(nVertices))
   allocate(yVertex(nVertices))
   allocate(zVertex(nVertices))

   do iRow = 1, ny
   do iCol = 1, nx
      iCell = cellIdx(iCol, iRow)
      if (mod(iRow,2) == 1) then
         xCell(iCell) = dc*real(iCol) - 0.5*dc
         yCell(iCell) = dc*real(iRow)*sqrt(THREE) / TWO
         zCell(iCell) = 0.0
      else
         xCell(iCell) = dc*real(iCol)
         yCell(iCell) = dc*real(iRow)*sqrt(THREE) / TWO
         zCell(iCell) = 0.0
      end if
      xVertex(verticesOnCell(1,iCell)) = xCell(iCell) - 0.5*dc
      yVertex(verticesOnCell(1,iCell)) = yCell(iCell) + dc * sqrt(THREE) / SIX
      zVertex(verticesOnCell(1,iCell)) = 0.0
      xVertex(verticesOnCell(2,iCell)) = xCell(iCell) - 0.5*dc
      yVertex(verticesOnCell(2,iCell)) = yCell(iCell) - dc * sqrt(THREE) / SIX
      zVertex(verticesOnCell(2,iCell)) = 0.0
   end do
   end do

   deallocate(verticesOnCell)

   call write_netcdf_fields( xCell=xCell, yCell=yCell, zCell=zCell )
   deallocate(xCell)
   deallocate(yCell)
   deallocate(zCell)

   call write_netcdf_fields( xVertex=xVertex, yVertex=yVertex, zVertex=zVertex ) 
   deallocate(xVertex)
   deallocate(yVertex)
   deallocate(zVertex)

   call write_netcdf_finalize()

end program hexagonal_periodic_grid
