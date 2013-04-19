module cell_indexing

! this subroutine provide index mapping for quad meshes dimensioned (nx, ny)

   integer, parameter :: maxEdges = 4

   integer :: nx, ny, nVertLevels, nTracers, nVertexDegree
   real (kind=8) :: dc
   integer, dimension(20) :: nproc


   contains


   subroutine cell_indexing_read_nl()

      implicit none

      namelist /periodic_grid/ nx, ny, dc, nVertLevels, nTracers, nproc, nVertexDegree

      nx = 200
      ny = 200
      dc = 100000.
      nVertLevels = 1
      nTracers = 2
      nproc(:) = -1
      nVertexDegree = 4

      open(20,file='namelist.input',status='old')
      read(20,periodic_grid)
      close(20)

   end subroutine cell_indexing_read_nl


   subroutine cellColRow(idx, iCol, iRow)

      implicit none

      integer, intent(in) :: idx
      integer, intent(out) :: iCol, iRow

      iRow = ((idx-1) / nx) + 1
      iCol = mod((idx-1), nx) + 1

   end subroutine cellColRow


   integer function cellIdx(iCol, iRow)

      implicit none

      integer, intent(in) :: iCol, iRow

      cellIdx = (iRow-1)*nx + iCol

   end function cellIdx


   integer function cellOnCell(iCol, iRow, neighborNumber)

      implicit none

      integer, intent(in) :: iCol, iRow, neighborNumber

      integer :: mx, px, my, py

      mx = iCol - 1
      if (mx == 0) mx = nx
      my = iRow - 1
      if (my == 0) my = ny
      px = iCol + 1
      if (px == nx + 1) px = 1
      py = iRow + 1
      if (py == ny + 1) py = 1

      if (neighborNumber == 1) then
         cellOnCell = cellIdx(mx, iRow)
      else if (neighborNumber == 2) then
         cellOnCell = cellIdx(iCol, my)
      else if (neighborNumber == 3) then
         cellOnCell = cellIdx(px, iRow)
      else if (neighborNumber == 4) then
         cellOnCell = cellIdx(iCol, py)
      endif

   end function cellOnCell


   integer function edgeOnCell(iCell, neighborNumber)

      implicit none

      integer, intent(in) :: iCell, neighborNumber

      integer :: myRow, myCol

      call cellColRow(iCell, myCol, myRow)
      
      if (neighborNumber == 1) then
         edgeOnCell = 2*(iCell - 1) + 1
      else if (neighborNumber == 2) then
         edgeOnCell = 2*(iCell - 1) + 2
      else if (neighborNumber == 3) then
         edgeOnCell = 2*(cellOnCell(myCol, myRow, 3) - 1) + 1
      else if (neighborNumber == 4) then
         edgeOnCell = 2*(cellOnCell(myCol, myRow, 4) - 1) + 2
      end if

   end function edgeOnCell


   integer function vertexOnCell(iCell, neighborNumber)

      implicit none

      integer, intent(in) :: iCell, neighborNumber

      integer :: myRow, myCol, mx, my, px, py

      call cellColRow(iCell, myCol, myRow)

      mx = myCol - 1
      if (mx == 0) mx = nx
      my = myRow - 1
      if (my == 0) my = ny
      px = myCol + 1
      if (px == nx + 1) px = 1
      py = myRow + 1
      if (py == ny + 1) py = 1

      if (neighborNumber == 1) then
         vertexOnCell = cellIdx(myCol,myRow)
      else if (neighborNumber == 2) then
         vertexOnCell = cellIdx(myCol,my)
      else if (neighborNumber == 3) then
         vertexOnCell = cellIdx(px,my)
      else if (neighborNumber == 4) then
         vertexOnCell = cellIdx(px,myRow)
      end if

   end function vertexOnCell


end module cell_indexing
