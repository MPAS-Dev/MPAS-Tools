module cullLoops

   public :: eliminateLoops

   contains

   subroutine eliminateLoops(nCells,nEdges,nVertices,maxEdges,vertexDegree, &
                nEdgesOnCell, cellsOnCell, verticesOnEdge, cellsOnVertex, edgesOnCell, lonCell, latCell, &
                xCell, yCell, zCell, xEdge, yEdge, zEdge, xVertex, yVertex, zVertex, &
                KMT)

       implicit none

       ! intent (in)
       integer :: nCells, nEdges, nVertices, maxEdges, vertexDegree
       integer :: nEdgesOnCell(nCells), cellsOnCell(maxEdges,nCells), verticesOnEdge(2,nEdges)
       integer :: cellsOnVertex(vertexDegree,nVertices), edgesOnCell(maxEdges,nCells)
       real :: lonCell(nCells), latCell(nCells)
       real :: xCell(nCells), yCell(nCells), zCell(nCells)
       real :: xEdge(nEdges), yEdge(nEdges), zEdge(nEdges)
       real :: xVertex(nVertices), yVertex(nVertices), zVertex(nVertices)
       integer :: edgeList(nEdges), iCellMask(nCells)

       ! intent(inout)
       integer, intent(inout) :: KMT(ncells)

       ! local workspace
       integer :: iCell, jCell, oCell, lCell, iEdge, i, kCell, iSharedEdge, iStartEdge, iSave, iSweep
       integer :: iEdgeCounter, nEdgesInLoop(nCells), iCellAhead, LeftTurns, RightTurns
       logical :: connected, atBoundary, moveSouth, moveEast, atGrenwich
       real :: lat, rlat, rlon, rCenter(3), s(3), t(3), q(3), rCross, mylon, mylat, pi

       integer, dimension(:), pointer :: cellStack
       integer, dimension(:), pointer :: oceanMask
       integer :: iCellStart, nStack, addedCells
       real :: latStart, lonStart

       write(6,*) 'Culling inland seas.....'

       allocate(cellStack(nCells/2))
       allocate(oceanMask(nCells))

       oceanMask = 0
       addedCells = 0

       iCellStart = maxloc(kmt, dim=1)

       write(6,*) 'Starting index. ', iCellStart
       write(6,*) 'lat, lon: ', latCell(iCellStart), lonCell(iCellStart)
       write(6,*) 'Starting kmt: ', kmt(iCellStart)

       nStack = 1
       cellStack(nStack) = iCellStart
       oceanMask(iCellStart) = 1
       addedCells = 1

       do while(nStack > 0)
         oCell = cellStack(nStack)
         nStack = nStack - 1
         !write(6,*) '   Working on cell ', oCell, addedCells, nStack

         do i = 1, nEdgesOnCell(oCell)
            iCell = cellsOnCell(i, oCell) 

            if(kmt(iCell) > 0 .and. oceanMask(iCell) == 0) then
              nStack = nStack + 1
              cellStack(nStack) = iCell
              oceanMask(iCell) = 1
              addedCells = addedCells + 1
            end if
         end do
       end do

       where(oceanMask == 0) kmt(:) = 0

       write(6,*) addedCells, ' total cells have been in the stack.'
       write(6,*) 'Done culling inland seas.....'

       deallocate(cellStack)
       deallocate(oceanMask)

   end subroutine eliminateLoops

end module cullLoops
