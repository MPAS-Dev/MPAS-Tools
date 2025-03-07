program project_hexes
  
  use precision
  use projection_setup
  use mpas_geometry_utilities
  use ph_utils
  use ph_write_mesh
  use projections

  implicit none

  integer, parameter :: bZoneLevels = 7
  integer :: nCellsX, nCellsY
  real(dp), allocatable, dimension(:,:) :: xCell2d, yCell2d, zCell2d
  real(dp), allocatable, dimension(:,:,:) :: xVertex2d, yVertex2d, zVertex2d
  real(dp), allocatable, dimension(:,:,:) :: xEdge2d, yEdge2d, zEdge2d
  
  integer, allocatable, dimension(:,:) :: cellMask2d
  integer, allocatable, dimension(:,:) :: cellIndex1d

  integer, allocatable, dimension(:,:,:) :: edgeMask2d
  integer, allocatable, dimension(:,:,:) :: edgeIndex1d

  integer, allocatable, dimension(:,:,:) :: vertexMask2D
  integer, allocatable, dimension(:,:,:) :: vertexIndex1D
 
  integer :: nxCells, nyCells

  integer :: nCells, nEdges, nVertices, icell, iEdge, iVertex
  integer :: nCellsLevel(0:bZoneLevels)
  integer :: nEdgesLevel(0:bZoneLevels)
  integer :: nVerticesLevel(0:bZoneLevels)

  real(dp), dimension(6) :: de,dv
  real(dp) :: dh, dhp, pii
  integer :: i,j,iCenter,jCenter,ib,ie,iv,icount

  logical, parameter :: check_lengths = .false.

  ! mpas arrays for mesh description

  real(dp), allocatable, dimension(:) :: xCell, yCell, zCell, latCell, lonCell, mapFactorCell, meshDensity
  real(dp), allocatable, dimension(:) :: xEdge, yEdge, zEdge, latEdge, lonEdge, mapFactorEdge
  real(dp), allocatable, dimension(:) :: xVertex, yVertex, zVertex, latVertex, lonVertex, mapFactorVertex
  real(dp), allocatable, dimension(:) :: areaCell
  real(dp), allocatable, dimension(:) :: areaTriangle
  real(dp), allocatable, dimension(:) :: angleEdge
  real(dp), allocatable, dimension(:,:) :: kiteAreasOnVertex
  real(dp), allocatable, dimension(:,:) :: weightsOnEdge
  integer, allocatable, dimension(:) :: nEdgesOnCell
  integer, allocatable, dimension(:) :: cellMask
  integer, allocatable, dimension(:) :: vertexMask
  integer, allocatable, dimension(:) :: edgeMask
  integer, allocatable, dimension(:) :: nEdgesOnEdge
  integer, allocatable, dimension(:,:) :: edgesOnEdge

  integer, allocatable, dimension(:) :: indexToCellID, indexToVertexID, indexToEdgeID

  real(dp), allocatable, dimension(:) :: dcEdge, dvEdge
  real(dp), parameter :: earth_radius = 6378.14_dp*1000.0_dp  ! meters
  real(dp) :: invRadius, invRadius2

  real(dp) :: nominalMinDc

  integer, allocatable, dimension(:,:) :: cellsOnCell
  integer, allocatable, dimension(:,:) :: verticesOnCell
  integer, allocatable, dimension(:,:) :: edgesOnCell
  integer, allocatable, dimension(:,:) :: cellsOnVertex
  integer, allocatable, dimension(:,:) :: cellsOnEdge
  integer, allocatable, dimension(:,:) :: verticesOnEdge
  integer, allocatable, dimension(:,:) :: edgesOnVertex

  integer, parameter :: vertexDegree = 3
  integer, parameter :: maxEdges = 10
  integer, parameter :: maxEdges2 = 20


  call read_namelist()

  if(trim(projection_type) == 'lambert_conformal') then
     write(6,*) " projection is ",trim(projection_type)
  else
     write(6,*) ' unidentified projection '
     stop
  end if

  pii = 2.*asin(1.0_dp)

  dh = cell_spacing_km * 1000._dp
  dhp = dh*0.5_dp*sqrt(3.0_dp)

  write(6,*) ' XL and YL ',mesh_length_x_km,mesh_length_y_km
  nCellsY = int(mesh_length_y_km/cell_spacing_km)
  nCellsX = int(mesh_length_x_km/cell_spacing_km*2.0_dp/sqrt(3.0_dp))

  nxCells = nCellsX + 20
  nyCells = nCellsY + 20
  
  iCenter = nxCells/2 + 1
  jCenter = nyCells/2 + 1

  write(6,*) " "
  write(6,*) " nCellsX, nCellsY ",nCellsX,nCellsY
  write(6,*) " nxCells, nyCells ",nxCells,nyCells
  write(6,*) " iCenter, jCenter ",iCenter,jCenter
  write(6,*) " dh, dhp ",dh,dhp
  write(6,*) " "
  
  allocate(xCell2d(nxCells,nyCells))
  allocate(yCell2d(nxCells,nyCells))
  allocate(zCell2d(nxCells,nyCells))

  allocate(xVertex2d(2,nxCells,nyCells))
  allocate(yVertex2d(2,nxCells,nyCells))
  allocate(zVertex2d(2,nxCells,nyCells))

  allocate(xEdge2d(3,nxCells,nyCells))
  allocate(yEdge2d(3,nxCells,nyCells))
  allocate(zEdge2d(3,nxCells,nyCells))

  allocate(cellMask2d(nxCells,nyCells))
  allocate(cellIndex1d(nxCells,nyCells))

  allocate(vertexMask2d(2,nxCells,nyCells))
  allocate(vertexIndex1d(2,nxCells,nyCells))

  allocate(edgeMask2d(3,nxCells,nyCells))
  allocate(edgeIndex1d(3,nxCells,nyCells))

  write(6,*) " setting initial values for masks and locations "

  cellMask2d(:,:) = -1
  cellIndex1d(:,:) = -1
  vertexMask2d(:,:,:) = -1
  vertexIndex1d(:,:,:) = -1
  edgeMask2d(:,:,:) = -1
  edgeIndex1d(:,:,:) = -1

  zCell2d(:,:) = 0.0_dp
  zVertex2d(:,:,:) = 0.0_dp
  zEdge2d(:,:,:) = 0.0_dp

  ! set cell, edge and vertex x, y, and z positions for perfect hexegons on a Cartesian plane

  ! cells

  write(6,*) " setting cell locations "
  do j=1,nyCells
     do i=1,nxCells
        xCell2d(i,j) = (i-iCenter)*dhp
        yCell2d(i,j) = (j-jCenter)*dh
     enddo
  enddo
  do i=2,nxCells-1,2
     yCell2d(i,:) = yCell2d(i,:) - dh*0.5_dp
  enddo

  ! vertices

  write(6,*) " setting vertex locations "
  do j=1,nyCells
     do i=1,nxCells
        xVertex2d(1,i,j) = xCell2d(i,j) - dhp/3.0_dp  ! = dh/sqrt(3)
        xVertex2d(2,i,j) = xCell2d(i,j) - (2.0_dp)*dhp/(3.0_dp)
        yVertex2d(1,i,j) = yCell2d(i,j) + 0.5_dp*dh
        yVertex2d(2,i,j) = yCell2d(i,j)
     enddo
  enddo

  ! edges

  write(6,*) " setting edge locations "

  do j=1,nyCells
     do i=1,nxCells
        xEdge2d(1,i,j) = xCell2d(i,j) - dhp*0.5_dp
        xEdge2d(2,i,j) = xCell2d(i,j) - dhp*0.5_dp
        xEdge2d(3,i,j) = xCell2d(i,j)
        yEdge2d(1,i,j) = yCell2d(i,j) + dh*0.25_dp
        yEdge2d(2,i,j) = yCell2d(i,j) - dh*0.25_dp
        yEdge2d(3,i,j) = yCell2d(i,j) - dh*0.5_dp
     enddo
  enddo

  ! set the interior cells on the mesh

  write(6,*) " setting interior cell mask values "

  cellMask2d(11:10+nCellsX-2,10) = 0
  do j=11,10+nCellsY-2
     cellMask2d(10:10+nCellsX-1,j) = 0
  enddo

  do i=12,10+nCellsX-2,2
     cellMask2d(i,10+nCellsY-1) = 0
  end do

  ! set boundary condition zone cells - brute force

  do ib=1,7
     do j=2,nyCells-1
        do i=2,nxCells-1,2   ! lower cells in row
           if(cellMask2d(i,j).eq.-1) then  ! should this cell be in bc region ib
              if( (cellMask2d(i+1,j)   .eq. ib-1) .or.    &
                  (cellMask2d(i,j+1)   .eq. ib-1) .or.    &
                  (cellMask2d(i-1,j)   .eq. ib-1) .or.    &
                  (cellMask2d(i-1,j-1) .eq. ib-1) .or.    &
                  (cellMask2d(i,j-1)   .eq. ib-1) .or.    &
                  (cellMask2d(i+1,j-1) .eq. ib-1)      ) then
                 cellMask2d(i,j) = ib
                 ! write(6,*) " cell ",i,j," set to ",ib
              end if
           end if
        enddo
        do i=3,nxCells-2,2   ! upper cells in row
           if(cellMask2d(i,j).eq.-1) then  ! should this cell be in bc region ib
              if( (cellMask2d(i+1,j)     .eq. ib-1) .or.    &
                  (cellMask2d(i+1,j+1)   .eq. ib-1) .or.    &
                  (cellMask2d(i,j+1)     .eq. ib-1) .or.    &
                  (cellMask2d(i-1,j+1)   .eq. ib-1) .or.    &
                  (cellMask2d(i-1,j)     .eq. ib-1) .or.    &
                  (cellMask2d(i,j-1)     .eq. ib-1)      ) then
                 cellMask2d(i,j) = ib
                 ! write(6,*) " cell ",i,j," set to ",ib
              end if
           end if
        enddo
     enddo
  enddo
      
  if(check_lengths) then

  !  length checks
     
  write(6,*) " checking edge distance from cell center for each edge on cell "
  ! this indicates ordering of edges around the cell 
  do j=2, nyCells-1
     do i=2, nxCells-1, 2
        de(1) = sqrt( (xEdge2d(1,i  ,j  )-xCell2d(i,j))**2 + (yEdge2d(1,i  ,j  )-yCell2d(i,j))**2 )/dh ! owned edge
        de(2) = sqrt( (xEdge2d(2,i  ,j  )-xCell2d(i,j))**2 + (yEdge2d(2,i  ,j  )-yCell2d(i,j))**2 )/dh ! owned edge
        de(3) = sqrt( (xEdge2d(3,i  ,j  )-xCell2d(i,j))**2 + (yEdge2d(3,i  ,j  )-yCell2d(i,j))**2 )/dh ! owned edge
        de(4) = sqrt( (xEdge2d(1,i+1,j-1)-xCell2d(i,j))**2 + (yEdge2d(1,i+1,j-1)-yCell2d(i,j))**2 )/dh 
        de(5) = sqrt( (xEdge2d(2,i+1,j  )-xCell2d(i,j))**2 + (yEdge2d(2,i+1,j  )-yCell2d(i,j))**2 )/dh
        de(6) = sqrt( (xEdge2d(3,i  ,j+1)-xCell2d(i,j))**2 + (yEdge2d(3,i  ,j+1)-yCell2d(i,j))**2 )/dh
        write(6,*) " cell ",i,j," de ",de(1),de(2),de(3),de(4),de(5),de(6)
     end do

     do i=3, nxCells-2, 2
        de(1) = sqrt( (xEdge2d(1,i  ,j  )-xCell2d(i,j))**2 + (yEdge2d(1,i  ,j  )-yCell2d(i,j))**2 )/dh ! owned edge
        de(2) = sqrt( (xEdge2d(2,i  ,j  )-xCell2d(i,j))**2 + (yEdge2d(2,i  ,j  )-yCell2d(i,j))**2 )/dh ! owned edge
        de(3) = sqrt( (xEdge2d(3,i  ,j  )-xCell2d(i,j))**2 + (yEdge2d(3,i  ,j  )-yCell2d(i,j))**2 )/dh ! owned edge
        de(4) = sqrt( (xEdge2d(1,i+1,j  )-xCell2d(i,j))**2 + (yEdge2d(1,i+1,j  )-yCell2d(i,j))**2 )/dh 
        de(5) = sqrt( (xEdge2d(2,i+1,j+1)-xCell2d(i,j))**2 + (yEdge2d(2,i+1,j+1)-yCell2d(i,j))**2 )/dh
        de(6) = sqrt( (xEdge2d(3,i  ,j+1)-xCell2d(i,j))**2 + (yEdge2d(3,i  ,j+1)-yCell2d(i,j))**2 )/dh
        write(6,*) " cell ",i,j," de ",de(1),de(2),de(3),de(4),de(5),de(6)
     end do
  end do

  write(6,*) " "
  write(6,*) " checking vertex distance from center for each vertex on cell "
  ! this indicates ordering around the cell of vertices
  do j=2, nyCells-1
     do i=2, nxCells-1, 2
        dv(1) = sqrt( (xVertex2d(1,i  ,j  )-xCell2d(i,j))**2 + (yVertex2d(1,i  ,j  )-yCell2d(i,j))**2 )/dh ! owned vertex
        dv(2) = sqrt( (xVertex2d(2,i  ,j  )-xCell2d(i,j))**2 + (yVertex2d(2,i  ,j  )-yCell2d(i,j))**2 )/dh ! owned vertex
        dv(3) = sqrt( (xVertex2d(1,i  ,j-1)-xCell2d(i,j))**2 + (yVertex2d(1,i  ,j-1)-yCell2d(i,j))**2 )/dh
        dv(4) = sqrt( (xVertex2d(2,i+1,j-1)-xCell2d(i,j))**2 + (yVertex2d(2,i+1,j-1)-yCell2d(i,j))**2 )/dh 
        dv(5) = sqrt( (xVertex2d(1,i+1,j-1)-xCell2d(i,j))**2 + (yVertex2d(1,i+1,j-1)-yCell2d(i,j))**2 )/dh
        dv(6) = sqrt( (xVertex2d(2,i+1,j  )-xCell2d(i,j))**2 + (yVertex2d(2,i+1,j  )-yCell2d(i,j))**2 )/dh
        write(6,*) " cell ",i,j," dv ",dv(1),dv(2),dv(3),dv(4),dv(5),dv(6)
     end do

     do i=3, nxCells-2, 2
        dv(1) = sqrt( (xVertex2d(1,i  ,j  )-xCell2d(i,j))**2 + (yVertex2d(1,i  ,j  )-yCell2d(i,j))**2 )/dh ! owned vertex
        dv(2) = sqrt( (xVertex2d(2,i  ,j  )-xCell2d(i,j))**2 + (yVertex2d(2,i  ,j  )-yCell2d(i,j))**2 )/dh ! owned vertex
        dv(3) = sqrt( (xVertex2d(1,i  ,j-1)-xCell2d(i,j))**2 + (yVertex2d(1,i  ,j-1)-yCell2d(i,j))**2 )/dh
        dv(4) = sqrt( (xVertex2d(2,i+1,j  )-xCell2d(i,j))**2 + (yVertex2d(2,i+1,j  )-yCell2d(i,j))**2 )/dh 
        dv(5) = sqrt( (xVertex2d(1,i+1,j  )-xCell2d(i,j))**2 + (yVertex2d(1,i+1,j  )-yCell2d(i,j))**2 )/dh
        dv(6) = sqrt( (xVertex2d(2,i+1,j+1)-xCell2d(i,j))**2 + (yVertex2d(2,i+1,j+1)-yCell2d(i,j))**2 )/dh
        write(6,*) " cell ",i,j," dv ",dv(1),dv(2),dv(3),dv(4),dv(5),dv(6)
     end do
     
  end do

  write(6,*) " "
  write(6,*) " checking edge lengths for each edge on cell i,j  "

  ! this indicates the vertices on each edge
  do j=2, nyCells-1
     do i=2, nxCells-1, 2
        de(1) = sqrt( (xVertex2d(2,i  ,j  ) - xVertex2d(1,i  ,j  ))**2 + (yVertex2d(2,i  ,j  ) - yVertex2d(1,i  ,j  ))**2 )/dh
        de(2) = sqrt( (xVertex2d(1,i  ,j-1) - xVertex2d(2,i  ,j  ))**2 + (yVertex2d(1,i  ,j-1) - yVertex2d(2,i  ,j  ))**2 )/dh
        de(3) = sqrt( (xVertex2d(2,i+1,j-1) - xVertex2d(1,i  ,j-1))**2 + (yVertex2d(2,i+1,j-1) - yVertex2d(1,i  ,j-1))**2 )/dh
        de(4) = sqrt( (xVertex2d(1,i+1,j-1) - xVertex2d(2,i+1,j-1))**2 + (yVertex2d(1,i+1,j-1) - yVertex2d(2,i+1,j-1))**2 )/dh
        de(5) = sqrt( (xVertex2d(2,i+1,j  ) - xVertex2d(1,i+1,j-1))**2 + (yVertex2d(2,i+1,j  ) - yVertex2d(1,i+1,j-1))**2 )/dh
        de(6) = sqrt( (xVertex2d(1,i  ,j  ) - xVertex2d(2,i+1,j  ))**2 + (yVertex2d(1,i  ,j  ) - yVertex2d(2,i+1,j  ))**2 )/dh
        write(6,*) " cell ",i,j," dvEdge ",dv(1),dv(2),dv(3),dv(4),dv(5),dv(6)
     enddo
  end do

  do j=2, nyCells-1
     do i=3, nxCells-2, 2
        de(1) = sqrt( (xVertex2d(2,i  ,j  ) - xVertex2d(1,i  ,j  ))**2 + (yVertex2d(2,i  ,j  ) - yVertex2d(1,i  ,j  ))**2 )/dh
        de(2) = sqrt( (xVertex2d(1,i  ,j-1) - xVertex2d(2,i  ,j  ))**2 + (yVertex2d(1,i  ,j-1) - yVertex2d(2,i  ,j  ))**2 )/dh
        de(3) = sqrt( (xVertex2d(2,i+1,j  ) - xVertex2d(1,i  ,j-1))**2 + (yVertex2d(2,i+1,j  ) - yVertex2d(1,i  ,j-1))**2 )/dh
        de(4) = sqrt( (xVertex2d(1,i+1,j  ) - xVertex2d(2,i+1,j  ))**2 + (yVertex2d(1,i+1,j  ) - yVertex2d(2,i+1,j  ))**2 )/dh
        de(5) = sqrt( (xVertex2d(2,i+1,j+1) - xVertex2d(1,i+1,j  ))**2 + (yVertex2d(2,i+1,j+1) - yVertex2d(1,i+1,j  ))**2 )/dh
        de(6) = sqrt( (xVertex2d(1,i  ,j  ) - xVertex2d(2,i+1,j+1))**2 + (yVertex2d(1,i  ,j  ) - yVertex2d(2,i+1,j+1))**2 )/dh
        write(6,*) " cell ",i,j," dvEdge ",dv(1),dv(2),dv(3),dv(4),dv(5),dv(6)
     enddo
  end do

  end if

  ! get cell count for 1D unstructured mesh storage

  nCells = 0
  do j=1,nyCells
     do i=1,nxCells
        if(cellMask2d(i,j) .ge. 0) nCells = nCells+1
     end do
  end do

  write(6,*) " "
  write(6,*) " number of cells including boundary zone ",nCells
  
  ! create the 1D ordering, starting with interior cells

  icell = 0
  do i=0,bZoneLevels
     call cell_ordering_1d( cellMask2d, cellIndex1d, nxCells, nyCells, icell, nCellsLevel(i), i )
     write(6,*) ' cell mask value, number of cells ',i,nCellsLevel(i),icell
  end do

  ! set vertex mask for all edges.
  ! We start at outermost cells and proceed inward, thus the vertices will
  ! have a bdyZone index equal to the lowest of the cells sharing that vertex
  
  do ib=bZoneLevels,0,-1
     do j=2,nyCells-1
        do i=2,nxCells-1,2   ! lower cells in row
           if(cellMask2d(i,j) .eq. ib) then
              vertexMask2d(:,  i,j  ) = ib
              vertexMask2d(1,i  ,j-1) = ib
              vertexMask2d(2,i+1,j-1) = ib
              vertexMask2d(1,i+1,j-1) = ib
              vertexMask2d(2,i+1,j  ) = ib
           end if
        end do
        do i=3,nxCells-3,2   ! lower cells in row
           if(cellMask2d(i,j) .eq. ib) then
              vertexMask2d(:,  i,j  ) = ib
              vertexMask2d(1,i  ,j-1) = ib
              vertexMask2d(2,i+1,j  ) = ib
              vertexMask2d(1,i+1,j  ) = ib
              vertexMask2d(2,i+1,j+1) = ib
           end if
        end do
     end do
  end do

  ! count total number of vertices on mesh

  nVertices = 0
  do j=2,nyCells-1
     do i=2,nxCells-1
        do ie=1,2
           if(vertexMask2d(ie,i,j) .ge. 0) nVertices = nVertices + 1
        end do
     end do
  end do
  write(6,*) " number of vertices on mesh ",nVertices

  ivertex = 0
  do i=0,bZoneLevels
     call vertex_ordering_1d( vertexMask2d, vertexIndex1d, nxCells, nyCells, ivertex, nVerticesLevel(i), i )
     write(6,*) ' vertex mask value, number of vertices ',i,nVerticesLevel(i),ivertex
  end do

  ! set edge mask for all edges.
  ! We start at outermost cells and proceed inward, thus the edges will
  ! have a bdyZone index equal to the lowest of the cells sharing that edge
  
  do ib=bZoneLevels,0,-1
     do j=2,nyCells-1
        do i=2,nxCells-1,2   ! lower cells in row
           if(cellMask2d(i,j) .eq. ib) then
              edgeMask2d(1,i  ,j  ) = ib
              edgeMask2d(2,i  ,j  ) = ib
              edgeMask2d(3,i  ,j  ) = ib
              edgeMask2d(1,i+1,j-1) = ib
              edgeMask2d(2,i+1,j  ) = ib
              edgeMask2d(3,i  ,j+1) = ib
           end if
        end do
        do i=3,nxCells-3,2   ! lower cells in row
           if(cellMask2d(i,j) .eq. ib) then
              edgeMask2d(1,i  ,j  ) = ib
              edgeMask2d(2,i  ,j  ) = ib
              edgeMask2d(3,i  ,j  ) = ib
              edgeMask2d(1,i+1,j  ) = ib
              edgeMask2d(2,i+1,j+1) = ib
              edgeMask2d(3,i  ,j+1) = ib
           end if
        end do
     end do
  end do

  ! count total number of edges on mesh

  nEdges = 0
  do j=2,nyCells-1
     do i=2,nxCells-1
        do ie=1,3
           if(edgeMask2d(ie,i,j) .ge. 0) nEdges = nEdges + 1
        end do
     end do
  end do
  write(6,*) " number of edges on mesh ",nEdges

  iedge = 0
  do i=0,bZoneLevels
     call edge_ordering_1d( edgeMask2d, edgeIndex1d, nxCells, nyCells, iedge, nEdgesLevel(i), i )
     write(6,*) ' edge mask value, number of edges ',i,nEdgesLevel(i),iedge
  end do

  allocate(xCell(nCells))
  allocate(yCell(nCells))
  allocate(zCell(nCells))
  allocate(cellMask(nCells))

  allocate(xVertex(nVertices))
  allocate(yVertex(nVertices))
  allocate(zVertex(nVertices))
  allocate(VertexMask(nVertices))

  allocate(xEdge(nEdges))
  allocate(yEdge(nEdges))
  allocate(zEdge(nEdges))
  allocate(edgeMask(nEdges))

  allocate(dcEdge(nEdges))
  allocate(dvEdge(nEdges))

  icount = 0
  do j=2,nyCells-1
     do i=2,nxCells-1
        if(cellMask2d(i,j) .ge. 0) then
           icount = icount+1
           xCell(CellIndex1d(i,j)) = xCell2d(i,j)
           yCell(CellIndex1d(i,j)) = yCell2d(i,j)
           zCell(CellIndex1d(i,j)) = zCell2d(i,j)
           cellMask(CellIndex1d(i,j)) = cellMask2d(i,j)
        end if
     end do
  end do
  write(6,*) " x,y set for cells ",icount

  icount = 0
  do j=2,nyCells-1
     do i=2,nxCells-1
        do iv = 1,2
           if(vertexMask2d(iv,i,j) .ge. 0) then
              icount = icount+1
              xVertex(vertexIndex1d(iv,i,j)) = xVertex2d(iv,i,j)
              yVertex(vertexIndex1d(iv,i,j)) = yVertex2d(iv,i,j)
              zVertex(vertexIndex1d(iv,i,j)) = zVertex2d(iv,i,j)
              vertexMask(vertexIndex1d(iv,i,j)) = vertexMask2d(iv,i,j)
           end if
        end do
     end do
  end do
  write(6,*) " x,y set for vertices ",icount

  icount = 0
  do j=2,nyCells-1
     do i=2,nxCells-1
        do ie = 1,3
           if(edgeMask2d(ie,i,j) .ge. 0) then
              icount = icount+1
              xEdge(edgeIndex1d(ie,i,j)) = xEdge2d(ie,i,j)
              yEdge(edgeIndex1d(ie,i,j)) = yEdge2d(ie,i,j)
              zEdge(edgeIndex1d(ie,i,j)) = zEdge2d(ie,i,j)
              edgeMask(edgeIndex1d(ie,i,j)) = edgeMask2d(ie,i,j)
           end if
        end do
     end do
  end do
  write(6,*) " x,y set for edges ",icount

  ! we are finished with the 2d arrays except for the masks

  deallocate(xCell2d)
  deallocate(yCell2d)
  deallocate(zCell2d)
  deallocate(xVertex2d)
  deallocate(yVertex2d)
  deallocate(zVertex2d)
  deallocate(xEdge2d)
  deallocate(yEdge2d)
  deallocate(zEdge2d)

  ! 1d ordering cell neighbors

  allocate (cellsOnCell(maxEdges,nCells))
  allocate (edgesOnCell(maxEdges,nCells))
  allocate (verticesOnCell(maxEdges,nCells))
  allocate (nEdgesOnCell(nCells))

  allocate (cellsOnVertex(3,nVertices))
  allocate (edgesOnVertex(3,nVertices))
  allocate (cellsOnEdge(2,nEdges))  
  allocate (verticesOnEdge(2,nEdges))

  cellsOnEdge(:,:) = 0.

  nEdgesOnCell(:) = 6 ! perfect hex mesh

  write(6,*) " setting cellsOnCell, edgesOnCell, verticesOnCell, "
  write(6,*) " cellsOnVertex, edgesOnVertex, verticesOnEdge "

  icount = 0
  do j=2,nyCells-1
     do i=2,nxCells-1,2   ! lower cells in row
        if(cellMask2d(i,j).ge.0) then  

           ! cells - clockwise ordering starting from upper left
           cellsOnCell(1,cellIndex1d(i,j)) = cellIndex1d(i-1,j  )
           cellsOnCell(2,cellIndex1d(i,j)) = cellIndex1d(i-1,j-1)
           cellsOnCell(3,cellIndex1d(i,j)) = cellIndex1d(i  ,j-1)
           cellsOnCell(4,cellIndex1d(i,j)) = cellIndex1d(i+1,j-1)
           cellsOnCell(5,cellIndex1d(i,j)) = cellIndex1d(i+1,j  )
           cellsOnCell(6,cellIndex1d(i,j)) = cellIndex1d(i  ,j+1)

           ! vertices - same ordering as cells
           verticesOnCell(1,cellIndex1d(i,j)) = vertexIndex1d(1,i  ,j  )
           verticesOnCell(2,cellIndex1d(i,j)) = vertexIndex1d(2,i  ,j  )
           verticesOnCell(3,cellIndex1d(i,j)) = vertexIndex1d(1,i  ,j-1)
           verticesOnCell(4,cellIndex1d(i,j)) = vertexIndex1d(2,i+1,j-1)
           verticesOnCell(5,cellIndex1d(i,j)) = vertexIndex1d(1,i+1,j-1)
           verticesOnCell(6,cellIndex1d(i,j)) = vertexIndex1d(2,i+1,j  )

           ! edges - same ordering as cells
           edgesOnCell(1,cellIndex1d(i,j)) = edgeIndex1d(1,i  ,j  )
           edgesOnCell(2,cellIndex1d(i,j)) = edgeIndex1d(2,i  ,j  )
           edgesOnCell(3,cellIndex1d(i,j)) = edgeIndex1d(3,i  ,j  )
           edgesOnCell(4,cellIndex1d(i,j)) = edgeIndex1d(1,i+1,j-1)
           edgesOnCell(5,cellIndex1d(i,j)) = edgeIndex1d(2,i+1,j  )
           edgesOnCell(6,cellIndex1d(i,j)) = edgeIndex1d(3,i  ,j+1)

           ! cellsOnEdge get set twice (from each cell owning the edge),
           ! but we won't miss any doing it this way.  Also, set so that cell1 < cell2
           call setCellsOnEdge( cellsOnEdge(1,edgeIndex1d(1,i  ,j  )), cellIndex1d(i,j), cellIndex1d(i-1,j  ) )
           call setCellsOnEdge( cellsOnEdge(1,edgeIndex1d(2,i  ,j  )), cellIndex1d(i,j), cellIndex1d(i-1,j-1) )
           call setCellsOnEdge( cellsOnEdge(1,edgeIndex1d(3,i  ,j  )), cellIndex1d(i,j), cellIndex1d(i  ,j-1) )
           call setCellsOnEdge( cellsOnEdge(1,edgeIndex1d(1,i+1,j-1)), cellIndex1d(i,j), cellIndex1d(i+1,j-1) )
           call setCellsOnEdge( cellsOnEdge(1,edgeIndex1d(2,i+1,j  )), cellIndex1d(i,j), cellIndex1d(i+1,j  ) )
           call setCellsOnEdge( cellsOnEdge(1,edgeIndex1d(3,i  ,j+1)), cellIndex1d(i,j), cellIndex1d(i  ,j+1) )

           ! vertices on edge get set twice, and we'll fix the ordering later to obey right-hand-rule 
           verticesOnEdge(1,edgeIndex1d(1,i  ,j  )) = vertexIndex1d(1,i  ,j  )
           verticesOnEdge(2,edgeIndex1d(1,i  ,j  )) = vertexIndex1d(2,i  ,j  )

           verticesOnEdge(1,edgeIndex1d(2,i  ,j  )) = vertexIndex1d(2,i  ,j  )
           verticesOnEdge(2,edgeIndex1d(2,i  ,j  )) = vertexIndex1d(1,i  ,j-1)

           verticesOnEdge(1,edgeIndex1d(3,i  ,j  )) = vertexIndex1d(1,i  ,j-1)
           verticesOnEdge(2,edgeIndex1d(3,i  ,j  )) = vertexIndex1d(2,i+1,j-1)

           verticesOnEdge(1,edgeIndex1d(1,i+1,j-1)) = vertexIndex1d(2,i+1,j-1)
           verticesOnEdge(2,edgeIndex1d(1,i+1,j-1)) = vertexIndex1d(1,i+1,j-1)

           verticesOnEdge(1,edgeIndex1d(2,i+1,j  )) = vertexIndex1d(1,i+1,j-1)
           verticesOnEdge(2,edgeIndex1d(2,i+1,j  )) = vertexIndex1d(2,i+1,j  )

           verticesOnEdge(1,edgeIndex1d(3,i  ,j+1)) = vertexIndex1d(2,i+1,j  )
           verticesOnEdge(2,edgeIndex1d(3,i  ,j+1)) = vertexIndex1d(1,i  ,j  )

           ! cellsOnVertex get set 3 times for each vertex (from each cell sharing the vertex),
           ! except for the outer vertices on the regional mesh

           cellsOnVertex(1,vertexIndex1d(1,i  ,j  )) = cellIndex1d(i  ,j  )
           cellsOnVertex(2,vertexIndex1d(1,i  ,j  )) = cellIndex1d(i  ,j+1)
           cellsOnVertex(3,vertexIndex1d(1,i  ,j  )) = cellIndex1d(i-1,j  )

           cellsOnVertex(1,vertexIndex1d(2,i  ,j  )) = cellIndex1d(i  ,j  )
           cellsOnVertex(2,vertexIndex1d(2,i  ,j  )) = cellIndex1d(i-1,j  )
           cellsOnVertex(3,vertexIndex1d(2,i  ,j  )) = cellIndex1d(i-1,j-1)

           cellsOnVertex(1,vertexIndex1d(1,i  ,j-1)) = cellIndex1d(i  ,j  )
           cellsOnVertex(2,vertexIndex1d(1,i  ,j-1)) = cellIndex1d(i-1,j-1)
           cellsOnVertex(3,vertexIndex1d(1,i  ,j-1)) = cellIndex1d(i  ,j-1)

           cellsOnVertex(1,vertexIndex1d(2,i+1,j-1)) = cellIndex1d(i  ,j  )
           cellsOnVertex(2,vertexIndex1d(2,i+1,j-1)) = cellIndex1d(i  ,j-1)
           cellsOnVertex(3,vertexIndex1d(2,i+1,j-1)) = cellIndex1d(i+1,j-1)

           cellsOnVertex(1,vertexIndex1d(1,i+1,j-1)) = cellIndex1d(i  ,j  )
           cellsOnVertex(2,vertexIndex1d(1,i+1,j-1)) = cellIndex1d(i+1,j-1)
           cellsOnVertex(3,vertexIndex1d(1,i+1,j-1)) = cellIndex1d(i+1,j  )
           
           cellsOnVertex(1,vertexIndex1d(2,i+1,j  )) = cellIndex1d(i  ,j  )
           cellsOnVertex(2,vertexIndex1d(2,i+1,j  )) = cellIndex1d(i+1,j  )
           cellsOnVertex(3,vertexIndex1d(2,i+1,j  )) = cellIndex1d(i  ,j+1)
           
           ! edgesOnVertex get set 3 times for each vertex (from each edge sharing the vertex),
           ! except for the outer vertices on the regional mesh

           edgesOnVertex(1,vertexIndex1d(1,i  ,j  )) = edgeIndex1d(1,i  ,j  )
           edgesOnVertex(2,vertexIndex1d(1,i  ,j  )) = edgeIndex1d(3,i  ,j+1)
           edgesOnVertex(3,vertexIndex1d(1,i  ,j  )) = edgeIndex1d(2,i  ,j+1)

           edgesOnVertex(1,vertexIndex1d(2,i  ,j  )) = edgeIndex1d(2,i  ,j  )
           edgesOnVertex(2,vertexIndex1d(2,i  ,j  )) = edgeIndex1d(1,i  ,j  )
           edgesOnVertex(3,vertexIndex1d(2,i  ,j  )) = edgeIndex1d(3,i-1,j  )

           edgesOnVertex(1,vertexIndex1d(1,i  ,j-1)) = edgeIndex1d(3,i  ,j  )
           edgesOnVertex(2,vertexIndex1d(1,i  ,j-1)) = edgeIndex1d(2,i  ,j  )
           edgesOnVertex(3,vertexIndex1d(1,i  ,j-1)) = edgeIndex1d(1,i  ,j-1)

           edgesOnVertex(1,vertexIndex1d(2,i+1,j-1)) = edgeIndex1d(3,i  ,j  )
           edgesOnVertex(2,vertexIndex1d(2,i+1,j-1)) = edgeIndex1d(1,i+1,j-1)
           edgesOnVertex(3,vertexIndex1d(2,i+1,j-1)) = edgeIndex1d(2,i+1,j-1)

           edgesOnVertex(1,vertexIndex1d(1,i+1,j-1)) = edgeIndex1d(1,i+1,j-1)
           edgesOnVertex(2,vertexIndex1d(1,i+1,j-1)) = edgeIndex1d(3,i+1,j  )
           edgesOnVertex(3,vertexIndex1d(1,i+1,j-1)) = edgeIndex1d(2,i+1,j  )
           
           edgesOnVertex(1,vertexIndex1d(2,i+1,j  )) = edgeIndex1d(2,i+1,j  )
           edgesOnVertex(2,vertexIndex1d(2,i+1,j  )) = edgeIndex1d(1,i+1,j  )
           edgesOnVertex(3,vertexIndex1d(2,i+1,j  )) = edgeIndex1d(3,i  ,j+1)
           
           icount = icount+1
        end if
     end do

     do i=3,nxCells-2,2   ! upper cells in row
        if(cellMask2d(i,j).ge.0) then  ! clockwise ordering starting from upper left

           ! cells - clockwise ordering starting from upper left
           cellsOnCell(1,cellIndex1d(i,j)) = cellIndex1d(i-1,j+1)
           cellsOnCell(2,cellIndex1d(i,j)) = cellIndex1d(i-1,j  )
           cellsOnCell(3,cellIndex1d(i,j)) = cellIndex1d(i  ,j-1)
           cellsOnCell(4,cellIndex1d(i,j)) = cellIndex1d(i+1,j  )
           cellsOnCell(5,cellIndex1d(i,j)) = cellIndex1d(i+1,j+1)
           cellsOnCell(6,cellIndex1d(i,j)) = cellIndex1d(i  ,j+1)

           ! vertices - same ordering as cells
           verticesOnCell(1,cellIndex1d(i,j)) = vertexIndex1d(1,i  ,j  )
           verticesOnCell(2,cellIndex1d(i,j)) = vertexIndex1d(2,i  ,j  )
           verticesOnCell(3,cellIndex1d(i,j)) = vertexIndex1d(1,i  ,j-1)
           verticesOnCell(4,cellIndex1d(i,j)) = vertexIndex1d(2,i+1,j  )
           verticesOnCell(5,cellIndex1d(i,j)) = vertexIndex1d(1,i+1,j  )
           verticesOnCell(6,cellIndex1d(i,j)) = vertexIndex1d(2,i+1,j+1)

           ! edges - same ordering as cells
           edgesOnCell(1,cellIndex1d(i,j)) = edgeIndex1d(1,i  ,j  )
           edgesOnCell(2,cellIndex1d(i,j)) = edgeIndex1d(2,i  ,j  )
           edgesOnCell(3,cellIndex1d(i,j)) = edgeIndex1d(3,i  ,j  )
           edgesOnCell(4,cellIndex1d(i,j)) = edgeIndex1d(1,i+1,j  )
           edgesOnCell(5,cellIndex1d(i,j)) = edgeIndex1d(2,i+1,j+1)
           edgesOnCell(6,cellIndex1d(i,j)) = edgeIndex1d(3,i  ,j+1)

           ! cellsOnEdge get set twice (from each cell owning the edge),
           ! but we won't miss any doing it this way.  Also, set so that cell1 < cell2
           call setCellsOnEdge( cellsOnEdge(1,edgeIndex1d(1,i  ,j  )), cellIndex1d(i,j), cellIndex1d(i-1,j+1) )
           call setCellsOnEdge( cellsOnEdge(1,edgeIndex1d(2,i  ,j  )), cellIndex1d(i,j), cellIndex1d(i-1,j  ) )
           call setCellsOnEdge( cellsOnEdge(1,edgeIndex1d(3,i  ,j  )), cellIndex1d(i,j), cellIndex1d(i  ,j-1) )
           call setCellsOnEdge( cellsOnEdge(1,edgeIndex1d(1,i+1,j  )), cellIndex1d(i,j), cellIndex1d(i+1,j  ) )
           call setCellsOnEdge( cellsOnEdge(1,edgeIndex1d(2,i+1,j+1)), cellIndex1d(i,j), cellIndex1d(i+1,j+1) )
           call setCellsOnEdge( cellsOnEdge(1,edgeIndex1d(3,i  ,j+1)), cellIndex1d(i,j), cellIndex1d(i  ,j+1) )

           ! vertices on edge get set twice, and we'll fix the ordering later to obey right-hand-rule 
           verticesOnEdge(1,edgeIndex1d(1,i  ,j  )) = vertexIndex1d(1,i  ,j  )
           verticesOnEdge(2,edgeIndex1d(1,i  ,j  )) = vertexIndex1d(2,i  ,j  )

           verticesOnEdge(1,edgeIndex1d(2,i  ,j  )) = vertexIndex1d(2,i  ,j  )
           verticesOnEdge(2,edgeIndex1d(2,i  ,j  )) = vertexIndex1d(1,i  ,j-1)

           verticesOnEdge(1,edgeIndex1d(3,i  ,j  )) = vertexIndex1d(1,i  ,j-1)
           verticesOnEdge(2,edgeIndex1d(3,i  ,j  )) = vertexIndex1d(2,i+1,j  )

           verticesOnEdge(1,edgeIndex1d(1,i+1,j  )) = vertexIndex1d(2,i+1,j  )
           verticesOnEdge(2,edgeIndex1d(1,i+1,j  )) = vertexIndex1d(1,i+1,j  )

           verticesOnEdge(1,edgeIndex1d(2,i+1,j+1)) = vertexIndex1d(1,i+1,j  )
           verticesOnEdge(2,edgeIndex1d(2,i+1,j+1)) = vertexIndex1d(2,i+1,j+1)

           verticesOnEdge(1,edgeIndex1d(3,i  ,j+1)) = vertexIndex1d(2,i+1,j+1)
           verticesOnEdge(2,edgeIndex1d(3,i  ,j+1)) = vertexIndex1d(1,i  ,j  )

           ! cellsOnVertex get set 3 times for each vertex (from each cell sharing the vertex),
           ! except for the outer vertices on the regional mesh

           cellsOnVertex(1,vertexIndex1d(1,i  ,j  )) = cellIndex1d(i  ,j  )
           cellsOnVertex(2,vertexIndex1d(1,i  ,j  )) = cellIndex1d(i  ,j+1)
           cellsOnVertex(3,vertexIndex1d(1,i  ,j  )) = cellIndex1d(i-1,j+1)

           cellsOnVertex(1,vertexIndex1d(2,i  ,j  )) = cellIndex1d(i  ,j  )
           cellsOnVertex(2,vertexIndex1d(2,i  ,j  )) = cellIndex1d(i-1,j+1)
           cellsOnVertex(3,vertexIndex1d(2,i  ,j  )) = cellIndex1d(i-1,j  )

           cellsOnVertex(1,vertexIndex1d(1,i  ,j-1)) = cellIndex1d(i  ,j  )
           cellsOnVertex(2,vertexIndex1d(1,i  ,j-1)) = cellIndex1d(i-1,j  )
           cellsOnVertex(3,vertexIndex1d(1,i  ,j-1)) = cellIndex1d(i  ,j-1)

           cellsOnVertex(1,vertexIndex1d(2,i+1,j  )) = cellIndex1d(i  ,j  )
           cellsOnVertex(2,vertexIndex1d(2,i+1,j  )) = cellIndex1d(i  ,j-1)
           cellsOnVertex(3,vertexIndex1d(2,i+1,j  )) = cellIndex1d(i+1,j  )

           cellsOnVertex(1,vertexIndex1d(1,i+1,j  )) = cellIndex1d(i  ,j  )
           cellsOnVertex(2,vertexIndex1d(1,i+1,j  )) = cellIndex1d(i+1,j  )
           cellsOnVertex(3,vertexIndex1d(1,i+1,j  )) = cellIndex1d(i+1,j+1)
           
           cellsOnVertex(1,vertexIndex1d(2,i+1,j+1)) = cellIndex1d(i  ,j  )
           cellsOnVertex(2,vertexIndex1d(2,i+1,j+1)) = cellIndex1d(i+1,j+1)
           cellsOnVertex(3,vertexIndex1d(2,i+1,j+1)) = cellIndex1d(i  ,j+1)
           
           ! edgesOnVertex get set 3 times for each vertex (from each edge sharing the vertex),
           ! except for the outer vertices on the regional mesh

           edgesOnVertex(1,vertexIndex1d(1,i  ,j  )) = edgeIndex1d(1,i  ,j  )
           edgesOnVertex(2,vertexIndex1d(1,i  ,j  )) = edgeIndex1d(3,i  ,j+1)
           edgesOnVertex(3,vertexIndex1d(1,i  ,j  )) = edgeIndex1d(2,i  ,j+1)

           edgesOnVertex(1,vertexIndex1d(2,i  ,j  )) = edgeIndex1d(2,i  ,j  )
           edgesOnVertex(2,vertexIndex1d(2,i  ,j  )) = edgeIndex1d(1,i  ,j  )
           edgesOnVertex(3,vertexIndex1d(2,i  ,j  )) = edgeIndex1d(3,i-1,j+1)

           edgesOnVertex(1,vertexIndex1d(1,i  ,j-1)) = edgeIndex1d(3,i  ,j  )
           edgesOnVertex(2,vertexIndex1d(1,i  ,j-1)) = edgeIndex1d(2,i  ,j  )
           edgesOnVertex(3,vertexIndex1d(1,i  ,j-1)) = edgeIndex1d(1,i  ,j-1)

           edgesOnVertex(1,vertexIndex1d(2,i+1,j  )) = edgeIndex1d(3,i  ,j  )
           edgesOnVertex(2,vertexIndex1d(2,i+1,j  )) = edgeIndex1d(2,i+1,j  )
           edgesOnVertex(3,vertexIndex1d(2,i+1,j  )) = edgeIndex1d(1,i+1,j  )

           edgesOnVertex(1,vertexIndex1d(1,i+1,j  )) = edgeIndex1d(1,i+1,j  )
           edgesOnVertex(2,vertexIndex1d(1,i+1,j  )) = edgeIndex1d(3,i+1,j+1)
           edgesOnVertex(3,vertexIndex1d(1,i+1,j  )) = edgeIndex1d(2,i+1,j+1)
           
           edgesOnVertex(1,vertexIndex1d(2,i+1,j+1)) = edgeIndex1d(2,i+1,j+1)
           edgesOnVertex(2,vertexIndex1d(2,i+1,j+1)) = edgeIndex1d(1,i+1,j+1)
           edgesOnVertex(3,vertexIndex1d(2,i+1,j+1)) = edgeIndex1d(3,i  ,j+1)
           
           icount = icount+1
        end if
     enddo
  enddo
  write(6,*) " cellsOnCell, edgesOnCell, verticesOnCell, "
  write(6,*) " cellsOnVertex, edgesOnVertex, verticesOnEdge set "

  call check_cell_spacing(cellsOnCell, nEdgesOnCell, xCell, yCell, dh, nCells, maxEdges)
  call check_cell_center_edge_spacing( edgesOnCell, nEdgesOnCell, xCell, yCell,    &
                                       xEdge, yEdge, dh, nCells, nEdges, maxEdges, &
                                       cellMask )
  call check_cell_center_vertex_spacing( verticesOnCell, nEdgesOnCell, xCell, yCell,        &
                                         xVertex, yVertex, dh, nCells, nVertices, maxEdges, &
                                         cellMask )
  call check_cells_on_edge_spacing( cellsOnEdge, xCell, yCell, dh, nCells, nEdges )

  call check_vertices_on_edge_spacing( verticesOnEdge, xVertex, yVertex, dh, nEdges, nVertices )
  
  call check_cells_on_vertex_spacing( cellsOnVertex, xVertex, yVertex, xCell, yCell, dh, nCells, nVertices )

  call check_edges_on_vertex_spacing( edgesOnVertex, xVertex, yVertex, xEdge, yEdge, dh, nVertices, nEdges )

  ! planar geometry in MPAS 1D arrays is set
  ! next - deallocate 2D representation

  deallocate(cellMask2d)
  deallocate(vertexMask2d)
  deallocate(edgeMask2d)
  deallocate(cellIndex1d)
  deallocate(vertexIndex1d)
  deallocate(edgeIndex1d)

  !  allocate latitude and longitude space for projection to the sphere

  allocate(latCell(nCells))
  allocate(lonCell(nCells))
  allocate(latEdge(nEdges))
  allocate(lonEdge(nEdges))
  allocate(latVertex(nVertices))
  allocate(lonVertex(nVertices))

  allocate(mapFactorCell(nCells))
  allocate(meshDensity(nCells))
  allocate(mapFactorVertex(nVertices))
  allocate(mapFactorEdge(nEdges))

  allocate(areaCell(nCells))
  allocate(areaTriangle(nVertices))
  allocate(kiteAreasOnVertex(vertexDegree, nVertices))
  allocate(angleEdge(nEdges))

  allocate(weightsOnEdge(maxEdges2, nEdges))
  allocate(nEdgesOnEdge(nEdges))
  allocate(edgesOnEdge(maxEdges2, nEdges ))

  allocate(indexToCellID(nCells))
  allocate(indexToVertexID(nVertices))
  allocate(indexToEdgeID(nEdges))

  write(6,*) " "
  write(6,*) " centering in x "
  call center_coordinate(xCell,xVertex,xEdge,nCells,nVertices,nEdges)
  write(6,*) " centering in y "
  call center_coordinate(yCell,yVertex,yEdge,nCells,nVertices,nEdges)

  ! here is where we project the location of the cell, vertex and edge points to the sphere
  
  if(trim(projection_type) == 'lambert_conformal') then

    write(6,*) " "
    write(6,*) " projecting cell center locations "
    call xy_to_lambert_conformal_projection( xCell, yCell, zCell, latCell, lonCell, mapFactorCell, nCells )
    write(6,*) " "
    write(6,*) " projecting vertex locations "
    call xy_to_lambert_conformal_projection( xVertex, yVertex, zVertex, latVertex, lonVertex, mapFactorVertex, nVertices )
    write(6,*) " "
    write(6,*) " projecting edge locations "
    call xy_to_lambert_conformal_projection( xEdge, yEdge, zEdge, latEdge, lonEdge, mapFactorEdge, nEdges )

 end if

 ! now we compute all the other quantities needed for the MPAS mesh based on the points on the sphere
 
  write(6,*) " setting meshDensity "
  call set_meshDensity(meshDensity, mapFactorCell, nCells)
  write(6,*) " finished setting meshDensity "

  write(6,*) " setting indexToID for cells, vertices and edges "
  call set_indexToID( indexToCellID, indexToVertexID, indexToEdgeID, nCells, nVertices, nEdges )
  write(6,*) " finished indexToID set "

  write(6,*) " setting verticeOnEdge ordering "
  call verticesOnEdge_ordering( verticesOnEdge, cellsOnEdge, &
       xCell, yCell, zCell, xEdge, yEdge, zEdge, xVertex, yVertex, zVertex, nCells, nEdges, nVertices )
  write(6,*) " finished vertex ordering "
  
  ! set lengths between cell centers sharing an edge, edges

  write(6,*) " setting dcEdge "
  call set_dcEdge( dcEdge, xCell, yCell, zCell, xEdge, yEdge, zEdge, nominalMinDc, cellsOnEdge, nEdges, nCells )
  write(6,*) " setting dvEdge "
  call set_dvEdge( dvEdge, xVertex, yVertex, zVertex, verticesOnEdge, nEdges, nVertices )
  write(6,*) " finished setting lengths "

  write(6,*) " setting areaCell "
  call set_areaCell( areaCell, verticesOnCell, nEdgesOnCell, xCell, yCell, zCell, &
       xVertex, yVertex, zVertex, nCells, nVertices, maxEdges )
  write(6,*) " finished setting areaCell "

  write(6,*) " setting kite areas "
  call set_kiteAreasOnVertex( kiteAreasOnVertex, xVertex, yVertex, zVertex, &
       xCell, yCell, zCell,       &
       xEdge, yEdge, zEdge, &
       edgesOnVertex, cellsOnEdge, cellsOnVertex, &
       nVertices, nCells, nEdges, vertexDegree )
  write(6,*) " finished setting kite areas "

  write(6,*) " setting triangleArea "
  call set_areaTriangle( areaTriangle, kiteAreasOnVertex, nVertices, vertexDegree )
  write(6,*) " finished setting triangleArea "

  write(6,*) " setting angleEdge "
  call set_angleEdge( angleEdge, xEdge, yEdge, zEdge, latEdge, lonEdge, &
       verticesOnEdge, xVertex, yVertex, zVertex, nEdges, nVertices )
  write(6,*) " finished setting angleEdge "

  write(6,*) " setting weights for tangential velocity reconstruction "
  call set_weightsOnEdge( weightsOnEdge, nEdgesOnEdge, edgesOnEdge, areaCell, angleEdge, &
       dcEdge, dvEdge, &
       kiteAreasOnVertex, edgesOnCell, cellsOnVertex, cellsOnEdge, verticesOnCell, verticesOnEdge, &
       nEdgesOnCell, nCells, nVertices, nEdges, maxEdges, maxEdges2, vertexDegree )
  write(6,*) " finished setting weights for tangential velocity reconstruction "

  ! set mesh on a unit sphere - this is what init_atmosphere expects - it will expand it to the full sphere.

  invRadius = 1.0_dp/earth_radius
  invRadius2 = invRadius*invRadius
  
  xCell(:) = xCell(:)*invRadius
  yCell(:)= yCell(:)*invRadius
  zCell(:) = zCell(:)*invRadius
  xVertex(:) = xVertex(:)*invRadius
  yVertex(:) = yVertex(:)*invRadius
  zVertex(:) = zVertex(:)*invRadius
  xEdge(:) = xEdge(:)*invRadius
  yEdge(:) = yEdge(:)*invRadius
  zEdge(:) = zEdge(:)*invRadius
  dcEdge(:) = dcEdge(:)*invRadius
  dvEdge(:) = dvEdge(:)*invRadius
  areaCell(:) = areaCell(:)*invRadius2
  areaTriangle(:) = areaTriangle(:)*invRadius2
  kiteAreasOnVertex(:,:) = kiteAreasOnVertex(:,:)*invRadius2
  nominalMinDc = nominalMinDc*invRadius

  write(6,*) " resetting exterior index "
  call set_outside_index( cellsOnCell, EdgesOnCell, verticesOnCell, cellsOnVertex, cellsOnEdge, &
                          edgesOnVertex, verticesOnEdge, nCells, nEdges, nVertices,             &
                          maxEdges, vertexDegree )

  write(6,*) " creating mesh file "

      call write_mpas_mesh( xCell, yCell, zCell, latCell, lonCell,                                 &
                            xVertex, yVertex, zVertex, latVertex, lonVertex,                       &
                            xEdge, yEdge, zEdge, latEdge, lonEdge,                                 &
                            cellMask, VertexMask, edgeMask,                                        &
                            cellsOnCell, edgesOnCell, verticesOnCell,                              &
                            cellsOnVertex, cellsOnEdge, edgesOnVertex,                             &
                            verticesOnEdge, nEdgesOnCell,                                          &
                            indexToCellID, indexToVertexID, indexToEdgeID,                         &
                            nEdgesOnEdge, edgesOnEdge, areaCell, areaTriangle, kiteAreasOnVertex,  &
                            dvEdge, dcEdge, angleEdge, weightsOnEdge, meshDensity, nominalMinDc,   &
                            nCells, nEdges, nVertices,                                             &
                            maxEdges, maxEdges2, vertexDegree )

      call write_graph_info( cellsOnCell, cellsOnEdge, nEdgesOnCell, nCells, nEdges, maxEdges )

end program project_hexes

