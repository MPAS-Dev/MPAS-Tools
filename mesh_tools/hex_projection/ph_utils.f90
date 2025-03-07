module ph_utils

  use precision
  use mpas_geometry_utilities
  
  contains

  !------------------

subroutine cell_ordering_1d( cellMask2d, cellIndex1d, nxCells, nyCells, icell, nLevel, level )

  implicit none
  integer :: nxCells, nyCells, icell, nLevel, level
  integer, dimension(nxCells,nyCells) :: cellMask2d, cellIndex1d

  integer :: i,j,icellStart

  icellStart = icell
  do j=1,nyCells
     do i=1,nxCells
        if(cellMask2d(i,j) .eq. level) then
           icell = icell+1
           cellIndex1d(i,j) = icell
        end if
     end do
  end do
  nLevel = icell - icellStart

end subroutine cell_ordering_1d

subroutine vertex_ordering_1d( vertexMask2d, vertexIndex1d, nxCells, nyCells, ivertex, nLevel, level )

  implicit none
  integer :: nxCells, nyCells, ivertex, nLevel, level
  integer, dimension(2,nxCells,nyCells) :: vertexMask2d, vertexIndex1d

  integer :: i,j,ivertexStart,iv

  ivertexStart = ivertex
  do j=1,nyCells
     do i=1,nxCells
        do iv=1,2
           if(vertexMask2d(iv,i,j) .eq. level) then
              ivertex = ivertex+1
              vertexIndex1d(iv,i,j) = ivertex
           end if
        end do
     end do
  end do
  nLevel = ivertex - ivertexStart

end subroutine vertex_ordering_1d

subroutine edge_ordering_1d( edgeMask2d, edgeIndex1d, nxCells, nyCells, iedge, nLevel, level )

  implicit none
  integer :: nxCells, nyCells, iedge, nLevel, level
  integer, dimension(3,nxCells,nyCells) :: edgeMask2d, edgeIndex1d

  integer :: i,j,iedgeStart,ie

  iedgeStart = iedge
  do j=1,nyCells
     do i=1,nxCells
        do ie=1,3
           if(edgeMask2d(ie,i,j) .eq. level) then
              iedge = iedge+1
              edgeIndex1d(ie,i,j) = iedge
           end if
        end do
     end do
  end do
  nLevel = iedge - iedgeStart

end subroutine edge_ordering_1d

subroutine check_cell_spacing(cellsOnCell, nEdgesOnCell, xCell, yCell, dh, nCells, maxEdges)
  
  implicit none

  integer :: nCells, MaxEdges
  real(dp) :: dh
  real(dp), dimension(nCells) :: xCell, yCell
  integer, dimension(nCells) :: nEdgesOnCell
  integer, dimension(maxEdges,nCells) :: cellsOnCell

  integer :: i,ie,icell,ibad
  real(dp) :: dce

  ibad = 0
  do i=1,nCells
     do ie=1,nEdgesOnCell(i)
        icell = cellsOnCell(ie,i)
        if(icell.gt.0) then
           dce = 1.0_dp - sqrt( (xCell(i)-xCell(icell))**2 + (yCell(i)-yCell(icell))**2 )/dh
           if(abs(dce) .gt. 1.0e-10) ibad = ibad+1
        end if
     end do
  end do

  write(6,*) " cell spacing check on 1D mesh perfect hexes, bad distances ",ibad

end subroutine check_cell_spacing

subroutine check_cell_Center_edge_spacing(edgesOnCell, nEdgesOnCell, xCell, yCell, &
     xEdge, yEdge, dh, nCells, nEdges, maxEdges,  &
  cellMask )
  
  implicit none

  integer :: nCells, nEdges,MaxEdges
  real(dp) :: dh
  real(dp), dimension(nCells) :: xCell, yCell
  real(dp), dimension(nEdges) :: xEdge, yEdge
  integer, dimension(nCells) :: nEdgesOnCell, cellMask
  integer, dimension(maxEdges,nCells) :: edgesOnCell

  integer :: i,ie,iedge,ibad
  real(dp) :: dce

  ibad = 0
  do i=1,nCells
     do ie=1,nEdgesOnCell(i)
        iedge = edgesOnCell(ie,i)
        dce = 0.5_dp - sqrt( (xCell(i)-xEdge(iedge))**2 + (yCell(i)-yEdge(iedge))**2 )/dh
        if(abs(dce) .gt. 1.0e-10) then
           ibad = ibad+1
           write(6,*) " bad edge distance ",i,cellmask(i),iedge,dce+0.5_dp
        end if
     end do
  end do

  write(6,*) " cell_center-edge spacing check on 1D mesh perfect hexes, bad distances ",ibad

end subroutine check_cell_center_edge_spacing

subroutine check_cell_center_vertex_spacing(verticesOnCell, nEdgesOnCell, xCell, yCell, &
                              xVertex, yVertex, dh, nCells, nVertices, maxEdges,  &
                             cellMask )
  
  implicit none

  integer :: nCells, nVertices, MaxEdges
  real(dp) :: dh
  real(dp), dimension(nCells) :: xCell, yCell
  real(dp), dimension(nVertices) :: xVertex, yVertex
  integer, dimension(nCells) :: nEdgesOnCell, cellMask
  integer, dimension(maxEdges,nCells) :: verticesOnCell

  integer :: i,iv,ivertex,ibad
  real(dp) :: dce, dexact

  dexact = 1.0_dp/sqrt(3.0_dp)
  ibad = 0
  !  do i=1,nCells
    do i=1,10
     do iv=1,nEdgesOnCell(i)
        ivertex = verticesOnCell(iv,i)
        dce = dexact - sqrt( (xCell(i)-xVertex(ivertex))**2 + (yCell(i)-yVertex(ivertex))**2 )/dh
        if(abs(dce) .gt. 1.0e-10) then
           ibad = ibad+1
           write(6,*) " bad vertex distance ",i,cellmask(i),ivertex,dce+dexact
        end if
     end do
  end do

  write(6,*) " cell_center-vertex spacing check on 1D mesh perfect hexes, bad distances ",ibad

end subroutine check_cell_center_vertex_spacing

!-------------------

subroutine setCellsOnEdge(cellsOnEdge,cell1,cell2)
    
  implicit none

  integer, dimension(2) :: cellsOnEdge
  integer :: cell1, cell2

  if(cell1 .lt. cell2) then
     cellsOnEdge(1) = cell1
     cellsOnEdge(2) = cell2
  else
     cellsOnEdge(1) = cell2
     cellsOnEdge(2) = cell1
  end if

end subroutine setCellsOnEdge

!-------------------

subroutine check_cells_on_edge_spacing( cellsOnEdge, xCell, yCell, dh, nCells, nEdges )

  implicit none

  integer :: nCells, nEdges
  integer, dimension(2,nEdges) :: cellsOnEdge
  real(dp), dimension(nCells) :: xCell, yCell
  real(dp) :: dh, dhe

  integer :: i, cell1, cell2, level7_count, bad_count

  level7_count = 0
  bad_count = 0

  do i=1,nEdges
     cell1 = cellsOnEdge(1,i)
     cell2 = cellsOnEdge(2,i)
     if( (cell1.gt.0) .and. (cell2.gt.0)) then
        dhe = sqrt( (xCell(cell1)-xCell(cell2))**2 + (yCell(cell1)-yCell(cell2))**2 )/dh
        if( abs(1.0-dhe) .gt. 1.e-10) bad_count = bad_count+1
     else if( (cell1.le.0) .and. (cell2.le.0)) then
        write(6,*) " error stop, cellsOnedge check, bot cells <= 0, edge ",i
        stop
     else
        level7_count = level7_count + 1
     end if
  end do

  write(6,*) " cellsOnEdge check, bad_distances and one-cell-on-edge count ",bad_count,level7_count

end subroutine check_cells_on_edge_spacing

!------------------

subroutine check_vertices_on_edge_spacing( verticesOnEdge, xVertex, yVertex, dh, nEdges, nVertices )

  implicit none

  integer :: nEdges, nVertices

  integer, dimension(2,nEdges) :: verticesOnEdge
  real(dp), dimension(nVertices) :: xVertex, yVertex
  real(dp) :: dh

  integer :: i, bad_count, vertex1, vertex2
  real(dp) :: dhe, dhe_exact

  bad_count = 0

  dhe_exact = sqrt(3.0_dp)/3.0_dp
  do i=1,nEdges
     vertex1 = verticesOnEdge(1,i)
     vertex2 = verticesOnEdge(2,i)
     dhe = sqrt( (xVertex(vertex1)-xVertex(vertex2))**2 + (yVertex(vertex1)-yVertex(vertex2))**2 )/dh
     if( abs(dhe-dhe_exact) .gt. 1.e-10 ) bad_count = bad_count+1
  end do

  write(6,*) " edge lengths from vertices on edge locations, bad lengths ",bad_count

end subroutine check_vertices_on_edge_spacing

!--------------------------

subroutine check_cells_on_vertex_spacing( cellsOnVertex, xVertex, yVertex, xCell, yCell, dh, nCells, nVertices )

  implicit none

  integer :: nCells, nVertices
  integer, dimension(3,nVertices) :: cellsOnVertex

  real(dp), dimension(nCells) :: xCell, yCell
  real(dp), dimension(nVertices) :: xVertex, yVertex
  real(dp) :: dh

  integer :: good_cells, i, bad_count_dc, bad_count_dcv
  integer :: vertex_one_cell_count, vertex_two_cell_count
  real(dp) :: dc12, dc23, dc31, dcv1, dexact
  integer :: cell1, cell2, cell3

  bad_count_dcv = 0
  bad_count_dc = 0
  vertex_one_cell_count = 0
  vertex_two_cell_count = 0
  dexact = 1.0_dp/sqrt(3.0_dp)
        
  do i=1,nVertices

     cell1 = cellsOnVertex(1,i)
     cell2 = cellsOnVertex(2,i)
     cell3 = cellsOnVertex(3,i)

     ! note: cells not contained in the mesh have a cell index of -1, hence the logic below
     good_cells = min(3,min(1,cell1) + min(1,cell2) + min(1, cell3))
     if(good_cells .gt. 3) good_cells = 3
     if(good_cells .eq. 1) good_cells = 2
     if(good_cells .eq. -1) good_cells = 1
     if(good_cells .eq. -3) then
        write(6,*) " no cells associated with vertex ",i," error stop "
        stop
     end if

     ! check cell separations - should be equal to dh
     
     if(good_cells .eq. 3) then  ! vertex not on outer edge 

        dc12 = sqrt( (xCell(cell1)-xCell(cell2))**2 +(yCell(cell1)-yCell(cell2))**2 )/dh
        dc23 = sqrt( (xCell(cell3)-xCell(cell2))**2 +(yCell(cell3)-yCell(cell2))**2 )/dh
        dc31 = sqrt( (xCell(cell1)-xCell(cell3))**2 +(yCell(cell1)-yCell(cell3))**2 )/dh
        if( abs(1.0_dp - dc12) .gt. 1.e-10 ) bad_count_dc = bad_count_dc + 1
        if( abs(1.0_dp - dc23) .gt. 1.e-10 ) bad_count_dc = bad_count_dc + 1
        if( abs(1.0_dp - dc31) .gt. 1.e-10 ) bad_count_dc = bad_count_dc + 1

     else if (good_cells .eq. 2) then  ! vertex on an outer edge
        if(cell1 .eq. -1) then
           dc23 = sqrt( (xCell(cell3)-xCell(cell2))**2 +(yCell(cell3)-yCell(cell2))**2 )/dh
           if( abs(1.0_dp - dc23) .gt. 1.e-10 ) bad_count_dc = bad_count_dc + 1
        else if(cell2 .eq. -1) then
           dc31 = sqrt( (xCell(cell1)-xCell(cell3))**2 +(yCell(cell1)-yCell(cell3))**2 )/dh
           if( abs(1.0_dp - dc31) .gt. 1.e-10 ) bad_count_dc = bad_count_dc + 1
        else
           dc12 = sqrt( (xCell(cell1)-xCell(cell2))**2 +(yCell(cell1)-yCell(cell2))**2 )/dh
           if( abs(1.0_dp - dc12) .gt. 1.e-10 ) bad_count_dc = bad_count_dc + 1
        end if
        vertex_two_cell_count = vertex_two_cell_count + 1
           
     else 
        vertex_one_cell_count = vertex_one_cell_count + 1
     end if

     ! check cell - vertex separation
     
     if(cell1 .gt. 0) then
        dcv1 = sqrt( (xVertex(i)-xCell(cell1))**2 + (yVertex(i)-yCell(cell1))**2 ) / dh
        if( abs(dexact - dcv1) .gt. 1.e-10 ) bad_count_dcv = bad_count_dcv + 1
     end if
     if(cell2 .gt. 0) then
        dcv1 = sqrt( (xVertex(i)-xCell(cell2))**2 + (yVertex(i)-yCell(cell2))**2 ) / dh
        if( abs(dexact - dcv1) .gt. 1.e-10 ) bad_count_dcv = bad_count_dcv + 1
     end if
     if(cell3 .gt. 0) then
        dcv1 = sqrt( (xVertex(i)-xCell(cell3))**2 + (yVertex(i)-yCell(cell3))**2 ) / dh
        if( abs(dexact - dcv1) .gt. 1.e-10 ) bad_count_dcv = bad_count_dcv + 1
     end if

  end do

  write(6,*) " cells on vertex, bad cell seperation distance count and one,two-cell count ",  &
              bad_count_dc,vertex_one_cell_count, vertex_two_cell_count
  write(6,*) " cells on vertex, bad cell-vertex seperation distance ", bad_count_dcv

end subroutine check_cells_on_vertex_spacing

!---------------

subroutine check_edges_on_vertex_spacing( edgesOnVertex, xVertex, yVertex, xEdge, yEdge, dh, nVertices, nEdges )

  implicit none

  integer :: nEdges, nVertices
  integer, dimension(3,nVertices) :: edgesOnVertex
  real(dp), dimension(nEdges) :: xEdge, yEdge
  real(dp), dimension(nVertices) :: xVertex, yVertex
  real(dp) :: dh

  integer :: i, nev, bad_lengths
  integer :: edge1, edge2, edge3
  real(dp) :: dve1, dve2, dve3, dve_exact

  bad_lengths = 0
  dve_exact = sqrt(3.0_dp)/6.0_dp

  do i=1,nVertices
     edge1 = max(0,edgesOnVertex(1,i))
     edge2 = max(0,edgesOnVertex(2,i))
     edge3 = max(0,edgesOnVertex(3,i))
     nev = min(1,edge1) + min(1,edge2) + min(1,edge3)

     if(nev .eq. 3) then ! 3 edges on vertex

        ! check distances from vertex to edge center
        dve1 = sqrt( (xVertex(i)-xEdge(edge1))**2 + (yVertex(i)-yEdge(edge1))**2 ) / dh
        dve2 = sqrt( (xVertex(i)-xEdge(edge2))**2 + (yVertex(i)-yEdge(edge2))**2 ) / dh
        dve3 = sqrt( (xVertex(i)-xEdge(edge3))**2 + (yVertex(i)-yEdge(edge3))**2 ) / dh
        if( abs(dve_exact - dve1) .gt. 1.e-10 ) bad_lengths = bad_lengths + 1
        if( abs(dve_exact - dve2) .gt. 1.e-10 ) bad_lengths = bad_lengths + 1
        if( abs(dve_exact - dve3) .gt. 1.e-10 ) bad_lengths = bad_lengths + 1

     else if(nev .eq. 2) then

        if(edge1 .eq. 0) then
           dve2 = sqrt( (xVertex(i)-xEdge(edge2))**2 + (yVertex(i)-yEdge(edge2))**2 ) / dh
           dve3 = sqrt( (xVertex(i)-xEdge(edge3))**2 + (yVertex(i)-yEdge(edge3))**2 ) / dh
           if( abs(dve_exact - dve2) .gt. 1.e-10 ) bad_lengths = bad_lengths + 1
           if( abs(dve_exact - dve3) .gt. 1.e-10 ) bad_lengths = bad_lengths + 1
        else if (edge2 .eq. 0) then
           dve1 = sqrt( (xVertex(i)-xEdge(edge1))**2 + (yVertex(i)-yEdge(edge1))**2 ) / dh
           dve3 = sqrt( (xVertex(i)-xEdge(edge3))**2 + (yVertex(i)-yEdge(edge3))**2 ) / dh
           if( abs(dve_exact - dve1) .gt. 1.e-10 ) bad_lengths = bad_lengths + 1
           if( abs(dve_exact - dve3) .gt. 1.e-10 ) bad_lengths = bad_lengths + 1
        else
           dve1 = sqrt( (xVertex(i)-xEdge(edge1))**2 + (yVertex(i)-yEdge(edge1))**2 ) / dh
           dve2 = sqrt( (xVertex(i)-xEdge(edge2))**2 + (yVertex(i)-yEdge(edge2))**2 ) / dh
           if( abs(dve_exact - dve1) .gt. 1.e-10 ) bad_lengths = bad_lengths + 1
           if( abs(dve_exact - dve2) .gt. 1.e-10 ) bad_lengths = bad_lengths + 1
        end if

     else

        write(6,*) " error stop, only one edge on vertex ",i
        write(6,*) edge1, edge2, edge3
        stop

     end if

  end do

  write(6,*) " bad lengths in edgesOnVertex ",bad_lengths

end subroutine check_edges_on_vertex_spacing

subroutine  verticesOnEdge_ordering( verticesOnEdge, cellsOnEdge, &
     xCell, yCell, zCell, xEdge, yEdge, zEdge, xVertex, yVertex, zVertex, nCells, nEdges, nVertices )

    implicit none

    integer :: nEdges, nCells, nVertices
    integer, dimension(2,nEdges), intent(in) :: cellsOnEdge
    integer, dimension(2,nEdges), intent(inout) :: verticesOnEdge
    real(dp), dimension(nCells) :: xCell, yCell, zCell
    real(dp), dimension(nEdges) :: xEdge, yEdge, zEdge
    real(dp), dimension(nVertices) :: xVertex, yVertex, zVertex
    
    integer :: iEdge, vertex_tmp, cell1, cell2, vertex1, vertex2
    real(dp), dimension(3) :: xc, xe, xv
    real(dp) :: areaSigned, radius, pii

   integer, parameter :: test_option = 1

    pii = 2.0_dp*asin(1.0_dp)

    if(test_option == 1) then

    do iEdge=1, nEdges

       cell1 = cellsOnEdge(1,iEdge)
       cell2 = cellsOnEdge(2,iEdge)
       vertex1 = verticesOnEdge(1,iEdge)
       vertex2 = verticesOnEdge(2,iEdge)

       xe(1) = xEdge(iEdge)
       xe(2) = yEdge(iEdge)
       xe(3) = zEdge(iEdge)

       xv(1) = xVertex(vertex2)
       xv(2) = yVertex(vertex2)
       xv(3) = zVertex(vertex2)

       radius = sqrt((xe(1)*xe(1) + xe(2)*xe(2) + xe(3)*xe(3)))

       if(cell2 .gt. cell1) then 

          xc(1) = xCell(cell2)
          xc(2) = yCell(cell2)
          xc(3) = zCell(cell2)

          xe(:) = xe(:)/radius
          xc(:) = xc(:)/radius
          xv(:) = xv(:)/radius

           areaSigned = mpas_sphere_angle(xe(1),xe(2),xe(3),xc(1),xc(2),xc(3),xv(1),xv(2),xv(3))

          if(areaSigned .le. 0.0_dp) then  ! need to switch vertices
             verticesOnEdge(1,iEdge) = vertex2
             verticesOnEdge(2,iEdge) = vertex1
          end if

       else if (cell1 .gt. cell2) then

          write(6,*) " cell 1 gt cell2 "

          xc(1) = xCell(cell1)
          xc(2) = yCell(cell1)
          xc(3) = zCell(cell1)

          areaSigned = mpas_triangle_signed_area_sphere(xe,xc,xv,radius)

          if(areaSigned .le. 0.0_dp) then  ! need to switch vertices
             verticesOnEdge(1,iEdge) = vertex2
             verticesOnEdge(2,iEdge) = vertex1
          end if

       else

          write(6,*) " cell1 and 2 are equal in verticesOnEdge ordering "
          stop

       end if

    end do

 else
       
       do iEdge = 1, nEdges

          if( cellsOnEdge(1,iEdge) .gt. cellsOnEdge(2,iEdge)) then
             vertex_tmp = verticesOnEdge(1,iEdge)
             verticesOnEdge(1,iEdge) = verticesOnEdge(2,iEdge)
             verticesOnEdge(2,iEdge) = vertex_tmp
          end if
       
       end do
       
 end if
    
contains

  real(dp) function mpas_sign_check(x1,x2,x3,y1,y2,y3,z1,z2,z3)

    implicit none
    real(dp) :: x1,x2,x3,y1,y2,y3,z1,z2,z3
    real(dp) :: a1,a2,a3,b1,b2,b3
    real(dp) :: v1,v2,v3

    real(dp) :: radius1, radius2

    radius1 = x1*x1 + x2*x2 + x3*x3

    a1 = y1-x1
    a2 = y2-x2
    a3 = y3-x3

    b1 = z1-x1
    b2 = z2-x2
    b3 = z3-x3

    v1 = a2*b3 - a3*b2
    v2 = a3*b1 - a1*b3
    v3 = a1*b2 - a2*b1

    v1 = v1 + x1
    v2 = v2 + x2
    v3 = v3 + x3

    radius2 = v1*v1 + v2*v2 + v3*v3

!    write(6,*) " sign check ",radius1,radius2
    mpas_sign_check = 1.0
    if( radius2 .lt. radius1 )  mpas_sign_check = -1.0
    
  end function mpas_sign_check
  
end subroutine verticesOnEdge_ordering

  !-----------------

  subroutine cellsOnEdge_ordering( cellsOnEdge, nEdges )

    implicit none

    integer, intent(in) :: nEdges
    integer, dimension(2,nEdges), intent(inout) :: cellsOnEdge

    integer :: cell1, cell2, iEdge, nOuterEdges

    nOuterEdges = 0

    do iEdge = 1, nEdges
       cell1 = cellsOnEdge(1,iEdge)
       cell2 = cellsOnEdge(2,iEdge)

       if((cell1 == 0) .and. (cell2 == 0)) then
          write(6,*) " cell1 and cell2 are zero for edge ",iEdge
          write(6,*) " error stop "
          stop
       end if

       if(cell1 .le. 0) then ! this is our convention - inflow is positive for all outer edges
          nOuterEdges = nOuterEdges + 1
          cellsOnEdge(1,iEdge) = 0
       end if
       if(cell2 .le. 0) then  ! switch to match our convention
          nOuterEdges = nOuterEdges + 1
          cellsOnEdge(2,iEdge) = cellsOnEdge(1,iEdge)
          cellsOnEdge(1,iEdge) = 0
       end if

    end do

    write(6,*) " number of regional mesh boundary edges ",nOuterEdges

  end subroutine cellsOnEdge_ordering

  !-----------------

  subroutine set_dcEdge( dcEdge, xCell, yCell, zCell, xEdge, yEdge, zEdge, nominalMinDc, cellsOnEdge, nEdges, nCells )

    implicit none

    integer, intent(in) :: nEdges, nCells
    integer, dimension(2,nEdges), intent(in) :: cellsOnEdge
    real(dp), dimension(nCells), intent(in) :: xCell, yCell, zCell
    real(dp), dimension(nEdges), intent(in) :: xEdge, yEdge, zEdge
    real(dp), dimension(nEdges), intent(out) :: dcEdge
    real(dp), intent(out) :: nominalMinDc

    integer :: iEdge, cell1, cell2
    real(dp) :: x1,y1,z1,x2,y2,z2

    do iEdge = 1, nEdges

       cell1 = cellsOnEdge(1,iEdge)
       cell2 = cellsOnEdge(2,iEdge)

       if(cell1 .le. 0) then  ! regional mesh

          x1 = xEdge(iEdge)
          y1 = yEdge(iEdge)
          z1 = zEdge(iEdge)

          x2 = xCell(cell2)
          y2 = yCell(cell2)
          z2 = zCell(cell2)

          dcEdge(iEdge) = 2.0_dp *  mpas_arc_length(x1,y1,z1,x2,y2,z2)

       else if(cell2 .le. 0) then

          x2 = xEdge(iEdge)
          y2 = yEdge(iEdge)
          z2 = zEdge(iEdge)

          x1 = xCell(cell1)
          y1 = yCell(cell1)
          z1 = zCell(cell1)

          dcEdge(iEdge) = 2.0_dp *  mpas_arc_length(x1,y1,z1,x2,y2,z2)

       else

          x1 = xCell(cell1)
          y1 = yCell(cell1)
          z1 = zCell(cell1)

          x2 = xCell(cell2)
          y2 = yCell(cell2)
          z2 = zCell(cell2)
          
          dcEdge(iEdge) = mpas_arc_length(x1,y1,z1,x2,y2,z2)

       end if

       ! if(iEdge .lt. 101) write(6,*) " edge, dcEdge ",iEdge,dcEdge(iEdge)

    end do

    nominalMinDc = minval(dcEdge(:))
    write(6,*) " maximum dcEdge ",maxval(dcEdge(:))
    write(6,*) " nominalMinDc   ",nominalMinDc
    
    
  end subroutine set_dcEdge

  !-----------------

  subroutine set_dvEdge( dvEdge, xVertex, yVertex, zVertex, verticesOnEdge, nEdges, nVertices )

    implicit none

    integer, intent(in) :: nEdges, nVertices
    integer, dimension(2,nEdges), intent(in) :: verticesOnEdge
    real(dp), dimension(nVertices), intent(in) :: xVertex, yVertex, zVertex
    real(dp), dimension(nEdges), intent(out) :: dvEdge

    integer :: iEdge, vertex1, vertex2
    real(dp) :: x1,y1,z1,x2,y2,z2

    do iEdge = 1, nEdges

       vertex1 = verticesOnEdge(1,iEdge)
       vertex2 = verticesOnEdge(2,iEdge)

       x1 = xVertex(vertex1)
       y1 = yVertex(vertex1)
       z1 = zVertex(vertex1)

       x2 = xVertex(vertex2)
       y2 = yVertex(vertex2)
       z2 = zVertex(vertex2)

       dvEdge(iEdge) = mpas_arc_length(x1,y1,z1,x2,y2,z2)

       ! if(iEdge .lt. 101) write(6,*) " edge, dvEdge ",iEdge,dvEdge(iEdge)

    end do
    
  end subroutine set_dvEdge

!----------------------------
     
  subroutine set_areaCell( areaCell, verticesOnCell, nEdgesOnCell, xCell, yCell, zCell, xVertex, &
                           yVertex, zVertex, nCells, nVertices, maxEdges )

    implicit none

    integer, intent(in) :: nCells, nVertices, maxEdges
    integer, dimension(maxEdges, nCells), intent(in) :: verticesOnCell
    integer, dimension(nCells), intent(in) :: nEdgesOnCell
    real(dp), dimension(nCells), intent(out) :: areaCell
    real(dp), dimension(nCells), intent(in) :: xCell, yCell, zCell
    real(dp), dimension(nVertices), intent(in) :: xVertex, yVertex, zVertex

    real(dp), dimension(3) :: p1, p2, p3
    real(dp) :: totalArea, triangleArea, sphere_radius
    integer :: vertex1, vertex2, iCell, iv, ivp1, negative_area

    sphere_radius = sqrt(xCell(1)*xCell(1) + yCell(1)*yCell(1) + zCell(1)*zCell(1))
    negative_area = 0

    do iCell = 1, nCells
       
       totalArea = 0.0_dp
       p1(1) = xCell(iCell)
       p1(2) = yCell(iCell)
       p1(3) = zCell(iCell)

       do iv = 1, nEdgesOnCell(iCell)

          ivp1 = iv+1
          if(ivp1 .gt. nEdgesOnCell(iCell)) ivp1 = 1

          vertex1 = verticesOnCell(iv,iCell)
          vertex2 = verticesOnCell(ivp1,iCell)

          p2(1) = xVertex(vertex1)
          p2(2) = yVertex(vertex1)
          p2(3) = zVertex(vertex1)

          p3(1) = xVertex(vertex2)
          p3(2) = yVertex(vertex2)
          p3(3) = zVertex(vertex2)

          triangleArea = mpas_triangle_signed_area_sphere( p1, p2, p3, sphere_radius )
          if(triangleArea .lt. 0.0) negative_area = negative_area + 1

          totalArea = totalArea + abs(triangleArea)

       end do

       areaCell(iCell) = TotalArea

       ! if(iCell .lt. 101) write(6,*) " iCell, areaCell ",iCell,areaCell(iCell)


    end do

    if(negative_area .eq. 0) then
       write(6,*) " all right-hand-rule positive areas for areaCell calculations "
    else
       write(6,*) negative_area, " negative right-hand-rule areas - problems with vertex ordering "
    end if

  end subroutine set_areaCell
  
  !--------------------------------

  subroutine set_kiteAreasOnVertex( kiteAreasOnVertex, xVertex, yVertex, zVertex, &
       xCell, yCell, zCell,       &
       xEdge, yEdge, zEdge, &
       edgesOnVertex, cellsOnEdge, cellsOnVertex, &
       nVertices, nCells, nEdges, vertexDegree )
    
    implicit none

    integer, intent(in) :: nVertices, nCells, nEdges, vertexDegree
    integer, dimension(2,nEdges), intent(in) :: cellsOnEdge
    integer, dimension(vertexDegree,nVertices), intent(in) :: cellsOnVertex
    integer, dimension(vertexDegree, nVertices), intent(in) :: edgesOnVertex

    real(dp), dimension(vertexDegree, nVertices), intent(out) :: kiteAreasOnVertex
    real(dp), dimension(nVertices), intent(in) :: xVertex, yVertex, zVertex
    real(dp), dimension(nCells), intent(in) :: xCell, yCell, zCell
    real(dp), dimension(nEdges), intent(in) :: xEdge, yEdge, zEdge

    integer :: iVertex, iCell, edge1, edge2, ie, iep1, icv
    integer :: cell1Edge1, cell1Edge2, cell2Edge1, cell2Edge2
    real(dp), dimension(3) :: pVertex, pEdge1, pEdge2, pCell
    real(dp) :: sphere_radius
    real(dp) :: kiteArea

    integer :: nKites

    sphere_radius = sqrt(xCell(1)*xCell(1) + yCell(1)*yCell(1) + zCell(1)*zCell(1))

    do iVertex=1,nVertices

       pVertex(1) = xVertex(iVertex)
       pVertex(2) = yVertex(iVertex)
       pVertex(3) = zVertex(iVertex)

       nKites = 0
       kiteArea = 0._dp

       do ie=1,vertexDegree

          iep1 = ie+1
          if(iep1 .gt. vertexDegree) iep1 = 1

          edge1 = edgesOnVertex(ie,iVertex)
          edge2 = edgesOnVertex(iep1,iVertex)

          if( (edge1 .gt. 0) .and. (edge2 .gt. 0) ) then  ! possible interior kite

             cell1Edge1 = cellsOnEdge(1,edge1)
             cell2Edge1 = cellsOnEdge(2,edge1)

             cell1Edge2 = cellsOnEdge(1,edge2)
             cell2Edge2 = cellsOnEdge(2,edge2)

             iCell = 0
             if(cell1Edge1 .eq. cell1Edge2) iCell = max(cell1Edge1, iCell)
             if(cell1Edge1 .eq. cell2Edge2) iCell = max(cell1Edge1, iCell)
             if(cell2Edge1 .eq. cell1Edge2) iCell = max(cell2Edge1, iCell)
             if(cell2Edge1 .eq. cell2Edge2) iCell = max(cell2Edge1, iCell)

             if(iCell .gt. 0) then ! interior kite, calculate area

                nKites = nKites + 1

                icv = 0
                if( iCell .eq. cellsOnVertex(1,iVertex)) icv = 1
                if( iCell .eq. cellsOnVertex(2,iVertex)) icv = 2
                if( iCell .eq. cellsOnVertex(3,iVertex)) icv = 3

                if( icv == 0 ) then
                   write(6,*) " error in finding cell on vertex ",iCell,iVertex
                   stop
                end if
                
                pEdge1(1) = xEdge(edge1)
                pEdge1(2) = yEdge(edge1)
                pEdge1(3) = zEdge(edge1)

                pEdge2(1) = xEdge(edge2)
                pEdge2(2) = yEdge(edge2)
                pEdge2(3) = zEdge(edge2)
                
                pCell(1) = xCell(iCell)
                pCell(2) = yCell(iCell)
                pCell(3) = zCell(iCell)
             
                kiteAreasOnVertex(icv,iVertex) = &
                    abs(mpas_triangle_signed_area_sphere( pVertex, pEdge1, pCell, sphere_radius )) &
                  + abs(mpas_triangle_signed_area_sphere( pVertex, pEdge2, pCell, sphere_radius ))
                kiteArea =  kiteAreasOnVertex(icv,iVertex)

             end if

          end if

       end do

       ! check for error condition (no kites - there should always be at least 1)
       ! fill in outside kite areas (currently 0) with an interior kite area

       if (nKites .eq. 0) then
          write(6,*) " no kites on vertex ",iVertex
          stop
       else
          do ie = 1,vertexDegree
             if( kiteAreasOnVertex(ie,iVertex) .le. 0._dp) kiteAreasOnVertex(ie,iVertex) = kiteArea
          end do
       end if

!       if(iVertex .lt. 10) write(6,*) iVertex, kiteAreasOnVertex(1,iVertex), &
!                      kiteAreasOnVertex(2,iVertex), kiteAreasOnVertex(3,iVertex)

    end do

  end subroutine set_kiteAreasOnVertex

  !----------------------------------

  subroutine set_areaTriangle( areaTriangle, kiteAreasOnVertex, nVertices, vertexDegree )
    
    implicit none

    integer, intent(in) :: nVertices, vertexDegree
    real(dp), dimension(nVertices), intent(out) :: areaTriangle
    real(dp), dimension(vertexDegree,nVertices), intent(in) :: kiteAreasOnVertex

    integer :: ikite, iVertex

    do iVertex = 1, nVertices

       areaTriangle(iVertex) = 0._dp

       do ikite = 1, 3
          areaTriangle(iVertex) = areaTriangle(iVertex) + kiteAreasOnVertex(ikite,iVertex)
       end do

       ! if(iVertex .lt. 50) write(6,*) " vertex, areaTriangle ",iVertex,areaTriangle(iVertex)

    end do

    
  end subroutine set_areaTriangle

  !---------------------

  subroutine  set_angleEdge( angleEdge, xEdge, yEdge, zEdge, latEdge, lonEdge, &
       verticesOnEdge, xVertex, yVertex, zVertex, nEdges, nVertices )

    implicit none

    integer, intent(in) :: nEdges, nVertices
    integer, dimension(2,nEdges), intent(in) :: verticesOnEdge
    real(dp), dimension(nEdges), intent(in) :: xEdge, yEdge, zEdge, latEdge, lonEdge
    real(dp), dimension(nVertices), intent(in) :: xVertex, yVertex, zVertex

    real(dp), dimension(nEdges), intent(out) :: angleEdge

    real(dp) :: xe, ye, ze 
    real(dp) :: xv, yv, zv
    real(dp) :: xv1, yv1, zv1, xv2, yv2, zv2
    real(dp) :: xn, yn, zn, latn, lonn
    integer :: iEdge, vertex1, vertex2

    real(dp) :: dn, dlat, radius, pii

    pii = 2.0_dp*asin(1.0_dp)
    radius = sqrt(xEdge(1)*xEdge(1) + yEdge(1)*yEdge(1) + zEdge(1)*zEdge(1))
    vertex1 = verticesOnEdge(1,1)
    vertex2 = verticesOnEdge(2,1)

    xv1 = xVertex(vertex1)
    yv1 = yVertex(vertex1)
    zv1 = zVertex(vertex1)
    xv2 = xVertex(vertex2)
    yv2 = yVertex(vertex2)
    zv2 = zVertex(vertex2)

    dn = mpas_arc_length( xv1, yv1, zv1, xv2, yv2, zv2)
    dlat = 0.5_dp*dn/radius
    ! write(6,*) " dlat = ",dlat

    do iEdge = 1, nEdges

       xe = xEdge(iEdge)/radius
       ye = yEdge(iEdge)/radius
       ze = zEdge(iEdge)/radius

       vertex1 = verticesOnEdge(1,iEdge)
       vertex2 = verticesOnEdge(2,iEdge)
       xv = xVertex(vertex2)/radius
       yv = yVertex(vertex2)/radius
       zv = zVertex(vertex2)/radius

       ! point due north of edge center, distance dn/2 away

       lonn = lonEdge(iEdge)
       latn = latEdge(iEdge) + dlat

       xn = cos(latn)*cos(lonn)
       yn = cos(latn)*sin(lonn)
       zn = sin(latn)

       angleEdge(iEdge) = mpas_sphere_angle(xe, ye, ze, xn, yn, zn, xv, yv, zv)

    end do

  end subroutine set_angleEdge
  
  !---------------------


  subroutine set_weightsOnEdge( weightsOnEdge, nEdgesOnEdge, edgesOnEdge, areaCell, angleEdge, &
       dcEdge, dvEdge,  &
       kiteAreasOnVertex, edgesOnCell, cellsOnVertex, cellsOnEdge, verticesOnCell, verticesOnEdge, &
       nEdgesOnCell, nCells, nVertices, nEdges, maxEdges, maxEdges2, vertexDegree )

    implicit none

    integer, intent(in) :: nCells, nVertices, nEdges, maxEdges, maxEdges2, vertexDegree

    integer, dimension(maxEdges2,nEdges), intent(out) :: edgesOnEdge
    integer, dimension(nEdges), intent(out) :: nEdgesOnEdge
    real(dp), dimension(maxEdges2,nEdges), intent(out) :: weightsOnEdge

    real(dp), dimension(nCells), intent(in) :: areaCell
    real(dp), dimension(vertexDegree,nVertices), intent(in) :: kiteAreasOnVertex
    real(dp), dimension(nEdges), intent(in) :: angleEdge, dcEdge, dvEdge
    integer, dimension(maxEdges, nCells), intent(in) :: edgesOnCell
    integer, dimension(nCells), intent(in) :: nEdgesOnCell
    integer, dimension(maxEdges, nCells), intent(in) :: verticesOnCell
    integer, dimension(2, nEdges) :: cellsOnEdge
    integer, dimension(vertexDegree, nVertices) :: cellsOnVertex
    integer, dimension(2, nEdges) :: VerticesOnEdge

    integer :: iEdge, iVertex, iCell, cell1, cell2, id, iv, ne, i, in, nw1, ie, ic
    integer :: vertexStart, edgeIndexCell
    integer :: vertexIndex1
    real(dp), dimension(maxEdges) :: RivCell
    real(dp), dimension(maxEdges2) :: RivWrap, weightsWrap
    integer, dimension(maxEdges2) :: edgeIndex

    real(dp), dimension(maxEdges2) :: u
    real(dp) :: vEdge, cosa, sina, pii, ve, error_max, error_abs_avg, tev2
    real(dp), parameter :: uZonal = 1.0, uMeridional = 0.0

    logical, parameter :: debug = .false.
    logical :: boundaryEdge
    integer :: nBoundaryEdges
    

    edgesOnEdge(1:maxEdges2,1:nEdges) = 0
    nEdgesOnEdge(1:nEdges) = 0
    weightsOnEdge(1:maxEdges2,1:nEdges) = 0.0_dp
    pii = 2.0_dp*asin(1.0_dp)
    nBoundaryEdges = 0

    error_max = 0.0_dp
    error_abs_avg = 0.0_dp

    do iEdge = 1,nEdges

       boundaryEdge = .false.
       cell1 = cellsOnEdge(1,iEdge)
       cell2 = cellsOnEdge(2,iEdge)
       nEdgesOnEdge(iEdge) = 0

       if( (cell1 .gt. 0) .and. (cell2 .gt. 0)) then

          if(debug) write(6,*) " setting weights for edge ",iEdge

          do ic = 1,2

          if(ic.eq.1) then
             iCell = cell1
             vertexStart = verticesOnEdge(2,iEdge)
             tev2 = -1.0_dp
             nw1 = 0
          else
             iCell = cell2
             vertexStart = verticesOnEdge(1,iEdge)
             tev2 = 1.0_dp
          end if

          ne = nEdgesOnCell(iCell)

          RivCell(1:maxEdges) = 0.0_dp
          do iv = 1, nEdgesOnCell(iCell)
             iVertex = verticesOnCell(iv,iCell)
             id = findIndex(iCell, CellsOnVertex(:,iVertex), vertexDegree)
             RivCell(iv) = kiteAreasOnVertex(id,iVertex)/areaCell(iCell)
          end do
          
          ! calculate weights for cell
          ! counterclockwise summation around cell
          
          iv = findIndex(vertexStart, verticesOnCell(:,iCell), ne)
          vertexIndex1 = iv
          if (debug) write(6,*) " vertexStart, iv, tev2 for cell ",vertexStart, iv, tev2,  ic

          weightsWrap(1:MaxEdges2) = 0.0_dp
          RivWrap(1:ne) = RivCell(1:ne)
          RivWrap(ne+1:ne+ne) = RivCell(1:ne)
          do id=iv,iv+ne-2
             weightsWrap(id) = (sum(RivWrap(iv:id)) - 0.5_dp)*tev2
             weightsOnEdge(nw1+id-iv+1,iEdge) = weightsWrap(id)
          end do
          
          ! find iEdge index in edgesOnCell for cell
          edgeIndexCell = findIndex(iEdge, edgesOnCell(:,iCell),nEdgesOnCell(iCell))
          ne = nEdgesOnCell(iCell)
          edgeIndex(1:ne) = edgesOnCell(1:ne,iCell)
          edgeIndex(ne+1:2*ne) = edgeIndex(1:ne)
          do i=1,ne-1
             in = i + nw1
             ie = edgeIndex(i+edgeIndexCell)
             EdgesOnEdge(in,iEdge) = ie
             weightsOnEdge(in,iEdge) = weightsOnEdge(in,iEdge)*dvEdge(ie)/dcEdge(iEdge)
             ! check for inflow - switch sign on weight accordingly
             if(cellsOnEdge(2,ie) .eq. iCell) then
                weightsOnEdge(in,iEdge) = -weightsOnEdge(in,iEdge)
             end if
          end do

          nw1 = ne - 1
          nEdgesOnEdge(iEdge) = nEdgesOnEdge(iEdge) + nw1

          end do
          
       else

          if(debug) write(6,*) " edge is on regional domain boundary "
          boundaryEdge = .true.
          nBoundaryEdges = nBoundaryEdges + 1
          nEdgesOnEdge(iEdge) = 0

       end if

       end do

       !  test to see if weights are correct

       do iEdge=1,nEdges
       
       if( nEdgesOnEdge(iEdge) .gt. 0) then
          
       vEdge = 0.0_dp
       do i=1, nEdgesOnEdge(iEdge)

          ! veolicities on edges contributing to iEdge.  U is norm component
          cosa = cos(angleEdge(EdgesOnEdge(i,iEdge)))
          sina = sin(angleEdge(EdgesOnEdge(i,iEdge)))
          u(i) =  uZonal*cosa + uMeridional*sina
          vEdge = vEdge + u(i)*weightsOnEdge(i,iEdge)
       end do

       cosa = cos(angleEdge(iEdge))
       sina = sin(angleEdge(iEdge))
       ve = -uZonal*sina + uMeridional*cosa
       error_max = max(error_max,abs(vEdge - ve))
       error_abs_avg = error_abs_avg + abs(vEdge - ve)

       end if

       end do
       
       error_abs_avg = error_abs_avg/(nEdges-nBoundaryEdges)

    write(6,*) " number of boundary edges ",nBoundaryEdges
    write(6,*) " max reconstruction error ",error_max
    write(6,*) " avg abs reconstruction error ",error_abs_avg

  contains

  integer function findIndex( index, indices, nIndices )

    ! find index in list of indices, return index location

    implicit none
    integer :: index, nIndices
    integer, dimension(nedges) :: Indices
    integer :: i

    findIndex = 0
    do i=1,nIndices
       if(indices(i) .eq. index) findIndex = i
    end do

  end function findIndex
  
  end subroutine set_weightsOnEdge

  !------------------------------------

  subroutine set_meshDensity( meshDensity, mapFactorCell, nCells )

    implicit none

    integer, intent(in) :: nCells
    real(dp), dimension(nCells), intent(in) :: mapFactorCell
    real(dp), dimension(nCells), intent(out) :: meshDensity

    integer iCell
    real(dp) :: dhmin

    dhmin = 1.0_dp/mapFactorCell(1)
    do iCell = 2, nCells
       dhmin = min(dhmin,1.0_dp/mapFactorCell(iCell))
    end do

    do iCell = 1, nCells
       meshDensity(iCell) = (dhmin*mapFactorCell(iCell))**4.0_dp
    end do

    write(6,*) " min and max meshDensity ",minval(meshDensity(:)),maxval(meshDensity(:))

  end subroutine set_meshDensity
  
  !-----------------------------------------------

  subroutine set_indexToID( indexToCellID, indexToVertexID, indexToEdgeID, nCells, nVertices, nEdges )

    implicit none

    integer, intent(in) :: nCells, nVertices, nEdges
    integer, dimension(nCells), intent(out) :: indexToCellID
    integer, dimension(nVertices), intent(out) :: indexToVertexID
    integer, dimension(nEdges), intent(out) :: indexToEdgeID

    integer :: i

    do i=1,nCells
       indexToCellID(i) = i
    end do
    do i=1,nVertices
       indexToVertexID(i) = i
    end do
    do i=1,nEdges
       indexToEdgeID(i) = i
    end do

  end subroutine set_indexToID

  !---------------

  subroutine center_coordinate(pCell,pVertex,pEdge,nCells,nVertices,nEdges)

    implicit none

    integer, intent(in) :: nCells, nVertices, nEdges
    real(dp), dimension(nCells) :: pCell
    real(dp), dimension(nVertices) :: pVertex
    real(dp), dimension(nEdges) :: pEdge

    real(dp) :: pshift

    pshift = 0.5_dp*(minval(pCell(:))+maxval(pCell(:)))
    write(6,*) " shift ",pshift


    pCell(:) = pCell(:)-pshift
    pVertex(:) = pVertex(:)-pshift
    pEdge(:) = pEdge(:)-pshift

  end subroutine center_coordinate
  
  !---------------

  subroutine  write_graph_info( cellsOnCell, cellsOnEdge, nEdgesOnCell, nCells, nEdges, maxEdges )

    implicit none

    integer, intent(in) :: nCells, nEdges, maxEdges
    integer, dimension(maxEdges,nCells), intent(in) :: cellsOnCell
    integer, dimension(nCells), intent(in) :: nEdgesOnCell
    integer, dimension(2,nEdges), intent(in) :: cellsOnEdge

    integer :: iCell, iEdge, nout, edgeCells, nEdgesInterior
    integer, dimension(maxEdges) :: cellsOut

    write(6,*) " creating graph.info file "

    open(unit=10,file="graph.info",form="formatted",status="replace")

    write(6,*) " min and max cellsOnEdge ",minval(cellsOnEdge(:,:)),maxval(cellsOnEdge(:,:))
    nEdgesInterior = 0
    do iEdge=1,nEdges
       if( (cellsOnEdge(1,iEdge).ne.0) .and. (cellsOnEdge(2,iEdge).ne.0)) nEdgesInterior = nEdgesInterior + 1
    end do
    write(6,*) " nEdges, nEdgesInterior ",nEdges,nEdgesInterior

    write(10,"(2i10)") nCells, nEdgesInterior

    edgeCells = 0
    
    do iCell=1,nCells
       nout = 0
       do iEdge=1,nEdgesOnCell(iCell)
          if(cellsOnCell(iEdge,iCell) .gt. 0) then
             nout = nout+1
             cellsOut(nout) = cellsOnCell(iEdge,iCell)
          end if
       end do
       write(10,"(10i10)") cellsOut(1:nout)
       if(nout .lt. nEdgesOnCell(iCell)) edgeCells = edgeCells+1
    end do

    close(unit=10)

    write(6,*) ' finished creating graph.info file, edgeCells = ',edgeCells

  end subroutine write_graph_info
  
  subroutine set_outside_index( cellsOnCell, EdgesOnCell, verticesOnCell, cellsOnVertex, cellsOnEdge, &
                                edgesOnVertex, verticesOnEdge, nCells, nEdges, nVertices,             &
                                maxEdges, vertexDegree )

    implicit none

    integer, intent(in) :: nCells, nEdges, nVertices, maxEdges, vertexDegree
    integer, dimension(maxEdges, nCells) :: cellsOnCell, edgesOnCell, verticesOnCell
    integer, dimension(vertexDegree, nVertices) :: cellsOnVertex, edgesOnVertex
    integer, dimension(2,nEdges) :: cellsOnEdge, verticesOnEdge

    integer :: iCell, iEdge, iVertex, j

    do iCell = 1, nCells
       do j=1,maxEdges
          cellsOnCell(j,iCell) = max(0,cellsOnCell(j,iCell))
          edgesOnCell(j,iCell) = max(0,edgesOnCell(j,iCell))
          verticesOnCell(j,iCell) = max(0,verticesOnCell(j,iCell))
       end do
    end do
    do iVertex = 1, nVertices
       do j=1,vertexDegree
          cellsOnVertex(j,iVertex) = max(0,cellsOnVertex(j,iVertex))
          edgesOnVertex(j,iVertex) = max(0,edgesOnVertex(j,iVertex))
       end do
    enddo
    do iEdge = 1, nEdges
       do j=1,2
          cellsOnEdge(j,iEdge) = max(0,cellsOnEdge(j,iEdge))
          verticesOnEdge(j,iEdge) = max(0,verticesOnEdge(j,iEdge))
       end do
    end do

  end subroutine set_outside_index
  
end module ph_utils
