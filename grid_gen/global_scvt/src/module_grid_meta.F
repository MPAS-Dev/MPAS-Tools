module grid_meta

   contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE WRITE_GRID
   !
   ! Create metadata for mesh and write out the complete grid information.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine write_grid(rlat, rlon, n, vclat, vclon, nvc, alist, clist)
   
      use data_types
      use grid_constants
      use sphere_utilities
      use write_netcdf
      use scvt
   
      implicit none
   
      integer, intent(in) :: n, nvc
      real(kind=RKIND), dimension(n), intent(inout) :: rlat, rlon
      real(kind=RKIND), dimension(nvc), intent(inout) :: vclat, vclon
      type (adjacency_list) :: alist, clist
   
      integer, parameter :: maxEdges = 10

      type tree_pointer_list
         integer :: stage
         type (binary_tree), pointer :: p
      end type tree_pointer_list
   
      integer :: i, ii, j, js, k, prev_edge, next_edge, vtx, nObtuse, jj, iEdge, iFlag
      integer :: nCells, nEdges, nVertices, nVertLevels, nTracers
      integer, dimension(3) :: cellsOnTriangle
      integer, allocatable, dimension(:) :: permutation, indexToCellID, indexToEdgeID, indexToVertexID
      integer, allocatable, dimension(:) :: nEdgesOnCell, nEdgesOnEdge, nEdgesOnVertex, nCellsOnVertex
      integer, allocatable, dimension(:) :: tempEdgesOnEdge
      integer, allocatable, dimension(:,:) :: cellsOnEdge, edgesOnCell, verticesOnCell, cellsOnCell, &
                                              verticesOnEdge, edgesOnEdge, edgesOnVertex, cellsOnVertex
      integer, allocatable, dimension(:) :: isObtuse
      real(kind=RKIND), allocatable, dimension(:,:) :: weightsOnEdge, kiteAreasOnVertex
      integer :: temp
      logical :: tdrtest = .true.
      real(kind=RKIND) :: sum_r, area, r, s, de, rtmp
      real(kind=RKIND), allocatable, dimension(:) :: latCell, lonCell, latEdge, lonEdge, angleEdge, latVertex, lonVertex, &
                                         lat1Edge, lon1Edge, lat2Edge, lon2Edge, dvEdge, dv1Edge, dv2Edge, dcEdge, &
                                         areaCell, areaTriangle, fEdge, fVertex, h_s, u_sbr
      real(kind=RKIND), allocatable, dimension(:,:,:) :: u, v, h, vh, circulation, vorticity, ke
      real(kind=RKIND), allocatable, dimension(:,:,:,:) :: tracers
      real(kind=RKIND), allocatable, dimension(:) :: xCell, yCell, zCell, xEdge, yEdge, zEdge, xVertex, yVertex, zVertex, meshDensity
      type (geo_point) :: vertex1GP, vertex2GP, cell1GP, cell2GP, cell3GP, edgeGP, edgeGP_prev, edgeGP_next, pCell
      type (geo_point) :: center
      type (geo_point), allocatable, dimension(:) :: points
      type (binary_tree), pointer :: treeRoot
      type (tree_pointer_list), allocatable, dimension(:) :: cursor
   

      !
      ! Compute number of cells
      !
      nCells = alist % nNodes
   
   
      !
      ! Compute number of edges
      !
      nEdges = alist % nNeighbors / 2

      
      !
      ! Compute number of vertices
      !
      nVertices = nEdges - nCells + 2      ! Euler's formula
   
   
      allocate(latCell(nCells)) 
      allocate(lonCell(nCells))
      allocate(meshDensity(nCells))
      allocate(xCell(nCells))
      allocate(yCell(nCells))
      allocate(zCell(nCells))
      allocate(latEdge(nEdges))
      allocate(lonEdge(nEdges))
      allocate(xEdge(nEdges))
      allocate(yEdge(nEdges))
      allocate(zEdge(nEdges))
      allocate(dvEdge(nEdges))
      allocate(dv1Edge(nEdges))
      allocate(dv2Edge(nEdges))
      allocate(dcEdge(nEdges))
      allocate(areaCell(nCells))
      allocate(areaTriangle(nVertices))
      allocate(angleEdge(nEdges))
      allocate(u_sbr(nEdges))
      allocate(latVertex(nVertices))
      allocate(lonVertex(nVertices))
      allocate(xVertex(nVertices))
      allocate(yVertex(nVertices))
      allocate(zVertex(nVertices))
      allocate(lat1Edge(nEdges))
      allocate(lon1Edge(nEdges))
      allocate(lat2Edge(nEdges))
      allocate(lon2Edge(nEdges))
      allocate(cellsOnEdge(2,nEdges))
      allocate(isObtuse(nVertices))
      allocate(edgesOnCell(maxEdges,nCells))
      allocate(verticesOnCell(maxEdges,nCells))
      allocate(verticesOnEdge(2,nEdges))
      allocate(edgesOnEdge(2*maxEdges,nEdges))
      allocate(edgesOnVertex(3,nVertices))
      allocate(cellsOnVertex(3,nVertices))
      allocate(tempEdgesOnEdge(2*maxEdges))
      allocate(weightsOnEdge(2*maxEdges,nEdges))
      allocate(kiteAreasOnVertex(3,nVertices))
      allocate(cellsOnCell(maxEdges,nCells))
      allocate(nEdgesOnCell(nCells))
      allocate(nEdgesOnEdge(nEdges))
      allocate(nEdgesOnVertex(nVertices))
      allocate(nCellsOnVertex(nVertices))
      allocate(cursor(nEdges))
   
   
      !
      ! Determine lat/lon for each cell
      !
      do i=1,nCells
         latCell(i) = rlat(i)
         lonCell(i) = rlon(i)
      end do


      !
      ! Determine lat/lon for each vertex
      !
      do i=1,nVertices
         latVertex(i) = vclat(i)
         lonVertex(i) = vclon(i)
      end do 

   
      !
      ! Construct a list of edges (actually, a BST)
      !
      nullify(treeRoot)
      do i=1,nCells
         do j = clist % start(i), clist % start(i) + clist % len(i) - 2 
            js = j+1
            call insert_edge_to_tree(i, &
                                     clist % neighbor(j), clist % neighbor(js), &
                                     vclat(clist % neighbor(j)), vclon(clist % neighbor(j)),   &
                                     vclat(clist % neighbor(js)), vclon(clist % neighbor(js)), &
                                     treeRoot)
         end do
   
         js = clist % start(i)
         call insert_edge_to_tree(i, &
                                  clist % neighbor(j), clist % neighbor(js), &
                                  vclat(clist % neighbor(j)), vclon(clist % neighbor(j)),   &
                                  vclat(clist % neighbor(js)), vclon(clist % neighbor(js)), &
                                  treeRoot)
      end do


      !
      ! Determine which cells are on each edge, and lat/lon of edge endpoints
      !
      i = 1
      j = 1
      cursor(i) % p => treeRoot
      cursor(i) % stage = 0
      do while (i > 0)
         if (associated(cursor(i) % p % left) .and. cursor(i) % stage < 1) then
            cursor(i) % stage = 1
            i = i + 1
            cursor(i) % p => cursor(i-1) % p % left
            cursor(i) % stage = 0
            cycle
         end if
         if (cursor(i) % stage < 2) then
            cellsOnEdge(1,j) = cursor(i) % p % node1
            cellsOnEdge(2,j) = cursor(i) % p % node2
            lat1Edge(j) = cursor(i) % p % lat1
            lon1Edge(j) = cursor(i) % p % lon1
            lat2Edge(j) = cursor(i) % p % lat2
            lon2Edge(j) = cursor(i) % p % lon2
            verticesOnEdge(1,j) = cursor(i) % p % vertex1
            verticesOnEdge(2,j) = cursor(i) % p % vertex2
            j = j + 1
            cursor(i) % stage = 2
         end if
         if (associated(cursor(i) % p % right) .and. cursor(i) % stage < 3) then
            cursor(i) % stage = 3
            i = i + 1
            cursor(i) % p => cursor(i-1) % p % right
            cursor(i) % stage = 0
            cycle
         end if
         i = i - 1
      end do

     
      !
      ! Compute distance between cell centroids for each edge
      !
      do i=1,nEdges
         cell1GP % lat = latCell(cellsOnEdge(1,i))
         cell1GP % lon = lonCell(cellsOnEdge(1,i))
         cell2GP % lat = latCell(cellsOnEdge(2,i))
         cell2GP % lon = lonCell(cellsOnEdge(2,i))
         dcEdge(i) = sphere_distance(cell1GP, cell2GP, 1.0_RKIND)
      end do

      
      !
      ! Determine which edges are on each cell
      !
      do i=1,nCells
         nEdgesOnCell(i) = 0
      end do

      do j=1,nEdges
         nEdgesOnCell(cellsOnEdge(1,j)) = nEdgesOnCell(cellsOnEdge(1,j)) + 1
         if (nEdgesOnCell(cellsOnEdge(1,j)) > maxEdges) then
            write(0,*) 'Houston, we have an edge problem.'
            stop
         end if
         edgesOnCell(nEdgesOnCell(cellsOnEdge(1,j)),cellsOnEdge(1,j)) = j

         nEdgesOnCell(cellsOnEdge(2,j)) = nEdgesOnCell(cellsOnEdge(2,j)) + 1
         if (nEdgesOnCell(cellsOnEdge(2,j)) > maxEdges) then
            write(0,*) 'Houston, we have an edge problem.'
            stop
         end if
         edgesOnCell(nEdgesOnCell(cellsOnEdge(2,j)),cellsOnEdge(2,j)) = j
      end do
      do i=1,nCells
         do j=nEdgesOnCell(i)+1,maxEdges
            edgesOnCell(j,i) = edgesOnCell(nEdgesOnCell(i),i)
         end do
      end do


      !
      ! Determine which cells are on each cell
      !
      do i=1,nCells
         k = 1
         do j = alist % start(i), alist % start(i) + alist % len(i) - 1 
            cellsOnCell(k,i) = alist % neighbor(j)
            k = k + 1
         end do
         do j=k,maxEdges
            cellsOnCell(j,i) = alist % neighbor(alist % start(i) + alist % len(i) - 1)
         end do
      end do

      ! Re-order edgesOnCell to be in the same order as cellsOnCell
      do i=1,nCells
         do j=1,nEdgesOnCell(i)
            if (cellsOnEdge(1,edgesOnCell(j,i)) /= cellsOnCell(j,i) .and. cellsOnEdge(2,edgesOnCell(j,i)) /= cellsOnCell(j,i)) then
               ii = 0
               do k=j+1,nEdgesOnCell(i)
                  if (cellsOnEdge(1,edgesOnCell(k,i)) == cellsOnCell(j,i) .or. cellsOnEdge(2,edgesOnCell(k,i)) == cellsOnCell(j,i)) then
                     ! Swap edgesOnCell(k,i) and edgesOnCell(j,i)
                     ii = edgesOnCell(j,i)
                     edgesOnCell(j,i) = edgesOnCell(k,i)
                     edgesOnCell(k,i) = ii
                     exit
                  end if
               end do
if (ii == 0) then
   write(0,*) 'We didn''t find an edge to match cellsOnCell'
   stop
end if
            end if
         end do
      end do
   
   
      !
      ! Compute lat/lon of each edge
      !
      do j=1,nEdges
         vertex1GP % lat = lat1Edge(j)
         vertex1GP % lon = lon1Edge(j)
         vertex2GP % lat = lat2Edge(j)
         vertex2GP % lon = lon2Edge(j)
         cell1GP % lat = latCell(cellsOnEdge(1,j))
         cell1GP % lon = lonCell(cellsOnEdge(1,j))
         cell2GP % lat = latCell(cellsOnEdge(2,j))
         cell2GP % lon = lonCell(cellsOnEdge(2,j))
         call compute_edge_latlon(cell1GP, cell2GP, vertex1GP, vertex2GP, edgeGP)

         dvEdge(j) = sphere_distance(vertex1GP, vertex2GP, 1.0_RKIND)
         rtmp = (vertex2GP % lat - vertex1GP % lat) / dvEdge(j)
         if (rtmp > 1.0_RKIND) rtmp = 1.0_RKIND
         if (rtmp < -1.0_RKIND) rtmp = -1.0_RKIND
         rtmp = acos(rtmp)
         angleEdge(j) = meridian_angle(edgeGP, vertex2GP) 
         angleEdge(j) = sign(rtmp, angleEdge(j))
         if (angleEdge(j) > pii) angleEdge(j) = angleEdge(j) - 2.0*pii
         if (angleEdge(j) < -pii) angleEdge(j) = angleEdge(j) + 2.0*pii
         latEdge(j) = edgeGP % lat 
         lonEdge(j) = edgeGP % lon 

         ! Compute u for solid body rotation
         u_sbr(j) = (sin(vertex2GP % lat) - sin(vertex1GP % lat)) / dvEdge(j)

         ! Make sure order of vertices on edge is correct
!         if (is_flipped_vertex_order2(edgeGP, cell2GP, vertex2GP)) then
         if (is_flipped_vertex_order(cell1GP % lat, cell1GP % lon, &
                                     cell2GP % lat, cell2GP % lon, &
                                     vertex1GP % lat, vertex1GP % lon, &
                                     vertex2GP % lat, vertex2GP % lon)) then
            temp = verticesOnEdge(1,j)
            verticesOnEdge(1,j) = verticesOnEdge(2,j)
            verticesOnEdge(2,j) = temp
            u_sbr(j) = -1.0_RKIND*u_sbr(j)
            angleEdge(j) = angleEdge(j) + pii
            if (angleEdge(j) > pii) angleEdge(j) = angleEdge(j) - 2.0*pii
            if (angleEdge(j) < -pii) angleEdge(j) = angleEdge(j) + 2.0*pii
         end if

      end do

 
      !
      ! Determine which vertices are on each cell
      !
      do i=1,nCells
         if (i == cellsOnEdge(1,edgesOnCell(1,i))) then
            verticesOnCell(1,i) = verticesOnEdge(1,edgesOnCell(1,i))
         else if (i == cellsOnEdge(2,edgesOnCell(1,i))) then
            verticesOnCell(1,i) = verticesOnEdge(2,edgesOnCell(1,i))
         else
write(0,*) 'THIS EDGE SHOULDN''T BELONG TO THIS CELL'
         end if
         do j=1,nEdgesOnCell(i)-1
            if (verticesOnEdge(1,edgesOnCell(j,i)) == verticesOnCell(j,i)) then
               verticesOnCell(j+1,i) = verticesOnEdge(2,edgesOnCell(j,i))
            else if (verticesOnEdge(2,edgesOnCell(j,i)) == verticesOnCell(j,i)) then
               verticesOnCell(j+1,i) = verticesOnEdge(1,edgesOnCell(j,i))
            else
write(0,*) 'Broken chain of vertex-edge-vertex.'
stop
            end if
         end do
      end do


      !
      ! Determine which edges are incident with each vertex
      !
      nEdgesOnVertex(:) = 0
      do j=1,nEdges
         nEdgesOnVertex(verticesOnEdge(1,j)) = nEdgesOnVertex(verticesOnEdge(1,j)) + 1
         if (nEdgesOnVertex(verticesOnEdge(1,j)) > 3) then
            write(0,*) 'We have too many edges incident with vertex ',verticesOnEdge(1,j)
            stop
         end if
         edgesOnVertex(nEdgesOnVertex(verticesOnEdge(1,j)),verticesOnEdge(1,j)) = j

         nEdgesOnVertex(verticesOnEdge(2,j)) = nEdgesOnVertex(verticesOnEdge(2,j)) + 1
         if (nEdgesOnVertex(verticesOnEdge(2,j)) > 3) then
            write(0,*) 'We have too many edges incident with vertex ',verticesOnEdge(2,j)
            stop
         end if
         edgesOnVertex(nEdgesOnVertex(verticesOnEdge(2,j)),verticesOnEdge(2,j)) = j
      end do


      !
      ! Determine which cells share each vertex
      !
      nCellsOnVertex(:) = 0
      do i=1,nCells
         do j=1,nEdgesOnCell(i)
            nCellsOnVertex(verticesOnCell(j,i)) = nCellsOnVertex(verticesOnCell(j,i)) + 1
            if (nCellsOnVertex(verticesOnCell(j,i)) > 3) then
               write(0,*) 'We have too many cells sharing vertex ', verticesOnCell(j,i)
               stop
            end if
            cellsOnVertex(nCellsOnVertex(verticesOnCell(j,i)),verticesOnCell(j,i)) = i
         end do
      end do


      !
      ! Determine which edges "neighbor" each edge
      !
      do j=1,nEdges
         allocate(points(nEdgesOnCell(cellsOnEdge(1,j))))
         allocate(permutation(nEdgesOnCell(cellsOnEdge(1,j))))
         js = 1
         points(js) % lat = latEdge(j)
         points(js) % lon = lonEdge(j)
         permutation(js) = j
         js = js + 1
         center % lat = latCell(cellsOnEdge(1,j))
         center % lon = lonCell(cellsOnEdge(1,j))

         do k=1,nEdgesOnCell(cellsOnEdge(1,j))
            if (edgesOnCell(k,cellsOnEdge(1,j)) /= j) then
               nEdgesOnEdge(j) = nEdgesOnEdge(j) + 1
               edgesOnEdge(nEdgesOnEdge(j),j) = edgesOnCell(k,cellsOnEdge(1,j))
               points(js) % lat = latEdge(edgesOnCell(k,cellsOnEdge(1,j)))
               points(js) % lon = lonEdge(edgesOnCell(k,cellsOnEdge(1,j)))
               permutation(js) = edgesOnCell(k,cellsOnEdge(1,j))
               js = js + 1
            end if
         end do
         call order_points_ccw(center, nEdgesOnCell(cellsOnEdge(1,j)), points, permutation)
         do k=2,nEdgesOnCell(cellsOnEdge(1,j))
            edgesOnEdge(k-1,j) = permutation(k)
         end do
         deallocate(points)
         deallocate(permutation)

         allocate(points(nEdgesOnCell(cellsOnEdge(2,j))))
         allocate(permutation(nEdgesOnCell(cellsOnEdge(2,j))))
         js = 1
         points(js) % lat = latEdge(j)
         points(js) % lon = lonEdge(j)
         permutation(js) = j
         js = js + 1
         center % lat = latCell(cellsOnEdge(2,j))
         center % lon = lonCell(cellsOnEdge(2,j))

         do k=1,nEdgesOnCell(cellsOnEdge(2,j))
            if (edgesOnCell(k,cellsOnEdge(2,j)) /= j) then
               nEdgesOnEdge(j) = nEdgesOnEdge(j) + 1
               edgesOnEdge(nEdgesOnEdge(j),j) = edgesOnCell(k,cellsOnEdge(2,j))
               points(js) % lat = latEdge(edgesOnCell(k,cellsOnEdge(2,j)))
               points(js) % lon = lonEdge(edgesOnCell(k,cellsOnEdge(2,j)))
               permutation(js) = edgesOnCell(k,cellsOnEdge(2,j))
               js = js + 1
            end if
         end do
         call order_points_ccw(center, nEdgesOnCell(cellsOnEdge(2,j)), points, permutation)
         do k=2,nEdgesOnCell(cellsOnEdge(2,j))
            edgesOnEdge(nEdgesOnCell(cellsOnEdge(1,j))+k-2,j) = permutation(k)
         end do
         deallocate(points)
         deallocate(permutation)
      end do


      !
      ! Compute area of each cell
      !
      do i=1,nCells
         cell1GP % lat = latCell(i) 
         cell1GP % lon = lonCell(i) 
         areaCell(i) = 0.0
         do j=1,nEdgesOnCell(i)-1
            vertex1GP % lat = latVertex(verticesOnCell(j,i)) 
            vertex1GP % lon = lonVertex(verticesOnCell(j,i)) 
            vertex2GP % lat = latVertex(verticesOnCell(j+1,i)) 
            vertex2GP % lon = lonVertex(verticesOnCell(j+1,i)) 
            areaCell(i) = areaCell(i) + triangle_area(cell1GP, vertex1GP, vertex2GP, 1.0_RKIND)
         end do
         vertex1GP % lat = latVertex(verticesOnCell(j,i)) 
         vertex1GP % lon = lonVertex(verticesOnCell(j,i)) 
         vertex2GP % lat = latVertex(verticesOnCell(1,i)) 
         vertex2GP % lon = lonVertex(verticesOnCell(1,i)) 
         areaCell(i) = areaCell(i) + triangle_area(cell1GP, vertex1GP, vertex2GP, 1.0_RKIND)
      end do


      !
      ! Compute area of triangles associated with each vertex
      !
      do i=1,nVertices
         cell1GP % lat = latCell(cellsOnVertex(1,i))
         cell1GP % lon = lonCell(cellsOnVertex(1,i))
         cell2GP % lat = latCell(cellsOnVertex(2,i))
         cell2GP % lon = lonCell(cellsOnVertex(2,i))
         cell3GP % lat = latCell(cellsOnVertex(3,i))
         cell3GP % lon = lonCell(cellsOnVertex(3,i))
         areaTriangle(i) = triangle_area(cell1GP, cell2GP, cell3GP, 1.0_RKIND)
      end do

      !
      ! Test to see if any of the triangles are obtuse
      !
      nObtuse = 0
      do i=1,nVertices
         cell1GP % lat = latCell(cellsOnVertex(1,i))
         cell1GP % lon = lonCell(cellsOnVertex(1,i))
         cell2GP % lat = latCell(cellsOnVertex(2,i))
         cell2GP % lon = lonCell(cellsOnVertex(2,i))
         cell3GP % lat = latCell(cellsOnVertex(3,i))
         cell3GP % lon = lonCell(cellsOnVertex(3,i))
         isObtuse(i) = obtuse(cell1GP, cell2GP, cell3GP)
         if(isObtuse(i).gt.0) nObtuse = nObtuse + 1
      end do
      write(6,*) ' number of obtuse triangles ', nObtuse


      kiteAreasOnVertex(:,:) = -1.0_RKIND

      !
      ! Compute weights used in tangential velocity reconstruction
      !
      do j=1,nEdges
         cell1GP % lat = latCell(cellsOnEdge(1,j))
         cell1GP % lon = lonCell(cellsOnEdge(1,j))
         cell2GP % lat = latCell(cellsOnEdge(2,j))
         cell2GP % lon = lonCell(cellsOnEdge(2,j))
         de = dcEdge(j)
         prev_edge = j
         sum_r = 0.0
         do i=1,nEdgesOnCell(cellsOnEdge(1,j))-1
            next_edge = edgesOnEdge(i,j)
            if ((verticesOnEdge(1,prev_edge) == verticesOnEdge(1,next_edge)) .or. &
                (verticesOnEdge(1,prev_edge) == verticesOnEdge(2,next_edge))) then
               vtx = verticesOnEdge(1,prev_edge)
            else if ((verticesOnEdge(2,prev_edge) == verticesOnEdge(1,next_edge)) .or. &
                     (verticesOnEdge(2,prev_edge) == verticesOnEdge(2,next_edge))) then
               vtx = verticesOnEdge(2,prev_edge)
            else
               write(0,*) 'Somehow these edges don''t share a vertex',j
               write(0,*) 'Cells are ',cellsOnEdge(1,j),cellsOnEdge(2,j)
               write(0,*) 'Edge ',prev_edge,' has vertices ',verticesOnEdge(1,prev_edge),verticesOnEdge(2,prev_edge)
               write(0,*) 'Edge ',next_edge,' has vertices ',verticesOnEdge(1,next_edge),verticesOnEdge(2,next_edge)
               write(0,*) 'edgesOnEdge()=',edgesOnEdge(:,j)
               write(0,*) 'verticesOnCell(1) = ',verticesOnCell(:,cellsOnEdge(1,j))
               write(0,*) 'edgesOnCell(1) = ',edgesOnCell(:,cellsOnEdge(1,j))
               write(0,*) 'verticesOnEdge(12821)=',verticesOnEdge(:,12821)
               write(0,*) 'verticesOnEdge(12823)=',verticesOnEdge(:,12823)
               write(0,*) 'verticesOnEdge(13243)=',verticesOnEdge(:,13243)
               write(0,*) 'verticesOnEdge(13245)=',verticesOnEdge(:,13245)
               write(0,*) 'verticesOnEdge(13448)=',verticesOnEdge(:,13448)
               write(0,*) 'verticesOnEdge(13449)=',verticesOnEdge(:,13449)
               stop
            end if
            cell1GP % lat = latCell(cellsOnEdge(1,j))
            cell1GP % lon = lonCell(cellsOnEdge(1,j))
            vertex1GP % lat = latVertex(vtx)
            vertex1GP % lon = lonVertex(vtx)
            edgeGP_prev % lat = latEdge(prev_edge)
            edgeGP_prev % lon = lonEdge(prev_edge)
            edgeGP_next % lat = latEdge(next_edge)
            edgeGP_next % lon = lonEdge(next_edge)

            if(isObtuse(vtx).eq.0) then
               area =        abs(triangle_area(cell1GP, vertex1GP, edgeGP_prev, 1.0_RKIND))
               area = area + abs(triangle_area(cell1GP, vertex1GP, edgeGP_next, 1.0_RKIND))
            else
               if(cellsOnVertex(isObtuse(vtx),vtx).eq.cellsOnEdge(1,j)) then
                 iFlag = 0
                 do ii=1,3
                  jj = edgesOnVertex(ii,vtx)
                  if(jj.ne.prev_edge.and.jj.ne.next_edge) then
                      write(6,*) jj, prev_edge, next_edge
                      if(iFlag.eq.1) then
                         write(6,*) ' can not be true '
                         stop
                      endif
                      iEdge = jj
                      iFlag = 1
                   endif
                  enddo
                  edgeGP % lat = latEdge(iEdge)
                  edgeGP % lon = lonEdge(iEdge)
                  area =        abs(triangle_area(cell1GP, edgeGP, edgeGP_prev, 1.0_RKIND))
                  area = area + abs(triangle_area(cell1GP, edgeGP, edgeGP_next, 1.0_RKIND))
               else
                  area = abs(triangle_area(cell1GP, edgeGP_prev, edgeGP_next, 1.0_RKIND))
               endif
            endif

            do ii=1,3
               if (cellsOnEdge(1,j) == cellsOnVertex(ii,vtx)) then
                  kiteAreasOnVertex(ii,vtx) = area
                  exit
               end if
            end do

            if(.not.tdrtest) then
              write(6,*) ' not tdrtest'
              r = area / areaCell(cellsOnEdge(1,j))
              sum_r = sum_r + r
              if (cellsOnEdge(1,j) == cellsOnEdge(1,edgesOnEdge(i,j))) then
                 s = 1.0_RKIND
              else
                 s = -1.0_RKIND
              end if
              weightsOnEdge(i,j) = s*(0.5-sum_r)*dvEdge(edgesOnEdge(i,j))/de
            endif

            prev_edge = next_edge
         end do

         prev_edge = j
         sum_r = 0.0
         do i=nEdgesOnCell(cellsOnEdge(1,j)),nEdgesOnEdge(j)
            next_edge = edgesOnEdge(i,j)
            if ((verticesOnEdge(1,prev_edge) == verticesOnEdge(1,next_edge)) .or. &
                (verticesOnEdge(1,prev_edge) == verticesOnEdge(2,next_edge))) then
               vtx = verticesOnEdge(1,prev_edge)
            else if ((verticesOnEdge(2,prev_edge) == verticesOnEdge(1,next_edge)) .or. &
                     (verticesOnEdge(2,prev_edge) == verticesOnEdge(2,next_edge))) then
               vtx = verticesOnEdge(2,prev_edge)
            else
               write(0,*) 'Somehow these edges don''t share a vertex'
               write(0,*) 'Edge ',prev_edge,' has vertices ',verticesOnEdge(1,prev_edge),verticesOnEdge(2,prev_edge)
               write(0,*) 'Edge ',next_edge,' has vertices ',verticesOnEdge(1,next_edge),verticesOnEdge(2,next_edge)
               write(0,*) 'VerticesOnCell(2) = ',verticesOnCell(:,cellsOnEdge(2,j))
               write(0,*) 'edgesOnCell(2) = ',edgesOnCell(:,cellsOnEdge(2,j))
               stop
            end if
            cell1GP % lat = latCell(cellsOnEdge(2,j))
            cell1GP % lon = lonCell(cellsOnEdge(2,j))
            vertex1GP % lat = latVertex(vtx)
            vertex1GP % lon = lonVertex(vtx)

            edgeGP_prev % lat = latEdge(prev_edge)
            edgeGP_prev % lon = lonEdge(prev_edge)
            edgeGP_next % lat = latEdge(next_edge)
            edgeGP_next % lon = lonEdge(next_edge)

            if(isObtuse(vtx).eq.0) then
               area =        abs(triangle_area(cell1GP, vertex1GP, edgeGP_prev, 1.0_RKIND))
               area = area + abs(triangle_area(cell1GP, vertex1GP, edgeGP_next, 1.0_RKIND))
            else
               if(cellsOnVertex(isObtuse(vtx),vtx).eq.cellsOnEdge(2,j)) then
                 iFlag = 0
                 do ii=1,3
                  jj = edgesOnVertex(ii,vtx)
                  if(jj.ne.prev_edge.and.jj.ne.next_edge) then
                      write(6,*) jj, prev_edge, next_edge       
                      if(iFlag.eq.1) then
                         write(6,*) ' can not be true '
                         stop
                      endif
                      iEdge = jj
                      iFlag = 1
                   endif
                  enddo
                  edgeGP % lat = latEdge(iEdge)
                  edgeGP % lon = lonEdge(iEdge)
                  area =        abs(triangle_area(cell1GP, edgeGP, edgeGP_prev, 1.0_RKIND))
                  area = area + abs(triangle_area(cell1GP, edgeGP, edgeGP_next, 1.0_RKIND))
               else
                  area = abs(triangle_area(cell1GP, edgeGP_prev, edgeGP_next, 1.0_RKIND))
               endif
            endif

            do ii=1,3
               if (cellsOnEdge(2,j) == cellsOnVertex(ii,vtx)) then
                  kiteAreasOnVertex(ii,vtx) = area
                  exit
               end if
            end do

            if(.not.tdrtest) then
              write(6,*) ' not tdrtest'
              r = area / areaCell(cellsOnEdge(2,j))
              sum_r = sum_r + r
              if (cellsOnEdge(2,j) == cellsOnEdge(1,edgesOnEdge(i,j))) then
                 s = -1.0_RKIND
              else
                 s = 1.0_RKIND
              end if
              weightsOnEdge(i,j) = s*(0.5-sum_r)*dvEdge(edgesOnEdge(i,j))/de
            endif

            prev_edge = next_edge
         end do
      end do

!----

      if(tdrtest) then
 
      write(6,*) ' testing tdr '

      areaTriangle = 0.0
      areaCell = 0.0
      do i=1,nVertices
        do ii=1,3
           jj = cellsOnVertex(ii,i)
           areaCell(jj) = areaCell(jj) + kiteAreasOnVertex(ii,i)
           areaTriangle(i) = areaTriangle(i) + kiteAreasOnVertex(ii,i)
         enddo
      end do

      do j=1,nEdges
         de = dcEdge(j)
         prev_edge = j
         sum_r = 0.0
         do i=1,nEdgesOnCell(cellsOnEdge(1,j))-1
            next_edge = edgesOnEdge(i,j)
            if ((verticesOnEdge(1,prev_edge) == verticesOnEdge(1,next_edge)) .or. &
                (verticesOnEdge(1,prev_edge) == verticesOnEdge(2,next_edge))) then
               vtx = verticesOnEdge(1,prev_edge)
            else if ((verticesOnEdge(2,prev_edge) == verticesOnEdge(1,next_edge)) .or. &
                     (verticesOnEdge(2,prev_edge) == verticesOnEdge(2,next_edge))) then
               vtx = verticesOnEdge(2,prev_edge)
            else
               stop
            end if

            do ii=1,3
             if(cellsOnVertex(ii,vtx).eq.cellsOnEdge(1,j)) then
               area = kiteAreasOnVertex(ii,vtx)
               exit
             endif
            enddo

            r = area / areaCell(cellsOnEdge(1,j))
            sum_r = sum_r + r
            if (cellsOnEdge(1,j) == cellsOnEdge(1,edgesOnEdge(i,j))) then
               s = 1.0_RKIND
            else
               s = -1.0_RKIND
            end if
            weightsOnEdge(i,j) = s*(0.5-sum_r)*dvEdge(edgesOnEdge(i,j))/de
            prev_edge = next_edge
         end do

         prev_edge = j
         sum_r = 0.0
         do i=nEdgesOnCell(cellsOnEdge(1,j)),nEdgesOnEdge(j)
            next_edge = edgesOnEdge(i,j)
            if ((verticesOnEdge(1,prev_edge) == verticesOnEdge(1,next_edge)) .or. &
                (verticesOnEdge(1,prev_edge) == verticesOnEdge(2,next_edge))) then
               vtx = verticesOnEdge(1,prev_edge)
            else if ((verticesOnEdge(2,prev_edge) == verticesOnEdge(1,next_edge)) .or. &
                     (verticesOnEdge(2,prev_edge) == verticesOnEdge(2,next_edge))) then
               vtx = verticesOnEdge(2,prev_edge)
            else
               stop
            end if

            do ii=1,3
             if(cellsOnVertex(ii,vtx).eq.cellsOnEdge(2,j)) then
               area = kiteAreasOnVertex(ii,vtx)
               exit
             endif
            enddo

            r = area / areaCell(cellsOnEdge(2,j))
            sum_r = sum_r + r
            if (cellsOnEdge(2,j) == cellsOnEdge(1,edgesOnEdge(i,j))) then
               s = -1.0_RKIND
            else
               s = 1.0_RKIND
            end if
            weightsOnEdge(i,j) = s*(0.5-sum_r)*dvEdge(edgesOnEdge(i,j))/de
            prev_edge = next_edge
         end do
      end do

      write(6,*) ' areas ', minval(areaCell), maxval(areaCell), maxval(areaCell)/minval(areaCell)
      write(6,*) ' dcEdge ', minval(dcEdge), maxval(dcEdge), maxval(dcEdge)/minval(dcEdge)

     endif  ! tdrtest

!-----

      do i=1,nEdges
         edgeGP % lat = latEdge(i)
         edgeGP % lon = lonEdge(i)
         vertex1GP % lat = latVertex(verticesOnEdge(1,i))
         vertex1GP % lon = lonVertex(verticesOnEdge(1,i))
         vertex2GP % lat = latVertex(verticesOnEdge(2,i))
         vertex2GP % lon = lonVertex(verticesOnEdge(2,i))
         dv1Edge(i) = sphere_distance(edgeGP, vertex1GP, 1.0_RKIND) 
         dv2Edge(i) = sphere_distance(edgeGP, vertex2GP, 1.0_RKIND) 
      end do


      write(0,*) 'There are ', nCells, ' cells and ', nEdges, ' edges' 
      write(0,*) ' '
      write(0,*) 'Total area of all cells is ',sum(areaCell(:))
      write(0,*) 'Total triangle area=',sum(areaTriangle(:))
      write(0,*) 'Total kite area=',sum(kiteAreasOnVertex(:,:))


      do i=1,nCells
         cell1GP % lat = latCell(i)
         cell1GP % lon = lonCell(i)
         call convert_lx(xCell(i), yCell(i), zCell(i), 1.0_RKIND, cell1GP)
      end do
      do i=1,nVertices
         vertex1GP % lat = latVertex(i)
         vertex1GP % lon = lonVertex(i)
         call convert_lx(xVertex(i), yVertex(i), zVertex(i), 1.0_RKIND, vertex1GP)
      end do
      do i=1,nEdges
         edgeGP % lat = latEdge(i)
         edgeGP % lon = lonEdge(i)
         call convert_lx(xEdge(i), yEdge(i), zEdge(i), 1.0_RKIND, edgeGP)
      end do


      nVertLevels = 1
      nTracers = 1

      allocate(indexToCellID(nCells))
      allocate(indexToEdgeID(nEdges))
      allocate(indexToVertexID(nVertices))
      allocate(fEdge(nEdges))
      allocate(fVertex(nVertices))
      allocate(h_s(nCells))
      allocate(u(nVertLevels,nEdges,1))
      allocate(v(nVertLevels,nEdges,1))
      allocate(h(nVertLevels,nCells,1))
      allocate(vh(nVertLevels,nEdges,1))
      allocate(circulation(nVertLevels,nVertices,1))
      allocate(vorticity(nVertLevels,nVertices,1))
      allocate(ke(nVertLevels,nCells,1))
      allocate(tracers(nTracers,nVertLevels,nCells,1))

      u(1,:,1) = u_sbr(:)
      h(:,:,1) = 1000.0
      v(:,:,1) = 0.0
      vh(:,:,1) = 0.0
      circulation(:,:,1) = 0.0
      vorticity(:,:,1) = 0.0
      ke(:,:,1) = 0.0
      tracers(:,:,:,1) = 0.0

      h_s(:) = 0.0
      do i=1,nCells
         indexToCellID(i) = i
      end do
      do i=1,nEdges
         fEdge(i) = 2.0 * 7.292E-5 * sin(latEdge(i))
         indexToEdgeID(i) = i
      end do
      do i=1,nVertices
         fVertex(i) = 2.0 * 7.292E-5 * sin(latVertex(i))
         indexToVertexID(i) = i
      end do

      do i=1,nCells
        pCell%lat = latCell(i)
        pCell%lon = lonCell(i)
        if(pCell%lon.gt. pii) pCell%lon=pCell%lon-2.0*pii
        if(pCell%lon.ge.-pii) pCell%lon=pCell%lon+2.0*pii
        meshDensity(i) = density_for_point(pCell)
      enddo

      call write_netcdf_init( &
                              nCells, &
                              nEdges, &
                              nVertices, &
                              maxEdges, &
                              nVertLevels, &
                              nTracers &
                              )

      call write_netcdf_fields( &
                                1, &
                                latCell, &
                                lonCell, &
                                meshDensity, &
                                xCell, &
                                yCell, &
                                zCell, &
                                indexToCellID, &
                                latEdge, &
                                lonEdge, &
                                xEdge, &
                                yEdge, &
                                zEdge, &
                                indexToEdgeID, &
                                latVertex, &
                                lonVertex, &
                                xVertex, &
                                yVertex, &
                                zVertex, &
                                indexToVertexID, &
                                cellsOnEdge, &
                                nEdgesOnCell, &
                                nEdgesOnEdge, &
                                edgesOnCell, &
                                edgesOnEdge, &
                                weightsOnEdge, &
                                dvEdge, &
                                dv1Edge, &
                                dv2Edge, &
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
      ! Write out a file compatible with metis for block decomposition
      !
      open(42,file='graph.info',form='formatted')
      write(42,*) nCells, nEdges
      do i=1,nCells
         do j=1,nEdgesOnCell(i)
            write(42,'(1x,i8)',advance='no') cellsOnCell(j,i)
         end do
         write(42,'(1x)')
      end do
      close(42)


      !
      ! The following information is written out to permit further refinement
      !    of a grid using the grid_gen code
      !
      open(22,file='locs.dat.out',form='formatted')
      write(22,*) nCells,nEdges
      do i=1,nCells
         write(22,'(10x,f22.10,f23.10,f23.10)') xCell(i), yCell(i), zCell(i)
      end do
      do i=1,nEdges
         write(22,'(10x,f22.10,f23.10,f23.10)') xEdge(i), yEdge(i), zEdge(i)
      end do
   !  do i=1,nCells
   !     write(22,'(f13.10,1x,f13.10)') latCell(i), lonCell(i)
   !  end do
   !  do i=1,nEdges
   !     write(22,'(f13.10,1x,f13.10)') latEdge(i), lonEdge(i)
   !  end do
      close(22)

      !
      ! Write out data for visualization in OpenDx
      !
      call write_OpenDX(       nCells, &
                               nVertices, &
                               xCell, &
                               yCell, &
                               zCell, &
                               xVertex, &
                               yVertex, &
                               zVertex, &
                               nEdgesOnCell, &
                               verticesOnCell, &
                               areaCell )


      deallocate(indexToCellID)
      deallocate(indexToEdgeID)
      deallocate(indexToVertexID)
      deallocate(fEdge)
      deallocate(fVertex)
      deallocate(h_s)
      deallocate(u)
      deallocate(v)
      deallocate(h)
      deallocate(vh)
      deallocate(circulation)
      deallocate(vorticity)
      deallocate(ke)
      deallocate(tracers)


      deallocate(cursor)
      deallocate(latCell)
      deallocate(lonCell)
      deallocate(xCell)
      deallocate(yCell)
      deallocate(zCell)
      deallocate(latEdge)
      deallocate(lonEdge)
      deallocate(xEdge)
      deallocate(yEdge)
      deallocate(zEdge)
      deallocate(dvEdge)
      deallocate(dv1Edge)
      deallocate(dv2Edge)
      deallocate(dcEdge)
      deallocate(areaCell)
      deallocate(areaTriangle)
      deallocate(angleEdge)
      deallocate(u_sbr)
      deallocate(latVertex)
      deallocate(lonVertex)
      deallocate(xVertex)
      deallocate(yVertex)
      deallocate(zVertex)
      deallocate(lat1Edge)
      deallocate(lon1Edge)
      deallocate(lat2Edge)
      deallocate(lon2Edge)
      deallocate(cellsOnEdge)
      deallocate(edgesOnCell)
      deallocate(verticesOnCell)
      deallocate(verticesOnEdge)
      deallocate(edgesOnEdge)
      deallocate(edgesOnVertex)
      deallocate(cellsOnVertex)
      deallocate(tempEdgesOnEdge)
      deallocate(weightsOnEdge)
      deallocate(kiteAreasOnVertex)
      deallocate(cellsOnCell)
      deallocate(nEdgesOnCell)
      deallocate(nEdgesOnEdge)
      deallocate(nEdgesOnVertex)
      deallocate(nCellsOnVertex)
   
   end subroutine write_grid
   
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE INSERT_EDGE_TO_TREE
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine insert_edge_to_tree(cellID, vertex1ID, vertex2ID, lat1, lon1, lat2, lon2, root)
   
      use data_types
   
      implicit none
   
      integer, intent(in) :: cellID
      integer, intent(in) :: vertex1ID, vertex2ID
      real(kind=RKIND), intent(in) :: lat1, lon1, lat2, lon2
      type (binary_tree), pointer :: root
   
      real(kind=RKIND) :: tLat1, tLat2, tLon1, tLon2 
      integer :: tID1, tID2
      logical :: found
      type (binary_tree), pointer :: pre_cursor, cursor
   
      if (point_compare(lat1, lon1, lat2, lon2) > 0) then
         tLat1 = lat1
         tLon1 = lon1
         tLat2 = lat2
         tLon2 = lon2
         tID1 = vertex1ID
         tID2 = vertex2ID
      else
         tLat1 = lat2
         tLon1 = lon2
         tLat2 = lat1
         tLon2 = lon1
         tID1 = vertex2ID
         tID2 = vertex1ID
      end if
   
      if (.not. associated(root)) then
         allocate(root)
         root % lat1 = tLat1
         root % lon1 = tLon1
         root % lat2 = tLat2
         root % lon2 = tLon2
         root % vertex1 = tID1
         root % vertex2 = tID2
         root % node1 = cellID
         nullify(root % left)
         nullify(root % right)
         nullify(root % parent)
   
      else
   
         found = .false.
         cursor => root
         do while (associated(cursor) .and. .not. found)
            if (edge_compare(tLat1, tLon1, tLat2, tLon2, cursor % lat1, cursor % lon1, cursor % lat2, cursor % lon2) == 0) then
               found = .true.
            else if (edge_compare(tLat1, tLon1, tLat2, tLon2, cursor % lat1, cursor % lon1, cursor % lat2, cursor % lon2) > 0) then
               pre_cursor => cursor
               cursor => cursor % right
            else 
               pre_cursor => cursor
               cursor => cursor % left
            end if
         end do
   
         if (.not. found) then
            if (edge_compare(tLat1, tLon1, tLat2, tLon2, pre_cursor % lat1, pre_cursor % lon1, pre_cursor % lat2, pre_cursor % lon2) > 0) then
               allocate(pre_cursor % right)
               cursor => pre_cursor % right
            else
               allocate(pre_cursor % left)
               cursor => pre_cursor % left
            end if
            cursor % lat1 = tLat1
            cursor % lon1 = tLon1
            cursor % lat2 = tLat2
            cursor % lon2 = tLon2
            cursor % vertex1 = tID1
            cursor % vertex2 = tID2
            cursor % node1 = cellID
            cursor % parent => pre_cursor
            nullify(cursor % left)
            nullify(cursor % right)
         else
            cursor % node2 = cellID
         end if
   
      end if
   
   end subroutine insert_edge_to_tree
   
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION POINT_COMPARE
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer function point_compare(lat1, lon1, lat2, lon2)
   
      implicit none
   
      real(kind=RKIND), intent(in) :: lat1, lon1, lat2, lon2
   
      point_compare = -1
   
      if (lat1 > lat2) then
         point_compare = 1
      else if (lat1 == lat2) then
         if (lon1 > lon2) then
            point_compare = 1
         else if (lon1 == lon2) then
            point_compare = 0
         end if
      end if
   
   end function point_compare
   
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION EDGE_COMPARE
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer function edge_compare(aLat1, aLon1, aLat2, aLon2, bLat1, bLon1, bLat2, bLon2)
   
      implicit none
   
      real(kind=RKIND), intent(in) :: aLat1, aLon1, aLat2, aLon2, bLat1, bLon1, bLat2, bLon2 
   
      edge_compare = -1
   
      if (point_compare(aLat1,aLon1,bLat1,bLon1) > 0) then
         edge_compare = 1
      else if (point_compare(aLat1,aLon1,bLat1,bLon1) == 0) then
         if (point_compare(aLat2,aLon2,bLat2,bLon2) > 0) then
            edge_compare = 1
         else if (point_compare(aLat2,aLon2,bLat2,bLon2) == 0) then
            edge_compare = 0
         end if
      end if
   
   end function edge_compare
   
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE COMPUTE_EDGE_LATLON
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine compute_edge_latlon(cell1, cell2, vertex1, vertex2, edge)
   
      use data_types
      use sphere_utilities
   
      implicit none
   
      type (geo_point), intent(in) :: cell1, cell2, vertex1, vertex2
      type (geo_point), intent(out) :: edge
    
      call gc_intersect(cell1, cell2, vertex1, vertex2, edge) 

   end subroutine compute_edge_latlon

 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE is_flipped_vertex_order
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   logical function is_flipped_vertex_order(latCell1, lonCell1, &
                                            latCell2, lonCell2, &
                                            latVertex1, lonVertex1, &
                                            latVertex2, lonVertex2)

      use data_types
      use sphere_utilities

      implicit none

      real(kind=RKIND), intent(in) :: latCell1, lonCell1, &
                          latCell2, lonCell2, &
                          latVertex1, lonVertex1, &
                          latVertex2, lonVertex2

      real(kind=RKIND) :: xCell1, yCell1, zCell1
      real(kind=RKIND) :: xCell2, yCell2, zCell2
      real(kind=RKIND) :: xVertex1, yVertex1, zVertex1
      real(kind=RKIND) :: xVertex2, yVertex2, zVertex2
      real(kind=RKIND) :: xV1, yV1, zV1
      real(kind=RKIND) :: xV2, yV2, zV2
      real(kind=RKIND) :: ci, cj, ck
      type (geo_point) :: cell1, cell2, vertex1, vertex2

      cell1 % lat = latCell1
      cell1 % lon = lonCell1
      cell2 % lat = latCell2
      cell2 % lon = lonCell2
      vertex1 % lat = latVertex1
      vertex1 % lon = lonVertex1
      vertex2 % lat = latVertex2
      vertex2 % lon = lonVertex2
      
      call convert_lx(xCell1, yCell1, zCell1, 1.0_RKIND, cell1) 
      call convert_lx(xCell2, yCell2, zCell2, 1.0_RKIND, cell2) 
      call convert_lx(xVertex1, yVertex1, zVertex1, 1.0_RKIND, vertex1) 
      call convert_lx(xVertex2, yVertex2, zVertex2, 1.0_RKIND, vertex2) 

      xV1 = xCell2 - xCell1
      yV1 = yCell2 - yCell1
      zV1 = zCell2 - zCell1
      xV2 = xVertex2 - xVertex1
      yV2 = yVertex2 - yVertex1
      zV2 = zVertex2 - zVertex1

      ci = yV1*zV2 - zV1*yV2
      cj = zV1*xV2 - xV1*zV2
      ck = xV1*yV2 - yV1*xV2 

      if ((ci*xCell1 + cj*yCell1 + ck*zCell1) >= 0.0) then
         is_flipped_vertex_order = .false.
      else
         is_flipped_vertex_order = .true.
      end if

   end function is_flipped_vertex_order


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE is_flipped_vertex_order2
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   logical function is_flipped_vertex_order2(edge, cell2, vertex2)

      use data_types
      use sphere_utilities
      use grid_constants

      implicit none

      type (geo_point), intent(in) :: edge, cell2, vertex2

      real(kind=RKIND) :: xEdge, yEdge, zEdge
      real(kind=RKIND) :: xCell2, yCell2, zCell2
      real(kind=RKIND) :: xVertex2, yVertex2, zVertex2
      real(kind=RKIND) :: angle

      call convert_lx(xEdge, yEdge, zEdge, 1.0_RKIND, edge) 
      call convert_lx(xCell2, yCell2, zCell2, 1.0_RKIND, cell2) 
      call convert_lx(xVertex2, yVertex2, zVertex2, 1.0_RKIND, vertex2) 

      angle = plane_angle(xEdge, yEdge, zEdge, &
                          xCell2, yCell2, zCell2, &
                          xVertex2, yVertex2, zVertex2, &
                          xEdge, yEdge, zEdge)

      if (angle > 0.0 .and. angle < pii) then
         is_flipped_vertex_order2 = .false.
      else
         is_flipped_vertex_order2 = .true.
      end if

   end function is_flipped_vertex_order2


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE ORDER_POINTS_CCW
   !
   ! Given a center around which the ordering should be done, the array of points
   !   is re-ordered in CCW order, taking the first point in the array to be the
   !   first point in the ordering, and the vector from the origin to center 
   !   as the normal vector of the suface containing the points at the center.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine order_points_ccw(center, npts, points, permutation)

      use data_types
      use sphere_utilities
      use grid_constants

      implicit none
      
      type (geo_point), intent(in) :: center
      integer, intent(in) :: npts
      integer, dimension(npts), intent(inout) :: permutation
      type (geo_point), dimension(npts), intent(inout) :: points

      integer :: i, j
      integer :: itemp
      real(kind=RKIND) :: rtemp
      real(kind=RKIND) :: nx, ny, nz
      real(kind=RKIND) :: px, py, pz
      real(kind=RKIND) :: p0x, p0y, p0z
      real(kind=RKIND), dimension(npts) :: angle
      type (geo_point) :: ptemp

      call convert_lx(nx, ny, nz, 1.0_RKIND, center)
      call convert_lx(p0x, p0y, p0z, 1.0_RKIND, points(1))

      angle(1) = 0.0

      do i=2,npts 
         call convert_lx(px, py, pz, 1.0_RKIND, points(i))
         angle(i) = plane_angle(nx, ny, nz, p0x, p0y, p0z, px, py, pz, nx, ny, nz)
         if (angle(i) < 0.0) angle(i) = angle(i) + 2.0*pii
         if (angle(i) > 2.0*pii) angle(i) = angle(i) - 2.0*pii
      end do

      do i=2,npts 
         do j=i+1,npts
            if (angle(j) < angle(i)) then
               rtemp = angle(i)
               angle(i) = angle(j)
               angle(j) = rtemp

               itemp = permutation(i)
               permutation(i) = permutation(j)
               permutation(j) = itemp
     
               ptemp = points(i)
               points(i) = points(j)
               points(j) = ptemp
            end if
         end do
      end do

   end subroutine order_points_ccw

   subroutine write_OpenDX(    nCells, &
                               nVertices, &
                               xCell, &
                               yCell, &
                               zCell, &
                               xVertex, &
                               yVertex, &
                               zVertex, &
                               nEdgesOnCell, &
                               verticesOnCell, &
                               areaCell )

      integer, intent(in) :: nCells
      integer, intent(in) :: nVertices
      real(kind=RKIND), dimension(:), intent(in) :: xCell
      real(kind=RKIND), dimension(:), intent(in) :: yCell
      real(kind=RKIND), dimension(:), intent(in) :: zCell
      real(kind=RKIND), dimension(:), intent(in) :: xVertex
      real(kind=RKIND), dimension(:), intent(in) :: yVertex
      real(kind=RKIND), dimension(:), intent(in) :: zVertex
      integer, dimension(:), intent(in) :: nEdgesOnCell
      integer, dimension(:,:), intent(in) :: verticesOnCell
      real(kind=RKIND), dimension(:), intent(in) :: areaCell

      character(len=80) :: a, b, c, d, e, f
      integer :: i, j, k, nVerticesTotal, iEdge, iLoop

      nVerticesTotal = 0
      do i=1,nCells
       nVerticesTotal = nVerticesTotal + nEdgesOnCell(i)
      enddo

      open(unit=1,file='dx/voronoi.dx',form='formatted',status='unknown')

      a = trim('object "positions list" class array type float rank 1 shape 3 items')
      b = trim('ascii data file vor.position.data')
      write(1,10) a, nVerticesTotal
      write(1,10) b
      write(1,*)
      10 format(a70,i10)
      
      a = trim('object "edge list" class array type int rank 0 items')
      b = trim('ascii data file vor.edge.data')
      c = trim('attribute "ref" string "positions"')
      write(1,10) a, nVerticesTotal
      write(1,10) b
      write(1,10) c
      write(1,*)
      
      a = trim('object "loops list" class array type int rank 0 items')
      b = trim('ascii data file vor.loop.data')
      c = trim('attribute "ref" string "edges"')
      write(1,10) a, nCells
      write(1,10) b
      write(1,10) c
      write(1,*)
      
      a = trim('object "face list" class array type int rank 0 items')
      b = trim('ascii data file vor.face.data')
      c = trim('attribute "ref" string "loops"')
      write(1,10) a, nCells
      write(1,10) b
      write(1,10) c
      write(1,*)
      
      a = trim('object 0  class array type float rank 0 items')
      b = trim('data file vor.area.data')
      c = trim('attribute "dep" string "faces"')
      write(1,10) a, nCells
      write(1,10) b
      write(1,10) c
      write(1,*)
      
      a = trim('object "area" class field')
      b = trim('component "positions"     "positions list"')
      c = trim('component "edges"         "edge list"')
      d = trim('component "loops"         "loops list"')
      e = trim('component "faces"         "face list"')
      f = trim('component "data"           0')
      write(1,10) a
      write(1,10) b
      write(1,10) c
      write(1,10) d
      write(1,10) e
      write(1,10) f

      close(1)

      open(unit=10,file='dx/vor.area.data',form='formatted',status='unknown')
      open(unit=11,file='dx/vor.face.data',form='formatted',status='unknown')
      open(unit=12,file='dx/vor.loop.data',form='formatted',status='unknown')
      open(unit=13,file='dx/vor.edge.data',form='formatted',status='unknown')
      open(unit=14,file='dx/vor.position.data',form='formatted',status='unknown')

      iLoop = 0
      iEdge = 0
      do i=1,nCells
       write(10,20) areaCell(i)
       write(11,21) i-1
       write(12,21) iLoop
       iLoop = iLoop + nEdgesOnCell(i)
       do j=1,nEdgesOnCell(i)
         write(13,21) iEdge
         iEdge = iEdge + 1
         k = verticesOnCell(j,i)
         write(14,22) xVertex(k), yVertex(k), zVertex(k)
       enddo
      enddo

 20   format(e20.10)
 21   format(i20)
 22   format(3e20.10)

      close(10)
      close(11)
      close(12)
      close(13)
      close(14)


   end subroutine write_OpenDX

   
end module grid_meta
