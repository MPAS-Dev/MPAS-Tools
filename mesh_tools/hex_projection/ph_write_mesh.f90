  module ph_write_mesh

  use precision

  contains


  subroutine write_mpas_mesh( xCell, yCell, zCell, latCell, lonCell,                                &
                            xVertex, yVertex, zVertex, latVertex, lonVertex,                        &
                            xEdge, yEdge, zEdge, latEdge, lonEdge,                                  &
                            cellMask, vertexMask, edgeMask,                                         &
                            cellsOnCell, edgesOnCell, verticesOnCell,                               &
                            cellsOnVertex, cellsOnEdge, edgesOnVertex,                              &
                            verticesOnEdge, nEdgesOnCell,                                           &
                            indexToCellID, indexToVertexID, indexToEdgeID,                          &
                            nEdgesOnEdge, edgesOnEdge, areaCell, areaTriangle, kiteAreasOnVertex,   &
                            dvEdge, dcEdge, angleEdge, weightsOnEdge, meshDensity, nominalMinDc,    &
                            nCells, nEdges, nVertices,                                              &
                            maxEdges, maxEdges2, vertexDegree )

  use netcdf

  implicit none

  integer :: nCells, nEdges, nVertices, maxEdges, maxEdges2, vertexDegree
  real(dp), dimension(nCells) :: xCell, yCell, zCell, latCell, lonCell
  real(dp), dimension(nEdges) :: xEdge, yEdge, zEdge, latEdge, lonEdge
  real(dp), dimension(nVertices) :: xVertex, yVertex, zVertex, latVertex, lonVertex
  integer, dimension(maxEdges, nCells) :: cellsOnCell, edgesOnCell, verticesOnCell
  integer, dimension(3,nVertices) :: cellsOnVertex, edgesOnVertex
  integer, dimension(2,nEdges) :: cellsOnEdge, verticesOnEdge
  integer, dimension(nCells) :: nEdgesOnCell
  integer, dimension(nCells) :: cellMask
  integer, dimension(nVertices) :: vertexMask
  integer, dimension(nEdges) :: edgeMask

  integer, dimension(nCells), intent(in) :: indexToCellID
  integer, dimension(nVertices), intent(in) :: indexToVertexID
  integer, dimension(nEdges), intent(in) :: indexToEdgeID
  integer, dimension(nEdges), intent(in) :: nEdgesOnEdge
  integer, dimension(maxEdges2, nEdges), intent(in) :: EdgesOnEdge

  real(dp), dimension(nCells), intent(in) :: areaCell, meshDensity
  real(dp), dimension(nVertices), intent(in) :: areaTriangle
  real(dp), dimension(vertexDegree, nVertices), intent(in) :: kiteAreasOnVertex
  real(dp), dimension(nEdges), intent(in) :: dcEdge, dvEdge, angleEdge
  real(dp), dimension(maxEdges2, nEdges), intent(in) :: weightsOnEdge
  real(dp), intent(in) :: nominalMinDc

  !--

  integer :: nCellsID, nEdgesID, nVerticesID, maxEdgesID
  integer :: maxEdges2ID
  integer :: xCellID, yCellID, zCellID, latCellID, lonCellID
  integer :: xEdgeID, yEdgeID, zEdgeID, latEdgeID, lonEdgeID
  integer :: xVertexID, yVertexID, zVertexID, latVertexID, lonVertexID
  integer :: cellsOnCellID, edgesOnCellID, verticesOnCellID
  integer :: cellsOnVertexID, edgesOnVertexID
  integer :: cellsOnEdgeID, verticesOnEdgeID
  integer :: nEdgesOnCellID
  integer :: cellMaskID, VertexMaskID, edgeMaskID
  integer :: indexToCellIDID, indexToVertexIDID, indexToEdgeIDID

  integer :: nEdgesOnEdgeID, edgesOnEdgeID, areaCellID, areaTriangleID, kiteAreasOnVertexID
  integer :: dvEdgeID, dcEdgeID, angleEdgeID, weightsOnEdgeID, meshDensityID, nominalMinDcID

  integer :: status, ncid, status1, status2
  integer :: twoID, threeID
  integer :: id1(1), id2(2)

  logical :: test_bounds = .true.
  real(dp) :: avglat, avglon, r2d

  r2d = 90./asin(1.0_dp)

  if(test_bounds) then
    avglat = r2d*sum(latVertex)/real(nVertices)
    avglon = r2d*sum(lonVertex)/real(nVertices)
    write(6,*) ' avg latVertex, lonVertex ',avglat,avglon
    avglat = r2d*sum(latEdge)/real(nEdges)
    avglon = r2d*sum(lonEdge)/real(nEdges)
    write(6,*) ' avg latEdge, lonEdge ',avglat,avglon
    avglat = r2d*sum(latCell)/real(nCells)
    avglon = r2d*sum(lonCell)/real(nCells)
    write(6,*) ' avg latEdge, lonEdge ',avglat,avglon
    !  return
 end if
 
  write(6,*) ' writing netcdf mesh file '
  write(6,*) ' netcdf create mesh '
  write(6,*) ' defining dimensions '

      status = nf90_create("mpas_hex_mesh.nc", NF90_64BIT_DATA, ncid)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_create for netcdf file"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      ! setting dimensions

      status = nf90_def_dim(ncid, 'nCells', nCells, nCellsID)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_dim for 'nCells'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_def_dim(ncid, 'nVertices', nVertices, nVerticesID)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_dim for 'nVertices'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_def_dim(ncid, 'nEdges', nEdges, nEdgesID)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_dim for 'nEdges'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_def_dim(ncid, 'maxEdges', maxEdges, maxEdgesID)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_dim for 'maxEdges'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      
      status = nf90_def_dim(ncid, 'maxEdges2', maxEdges2, maxEdges2ID)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_dim for 'maxEdges2'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      
      status = nf90_def_dim(ncid, 'TWO', 2, twoID)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_dim for 'TWO'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      
      status = nf90_def_dim(ncid, 'vertexDegree', vertexDegree, threeID)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_dim for 'vertexDegree'"
          write(0,*) trim(nf90_strerror(status))
          stop 
       end if
       
       write(6,*) ' defining data '
      ! defining arrays

      id1(1) = nCellsID
      status = nf90_def_var(ncid, 'nEdgesOnCell', NF90_INT, id1, nEdgesOnCellID)
      status1 = nf90_put_att(ncid, nEdgesOnCellID, 'units' , '-')
      status2 = nf90_put_att(ncid, nEdgesOnCellID, 'long_name', 'Number of edges forming the boundary of a cell')
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_var for 'nEdgesOnCell'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_def_var(ncid, 'indexToCellID', NF90_INT, id1, indexToCellIDID)
      status1 = nf90_put_att(ncid, indexToCellIDID, 'units' , '-')
      status2 = nf90_put_att(ncid, indexToCellIDID, 'long_name', 'Mapping from local array index to global cell ID')
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_var for 'indexToCellID'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_def_var(ncid, 'bdyMaskCell', NF90_INT, id1, cellMaskID)
      status1 = nf90_put_att(ncid, cellMaskID, 'units' , '-')
      status2 = nf90_put_att(ncid, cellMaskID, 'long_name', 'Limited-area specified/relaxation zone index for cells')
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_var for 'cellMask'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      id1 = nVerticesID
      status = nf90_def_var(ncid, 'indexToVertexID', NF90_INT, id1, indexToVertexIDID)
      status1 = nf90_put_att(ncid, indexToVertexIDID, 'units' , '-')
      status2 = nf90_put_att(ncid, indexToVertexIDID, 'long_name', 'Mapping from local array index to global vertex ID')
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_var for 'indexToVertexID'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_def_var(ncid, 'bdyMaskVertex', NF90_INT, id1, vertexMaskID)
      status1 = nf90_put_att(ncid, vertexMaskID, 'units' , '-')
      status2 = nf90_put_att(ncid, vertexMaskID, 'long_name', 'Limited-area specified/relaxation zone index for vertices')
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_var for 'cellMask'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      
      id1 = nEdgesID
      status = nf90_def_var(ncid, 'indexToEdgeID', NF90_INT, id1, indexToEdgeIDID)
      status1 = nf90_put_att(ncid, indexToEdgeIDID, 'units' , '-')
      status2 = nf90_put_att(ncid, indexToEdgeIDID, 'long_name', 'Mapping from local array index to global edge ID')
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_var for 'indexToEdgeID'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      
      status = nf90_def_var(ncid, 'bdyMaskEdge', NF90_INT, id1, edgeMaskID)
      status1 = nf90_put_att(ncid, edgeMaskID, 'units' , '-')
      status2 = nf90_put_att(ncid, edgeMaskID, 'long_name', 'Limited-area specified/relaxation zone index for vertices')
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_var for 'cellMask'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      
      id2(1) = maxEdgesID
      id2(2) = nCellsID
      status = nf90_def_var(ncid, 'cellsOnCell', NF90_INT, id2, cellsOnCellID)
      status1 = nf90_put_att(ncid, cellsOnCellID, 'units' , '-')
      status2 = nf90_put_att(ncid, cellsOnCellID, 'long_name', 'IDs of cells sharing an edge')
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_var for 'cellsOnCell'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_def_var(ncid, 'edgesOnCell', NF90_INT, id2, edgesOnCellID)
      status1 = nf90_put_att(ncid, edgesOnCellID, 'units' , '-')
      status2 = nf90_put_att(ncid, edgesOnCellID, 'long_name', 'IDs of edges forming the boundary of a cell')
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_var for 'edgesOnCell'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_def_var(ncid, 'verticesOnCell', NF90_INT, id2, verticesOnCellID)
      status1 = nf90_put_att(ncid, verticesOnCellID, 'units' , '-')
      status2 = nf90_put_att(ncid, verticesOnCellID, 'long_name', 'IDs of vertices (corner points) of a cell')
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_var for 'verticesOnCell'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      id2(1) = threeID
      id2(2) = nVerticesID
       
      status = nf90_def_var(ncid, 'cellsOnVertex', NF90_INT, id2, cellsOnVertexID)
      status1 = nf90_put_att(ncid, cellsOnVertexID, 'units' , '-')
      status2 = nf90_put_att(ncid, cellsOnVertexID, 'long_name', 'IDs of the cells that meet at a vertex')
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_var for 'cellsOnVertex'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_def_var(ncid, 'edgesOnVertex', NF90_INT, id2, edgesOnVertexID)
      status1 = nf90_put_att(ncid, edgesOnVertexID, 'units' , '-')
      status2 = nf90_put_att(ncid, edgesOnVertexID, 'long_name', 'IDs of the edges that meet at a vertex')
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_var for 'edgesOnVertex'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
       
      id1(1) = nEdgesID
      status = nf90_def_var(ncid, 'nEdgesOnEdge', NF90_INT, id1, nEdgesOnEdgeID)
      status1 = nf90_put_att(ncid, nEdgesOnEdgeID, 'units' , '-')
      status2 = nf90_put_att(ncid, nEdgesOnEdgeID, 'long_name', &
           'Number of edges involved in reconstruction of tangential velocity for an edge')
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_var for 'nEdgesOnEdge'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      
      id2(1) = maxEdges2ID
      id2(2) = nEdgesID
      status = nf90_def_var(ncid, 'edgesOnEdge', NF90_INT, id2, edgesOnEdgeID)
      status1 = nf90_put_att(ncid, edgesOnEdgeID, 'units' , '-')
      status2 = nf90_put_att(ncid, edgesOnEdgeID, 'long_name', &
                 'IDs of edges involved in reconstruction of tangential velocity for an edge')
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_var for 'edgesOnedge'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      id1(1) = nCellsID 
      status = nf90_def_var(ncid, 'areaCell', NF90_DOUBLE, id1, areaCellID)
      status1 = nf90_put_att(ncid, areaCellID, 'units' , 'm^2')
      status2 = nf90_put_att(ncid, areaCellID, 'long_name', 'Spherical area of a Voronoi cell')
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_var for 'areaCell'"
          write(0,*) trim(nf90_strerror(status))
          stop 
       end if

      status = nf90_def_var(ncid, 'meshDensity', NF90_DOUBLE, id1, meshDensityID)
      status1 = nf90_put_att(ncid, meshDensityID, 'units' , 'unitless')
      status2 = nf90_put_att(ncid, meshDensityID, 'long_name', &
                'Mesh density function (used when generating the mesh) evaluated at a cell')
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_var for 'meshDensity'"
          write(0,*) trim(nf90_strerror(status))
          stop 
       end if
       
      id1(1) = nVerticesID 
      status = nf90_def_var(ncid, 'areaTriangle', NF90_DOUBLE, id1, areaTriangleID)
      status1 = nf90_put_att(ncid, areaTriangleID, 'units' , 'm^2')
      status2 = nf90_put_att(ncid, areaTriangleID, 'long_name', 'Spherical area of a Delaunay triangle')
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_var for 'areaTriangle'"
          write(0,*) trim(nf90_strerror(status))
          stop 
       end if

      id2(1) = ThreeID 
      id2(2) = nVerticesID 
      status = nf90_def_var(ncid, 'kiteAreasOnVertex', NF90_DOUBLE, id2, kiteAreasOnVertexID)
      status1 = nf90_put_att(ncid, kiteAreasOnVertexID, 'units' , 'm^2')
      status2 = nf90_put_att(ncid, kiteAreasOnVertexID, 'long_name', &
           'Intersection areas between primal (Voronoi) and dual (triangular) mesh cells')
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_var for 'kiteAreasOnVertex'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
       
      id1(1) = nEdgesID 
      status = nf90_def_var(ncid, 'dcEdge', NF90_DOUBLE, id1, dcEdgeID)
      status1 = nf90_put_att(ncid, dcEdgeID, 'units' , 'm')
      status2 = nf90_put_att(ncid, dcEdgeID, 'long_name', 'Spherical distance between cells separated by an edge')
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_var for 'dcEdge'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_def_var(ncid, 'dvEdge', NF90_DOUBLE, id1, dvEdgeID)
      status1 = nf90_put_att(ncid, dvEdgeID, 'units' , 'm')
      status2 = nf90_put_att(ncid, dvEdgeID, 'long_name', 'Spherical distance between vertex endpoints of an edge')
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_var for 'dvEdge'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_def_var(ncid, 'angleEdge', NF90_DOUBLE, id1, angleEdgeID)
      status1 = nf90_put_att(ncid, angleEdgeID, 'units' , 'rad')
      status2 = nf90_put_att(ncid, angleEdgeID, 'long_name', &
              'Angle between local north and the positive tangential direction of an edge')
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_var for 'angleEdge'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      id2(1) = maxEdges2ID
      id2(2) = nEdgesID
      status = nf90_def_var(ncid, 'weightsOnEdge', NF90_DOUBLE, id2, weightsOnEdgeID)
      status1 = nf90_put_att(ncid, weightsOnEdgeID, 'units' , '-')
      status2 = nf90_put_att(ncid, weightsOnEdgeID, 'long_name', &
                             'Weights used in reconstruction of tangential velocity for an edge')
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_var for 'weightsOnEdge'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      id2(1) = twoID
      id2(2) = nEdgesID
       
      status = nf90_def_var(ncid, 'cellsOnEdge', NF90_INT, id2, cellsOnEdgeID)
      status1 = nf90_put_att(ncid, cellsOnEdgeID, 'units' , '-')
      status2 = nf90_put_att(ncid, cellsOnEdgeID, 'long_name', 'IDs of cells sharing by an edge')
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_var for 'cellsOnEdge'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_def_var(ncid, 'verticesOnEdge', NF90_INT, id2, verticesOnEdgeID)
      status1 = nf90_put_att(ncid, verticesOnEdgeID, 'units' , '-')
      status2 = nf90_put_att(ncid, verticesOnEdgeID, 'long_name', 'IDs of the two vertex endpoints of an edge')
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_var for 'verticesOnEdge'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      
      id1(1) = nCellsID
      status = nf90_def_var(ncid, 'xCell', NF90_DOUBLE, id1, xCellID)
      status1 = nf90_put_att(ncid, xCellID, 'units' , 'm')
      status2 = nf90_put_att(ncid, xCellID, 'long_name', 'Cartesian x-coordinate of cells')
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_var for 'xCell'"
          write(0,*) trim(nf90_strerror(status))
          stop 
       end if

      id1(1) = nCellsID
      status = nf90_def_var(ncid, 'yCell', NF90_DOUBLE, id1, yCellID)
      status1 = nf90_put_att(ncid, yCellID, 'units' , 'm')
      status2 = nf90_put_att(ncid, yCellID, 'long_name', 'Cartesian y-coordinate of cells')
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_var for 'yCell'"
          write(0,*) trim(nf90_strerror(status))
          stop 
       end if

      id1(1) = nCellsID
      status = nf90_def_var(ncid, 'zCell', NF90_DOUBLE, id1, zCellID)
      status1 = nf90_put_att(ncid, zCellID, 'units' , 'm')
      status2 = nf90_put_att(ncid, zCellID, 'long_name', 'Cartesian z-coordinate of cells')
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_var for 'zCell'"
          write(0,*) trim(nf90_strerror(status))
          stop 
       end if

      id1(1) = nCellsID
      status = nf90_def_var(ncid, 'latCell', NF90_DOUBLE, id1, latCellID)
      status1 = nf90_put_att(ncid, latCellID, 'units' , 'rad')
      status2 = nf90_put_att(ncid, latCellID, 'long_name', 'Latitude of cells')
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_var for 'latCell'"
          write(0,*) trim(nf90_strerror(status))
          stop 
       end if

      id1(1) = nCellsID
      status = nf90_def_var(ncid, 'lonCell', NF90_DOUBLE, id1, lonCellID)
      status1 = nf90_put_att(ncid, lonCellID, 'units' , 'rad')
      status2 = nf90_put_att(ncid, lonCellID, 'long_name', 'Longitude of cells')
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_var for 'lonCell'"
          write(0,*) trim(nf90_strerror(status))
          stop 
       end if
       
      id1(1) = nEdgesID
      status = nf90_def_var(ncid, 'xEdge', NF90_DOUBLE, id1, xEdgeID)
      status1 = nf90_put_att(ncid, xEdgeID, 'units' , 'm')
      status2 = nf90_put_att(ncid, xEdgeID, 'long_name', 'Cartesian x-coordinate of edges')
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_var for 'xEdge'"
          write(0,*) trim(nf90_strerror(status))
          stop 
       end if

      id1(1) = nEdgesID
      status = nf90_def_var(ncid, 'yEdge', NF90_DOUBLE, id1, yEdgeID)
      status1 = nf90_put_att(ncid, yEdgeID, 'units' , 'm')
      status2 = nf90_put_att(ncid, yEdgeID, 'long_name', 'Cartesian y-coordinate of edges')
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_var for 'yEdge'"
          write(0,*) trim(nf90_strerror(status))
          stop 
       end if

      id1(1) = nEdgesID
      status = nf90_def_var(ncid, 'zEdge', NF90_DOUBLE, id1, zEdgeID)
      status1 = nf90_put_att(ncid, zEdgeID, 'units' , 'm')
      status2 = nf90_put_att(ncid, zEdgeID, 'long_name', 'Cartesian z-coordinate of edges')
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_var for 'zEdge'"
          write(0,*) trim(nf90_strerror(status))
          stop 
       end if

      id1(1) = nEdgesID
      status = nf90_def_var(ncid, 'latEdge', NF90_DOUBLE, id1, latEdgeID)
      status1 = nf90_put_att(ncid, latEdgeID, 'units' , 'rad')
      status2 = nf90_put_att(ncid, latEdgeID, 'long_name', 'Latitude of edges')
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_var for 'latEdge'"
          write(0,*) trim(nf90_strerror(status))
          stop 
       end if

      id1(1) = nEdgesID
      status = nf90_def_var(ncid, 'lonEdge', NF90_DOUBLE, id1, lonEdgeID)
      status1 = nf90_put_att(ncid, lonEdgeID, 'units' , 'rad')
      status2 = nf90_put_att(ncid, lonEdgeID, 'long_name', 'Longitude of edges')
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_var for 'lonEdge'"
          write(0,*) trim(nf90_strerror(status))
          stop 
       end if
       
      id1(1) = nVerticesID
      status = nf90_def_var(ncid, 'xVertex', NF90_DOUBLE, id1, xVertexID)
      status1 = nf90_put_att(ncid, xVertexID, 'units' , 'm')
      status2 = nf90_put_att(ncid, xVertexID, 'long_name', 'Cartesian x-coordinate of vertices')
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_var for 'xVertex'"
          write(0,*) trim(nf90_strerror(status))
          stop 
       end if

      id1(1) = nVerticesID
      status = nf90_def_var(ncid, 'yVertex', NF90_DOUBLE, id1, yVertexID)
      status1 = nf90_put_att(ncid, yVertexID, 'units' , 'm')
      status2 = nf90_put_att(ncid, yVertexID, 'long_name', 'Cartesian y-coordinate of vertices')
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_var for 'yVertex'"
          write(0,*) trim(nf90_strerror(status))
          stop 
       end if

      id1(1) = nVerticesID
      status = nf90_def_var(ncid, 'zVertex', NF90_DOUBLE, id1, zVertexID)
      status1 = nf90_put_att(ncid, zVertexID, 'units' , 'm')
      status2 = nf90_put_att(ncid, zVertexID, 'long_name', 'Cartesian z-coordinate of vertices')
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_var for 'zVertex'"
          write(0,*) trim(nf90_strerror(status))
          stop 
       end if

      id1(1) = nVerticesID
      status = nf90_def_var(ncid, 'latVertex', NF90_DOUBLE, id1, latVertexID)
      status1 = nf90_put_att(ncid, latVertexID, 'units' , 'rad')
      status2 = nf90_put_att(ncid, latVertexID, 'long_name', 'Latitude of vertices')
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_var for 'latVertex'"
          write(0,*) trim(nf90_strerror(status))
          stop 
       end if

      id1(1) = nVerticesID
      status = nf90_def_var(ncid, 'lonVertex', NF90_DOUBLE, id1, lonVertexID)
      status1 = nf90_put_att(ncid, lonVertexID, 'units' , 'rad')
      status2 = nf90_put_att(ncid, lonVertexID, 'long_name', 'Longitude of vertices')
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_var for 'lonVertex'"
          write(0,*) trim(nf90_strerror(status))
          stop 
       end if

      status = nf90_def_var(ncid, 'nominalMinDc', NF90_DOUBLE, nominalMinDcID)
      status1 = nf90_put_att(ncid, nominalMinDcID, 'units' , 'm')
      status2 = nf90_put_att(ncid, nominalMinDcID, 'long_name', 'Nominal minimum dcEdge value where meshDensity == 1.0')
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_def_var for 'nominalMinDc'"
          write(0,*) trim(nf90_strerror(status))
          stop 
       end if

       ! global attribute

       write(6,*) ' setting global attributes '
       
       status = nf90_put_att(ncid, NF90_GLOBAL, 'mesh_spec', '1.1')       
       status = nf90_put_att(ncid, NF90_GLOBAL, 'on_a_sphere', 'YES')       
       status = nf90_put_att(ncid, NF90_GLOBAL, 'global', 'NO')       
       status = nf90_put_att(ncid, NF90_GLOBAL, 'sphere_radius', real((/1.0/)))       
       status = nf90_put_att(ncid, NF90_GLOBAL, 'is_periodic', 'NO')       
       status = nf90_put_att(ncid, NF90_GLOBAL, 'x_period', '0.')       
       status = nf90_put_att(ncid, NF90_GLOBAL, 'y_period', '0.')       
       status = nf90_put_att(ncid, NF90_GLOBAL, 'file_id', '-')       
       
! end file definition

       write(6,*) ' closing define mode '
       
      status = nf90_enddef(ncid)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_enddef file"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

!  write data
!  first we output connecticity
      
      write(6,*) ' writing data '
      
      status = nf90_put_var(ncid, indexToCellIDID, indexToCellID)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_put_var for 'indexToCellID'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_put_var(ncid, indexToVertexIDID, indexToVertexID)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_put_var for 'indexToVertexID'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_put_var(ncid, indexToEdgeIDID, indexToEdgeID)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_put_var for 'indexToEdgeID'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      write(6,*) ' index IDs written '
      
      status = nf90_put_var(ncid, nEdgesOnCellID, nEdgesOnCell)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_put_var for 'nEdgesOnCell'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_put_var(ncid, cellsOnCellID, cellsOnCell)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_put_var for 'cellsOnCell'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_put_var(ncid, edgesOnCellID, edgesOnCell)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_put_var for 'edgesOnCell'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_put_var(ncid, verticesOnCellID, verticesOnCell)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_put_var for 'verticesOnCell'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      write(6,*) ' items on cells written ' 
      
      status = nf90_put_var(ncid, cellsOnVertexID, cellsOnVertex)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_put_var for 'cellsOnVertex'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_put_var(ncid, edgesOnVertexID, edgesOnVertex)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_put_var for 'edgesOnVertex'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      
      status = nf90_put_var(ncid, cellsOnEdgeID, cellsOnEdge)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_put_var for 'cellsOnEdge'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      
      status = nf90_put_var(ncid, verticesOnEdgeID, verticesOnEdge)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_put_var for 'verticesOnEdge'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      
      status = nf90_put_var(ncid, cellMaskID, cellMask)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_put_var for 'cellMask'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_put_var(ncid, vertexMaskID, vertexMask)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_put_var for 'vertexMask'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_put_var(ncid, edgeMaskID, edgeMask)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_put_var for 'edgeMask'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      
! next - postions on the sphere

      status = nf90_put_var(ncid, xCellID, xCell)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_put_var for 'xCell'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_put_var(ncid, yCellID, yCell)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_put_var for 'yCell'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_put_var(ncid, zCellID, zCell)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_put_var for 'zCell'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_put_var(ncid, latCellID, latCell)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_put_var for 'latCell'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      
      status = nf90_put_var(ncid, lonCellID, lonCell)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_put_var for 'lonCell'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      
      status = nf90_put_var(ncid, xVertexID, xVertex)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_put_var for 'xVertex'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_put_var(ncid, yVertexID, yVertex)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_put_var for 'yVertex'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_put_var(ncid, zVertexID, zVertex)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_put_var for 'zVertex'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_put_var(ncid, latVertexID, latVertex)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_put_var for 'latVertex'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      
      status = nf90_put_var(ncid, lonVertexID, lonVertex)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_put_var for 'lonVertex'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_put_var(ncid, xEdgeID, xEdge)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_put_var for 'xEdge'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_put_var(ncid, yEdgeID, yEdge)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_put_var for 'yEdge'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_put_var(ncid, zEdgeID, zEdge)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_put_var for 'zEdge'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_put_var(ncid, latEdgeID, latEdge)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_put_var for 'latEdge'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      
      status = nf90_put_var(ncid, lonEdgeID, lonEdge)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_put_var for 'lonEdge'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_put_var(ncid, nEdgesOnEdgeID, nEdgesOnEdge)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_put_var for 'nEdgesOnEdge'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_put_var(ncid, edgesOnEdgeID, edgesOnEdge)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_put_var for 'edgesOnEdge'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      
      status = nf90_put_var(ncid, areaCellID, areaCell)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_put_var for 'areaCell'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_put_var(ncid, areaTriangleID, areaTriangle)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_put_var for 'areaTriangle'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_put_var(ncid, kiteAreasOnVertexID, kiteAreasOnVertex)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_put_var for 'kiteAreasOnVertex'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      
      status = nf90_put_var(ncid, dcEdgeID, dcEdge)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_put_var for 'dcEdge'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_put_var(ncid, dvEdgeID, dvEdge)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_put_var for 'dvEdge'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_put_var(ncid, angleEdgeID, angleEdge)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_put_var for 'angleEdge'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_put_var(ncid, weightsOnEdgeID, weightsOnEdge)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_put_var for 'weightsOnEdge'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_put_var(ncid, meshDensityID, meshDensity)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_put_var for 'meshDensity'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_put_var(ncid, nominalMinDcID, nominalMinDc)
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_put_var for 'meshDensity'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

! write complete
! close file

      write(6,*) ' closing file '
      
      status = nf90_close(ncid) 
      if (status /= nf90_noerr) then
          write(0,*) "project_hexes: Error occured in nf90_close"
          write(0,*) trim(nf90_strerror(status))
          stop
       end if

       write(6,*) ' finished file output '
   
    end subroutine write_mpas_mesh

  end module ph_write_mesh
  
