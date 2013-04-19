module read_netcdf
 
   integer :: rd_ncid
   integer :: rdDimIDTime
   integer :: rdDimIDnCells
   integer :: rdDimIDnEdges
   integer :: rdDimIDnVertices
   integer :: rdDimIDmaxEdges
   integer :: rdDimIDmaxEdges2
   integer :: rdDimIDnVertLevels
   integer :: rdDimIDTWO
   integer :: rdDimIDvertexDegree
   integer :: rdVarIDlatCell
   integer :: rdVarIDlonCell
   integer :: rdVarIDmeshDensity
   integer :: rdVarIDxCell
   integer :: rdVarIDyCell
   integer :: rdVarIDzCell
   integer :: rdVarIDindexToCellID
   integer :: rdVarIDlatEdge
   integer :: rdVarIDlonEdge
   integer :: rdVarIDxEdge
   integer :: rdVarIDyEdge
   integer :: rdVarIDzEdge
   integer :: rdVarIDindexToEdgeID
   integer :: rdVarIDlatVertex
   integer :: rdVarIDlonVertex
   integer :: rdVarIDxVertex
   integer :: rdVarIDyVertex
   integer :: rdVarIDzVertex
   integer :: rdVarIDindexToVertexID
   integer :: rdVarIDcellsOnEdge
   integer :: rdVarIDnEdgesOnCell
   integer :: rdVarIDnEdgesOnEdge
   integer :: rdVarIDedgesOnCell
   integer :: rdVarIDedgesOnEdge
   integer :: rdVarIDweightsOnEdge
   integer :: rdVarIDdvEdge
   integer :: rdVarIDdcEdge
   integer :: rdVarIDangleEdge
   integer :: rdVarIDareaCell
   integer :: rdVarIDareaTriangle
   integer :: rdVarIDcellsOnCell
   integer :: rdVarIDverticesOnCell
   integer :: rdVarIDverticesOnEdge
   integer :: rdVarIDedgesOnVertex
   integer :: rdVarIDcellsOnVertex
   integer :: rdVarIDkiteAreasOnVertex
   integer :: rdVarIDfEdge
   integer :: rdVarIDfVertex
   integer :: rdVarIDh_s
   integer :: rdVarIDu
   integer :: rdVarIDv
   integer :: rdVarIDh
 
   integer :: rdLocalnCells
   integer :: rdLocalnEdges
   integer :: rdLocalnVertices
   integer :: rdLocalmaxEdges
   integer :: rdLocalmaxEdges2
   integer :: rdLocalnVertLevels
   integer :: rdLocalTWO
   integer :: rdLocalvertexDegree
 
   contains
 
   subroutine read_netcdf_init( &
                               nCells, &
                               nEdges, &
                               nVertices, &
                               maxEdges, &
                               maxEdges2, &
                               nVertLevels, &
                               TWO, &
                               vertexDegree &
                               )
 
      implicit none
 
      include 'netcdf.inc'
 
      integer, intent(out) :: nCells
      integer, intent(out) :: nEdges
      integer, intent(out) :: nVertices
      integer, intent(out) :: maxEdges
      integer, intent(out) :: maxEdges2
      integer, intent(out) :: nVertLevels
      integer, intent(out) :: TWO
      integer, intent(out) :: vertexDegree
 
      integer :: nferr
 
 
      nferr = nf_open('grid.nc', NF_SHARE, rd_ncid)
 
      !
      ! Get IDs for variable dimensions
      !
      nferr = nf_inq_unlimdim(rd_ncid, rdDimIDTime)
      nferr = nf_inq_dimid(rd_ncid, 'nCells', rdDimIDnCells)
      nferr = nf_inq_dimlen(rd_ncid, rdDimIDnCells, rdLocalnCells)
      nferr = nf_inq_dimid(rd_ncid, 'nEdges', rdDimIDnEdges)
      nferr = nf_inq_dimlen(rd_ncid, rdDimIDnEdges, rdLocalnEdges)
      nferr = nf_inq_dimid(rd_ncid, 'nVertices', rdDimIDnVertices)
      nferr = nf_inq_dimlen(rd_ncid, rdDimIDnVertices, rdLocalnVertices)
      nferr = nf_inq_dimid(rd_ncid, 'maxEdges', rdDimIDmaxEdges)
      nferr = nf_inq_dimlen(rd_ncid, rdDimIDmaxEdges, rdLocalmaxEdges)
      nferr = nf_inq_dimid(rd_ncid, 'maxEdges2', rdDimIDmaxEdges2)
      nferr = nf_inq_dimlen(rd_ncid, rdDimIDmaxEdges2, rdLocalmaxEdges2)
      nferr = nf_inq_dimid(rd_ncid, 'nVertLevels', rdDimIDnVertLevels)
      nferr = nf_inq_dimlen(rd_ncid, rdDimIDnVertLevels, rdLocalnVertLevels)
      nferr = nf_inq_dimid(rd_ncid, 'vertexDegree', rdDimIDvertexDegree)
      nferr = nf_inq_dimlen(rd_ncid, rdDimIDvertexDegree, rdLocalvertexDegree)
      nferr = nf_inq_dimid(rd_ncid, 'TWO', rdDimIDTWO)
      nferr = nf_inq_dimlen(rd_ncid, rdDimIDTWO, rdLocalTWO)


      nCells = rdLocalnCells
      nEdges = rdLocalnEdges
      nVertices = rdLocalnVertices
      maxEdges = rdLocalmaxEdges
      maxEdges2 = rdLocalmaxEdges2
      nVertLevels = rdLocalnVertLevels
      vertexDegree = rdLocalvertexDegree
      TWO = rdLocalTWO
 
      !
      ! Get IDs for variables
      !
      nferr = nf_inq_varid(rd_ncid, 'latCell', rdVarIDlatCell)
      nferr = nf_inq_varid(rd_ncid, 'lonCell', rdVarIDlonCell)
      nferr = nf_inq_varid(rd_ncid, 'meshDensity', rdVarIDmeshDensity)
      nferr = nf_inq_varid(rd_ncid, 'xCell', rdVarIDxCell)
      nferr = nf_inq_varid(rd_ncid, 'yCell', rdVarIDyCell)
      nferr = nf_inq_varid(rd_ncid, 'zCell', rdVarIDzCell)
      nferr = nf_inq_varid(rd_ncid, 'indexToCellID', rdVarIDindexToCellID)
      nferr = nf_inq_varid(rd_ncid, 'latEdge', rdVarIDlatEdge)
      nferr = nf_inq_varid(rd_ncid, 'lonEdge', rdVarIDlonEdge)
      nferr = nf_inq_varid(rd_ncid, 'xEdge', rdVarIDxEdge)
      nferr = nf_inq_varid(rd_ncid, 'yEdge', rdVarIDyEdge)
      nferr = nf_inq_varid(rd_ncid, 'zEdge', rdVarIDzEdge)
      nferr = nf_inq_varid(rd_ncid, 'indexToEdgeID', rdVarIDindexToEdgeID)
      nferr = nf_inq_varid(rd_ncid, 'latVertex', rdVarIDlatVertex)
      nferr = nf_inq_varid(rd_ncid, 'lonVertex', rdVarIDlonVertex)
      nferr = nf_inq_varid(rd_ncid, 'xVertex', rdVarIDxVertex)
      nferr = nf_inq_varid(rd_ncid, 'yVertex', rdVarIDyVertex)
      nferr = nf_inq_varid(rd_ncid, 'zVertex', rdVarIDzVertex)
      nferr = nf_inq_varid(rd_ncid, 'indexToVertexID', rdVarIDindexToVertexID)
      nferr = nf_inq_varid(rd_ncid, 'cellsOnEdge', rdVarIDcellsOnEdge)
      nferr = nf_inq_varid(rd_ncid, 'nEdgesOnCell', rdVarIDnEdgesOnCell)
      nferr = nf_inq_varid(rd_ncid, 'nEdgesOnEdge', rdVarIDnEdgesOnEdge)
      nferr = nf_inq_varid(rd_ncid, 'edgesOnCell', rdVarIDedgesOnCell)
      nferr = nf_inq_varid(rd_ncid, 'edgesOnEdge', rdVarIDedgesOnEdge)
      nferr = nf_inq_varid(rd_ncid, 'weightsOnEdge', rdVarIDweightsOnEdge)
      nferr = nf_inq_varid(rd_ncid, 'dvEdge', rdVarIDdvEdge)
      nferr = nf_inq_varid(rd_ncid, 'dcEdge', rdVarIDdcEdge)
      nferr = nf_inq_varid(rd_ncid, 'angleEdge', rdVarIDangleEdge)
      nferr = nf_inq_varid(rd_ncid, 'areaCell', rdVarIDareaCell)
      nferr = nf_inq_varid(rd_ncid, 'areaTriangle', rdVarIDareaTriangle)
      nferr = nf_inq_varid(rd_ncid, 'cellsOnCell', rdVarIDcellsOnCell)
      nferr = nf_inq_varid(rd_ncid, 'verticesOnCell', rdVarIDverticesOnCell)
      nferr = nf_inq_varid(rd_ncid, 'verticesOnEdge', rdVarIDverticesOnEdge)
      nferr = nf_inq_varid(rd_ncid, 'edgesOnVertex', rdVarIDedgesOnVertex)
      nferr = nf_inq_varid(rd_ncid, 'cellsOnVertex', rdVarIDcellsOnVertex)
      nferr = nf_inq_varid(rd_ncid, 'kiteAreasOnVertex', rdVarIDkiteAreasOnVertex)
      nferr = nf_inq_varid(rd_ncid, 'fEdge', rdVarIDfEdge)
      nferr = nf_inq_varid(rd_ncid, 'fVertex', rdVarIDfVertex)
      nferr = nf_inq_varid(rd_ncid, 'h_s', rdVarIDh_s)
      nferr = nf_inq_varid(rd_ncid, 'u', rdVarIDu)
      nferr = nf_inq_varid(rd_ncid, 'v', rdVarIDv)
      nferr = nf_inq_varid(rd_ncid, 'h', rdVarIDh)
 
   end subroutine read_netcdf_init
 
 
   subroutine read_netcdf_fields( &
                                  time, &
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
                                  h &
                                 )
 
      implicit none
 
      include 'netcdf.inc'
 
      integer, intent(in) :: time
      real (kind=8), dimension(:), intent(out) :: latCell
      real (kind=8), dimension(:), intent(out) :: lonCell
      real (kind=8), dimension(:), intent(out) :: meshDensity
      real (kind=8), dimension(:), intent(out) :: xCell
      real (kind=8), dimension(:), intent(out) :: yCell
      real (kind=8), dimension(:), intent(out) :: zCell
      integer, dimension(:), intent(out) :: indexToCellID
      real (kind=8), dimension(:), intent(out) :: latEdge
      real (kind=8), dimension(:), intent(out) :: lonEdge
      real (kind=8), dimension(:), intent(out) :: xEdge
      real (kind=8), dimension(:), intent(out) :: yEdge
      real (kind=8), dimension(:), intent(out) :: zEdge
      integer, dimension(:), intent(out) :: indexToEdgeID
      real (kind=8), dimension(:), intent(out) :: latVertex
      real (kind=8), dimension(:), intent(out) :: lonVertex
      real (kind=8), dimension(:), intent(out) :: xVertex
      real (kind=8), dimension(:), intent(out) :: yVertex
      real (kind=8), dimension(:), intent(out) :: zVertex
      integer, dimension(:), intent(out) :: indexToVertexID
      integer, dimension(:,:), intent(out) :: cellsOnEdge
      integer, dimension(:), intent(out) :: nEdgesOnCell
      integer, dimension(:), intent(out) :: nEdgesOnEdge
      integer, dimension(:,:), intent(out) :: edgesOnCell
      integer, dimension(:,:), intent(out) :: edgesOnEdge
      real (kind=8), dimension(:,:), intent(out) :: weightsOnEdge
      real (kind=8), dimension(:), intent(out) :: dvEdge
      real (kind=8), dimension(:), intent(out) :: dcEdge
      real (kind=8), dimension(:), intent(out) :: angleEdge
      real (kind=8), dimension(:), intent(out) :: areaCell
      real (kind=8), dimension(:), intent(out) :: areaTriangle
      integer, dimension(:,:), intent(out) :: cellsOnCell
      integer, dimension(:,:), intent(out) :: verticesOnCell
      integer, dimension(:,:), intent(out) :: verticesOnEdge
      integer, dimension(:,:), intent(out) :: edgesOnVertex
      integer, dimension(:,:), intent(out) :: cellsOnVertex
      real (kind=8), dimension(:,:), intent(out) :: kiteAreasOnVertex
      real (kind=8), dimension(:), intent(out) :: fEdge
      real (kind=8), dimension(:), intent(out) :: fVertex
      real (kind=8), dimension(:), intent(out) :: h_s
      real (kind=8), dimension(:,:,:), intent(out) :: u
      real (kind=8), dimension(:,:,:), intent(out) :: v
      real (kind=8), dimension(:,:,:), intent(out) :: h

      logical :: meshDensityPresent
 
      integer :: nferr
      integer, dimension(1) :: start1, count1
      integer, dimension(2) :: start2, count2
      integer, dimension(3) :: start3, count3
      integer, dimension(4) :: start4, count4

      meshDensityPresent = .false.
 
      start1(1) = 1
 
      start2(1) = 1
      start2(2) = 1
 
      start3(1) = 1
      start3(2) = 1
      start3(3) = 1
 
      start1(1) = 1
      count1( 1) = rdLocalnCells
      count1( 1) = rdLocalnCells
      nferr = nf_get_vara_double(rd_ncid, rdVarIDlatCell, start1, count1, latCell)
 
      start1(1) = 1
      count1( 1) = rdLocalnCells
      count1( 1) = rdLocalnCells
      nferr = nf_get_vara_double(rd_ncid, rdVarIDlonCell, start1, count1, lonCell)

      start1(1) = 1
      count1( 1) = rdLocalnCells
      count1( 1) = rdLocalnCells
      nferr = nf_inq_varid(rd_ncid, 'meshDensity', rdVarIDmeshDensity)
      if(nferr.eq.0) then
         nferr = nf_get_vara_double(rd_ncid, rdVarIDmeshDensity, start1, count1, meshDensity)
      else
         meshDensity=1.0
         write(6,*) ' mesh density not present ', nferr, rdVarIDmeshDensity
      endif
 
      start1(1) = 1
      count1( 1) = rdLocalnCells
      count1( 1) = rdLocalnCells
      nferr = nf_get_vara_double(rd_ncid, rdVarIDxCell, start1, count1, xCell)
 
      start1(1) = 1
      count1( 1) = rdLocalnCells
      count1( 1) = rdLocalnCells
      nferr = nf_get_vara_double(rd_ncid, rdVarIDyCell, start1, count1, yCell)
 
      start1(1) = 1
      count1( 1) = rdLocalnCells
      count1( 1) = rdLocalnCells
      nferr = nf_get_vara_double(rd_ncid, rdVarIDzCell, start1, count1, zCell)
 
      start1(1) = 1
      count1( 1) = rdLocalnCells
      count1( 1) = rdLocalnCells
      nferr = nf_get_vara_int(rd_ncid, rdVarIDindexToCellID, start1, count1, indexToCellID)

      start1(1) = 1
      count1( 1) = rdLocalnEdges
      count1( 1) = rdLocalnEdges
      nferr = nf_get_vara_double(rd_ncid, rdVarIDlatEdge, start1, count1, latEdge)
 
      start1(1) = 1
      count1( 1) = rdLocalnEdges
      count1( 1) = rdLocalnEdges
      nferr = nf_get_vara_double(rd_ncid, rdVarIDlonEdge, start1, count1, lonEdge)
 
      start1(1) = 1
      count1( 1) = rdLocalnEdges
      count1( 1) = rdLocalnEdges
      nferr = nf_get_vara_double(rd_ncid, rdVarIDxEdge, start1, count1, xEdge)
 
      start1(1) = 1
      count1( 1) = rdLocalnEdges
      count1( 1) = rdLocalnEdges
      nferr = nf_get_vara_double(rd_ncid, rdVarIDyEdge, start1, count1, yEdge)
 
      start1(1) = 1
      count1( 1) = rdLocalnEdges
      count1( 1) = rdLocalnEdges
      nferr = nf_get_vara_double(rd_ncid, rdVarIDzEdge, start1, count1, zEdge)
 
      start1(1) = 1
      count1( 1) = rdLocalnEdges
      count1( 1) = rdLocalnEdges
      nferr = nf_get_vara_int(rd_ncid, rdVarIDindexToEdgeID, start1, count1, indexToEdgeID)
 
      start1(1) = 1
      count1( 1) = rdLocalnVertices
      count1( 1) = rdLocalnVertices
      nferr = nf_get_vara_double(rd_ncid, rdVarIDlatVertex, start1, count1, latVertex)
 
      start1(1) = 1
      count1( 1) = rdLocalnVertices
      count1( 1) = rdLocalnVertices
      nferr = nf_get_vara_double(rd_ncid, rdVarIDlonVertex, start1, count1, lonVertex)
 
      start1(1) = 1
      count1( 1) = rdLocalnVertices
      count1( 1) = rdLocalnVertices
      nferr = nf_get_vara_double(rd_ncid, rdVarIDxVertex, start1, count1, xVertex)
 
      start1(1) = 1
      count1( 1) = rdLocalnVertices
      count1( 1) = rdLocalnVertices
      nferr = nf_get_vara_double(rd_ncid, rdVarIDyVertex, start1, count1, yVertex)
 
      start1(1) = 1
      count1( 1) = rdLocalnVertices
      count1( 1) = rdLocalnVertices
      nferr = nf_get_vara_double(rd_ncid, rdVarIDzVertex, start1, count1, zVertex)
 
      start1(1) = 1
      count1( 1) = rdLocalnVertices
      count1( 1) = rdLocalnVertices
      nferr = nf_get_vara_int(rd_ncid, rdVarIDindexToVertexID, start1, count1, indexToVertexID)
 
      start2(2) = 1
      count2( 1) = rdLocalTWO
      count2( 2) = rdLocalnEdges
      nferr = nf_get_vara_int(rd_ncid, rdVarIDcellsOnEdge, start2, count2, cellsOnEdge)
 
      start1(1) = 1
      count1( 1) = rdLocalnCells
      count1( 1) = rdLocalnCells
      nferr = nf_get_vara_int(rd_ncid, rdVarIDnEdgesOnCell, start1, count1, nEdgesOnCell)
 
      start1(1) = 1
      count1( 1) = rdLocalnEdges
      count1( 1) = rdLocalnEdges
      nferr = nf_get_vara_int(rd_ncid, rdVarIDnEdgesOnEdge, start1, count1, nEdgesOnEdge)
 
      start2(2) = 1
      count2( 1) = rdLocalmaxEdges
      count2( 2) = rdLocalnCells
      nferr = nf_get_vara_int(rd_ncid, rdVarIDedgesOnCell, start2, count2, edgesOnCell)
 
      start2(2) = 1
      count2( 1) = rdLocalmaxEdges2
      count2( 2) = rdLocalnEdges
      nferr = nf_get_vara_int(rd_ncid, rdVarIDedgesOnEdge, start2, count2, edgesOnEdge)
 
      start2(2) = 1
      count2( 1) = rdLocalmaxEdges2
      count2( 2) = rdLocalnEdges
      nferr = nf_get_vara_double(rd_ncid, rdVarIDweightsOnEdge, start2, count2, weightsOnEdge)
 
      start1(1) = 1
      count1( 1) = rdLocalnEdges
      count1( 1) = rdLocalnEdges
      nferr = nf_get_vara_double(rd_ncid, rdVarIDdvEdge, start1, count1, dvEdge)
 
      start1(1) = 1
      count1( 1) = rdLocalnEdges
      count1( 1) = rdLocalnEdges
      nferr = nf_get_vara_double(rd_ncid, rdVarIDdcEdge, start1, count1, dcEdge)
 
      start1(1) = 1
      count1( 1) = rdLocalnEdges
      count1( 1) = rdLocalnEdges
      nferr = nf_get_vara_double(rd_ncid, rdVarIDangleEdge, start1, count1, angleEdge)
 
      start1(1) = 1
      count1( 1) = rdLocalnCells
      count1( 1) = rdLocalnCells
      nferr = nf_get_vara_double(rd_ncid, rdVarIDareaCell, start1, count1, areaCell)
 
      start1(1) = 1
      count1( 1) = rdLocalnVertices
      count1( 1) = rdLocalnVertices
      nferr = nf_get_vara_double(rd_ncid, rdVarIDareaTriangle, start1, count1, areaTriangle)
 
      start2(2) = 1
      count2( 1) = rdLocalmaxEdges
      count2( 2) = rdLocalnCells
      nferr = nf_get_vara_int(rd_ncid, rdVarIDcellsOnCell, start2, count2, cellsOnCell)
 
      start2(2) = 1
      count2( 1) = rdLocalmaxEdges
      count2( 2) = rdLocalnCells
      nferr = nf_get_vara_int(rd_ncid, rdVarIDverticesOnCell, start2, count2, verticesOnCell)
 
      start2(2) = 1
      count2( 1) = rdLocalTWO
      count2( 2) = rdLocalnEdges
      nferr = nf_get_vara_int(rd_ncid, rdVarIDverticesOnEdge, start2, count2, verticesOnEdge)
 
      start2(2) = 1
      count2( 1) = rdLocalvertexDegree
      count2( 2) = rdLocalnVertices
      nferr = nf_get_vara_int(rd_ncid, rdVarIDedgesOnVertex, start2, count2, edgesOnVertex)
 
      start2(2) = 1
      count2( 1) = rdLocalvertexDegree
      count2( 2) = rdLocalnVertices
      nferr = nf_get_vara_int(rd_ncid, rdVarIDcellsOnVertex, start2, count2, cellsOnVertex)
 
      start2(2) = 1
      count2( 1) = rdLocalvertexDegree
      count2( 2) = rdLocalnVertices
      nferr = nf_get_vara_double(rd_ncid, rdVarIDkiteAreasOnVertex, start2, count2, kiteAreasOnVertex)
 
      start1(1) = 1
      count1( 1) = rdLocalnEdges
      count1( 1) = rdLocalnEdges
      nferr = nf_get_vara_double(rd_ncid, rdVarIDfEdge, start1, count1, fEdge)
 
      start1(1) = 1
      count1( 1) = rdLocalnVertices
      count1( 1) = rdLocalnVertices
      nferr = nf_get_vara_double(rd_ncid, rdVarIDfVertex, start1, count1, fVertex)
 
      start1(1) = 1
      count1( 1) = rdLocalnCells
      count1( 1) = rdLocalnCells
      nferr = nf_get_vara_double(rd_ncid, rdVarIDh_s, start1, count1, h_s)
 
      start3(3) = time
      count3( 1) = rdLocalnVertLevels
      count3( 2) = rdLocalnEdges
      count3( 3) = 1
      nferr = nf_get_vara_double(rd_ncid, rdVarIDu, start3, count3, u)
 
      start3(3) = time
      count3( 1) = rdLocalnVertLevels
      count3( 2) = rdLocalnEdges
      count3( 3) = 1
      nferr = nf_get_vara_double(rd_ncid, rdVarIDv, start3, count3, v)
 
      start3(3) = time
      count3( 1) = rdLocalnVertLevels
      count3( 2) = rdLocalnCells
      count3( 3) = 1
      nferr = nf_get_vara_double(rd_ncid, rdVarIDh, start3, count3, h)
 
   end subroutine read_netcdf_fields
 
 
   subroutine read_netcdf_finalize()
 
      implicit none
 
      include 'netcdf.inc'
 
      integer :: nferr
 
      nferr = nf_close(rd_ncid)
 
   end subroutine read_netcdf_finalize
 
end module read_netcdf
