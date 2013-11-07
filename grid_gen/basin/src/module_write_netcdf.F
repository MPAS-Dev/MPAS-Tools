module write_netcdf
 
   integer :: wr_ncid
   integer :: wrDimIDTime
   integer :: wrDimIDnCells
   integer :: wrDimIDnEdges
   integer :: wrDimIDnVertices
   integer :: wrDimIDmaxEdges
   integer :: wrDimIDmaxEdges2
   integer :: wrDimIDTWO
   integer :: wrDimIDvertexDegree
   integer :: wrDimIDnVertLevels
   integer :: wrVarIDlatCell
   integer :: wrVarIDlonCell
   integer :: wrVarIDmeshDensity
   integer :: wrVarIDxCell
   integer :: wrVarIDyCell
   integer :: wrVarIDzCell
   integer :: wrVarIDindexToCellID
   integer :: wrVarIDlatEdge
   integer :: wrVarIDlonEdge
   integer :: wrVarIDxEdge
   integer :: wrVarIDyEdge
   integer :: wrVarIDzEdge
   integer :: wrVarIDindexToEdgeID
   integer :: wrVarIDlatVertex
   integer :: wrVarIDlonVertex
   integer :: wrVarIDxVertex
   integer :: wrVarIDyVertex
   integer :: wrVarIDzVertex
   integer :: wrVarIDindexToVertexID
   integer :: wrVarIDmaxLevelCell
   integer :: wrVarIDcellsOnEdge
   integer :: wrVarIDnEdgesOnCell
   integer :: wrVarIDnEdgesOnEdge
   integer :: wrVarIDedgesOnCell
   integer :: wrVarIDedgesOnEdge
   integer :: wrVarIDweightsOnEdge
   integer :: wrVarIDdvEdge
   integer :: wrVarIDdcEdge
   integer :: wrVarIDangleEdge
   integer :: wrVarIDareaCell
   integer :: wrVarIDareaTriangle
   integer :: wrVarIDcellsOnCell
   integer :: wrVarIDverticesOnCell
   integer :: wrVarIDverticesOnEdge
   integer :: wrVarIDedgesOnVertex
   integer :: wrVarIDcellsOnVertex
   integer :: wrVarIDkiteAreasOnVertex
   integer :: wrVarIDfEdge
   integer :: wrVarIDfVertex
   integer :: wrVarIDfCell
   integer :: wrVarIDbottomDepth
   integer :: wrVarIDnormalVelocity
   integer :: wrVarIDboundaryEdge
   integer :: wrVarIDboundaryVertex
   integer :: wrVarIDsurfaceWindStress
   integer :: wrVarIDlayerThickness
   integer :: wrVarIDdensity
   integer :: wrVarIDtemperature
   integer :: wrVarIDsalinity
   integer :: wrVarIDtracer1
   integer :: wrVarIDtemperatureRestore
   integer :: wrVarIDsalinityRestore
   integer :: wrVarIDrefBottomDepth
 
   integer :: wrLocalnCells
   integer :: wrLocalnEdges
   integer :: wrLocalnVertices
   integer :: wrLocalmaxEdges
   integer :: wrLocalnVertLevels
   integer :: wrLocalvertexDegree
 
   contains
 
   subroutine write_netcdf_init( &
                               nCells, &
                               nEdges, &
                               nVertices, &
                               maxEdges, &
                               nVertLevels, &
                               vertexDegree, &
                               sphere_radius, & 
                               on_a_sphere &
                               )
 
      implicit none
 
      include 'netcdf.inc'
 
      integer, intent(in) :: nCells
      integer, intent(in) :: nEdges
      integer, intent(in) :: nVertices
      integer, intent(in) :: maxEdges
      integer, intent(in) :: nVertLevels
      integer, intent(in) :: vertexDegree
      character (len=16) :: on_a_sphere
      real*8 :: sphere_radius

 
      integer :: nferr
      integer, dimension(10) :: dimlist
 
 
      wrLocalnCells = nCells
      wrLocalnEdges = nEdges
      wrLocalnVertices = nVertices
      wrLocalmaxEdges = maxEdges
      wrLocalnVertLevels = nVertLevels
      wrLocalvertexDegree = vertexDegree
 
      nferr = nf_create('ocean.nc', IOR(NF_CLOBBER,NF_64BIT_OFFSET), wr_ncid)
 
      !
      ! Define dimensions
      !
      nferr = nf_def_dim(wr_ncid, 'nCells', nCells, wrDimIDnCells)
      nferr = nf_def_dim(wr_ncid, 'nEdges', nEdges, wrDimIDnEdges)
      nferr = nf_def_dim(wr_ncid, 'nVertices', nVertices, wrDimIDnVertices)
      nferr = nf_def_dim(wr_ncid, 'maxEdges', maxEdges, wrDimIDmaxEdges)
      nferr = nf_def_dim(wr_ncid, 'maxEdges2', 2*maxEdges, wrDimIDmaxEdges2)
      nferr = nf_def_dim(wr_ncid, 'TWO', 2, wrDimIDTWO)
      nferr = nf_def_dim(wr_ncid, 'vertexDegree', vertexDegree, wrDimIDvertexDegree)
      nferr = nf_def_dim(wr_ncid, 'nVertLevels', nVertLevels, wrDimIDnVertLevels)
      nferr = nf_def_dim(wr_ncid, 'Time', NF_UNLIMITED, wrDimIDTime)
 
      !
      ! Define variables
      !
      dimlist( 1) = wrDimIDnCells
      nferr = nf_def_var(wr_ncid, 'latCell', NF_DOUBLE,  1, dimlist, wrVarIDlatCell)
      dimlist( 1) = wrDimIDnCells
      nferr = nf_def_var(wr_ncid, 'lonCell', NF_DOUBLE,  1, dimlist, wrVarIDlonCell)
      dimlist( 1) = wrDimIDnCells
      nferr = nf_def_var(wr_ncid, 'meshDensity', NF_DOUBLE,  1, dimlist, wrVarIDmeshDensity)
      dimlist( 1) = wrDimIDnCells
      nferr = nf_def_var(wr_ncid, 'xCell', NF_DOUBLE,  1, dimlist, wrVarIDxCell)
      dimlist( 1) = wrDimIDnCells
      nferr = nf_def_var(wr_ncid, 'yCell', NF_DOUBLE,  1, dimlist, wrVarIDyCell)
      dimlist( 1) = wrDimIDnCells
      nferr = nf_def_var(wr_ncid, 'zCell', NF_DOUBLE,  1, dimlist, wrVarIDzCell)
      dimlist( 1) = wrDimIDnCells
      nferr = nf_def_var(wr_ncid, 'indexToCellID', NF_INT,  1, dimlist, wrVarIDindexToCellID)
      dimlist( 1) = wrDimIDnEdges
      nferr = nf_def_var(wr_ncid, 'latEdge', NF_DOUBLE,  1, dimlist, wrVarIDlatEdge)
      dimlist( 1) = wrDimIDnEdges
      nferr = nf_def_var(wr_ncid, 'lonEdge', NF_DOUBLE,  1, dimlist, wrVarIDlonEdge)
      dimlist( 1) = wrDimIDnEdges
      nferr = nf_def_var(wr_ncid, 'xEdge', NF_DOUBLE,  1, dimlist, wrVarIDxEdge)
      dimlist( 1) = wrDimIDnEdges
      nferr = nf_def_var(wr_ncid, 'yEdge', NF_DOUBLE,  1, dimlist, wrVarIDyEdge)
      dimlist( 1) = wrDimIDnEdges
      nferr = nf_def_var(wr_ncid, 'zEdge', NF_DOUBLE,  1, dimlist, wrVarIDzEdge)
      dimlist( 1) = wrDimIDnEdges
      nferr = nf_def_var(wr_ncid, 'indexToEdgeID', NF_INT,  1, dimlist, wrVarIDindexToEdgeID)
      dimlist( 1) = wrDimIDnVertices
      nferr = nf_def_var(wr_ncid, 'latVertex', NF_DOUBLE,  1, dimlist, wrVarIDlatVertex)
      dimlist( 1) = wrDimIDnVertices
      nferr = nf_def_var(wr_ncid, 'lonVertex', NF_DOUBLE,  1, dimlist, wrVarIDlonVertex)
      dimlist( 1) = wrDimIDnVertices
      nferr = nf_def_var(wr_ncid, 'xVertex', NF_DOUBLE,  1, dimlist, wrVarIDxVertex)
      dimlist( 1) = wrDimIDnVertices
      nferr = nf_def_var(wr_ncid, 'yVertex', NF_DOUBLE,  1, dimlist, wrVarIDyVertex)
      dimlist( 1) = wrDimIDnVertices
      nferr = nf_def_var(wr_ncid, 'zVertex', NF_DOUBLE,  1, dimlist, wrVarIDzVertex)
      dimlist( 1) = wrDimIDnVertices
      nferr = nf_def_var(wr_ncid, 'indexToVertexID', NF_INT,  1, dimlist, wrVarIDindexToVertexID)
      dimlist( 1) = wrDimIDnCells
      nferr = nf_def_var(wr_ncid, 'maxLevelCell', NF_INT,  1, dimlist, wrVarIDmaxLevelCell)
      dimlist( 1) = wrDimIDTWO
      dimlist( 2) = wrDimIDnEdges
      nferr = nf_def_var(wr_ncid, 'cellsOnEdge', NF_INT,  2, dimlist, wrVarIDcellsOnEdge)
      dimlist( 1) = wrDimIDnCells
      nferr = nf_def_var(wr_ncid, 'nEdgesOnCell', NF_INT,  1, dimlist, wrVarIDnEdgesOnCell)
      dimlist( 1) = wrDimIDnEdges
      nferr = nf_def_var(wr_ncid, 'nEdgesOnEdge', NF_INT,  1, dimlist, wrVarIDnEdgesOnEdge)
      dimlist( 1) = wrDimIDmaxEdges
      dimlist( 2) = wrDimIDnCells
      nferr = nf_def_var(wr_ncid, 'edgesOnCell', NF_INT,  2, dimlist, wrVarIDedgesOnCell)
      dimlist( 1) = wrDimIDmaxEdges2
      dimlist( 2) = wrDimIDnEdges
      nferr = nf_def_var(wr_ncid, 'edgesOnEdge', NF_INT,  2, dimlist, wrVarIDedgesOnEdge)
      dimlist( 1) = wrDimIDmaxEdges2
      dimlist( 2) = wrDimIDnEdges
      nferr = nf_def_var(wr_ncid, 'weightsOnEdge', NF_DOUBLE,  2, dimlist, wrVarIDweightsOnEdge)
      dimlist( 1) = wrDimIDnEdges
      nferr = nf_def_var(wr_ncid, 'dvEdge', NF_DOUBLE,  1, dimlist, wrVarIDdvEdge)
      dimlist( 1) = wrDimIDnEdges
      nferr = nf_def_var(wr_ncid, 'dcEdge', NF_DOUBLE,  1, dimlist, wrVarIDdcEdge)
      dimlist( 1) = wrDimIDnEdges
      nferr = nf_def_var(wr_ncid, 'angleEdge', NF_DOUBLE,  1, dimlist, wrVarIDangleEdge)
      dimlist( 1) = wrDimIDnCells
      nferr = nf_def_var(wr_ncid, 'areaCell', NF_DOUBLE,  1, dimlist, wrVarIDareaCell)
      dimlist( 1) = wrDimIDnVertices
      nferr = nf_def_var(wr_ncid, 'areaTriangle', NF_DOUBLE,  1, dimlist, wrVarIDareaTriangle)
      dimlist( 1) = wrDimIDmaxEdges
      dimlist( 2) = wrDimIDnCells
      nferr = nf_def_var(wr_ncid, 'cellsOnCell', NF_INT,  2, dimlist, wrVarIDcellsOnCell)
      dimlist( 1) = wrDimIDmaxEdges
      dimlist( 2) = wrDimIDnCells
      nferr = nf_def_var(wr_ncid, 'verticesOnCell', NF_INT,  2, dimlist, wrVarIDverticesOnCell)
      dimlist( 1) = wrDimIDTWO
      dimlist( 2) = wrDimIDnEdges
      nferr = nf_def_var(wr_ncid, 'verticesOnEdge', NF_INT,  2, dimlist, wrVarIDverticesOnEdge)
      dimlist( 1) = wrDimIDvertexDegree
      dimlist( 2) = wrDimIDnVertices
      nferr = nf_def_var(wr_ncid, 'edgesOnVertex', NF_INT,  2, dimlist, wrVarIDedgesOnVertex)
      dimlist( 1) = wrDimIDvertexDegree
      dimlist( 2) = wrDimIDnVertices
      nferr = nf_def_var(wr_ncid, 'cellsOnVertex', NF_INT,  2, dimlist, wrVarIDcellsOnVertex)
      dimlist( 1) = wrDimIDvertexDegree
      dimlist( 2) = wrDimIDnVertices
      nferr = nf_def_var(wr_ncid, 'kiteAreasOnVertex', NF_DOUBLE,  2, dimlist, wrVarIDkiteAreasOnVertex)
      dimlist( 1) = wrDimIDnEdges
      nferr = nf_def_var(wr_ncid, 'fEdge', NF_DOUBLE,  1, dimlist, wrVarIDfEdge)
      dimlist( 1) = wrDimIDnVertices
      nferr = nf_def_var(wr_ncid, 'fVertex', NF_DOUBLE,  1, dimlist, wrVarIDfVertex)
      dimlist( 1) = wrDimIDnCells
      nferr = nf_def_var(wr_ncid, 'fCell', NF_DOUBLE,  1, dimlist, wrVarIDfCell)
      dimlist( 1) = wrDimIDnCells
      nferr = nf_def_var(wr_ncid, 'bottomDepth', NF_DOUBLE,  1, dimlist, wrVarIDbottomDepth)
      dimlist( 1) = wrDimIDnCells
      nferr = nf_def_var(wr_ncid, 'temperatureRestore', NF_DOUBLE,  1, dimlist, wrVarIDtemperatureRestore)
      dimlist( 1) = wrDimIDnCells
      nferr = nf_def_var(wr_ncid, 'salinityRestore', NF_DOUBLE,  1, dimlist, wrVarIDsalinityRestore)

      dimlist( 1) = wrDimIDnVertLevels
      nferr = nf_def_var(wr_ncid, 'refBottomDepth', NF_DOUBLE,  1, dimlist, wrVarIDrefBottomDepth)
      dimlist( 1) = wrDimIDnVertLevels
      dimlist( 2) = wrDimIDnEdges
      dimlist( 3) = wrDimIDTime
      nferr = nf_def_var(wr_ncid, 'normalVelocity', NF_DOUBLE,  3, dimlist, wrVarIDnormalVelocity)
      dimlist( 1) = wrDimIDnVertLevels
      dimlist( 2) = wrDimIDnEdges
      nferr = nf_def_var(wr_ncid, 'boundaryEdge', NF_INT,  2, dimlist, wrVarIDboundaryEdge)
      dimlist( 1) = wrDimIDnVertLevels
      dimlist( 2) = wrDimIDnVertices
      nferr = nf_def_var(wr_ncid, 'boundaryVertex', NF_INT,  2, dimlist, wrVarIDboundaryVertex)

      dimlist( 1) = wrDimIDnEdges
      nferr = nf_def_var(wr_ncid, 'surfaceWindStress', NF_DOUBLE,  1, dimlist, wrVarIDsurfaceWindStress)

      dimlist( 1) = wrDimIDnVertLevels
      dimlist( 2) = wrDimIDnCells
      dimlist( 3) = wrDimIDTime
      nferr = nf_def_var(wr_ncid, 'layerThickness', NF_DOUBLE,  3, dimlist, &
			wrVarIDlayerThickness)
   !  dimlist( 1) = wrDimIDnVertLevels
   !  dimlist( 2) = wrDimIDnCells
   !  dimlist( 3) = wrDimIDTime
   !  nferr = nf_def_var(wr_ncid, 'density', NF_DOUBLE,  3, dimlist, wrVarIDdensity) 
      dimlist( 1) = wrDimIDnVertLevels
      dimlist( 2) = wrDimIDnCells
      dimlist( 3) = wrDimIDTime
      nferr = nf_def_var(wr_ncid, 'temperature', NF_DOUBLE,  3, dimlist, wrVarIDtemperature)
      dimlist( 1) = wrDimIDnVertLevels
      dimlist( 2) = wrDimIDnCells
      dimlist( 3) = wrDimIDTime
      nferr = nf_def_var(wr_ncid, 'salinity', NF_DOUBLE,  3, dimlist, wrVarIDsalinity)
      dimlist( 1) = wrDimIDnVertLevels
      dimlist( 2) = wrDimIDnCells
      dimlist( 3) = wrDimIDTime
      ! If you do not want tracer1 in your input file, simply comment out these two lines (one of two)
      nferr = nf_def_var(wr_ncid, 'tracer1', NF_DOUBLE,  3, dimlist, wrVarIDtracer1)

 
      nferr = nf_put_att_text(wr_ncid, NF_GLOBAL, 'on_a_sphere', 16, on_a_sphere)
      nferr = nf_put_att_double(wr_ncid, NF_GLOBAL, 'sphere_radius', NF_DOUBLE, 1, sphere_radius)

      nferr = nf_enddef(wr_ncid)

   end subroutine write_netcdf_init
 
 
   subroutine write_netcdf_fields( &
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
                                  maxLevelCell, &
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
                                  fCell, &
                                  bottomDepth, &
                                  boundaryEdge, &
                                  boundaryVertex, &
                                  surfaceWindStress, &
                                  normalVelocity, &
                                  layerThickness, &
                                  density, &
                                  temperature, &
                                  salinity, &
                                  tracer1, &
                                  temperatureRestore, &
                                  salinityRestore, &
                                  refBottomDepth &
                                 )
 
      implicit none
 
      include 'netcdf.inc'
 
      integer, intent(in) :: time
      real (kind=8), dimension(:), intent(in) :: latCell
      real (kind=8), dimension(:), intent(in) :: lonCell
      real (kind=8), dimension(:), intent(in) :: meshDensity
      real (kind=8), dimension(:), intent(in) :: xCell
      real (kind=8), dimension(:), intent(in) :: yCell
      real (kind=8), dimension(:), intent(in) :: zCell
      integer, dimension(:), intent(in) :: indexToCellID
      real (kind=8), dimension(:), intent(in) :: latEdge
      real (kind=8), dimension(:), intent(in) :: lonEdge
      real (kind=8), dimension(:), intent(in) :: xEdge
      real (kind=8), dimension(:), intent(in) :: yEdge
      real (kind=8), dimension(:), intent(in) :: zEdge
      integer, dimension(:), intent(in) :: indexToEdgeID
      real (kind=8), dimension(:), intent(in) :: latVertex
      real (kind=8), dimension(:), intent(in) :: lonVertex
      real (kind=8), dimension(:), intent(in) :: xVertex
      real (kind=8), dimension(:), intent(in) :: yVertex
      real (kind=8), dimension(:), intent(in) :: zVertex
      integer, dimension(:), intent(in) :: indexToVertexID
      integer, dimension(:), intent(in) :: maxLevelCell
      integer, dimension(:,:), intent(in) :: cellsOnEdge
      integer, dimension(:), intent(in) :: nEdgesOnCell
      integer, dimension(:), intent(in) :: nEdgesOnEdge
      integer, dimension(:,:), intent(in) :: edgesOnCell
      integer, dimension(:,:), intent(in) :: edgesOnEdge
      real (kind=8), dimension(:,:), intent(in) :: weightsOnEdge
      real (kind=8), dimension(:), intent(in) :: dvEdge
      real (kind=8), dimension(:), intent(in) :: dcEdge
      real (kind=8), dimension(:), intent(in) :: angleEdge
      real (kind=8), dimension(:), intent(in) :: areaCell
      real (kind=8), dimension(:), intent(in) :: areaTriangle
      integer, dimension(:,:), intent(in) :: cellsOnCell
      integer, dimension(:,:), intent(in) :: verticesOnCell
      integer, dimension(:,:), intent(in) :: verticesOnEdge
      integer, dimension(:,:), intent(in) :: edgesOnVertex
      integer, dimension(:,:), intent(in) :: cellsOnVertex
      real (kind=8), dimension(:,:), intent(in) :: kiteAreasOnVertex
      real (kind=8), dimension(:), intent(in) :: fEdge
      real (kind=8), dimension(:), intent(in) :: fVertex
      real (kind=8), dimension(:), intent(in) :: fCell
      real (kind=8), dimension(:), intent(in) :: bottomDepth
      integer, dimension(:,:), intent(in) :: boundaryEdge
      integer, dimension(:,:), intent(in) :: boundaryVertex
      real (kind=8), dimension(:), intent(in) :: surfaceWindStress
      real (kind=8), dimension(:,:,:), intent(in) :: normalVelocity
      real (kind=8), dimension(:,:,:), intent(in) :: layerThickness
      real (kind=8), dimension(:,:,:), intent(in) :: density
      real (kind=8), dimension(:,:,:), intent(in) :: temperature
      real (kind=8), dimension(:,:,:), intent(in) :: salinity
      real (kind=8), dimension(:,:,:), intent(in) :: tracer1
      real (kind=8), dimension(:), intent(in) :: temperatureRestore
      real (kind=8), dimension(:), intent(in) :: salinityRestore
      real (kind=8), dimension(:), intent(in) :: refBottomDepth

 
      integer :: nferr
      integer, dimension(1) :: start1, count1
      integer, dimension(2) :: start2, count2
      integer, dimension(3) :: start3, count3
      integer, dimension(4) :: start4, count4
 
      start1(1) = 1
 
      start2(1) = 1
      start2(2) = 1
 
      start3(1) = 1
      start3(2) = 1
      start3(3) = 1
 
      start4(1) = 1
      start4(2) = 1
      start4(3) = 1
      start4(4) = 1
 
      start1(1) = 1
      count1( 1) = wrLocalnCells
      nferr = nf_put_vara_double(wr_ncid, wrVarIDlatCell, start1, count1, latCell)
 
      start1(1) = 1
      count1( 1) = wrLocalnCells
      nferr = nf_put_vara_double(wr_ncid, wrVarIDlonCell, start1, count1, lonCell)

      start1(1) = 1
      count1( 1) = wrLocalnCells
      nferr = nf_put_vara_double(wr_ncid, wrVarIDmeshDensity, start1, count1, meshDensity)
 
      start1(1) = 1
      count1( 1) = wrLocalnCells
      nferr = nf_put_vara_double(wr_ncid, wrVarIDxCell, start1, count1, xCell)
 
      start1(1) = 1
      count1( 1) = wrLocalnCells
      nferr = nf_put_vara_double(wr_ncid, wrVarIDyCell, start1, count1, yCell)
 
      start1(1) = 1
      count1( 1) = wrLocalnCells
      nferr = nf_put_vara_double(wr_ncid, wrVarIDzCell, start1, count1, zCell)
 
      start1(1) = 1
      count1( 1) = wrLocalnCells
      nferr = nf_put_vara_int(wr_ncid, wrVarIDindexToCellID, start1, count1, indexToCellID)
 
      start1(1) = 1
      count1( 1) = wrLocalnEdges
      nferr = nf_put_vara_double(wr_ncid, wrVarIDlatEdge, start1, count1, latEdge)
 
      start1(1) = 1
      count1( 1) = wrLocalnEdges
      nferr = nf_put_vara_double(wr_ncid, wrVarIDlonEdge, start1, count1, lonEdge)
 
      start1(1) = 1
      count1( 1) = wrLocalnEdges
      nferr = nf_put_vara_double(wr_ncid, wrVarIDxEdge, start1, count1, xEdge)
 
      start1(1) = 1
      count1( 1) = wrLocalnEdges
      nferr = nf_put_vara_double(wr_ncid, wrVarIDyEdge, start1, count1, yEdge)
 
      start1(1) = 1
      count1( 1) = wrLocalnEdges
      nferr = nf_put_vara_double(wr_ncid, wrVarIDzEdge, start1, count1, zEdge)
 
      start1(1) = 1
      count1( 1) = wrLocalnEdges
      nferr = nf_put_vara_int(wr_ncid, wrVarIDindexToEdgeID, start1, count1, indexToEdgeID)
 
      start1(1) = 1
      count1( 1) = wrLocalnVertices
      nferr = nf_put_vara_double(wr_ncid, wrVarIDlatVertex, start1, count1, latVertex)
 
      start1(1) = 1
      count1( 1) = wrLocalnVertices
      nferr = nf_put_vara_double(wr_ncid, wrVarIDlonVertex, start1, count1, lonVertex)
 
      start1(1) = 1
      count1( 1) = wrLocalnVertices
      nferr = nf_put_vara_double(wr_ncid, wrVarIDxVertex, start1, count1, xVertex)
 
      start1(1) = 1
      count1( 1) = wrLocalnVertices
      nferr = nf_put_vara_double(wr_ncid, wrVarIDyVertex, start1, count1, yVertex)
 
      start1(1) = 1
      count1( 1) = wrLocalnVertices
      nferr = nf_put_vara_double(wr_ncid, wrVarIDzVertex, start1, count1, zVertex)
 
      start1(1) = 1
      count1( 1) = wrLocalnVertices
      nferr = nf_put_vara_int(wr_ncid, wrVarIDindexToVertexID, start1, count1, indexToVertexID)
 
      start1(1) = 1
      count1( 1) = wrLocalnCells
      nferr = nf_put_vara_int(wr_ncid, wrVarIDmaxLevelCell, start1, count1, maxLevelCell)

      start2(2) = 1
      count2( 1) = 2
      count2( 2) = wrLocalnEdges
      nferr = nf_put_vara_int(wr_ncid, wrVarIDcellsOnEdge, start2, count2, cellsOnEdge)
 
      start1(1) = 1
      count1( 1) = wrLocalnCells
      nferr = nf_put_vara_int(wr_ncid, wrVarIDnEdgesOnCell, start1, count1, nEdgesOnCell)
 
      start1(1) = 1
      count1( 1) = wrLocalnEdges
      nferr = nf_put_vara_int(wr_ncid, wrVarIDnEdgesOnEdge, start1, count1, nEdgesOnEdge)
 
      start2(2) = 1
      count2( 1) = wrLocalmaxEdges
      count2( 2) = wrLocalnCells
      nferr = nf_put_vara_int(wr_ncid, wrVarIDedgesOnCell, start2, count2, edgesOnCell)
 
      start2(2) = 1
      count2( 1) = 2*wrLocalmaxEdges
      count2( 2) = wrLocalnEdges
      nferr = nf_put_vara_int(wr_ncid, wrVarIDedgesOnEdge, start2, count2, edgesOnEdge)
 
      start2(2) = 1
      count2( 1) = 2*wrLocalmaxEdges
      count2( 2) = wrLocalnEdges
      nferr = nf_put_vara_double(wr_ncid, wrVarIDweightsOnEdge, start2, count2, weightsOnEdge)
 
      start1(1) = 1
      count1( 1) = wrLocalnEdges
      nferr = nf_put_vara_double(wr_ncid, wrVarIDdvEdge, start1, count1, dvEdge)
 
      start1(1) = 1
      count1( 1) = wrLocalnEdges
      nferr = nf_put_vara_double(wr_ncid, wrVarIDdcEdge, start1, count1, dcEdge)
 
      start1(1) = 1
      count1( 1) = wrLocalnEdges
      nferr = nf_put_vara_double(wr_ncid, wrVarIDangleEdge, start1, count1, angleEdge)
 
      start1(1) = 1
      count1( 1) = wrLocalnCells
      nferr = nf_put_vara_double(wr_ncid, wrVarIDareaCell, start1, count1, areaCell)
 
      start1(1) = 1
      count1( 1) = wrLocalnVertices
      nferr = nf_put_vara_double(wr_ncid, wrVarIDareaTriangle, start1, count1, areaTriangle)
 
      start2(2) = 1
      count2( 1) = wrLocalmaxEdges
      count2( 2) = wrLocalnCells
      nferr = nf_put_vara_int(wr_ncid, wrVarIDcellsOnCell, start2, count2, cellsOnCell)
 
      start2(2) = 1
      count2( 1) = wrLocalmaxEdges
      count2( 2) = wrLocalnCells
      nferr = nf_put_vara_int(wr_ncid, wrVarIDverticesOnCell, start2, count2, verticesOnCell)
 
      start2(2) = 1
      count2( 1) = 2
      count2( 2) = wrLocalnEdges
      nferr = nf_put_vara_int(wr_ncid, wrVarIDverticesOnEdge, start2, count2, verticesOnEdge)
 
      start2(2) = 1
      count2( 1) = wrLocalvertexDegree
      count2( 2) = wrLocalnVertices
      nferr = nf_put_vara_int(wr_ncid, wrVarIDedgesOnVertex, start2, count2, edgesOnVertex)
 
      start2(2) = 1
      count2( 1) = wrLocalvertexDegree
      count2( 2) = wrLocalnVertices
      nferr = nf_put_vara_int(wr_ncid, wrVarIDcellsOnVertex, start2, count2, cellsOnVertex)
 
      start2(2) = 1
      count2( 1) = wrLocalvertexDegree
      count2( 2) = wrLocalnVertices
      nferr = nf_put_vara_double(wr_ncid, wrVarIDkiteAreasOnVertex, start2, count2, kiteAreasOnVertex)
 
      start1(1) = 1
      count1( 1) = wrLocalnEdges
      nferr = nf_put_vara_double(wr_ncid, wrVarIDfEdge, start1, count1, fEdge)
 
      start1(1) = 1
      count1( 1) = wrLocalnVertices
      nferr = nf_put_vara_double(wr_ncid, wrVarIDfVertex, start1, count1, fVertex)

      start1(1) = 1
      count1( 1) = wrLocalnCells
      nferr = nf_put_vara_double(wr_ncid, wrVarIDfCell, start1, count1, fCell)
 
      start1(1) = 1
      count1( 1) = wrLocalnCells
      nferr = nf_put_vara_double(wr_ncid, wrVarIDbottomDepth, start1, count1, bottomDepth)

      start2(2) = 1
      count2( 1) = wrLocalnVertLevels
      count2( 2) = wrLocalnEdges
      nferr = nf_put_vara_int(wr_ncid, wrVarIDboundaryEdge, start2, count2, boundaryEdge)

      start2(2) = 1
      count2( 1) = wrLocalnVertLevels
      count2( 2) = wrLocalnVertices
      nferr = nf_put_vara_int(wr_ncid, wrVarIDboundaryVertex, start2, count2, boundaryVertex)

      start1(1) = 1
      count1( 1) = wrLocalnEdges
      nferr = nf_put_vara_double(wr_ncid, wrVarIDsurfaceWindStress, start1, count1, surfaceWindStress)

      start1(1) = 1
      count1( 1) = wrLocalnCells
      nferr = nf_put_vara_double(wr_ncid, wrVarIDtemperatureRestore, start1, count1, temperatureRestore)

      start1(1) = 1
      count1( 1) = wrLocalnCells
      nferr = nf_put_vara_double(wr_ncid, wrVarIDsalinityRestore, start1, count1, salinityRestore)

      start1(1) = 1
      count1( 1) = wrLocalnVertLevels
      nferr = nf_put_vara_double(wr_ncid, wrVarIDrefBottomDepth, start1, count1, refBottomDepth)
 
      start3(3) = time
      count3( 1) = wrLocalnVertLevels
      count3( 2) = wrLocalnEdges
      count3( 3) = 1
      nferr = nf_put_vara_double(wr_ncid, wrVarIDnormalVelocity, start3, count3, normalVelocity) 
 
      start3(3) = time
      count3( 1) = wrLocalnVertLevels
      count3( 2) = wrLocalnCells
      count3( 3) = 1
      nferr = nf_put_vara_double(wr_ncid, wrVarIDlayerThickness, start3, count3, layerThickness) 
 
   !  start3(3) = time
   !  count3( 1) = wrLocalnVertLevels
   !  count3( 2) = wrLocalnCells
   !  count3( 3) = 1
   !  nferr = nf_put_vara_double(wr_ncid, wrVarIDdensity, start3, count3, density) 
 
      start3(3) = time
      count3( 1) = wrLocalnVertLevels
      count3( 2) = wrLocalnCells
      count3( 3) = 1
      nferr = nf_put_vara_double(wr_ncid, wrVarIDtemperature, start3, count3, temperature)

      start3(3) = time
      count3( 1) = wrLocalnVertLevels
      count3( 2) = wrLocalnCells
      count3( 3) = 1
      nferr = nf_put_vara_double(wr_ncid, wrVarIDsalinity, start3, count3, salinity)
 
      start3(3) = time
      count3( 1) = wrLocalnVertLevels
      count3( 2) = wrLocalnCells
      count3( 3) = 1
      ! If you do not want tracer1 in your input file, simply comment out these two lines (two of two)
      nferr = nf_put_vara_double(wr_ncid, wrVarIDtracer1, start3, count3, tracer1)
 
   end subroutine write_netcdf_fields
 
 
   subroutine write_netcdf_finalize()
 
      implicit none
 
      include 'netcdf.inc'
 
      integer :: nferr
 
      nferr = nf_close(wr_ncid)
 
   end subroutine write_netcdf_finalize
 
end module write_netcdf
