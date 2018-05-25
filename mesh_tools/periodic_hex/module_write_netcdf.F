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
   integer :: wrDimIDnTracers
   integer :: wrVarIDlatCell
   integer :: wrVarIDlonCell
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
   integer :: wrVarIDh_s
   integer :: wrVarIDu
   integer :: wrVarIDuBC
   integer :: wrVarIDv
   integer :: wrVarIDh
   integer :: wrVarIDvh
   integer :: wrVarIDcirculation
   integer :: wrVarIDvorticity
   integer :: wrVarIDke
   integer :: wrVarIDtracers
 
   integer :: wrLocalnCells
   integer :: wrLocalnEdges
   integer :: wrLocalnVertices
   integer :: wrLocalmaxEdges
   integer :: wrLocalnVertLevels
   integer :: wrLocalnTracers
 
   contains
 
   subroutine write_netcdf_init( &
                               nCells, &
                               nEdges, &
                               nVertices, &
                               maxEdges, &
                               nVertLevels, &
                               nTracers, &
                               vertexDegree, &
                               dc, &
                               nx, &
                               ny &
                               )
 
      implicit none
 
      include 'netcdf.inc'
 
      integer, intent(in) :: nCells
      integer, intent(in) :: nEdges
      integer, intent(in) :: nVertices
      integer, intent(in) :: maxEdges
      integer, intent(in) :: nVertLevels
      integer, intent(in) :: nTracers
      integer, intent(in) :: vertexDegree
      real (kind=8), intent(in) :: dc
      integer, intent(in) :: nx
      integer, intent(in) :: ny
 
      integer :: nferr
      integer, dimension(10) :: dimlist
      character (len=16) :: on_a_sphere
      character (len=16) :: is_periodic
      real (kind=8) :: sphere_radius
      real (kind=8) :: x_period, y_period
 
 
      wrLocalnCells = nCells
      wrLocalnEdges = nEdges
      wrLocalnVertices = nVertices
      wrLocalmaxEdges = maxEdges
      wrLocalnVertLevels = nVertLevels
      wrLocalnTracers = nTracers

      on_a_sphere = 'NO'
      is_periodic = 'YES'
      sphere_radius = 0.0
      x_period = (nx) * dc
      y_period = (ny) * (dc * sqrt(3.0)) / 2.0
 
      nferr = nf_create('grid.nc', IOR(NF_CLOBBER,NF_64BIT_OFFSET), wr_ncid)
 
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
      nferr = nf_def_dim(wr_ncid, 'nTracers', nTracers, wrDimIDnTracers)
      nferr = nf_def_dim(wr_ncid, 'Time', NF_UNLIMITED, wrDimIDTime)


      !
      ! Define attributes
      !
      nferr = nf_put_att_text(wr_ncid, NF_GLOBAL, 'on_a_sphere', 16, on_a_sphere)
      nferr = nf_put_att_text(wr_ncid, NF_GLOBAL, 'is_periodic', 16, is_periodic)
      nferr = nf_put_att_double(wr_ncid, NF_GLOBAL, 'sphere_radius', NF_DOUBLE, 1, sphere_radius)
      nferr = nf_put_att_double(wr_ncid, NF_GLOBAL, 'x_period', NF_DOUBLE, 1, x_period)
      nferr = nf_put_att_double(wr_ncid, NF_GLOBAL, 'y_period', NF_DOUBLE, 1, y_period)

 
      !
      ! Define variables
      !
      dimlist( 1) = wrDimIDnCells
      nferr = nf_def_var(wr_ncid, 'latCell', NF_DOUBLE,  1, dimlist, wrVarIDlatCell)
      dimlist( 1) = wrDimIDnCells
      nferr = nf_def_var(wr_ncid, 'lonCell', NF_DOUBLE,  1, dimlist, wrVarIDlonCell)
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
      nferr = nf_def_var(wr_ncid, 'h_s', NF_DOUBLE,  1, dimlist, wrVarIDh_s)
      dimlist( 1) = wrDimIDnVertLevels
      dimlist( 2) = wrDimIDnEdges
      dimlist( 3) = wrDimIDTime
      nferr = nf_def_var(wr_ncid, 'u', NF_DOUBLE,  3, dimlist, wrVarIDu)
      dimlist( 1) = wrDimIDnVertLevels
      dimlist( 2) = wrDimIDnEdges
      nferr = nf_def_var(wr_ncid, 'uBC', NF_INT,  2, dimlist, wrVarIDuBC)
      dimlist( 1) = wrDimIDnVertLevels
      dimlist( 2) = wrDimIDnEdges
      dimlist( 3) = wrDimIDTime
      nferr = nf_def_var(wr_ncid, 'v', NF_DOUBLE,  3, dimlist, wrVarIDv)
      dimlist( 1) = wrDimIDnVertLevels
      dimlist( 2) = wrDimIDnCells
      dimlist( 3) = wrDimIDTime
      nferr = nf_def_var(wr_ncid, 'h', NF_DOUBLE,  3, dimlist, wrVarIDh)
      dimlist( 1) = wrDimIDnVertLevels
      dimlist( 2) = wrDimIDnEdges
      dimlist( 3) = wrDimIDTime
      nferr = nf_def_var(wr_ncid, 'vh', NF_DOUBLE,  3, dimlist, wrVarIDvh)
      dimlist( 1) = wrDimIDnVertLevels
      dimlist( 2) = wrDimIDnVertices
      dimlist( 3) = wrDimIDTime
      nferr = nf_def_var(wr_ncid, 'circulation', NF_DOUBLE,  3, dimlist, wrVarIDcirculation)
      dimlist( 1) = wrDimIDnVertLevels
      dimlist( 2) = wrDimIDnVertices
      dimlist( 3) = wrDimIDTime
      nferr = nf_def_var(wr_ncid, 'vorticity', NF_DOUBLE,  3, dimlist, wrVarIDvorticity)
      dimlist( 1) = wrDimIDnVertLevels
      dimlist( 2) = wrDimIDnCells
      dimlist( 3) = wrDimIDTime
      nferr = nf_def_var(wr_ncid, 'ke', NF_DOUBLE,  3, dimlist, wrVarIDke)
      dimlist( 1) = wrDimIDnTracers
      dimlist( 2) = wrDimIDnVertLevels
      dimlist( 3) = wrDimIDnCells
      dimlist( 4) = wrDimIDTime
      nferr = nf_def_var(wr_ncid, 'tracers', NF_DOUBLE,  4, dimlist, wrVarIDtracers)
 
      nferr = nf_enddef(wr_ncid)
 
   end subroutine write_netcdf_init
 
 
   subroutine write_netcdf_fields( &
                                  time, &
                                  latCell, &
                                  lonCell, &
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
 
      implicit none
 
      include 'netcdf.inc'
 
      integer, intent(in) :: time
      real (kind=8), dimension(:), intent(in) :: latCell
      real (kind=8), dimension(:), intent(in) :: lonCell
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
      real (kind=8), dimension(:), intent(in) :: h_s
      integer, dimension(:,:), intent(in) :: uBC
      real (kind=8), dimension(:,:,:), intent(in) :: u
      real (kind=8), dimension(:,:,:), intent(in) :: v
      real (kind=8), dimension(:,:,:), intent(in) :: h
      real (kind=8), dimension(:,:,:), intent(in) :: vh
      real (kind=8), dimension(:,:,:), intent(in) :: circulation
      real (kind=8), dimension(:,:,:), intent(in) :: vorticity
      real (kind=8), dimension(:,:,:), intent(in) :: ke
      real (kind=8), dimension(:,:,:,:), intent(in) :: tracers
 
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
      count2( 1) = 3
      count2( 2) = wrLocalnVertices
      nferr = nf_put_vara_int(wr_ncid, wrVarIDedgesOnVertex, start2, count2, edgesOnVertex)
 
      start2(2) = 1
      count2( 1) = 3
      count2( 2) = wrLocalnVertices
      nferr = nf_put_vara_int(wr_ncid, wrVarIDcellsOnVertex, start2, count2, cellsOnVertex)
 
      start2(2) = 1
      count2( 1) = 3
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
      nferr = nf_put_vara_double(wr_ncid, wrVarIDh_s, start1, count1, h_s)

      start2(2) = 1
      count2( 1) = wrLocalnVertLevels
      count2( 2) = wrLocalnEdges
      nferr = nf_put_vara_int(wr_ncid, wrVarIDuBC, start2, count2, u)
 
      start3(3) = time
      count3( 1) = wrLocalnVertLevels
      count3( 2) = wrLocalnEdges
      count3( 3) = 1
      nferr = nf_put_vara_double(wr_ncid, wrVarIDu, start3, count3, u)
 
      start3(3) = time
      count3( 1) = wrLocalnVertLevels
      count3( 2) = wrLocalnEdges
      count3( 3) = 1
      nferr = nf_put_vara_double(wr_ncid, wrVarIDv, start3, count3, v)
 
      start3(3) = time
      count3( 1) = wrLocalnVertLevels
      count3( 2) = wrLocalnCells
      count3( 3) = 1
      nferr = nf_put_vara_double(wr_ncid, wrVarIDh, start3, count3, h)
 
      start3(3) = time
      count3( 1) = wrLocalnVertLevels
      count3( 2) = wrLocalnEdges
      count3( 3) = 1
      nferr = nf_put_vara_double(wr_ncid, wrVarIDvh, start3, count3, vh)
 
      start3(3) = time
      count3( 1) = wrLocalnVertLevels
      count3( 2) = wrLocalnVertices
      count3( 3) = 1
      nferr = nf_put_vara_double(wr_ncid, wrVarIDcirculation, start3, count3, circulation)
 
      start3(3) = time
      count3( 1) = wrLocalnVertLevels
      count3( 2) = wrLocalnVertices
      count3( 3) = 1
      nferr = nf_put_vara_double(wr_ncid, wrVarIDvorticity, start3, count3, vorticity)
 
      start3(3) = time
      count3( 1) = wrLocalnVertLevels
      count3( 2) = wrLocalnCells
      count3( 3) = 1
      nferr = nf_put_vara_double(wr_ncid, wrVarIDke, start3, count3, ke)
 
      start4(4) = time
      count4( 1) = wrLocalnTracers
      count4( 2) = wrLocalnVertLevels
      count4( 3) = wrLocalnCells
      count4( 4) = 1
      nferr = nf_put_vara_double(wr_ncid, wrVarIDtracers, start4, count4, tracers)
 
 
   end subroutine write_netcdf_fields
 
 
   subroutine write_netcdf_finalize()
 
      implicit none
 
      include 'netcdf.inc'
 
      integer :: nferr
 
      nferr = nf_close(wr_ncid)
 
   end subroutine write_netcdf_finalize
 
end module write_netcdf
