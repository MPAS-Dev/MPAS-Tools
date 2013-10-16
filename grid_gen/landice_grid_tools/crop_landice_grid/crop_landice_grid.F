program crop_landice_grid

use read_netcdf
use write_netcdf
use utilities
use cullLoops

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Program: crop_landice_grid.F
!
! This program culls cells from a mesh.  Presumably you will use it to remove non-ice cells.
! You define what cells to keep with a 1 value in an integer field called keepCellMask that you must add to a grid file named land_ice_grid.nc.
! It generates a culled version called land_ice_grid_culled.nc.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none

! original grid variables
integer :: time, nCells, nEdges, nVertices
integer :: maxEdges, maxEdges2, TWO, vertexDegree, nVertLevels
integer, allocatable, dimension(:) :: indexToCellID, indexToEdgeID, indexToVertexID
real, allocatable, dimension(:) :: xCell, yCell, zCell, latCell, lonCell, meshDensity
real, allocatable, dimension(:) :: xEdge, yEdge, zEdge, latEdge, lonEdge
real, allocatable, dimension(:) :: xVertex, yVertex, zVertex, latVertex, lonVertex
integer, allocatable, dimension(:) :: nEdgesOnCell, nEdgesOnEdge
integer, allocatable, dimension(:,:) :: cellsOnCell, edgesOnCell, verticesOnCell
integer, allocatable, dimension(:,:) :: cellsOnEdge, verticesOnEdge, edgesOnEdge
integer, allocatable, dimension(:,:) :: cellsOnVertex, edgesOnVertex
real, allocatable, dimension(:) :: areaCell, areaTriangle, dcEdge, dvEdge, angleEdge
real, allocatable, dimension(:,:) :: kiteAreasOnVertex, weightsOnEdge

integer nlon, nlat, ndepth
real(kind=4), allocatable, dimension(:) :: t_lon, t_lat

real, dimension(40) :: dz

   real (kind=8) :: ymid, ytmp, ymax, xmid, xloc, yloc, pert, ymin, distance, r, c1(3), c2(3)
   real (kind=8) :: latmid, lattmp, latmax, latmin
   integer :: cell1, cell2

! new grid variables
integer :: nCellsNew, nEdgesNew, nVerticesNew
integer :: maxEdgesNew, maxEdges2New, TWONew, vertexDegreeNew, nVertLevelsNew
integer, allocatable, dimension(:) :: indexToCellIDNew, indexToEdgeIDNew, indexToVertexIDNew
real, allocatable, dimension(:) :: xCellNew, yCellNew, zCellNew, latCellNew, lonCellNew, meshDensityNew, meshSpacingNew
real, allocatable, dimension(:) :: xEdgeNew, yEdgeNew, zEdgeNew, latEdgeNew, lonEdgeNew
real, allocatable, dimension(:) :: xVertexNew, yVertexNew, zVertexNew, latVertexNew, lonVertexNew
integer, allocatable, dimension(:) :: nEdgesOnCellNew, nEdgesOnEdgeNew, flipVerticesOnEdgeOrdering
integer, allocatable, dimension(:,:) :: cellsOnCellNew, edgesOnCellNew, verticesOnCellNew
integer, allocatable, dimension(:,:) :: cellsOnEdgeNew, verticesOnEdgeNew, edgesOnEdgeNew
integer, allocatable, dimension(:,:) :: cellsOnVertexNew, edgesOnVertexNew
integer, allocatable, dimension(:,:) :: boundaryEdgeNew, boundaryVertexNew
real, allocatable, dimension(:) :: areaCellNew, areaTriangleNew, dcEdgeNew, dvEdgeNew, angleEdgeNew
real, allocatable, dimension(:,:) :: kiteAreasOnVertexNew, weightsOnEdgeNew, normalsNew

integer, allocatable, dimension(:) :: keepCellMask

! mapping variables
integer, allocatable, dimension(:) :: kmt
integer, allocatable, dimension(:) :: cellMap, edgeMap, vertexMap

! work variables
integer :: i,j,jNew,k,jEdge,jEdgeNew,iVertex1New,iVertex2New,iCell1New,iCell2New
integer :: iCell, iCell1, iCell2, iCell3, iEdge, iVertex, iVertex1, iVertex2
integer :: iCellNew, iEdgeNew, iVertexNew, ndata, jCell1, jCell2, jCell, iter
real :: xin, yin, zin, ulon, ulat, ux, uy, uz, rlon, rlat, temp_t, temp_s


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Namelist variables
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Variables in namelist file
character (len=32) :: on_a_sphere, zLevel_thickness,bottom_topography, initial_conditions
logical :: expand_from_unit_sphere, eliminate_inland_seas, &
   check_mesh, &
   cut_domain_from_sphere, solid_boundary_in_y, solid_boundary_in_x

integer :: nVertLevelsMOD
real (kind=8) :: sphere_radius

!! specify namelist
!namelist /landicecrop/ on_a_sphere, sphere_radius, &
!   expand_from_unit_sphere, &
!   eliminate_inland_seas, check_mesh, &
!   cut_domain_from_sphere, solid_boundary_in_y, solid_boundary_in_x
!!   Lx

! Default namelist values.  Default set for realistic global IC.
on_a_sphere = 'NO'
sphere_radius = 0.0
expand_from_unit_sphere = .false.

eliminate_inland_seas=.false.
solid_boundary_in_y = .false.
solid_boundary_in_x = .false.

! This needs to be changed for correct periodic boundaries
! Lx is the TOTAL domain width, and needs to be exact for correct periodic
! boundaries in x.
!Lx = 3200.0e3  ! 40x80km=3200km


!! Read in namelist
!   open(20,file='namelist.landicecrop',status='old')
!   read(20,landicecrop)
!   close(20)


! get to work
write(6,*) ' starting'
write(6,*)

! get grid
write(6,*) ' calling read_grid'
write(6,*)
call read_grid
write(6,*) ' xCell 1: ',minval(xCell), maxval(xCell)

nVertLevelsMOD = nVertLevels

! copy dimensions
write(6,*) ' copy dimensions'
write(6,*)
call copy_dimensions
write(6,*) ' xCell 1: ',minval(xCell), maxval(xCell)

! define the kmt array
write(6,*) ' calling define_kmt'
write(6,*)
call define_kmt

! define the mapping between original and new cells, edges and vertices
write(6,*) ' calling define_mapping'
write(6,*)
call define_mapping

! copy the vector arrays form the original to new arrays
write(6,*) ' calling map_vectors'
write(6,*)
call map_vectors

! define the new connectivity variables
write(6,*) ' calling map_connectivity'
write(6,*)
call map_connectivity

! check the mesh
if (check_mesh) then
   call error_checking
endif

! dump new grid to netCDF
write(6,*) ' calling write_grid'
write(6,*)
call write_grid

! dump graph for partioning
write(6,*) ' call write_graph'
write(6,*)
call write_graph

!! write OpenDx file
!if (write_OpenDX_flag) then
!   write(6,*) ' calling write_OpenDX'
!   write(6,*)
!   call write_OpenDX(        on_a_sphere, &
!                             nCellsNew, &
!                             nVerticesNew, &
!                             nEdgesNew, &
!                             vertexDegreeNew, &
!                             maxEdgesNew, &
!                             xCellNew, &
!                             yCellNew, &
!                             zCellNew, &
!                             xVertexNew, &
!                             yVertexNew, &
!                             zVertexNew, &
!                             xEdgeNew, &
!                             yEdgeNew, &
!                             zEdgeNew, &
!                             nEdgesOnCellNew, &
!                             verticesOnCellNew, &
!                             verticesOnEdgeNew, &
!                             cellsOnVertexNew, &
!                             edgesOnCellNew, &
!                             areaCellNew, &
!                             maxLevelCellNew, &
!                             meshDensityNew, &
!                             bottomDepthNew, &
!                             temperatureNew(1,1,:), &
!                             kiteAreasOnVertexNew )
!endif

!!do iCell=1,nCellsNew
!  !ulon = 1.0; ulat = 0.0
!  !xin = xCellNew(iCell); yin = yCellNew(iCell); zin = zCellNew(iCell)
!  !call transform_from_lonlat_to_xyz(xin, yin, zin, ulon, ulat, ux, uy, uz)
!  !if(abs(ux).lt.1.0e-10) ux=0.0
!  !if(abs(uy).lt.1.0e-10) uy=0.0
!  !if(abs(uz).lt.1.0e-10) uz=0.0
!  !write(20,10) ux, uy, uz
!  !10 format(3e25.10)
!!enddo
  
write(6,*) ' finished'

contains

subroutine write_graph
implicit none
integer :: m,itmp(maxEdgesNew),k

      m=nEdgesNew
      do i=1,nCellsNew
      do j=1,nEdgesOnCellNew(i)
         if(cellsOnCellNew(j,i).eq.0) m=m-1
      enddo
      enddo

      open(42,file='graph.info',form='formatted')
      write(42,*) nCellsNew, m
      do i=1,nCellsNew
         itmp = 0; k = 0;
         do j=1,nEdgesOnCellNew(i)
            if(cellsOnCellNew(j,i).gt.0) then
              k=k+1; itmp(k)=cellsOnCellNew(j,i)
            endif
         enddo
         write(42,'(1x,12i8)',advance='no') (itmp(m),m=1,k)
         write(42,'(1x)')
      end do
      close(42)
end subroutine write_graph


subroutine error_checking
real :: p(3), q(3), r(3), angle, s(3), t(3), dot, mindot, maxdot, b(vertexDegree)
real :: work(nCellsNew)


! write
write(6,*)
write(6,*) ' error checking '
write(6,*)

! check to see if every edge is normal to associated cells
mindot =  2
maxdot = -2
do iEdge=1,nEdgesNew
  if(boundaryEdgeNew(1,iEdge).eq.1) cycle
  iCell1 = cellsOnEdgeNew(1,iEdge)
  iCell2 = cellsOnEdgeNew(2,iEdge)
  p(1)=xCellNew(iCell1); p(2)=yCellNew(iCell1); p(3)=zCellNew(iCell1)
  q(1)=xCellNew(iCell2); q(2)=yCellNew(iCell2); q(3)=zCellNew(iCell2)
  r(1)=xEdgeNew(iEdge);  r(2)=yEdgeNew(iEdge);  r(3)=zEdgeNew(iEdge)
  call unit_vector_in_3space(p)
  call unit_vector_in_3space(q)
  call unit_vector_in_3space(r)
  t = q - p
  s = r - p
  call unit_vector_in_3space(t)
  call unit_vector_in_3space(s)
  dot = s(1)*t(1)+s(2)*t(2)+s(3)*t(3)
  if(dot.lt.mindot) mindot=dot
  if(dot.gt.maxdot) maxdot=dot
enddo
write(6,10) 'alignment of edges and cells (should be ones)', mindot, maxdot
10 format(a60,5x,2e15.5)

! check to see if every segments connecting cells and vertices are orothogonal'
mindot =  2
maxdot = -2
do iEdge=1,nEdgesNew
  if(boundaryEdgeNew(1,iEdge).eq.1) cycle
  iCell1 = cellsOnEdgeNew(1,iEdge)
  iCell2 = cellsOnEdgeNew(2,iEdge)
  iVertex1 = verticesOnEdgeNew(1,iEdge)
  iVertex2 = verticesOnEdgeNew(2,iEdge)
  p(1)=xCellNew(iCell1); p(2)=yCellNew(iCell1); p(3)=zCellNew(iCell1)
  q(1)=xCellNew(iCell2); q(2)=yCellNew(iCell2); q(3)=zCellNew(iCell2)
  r(1)=xVertexNew(iVertex1); r(2)=yVertexNew(iVertex1); r(3)=zVertexNew(iVertex1)
  s(1)=xVertexNew(iVertex2); s(2)=yVertexNew(iVertex2); s(3)=zVertexNew(iVertex2)
  call unit_vector_in_3space(p)
  call unit_vector_in_3space(q)
  call unit_vector_in_3space(r)
  call unit_vector_in_3space(s)
  t = q - p
  s = s - r
  call unit_vector_in_3space(t)
  call unit_vector_in_3space(s)
  dot = s(1)*t(1)+s(2)*t(2)+s(3)*t(3)
  if(dot.lt.mindot) mindot=dot
  if(dot.gt.maxdot) maxdot=dot
enddo
write(6,10) 'orthogonality of cell and vertex edges (should be zeros)', mindot, maxdot

! check that the kiteareas sum to the areatriangle
mindot =  2
maxdot = -2
do iVertex=1,nVerticesNew
  b = 0
  do i=1,vertexDegree
    b(i) = kiteAreasOnVertexNew(i,iVertex)
  enddo
  angle = sum(b)
  if(angle - areaTriangleNew(iVertex).lt.mindot) mindot = angle - areaTriangleNew(iVertex)
  if(angle - areaTriangleNew(iVertex).gt.maxdot) maxdot = angle - areaTriangleNew(iVertex)
enddo
write(6,10) ' error in sum of kites and triangles (should be zeroes)', mindot, maxdot

! check that the kiteareas sum to the areaCell
mindot =  2
maxdot = -2
work = 0
do iVertex=1,nVerticesNew
  iCell1 = cellsOnVertexNew(1,iVertex)
  iCell2 = cellsOnVertexNew(2,iVertex)
  iCell3 = cellsOnVertexNew(3,iVertex)
  if(iCell1.ne.0) work(iCell1) = work(iCell1) + kiteAreasOnVertexNew(1,iVertex)
  if(iCell2.ne.0) work(iCell2) = work(iCell2) + kiteAreasOnVertexNew(2,iVertex)
  if(iCell3.ne.0) work(iCell3) = work(iCell3) + kiteAreasOnVertexNew(3,iVertex)
enddo
mindot = minval(areaCellNew - work)
maxdot = maxval(areaCellNew - work)
write(6,10) ' error in sum of kites and cells (should be zeroes)', mindot, maxdot

!check for connectivity inverses for cells/edges
do iCell=1,nCellsNew
  do i=1,nEdgesOnCellNew(iCell)
    iEdge=edgesOnCellNew(i,iCell)
    if(iEdge.le.0) stop ' iEdge le 0'
    iCell1 = cellsOnEdgeNew(1,iEdge)
    iCell2 = cellsOnEdgeNew(2,iEdge)
    if(iCell1.ne.iCell.and.iCell2.ne.iCell) stop ' cells/edges inverse failed'
  enddo
enddo
write(6,*) ' cellsOnEdge and edgesOnCell are duals for every cell/edge combination'

!check for connectivity inverses for cells/vertices
do iCell=1,nCellsNew
  do i=1,nEdgesOnCellNew(iCell)
    iVertex = verticesOnCellNew(i,iCell)
    if(iVertex.le.0) stop ' iVertex le 0'
    iCell1 = cellsOnVertexNew(1,iVertex)
    iCell2 = cellsOnVertexNew(2,iVertex)
    iCell3 = cellsOnVertexNew(3,iVertex)
    if(iCell1.ne.iCell.and.iCell2.ne.iCell.and.iCell3.ne.iCell) stop ' cells/vertices inverse failed'
  enddo
enddo
write(6,*) ' cellsOnVertex and verticesOnCell are duals for every cell/vertex combination'

!check edgesOnEdge
do iEdge=1,nEdgesNew
  iCell1 = cellsOnEdgeNew(1,iEdge)
  iCell2 = cellsOnEdgeNew(2,iEdge)
  if(nEdgesOnEdgeNew(iEdge).eq.0) then
    if(boundaryEdgeNew(1,iEdge).ne.1) stop ' stopping boundaryEdgeNew'
  endif
  do i=1,nEdgesOnEdgeNew(iEdge)
    jEdge = edgesOnEdgeNew(i,iEdge)
    jCell1 = cellsOnEdgeNew(1,jEdge)
    jCell2 = cellsOnEdgeNew(2,jEdge)
    if(jCell1.ne.iCell1.and.jCell1.ne.iCell2) then
    if(jCell2.ne.iCell1.and.jCell2.ne.iCell2) then
          write(6,*) 'error in edgesOnEdge'
          write(6,*) iCell1, iCell2, jCell1, jCell2
          stop
    endif
    endif
  enddo
enddo
write(6,*) ' edgesOnEdge is consistent with cellsOnEdge'

end subroutine error_checking


subroutine copy_dimensions

maxEdgesNew = maxEdges
maxEdges2New = maxEdges2
TWONew = TWO
vertexDegreeNew = vertexDegree
nVertLevelsNew = nVertLevelsMod

write(6,*)
write(6,*) ' new dimensions '
write(6,*) ' maxEdgesNew     : ', maxEdgesNew
write(6,*) ' maxEdges2New    : ', maxEdges2New
write(6,*) ' TWONew          : ', TWONew
write(6,*) ' vertexDegreeNew : ', vertexDegreeNew
write(6,*) ' nVertLevelsNew  : ', nVertLevelsNew

end subroutine copy_dimensions



subroutine read_grid
implicit none

call read_netcdf_init(nCells, nEdges, nVertices, maxEdges,maxEdges2,&
                       nVertLevels,TWO,vertexDegree)

write(6,*) ' init from grid '
write(6,*) 'nCells        :', nCells
write(6,*) 'nEdges        :', nEdges
write(6,*) 'nVertices     :', nVertices
write(6,*) 'maxEdges      :', maxEdges
write(6,*) 'maxEdges2     :', maxEdges2
write(6,*) 'nVertLevels   :', nVertLevels
write(6,*) 'vertexDegree  :', vertexDegree
write(6,*) 'TWO           :', TWO

allocate(xCell(nCells))
allocate(yCell(nCells))
allocate(zCell(nCells))
allocate(latCell(nCells))
allocate(lonCell(nCells))
allocate(meshDensity(nCells))
allocate(xEdge(nEdges))
allocate(yEdge(nEdges))
allocate(zEdge(nEdges))
allocate(latEdge(nEdges))
allocate(lonEdge(nEdges))
allocate(xVertex(nVertices))
allocate(yVertex(nVertices))
allocate(zVertex(nVertices))
allocate(latVertex(nVertices))
allocate(lonVertex(nVertices))
allocate(dcEdge(nEdges))
allocate(dvEdge(nEdges))

allocate(indexToCellID(nCells))
allocate(indexToEdgeID(nEdges))
allocate(indexToVertexID(nVertices))

allocate(cellsOnEdge(TWO,nEdges))
allocate(nEdgesOnCell(nCells))
allocate(nEdgesOnEdge(nEdges))
allocate(edgesOnCell(maxEdges,nCells))
allocate(edgesOnEdge(maxEdges2,nEdges))
allocate(weightsOnEdge(maxEdges2,nEdges))

allocate(angleEdge(nEdges))
allocate(areaCell(nCells))
allocate(areaTriangle(nVertices))
allocate(cellsOnCell(maxEdges,nCells))
allocate(verticesOnCell(maxEdges,nCells))
allocate(verticesOnEdge(TWO,nEdges))
allocate(edgesOnVertex(vertexDegree,nVertices))
allocate(cellsOnVertex(vertexDegree,nVertices))
allocate(kiteAreasOnVertex(vertexDegree,nVertices))

allocate(keepCellMask(nCells))


xCell=0; yCell=0; zCell=0; latCell=0; lonCell=0; meshDensity=1.0
xEdge=0; yEdge=0; zEdge=0; latEdge=0; lonEdge=0
xVertex=0; yVertex=0; zVertex=0; latVertex=0; lonVertex=0

indexToCellID=0; indexToEdgeID=0; indexToVertexID=0
cellsOnEdge=0; nEdgesOnCell=0; edgesOnCell=0
edgesOnEdge=0; weightsOnEdge=0
angleEdge=0; areaCell=0; areaTriangle=0
cellsOnCell=0; verticesOnCell=0; verticesOnEdge=0
edgesOnVertex=0; cellsOnVertex=0; kiteAreasOnVertex=0


call  read_netcdf_fields( &
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
                    keepCellMask &
                   )

write(6,*) ' values from read grid, min/max'
write(6,*) ' latCell : ', minval(latCell), maxval(latCell)
write(6,*) ' lonCell : ', minval(lonCell), maxval(lonCell)
write(6,*) ' meshDensity : ', minval(meshDensity),maxval(meshDensity)
write(6,*) ' xCell : ', minval(xCell), maxval(xCell)
write(6,*) ' yCell : ', minval(yCell), maxval(yCell)
write(6,*) ' zCell : ', minval(zCell), maxval(zCell)
write(6,*) ' indexToCellID : ', minval(indexToCellID), maxval(indexToCellID)
write(6,*) ' latEdge : ', minval(latEdge), maxval(latEdge)
write(6,*) ' lonEdge : ', minval(lonEdge), maxval(lonEdge)
write(6,*) ' xEdge : ', minval(xEdge), maxval(xEdge)
write(6,*) ' yEdge : ', minval(yEdge), maxval(yEdge)
write(6,*) ' zEdge : ', minval(zEdge), maxval(zEdge)
write(6,*) ' indexToEdgeID : ', minval(indexToEdgeID), maxval(indexToEdgeID)
write(6,*) ' latVertex : ', minval(latVertex), maxval(latVertex)
write(6,*) ' lonVertex : ', minval(lonVertex), maxval(lonVertex)
write(6,*) ' xVertex : ', minval(xVertex), maxval(xVertex)
write(6,*) ' yVertex : ', minval(yVertex), maxval(yVertex)
write(6,*) ' zVertex : ', minval(zVertex), maxval(zVertex)
write(6,*) ' indexToVertexID : ', minval(indexToVertexID), maxval(indexToVertexID)
write(6,*) ' cellsOnEdge : ', minval(cellsOnEdge), maxval(cellsOnEdge)
write(6,*) ' nEdgesOnCell : ', minval(nEdgesOnCell), maxval(nEdgesOnCell)
write(6,*) ' nEdgesOnEdge : ', minval(nEdgesOnEdge), maxval(nEdgesOnEdge)
write(6,*) ' edgesOnCell : ', minval(edgesOnCell), maxval(edgesOnCell)
write(6,*) ' edgesOnEdge : ', minval(edgesOnEdge), maxval(edgesOnEdge)
write(6,*) ' weightsOnEdge : ', minval(weightsOnEdge), maxval(weightsOnEdge)
write(6,*) ' dvEdge : ', minval(dvEdge), maxval(dvEdge)
write(6,*) ' dcEdge : ', minval(dcEdge), maxval(dcEdge)
write(6,*) ' angleEdge : ', minval(angleEdge), maxval(angleEdge)
write(6,*) ' areaCell : ', minval(areaCell), maxval(areaCell)
write(6,*) ' areaTriangle : ', minval(areaTriangle), maxval(areaTriangle)
write(6,*) ' cellsOnCell : ', minval(cellsOnCell), maxval(cellsOnCell)
write(6,*) ' verticesOnCell : ', minval(verticesOnCell), maxval(verticesOnCell)
write(6,*) ' verticesOnEdge : ', minval(verticesOnEdge), maxval(verticesOnEdge)
write(6,*) ' edgesOnVertex : ', minval(edgesOnVertex), maxval(edgesOnVertex)
write(6,*) ' cellsOnVertex : ', minval(cellsOnVertex), maxval(cellsOnVertex)
write(6,*) ' kiteAreasOnVertex : ', minval(kiteAreasOnVertex), maxval(kiteAreasOnVertex)
write(6,*) ' keepCellMask : sum (should be number of 1 values) ', sum(keepCellMask)


end subroutine read_grid


subroutine write_grid
implicit none

if (expand_from_unit_sphere) then
   xCellNew = xCellNew * sphere_radius
   yCellNew = yCellNew * sphere_radius
   zCellNew = zCellNew * sphere_radius
   xEdgeNew = xEdgeNew * sphere_radius
   yEdgeNew = yEdgeNew * sphere_radius
   zEdgeNew = zEdgeNew * sphere_radius
   xVertexNew = xVertexNew * sphere_radius
   yVertexNew = yVertexNew * sphere_radius
   zVertexNew = zVertexNew * sphere_radius
   dcEdgeNew = dcEdgeNew * sphere_radius
   dvEdgeNew = dvEdgeNew * sphere_radius
   areaCellNew = areaCellNew * (sphere_radius)**2
   areaTriangleNew = areaTriangleNew * (sphere_radius)**2
   kiteAreasOnVertexNew = kiteAreasOnVertexNew * (sphere_radius)**2
endif

call write_netcdf_init( &
                nCellsNew, &
                nEdgesNew, &
                nVerticesNew, &
                maxEdgesNew, &
                nVertLevelsNew, &
                vertexDegreeNew, &
                sphere_radius, &
                on_a_sphere &
                )

call write_netcdf_fields( &
                    1, &
                    latCellNew, &
                    lonCellNew, &
                    meshDensityNew, &
                    xCellNew, &
                    yCellNew, &
                    zCellNew, &
                    indexToCellIDNew, &
                    latEdgeNew, &
                    lonEdgeNew, &
                    xEdgeNew, &
                    yEdgeNew, &
                    zEdgeNew, &
                    indexToEdgeIDNew, &
                    latVertexNew, &
                    lonVertexNew, &
                    xVertexNew, &
                    yVertexNew, &
                    zVertexNew, &
                    indexToVertexIDNew, &
                    cellsOnEdgeNew, &
                    nEdgesOnCellNew, &
                    nEdgesOnEdgeNew, &
                    edgesOnCellNew, &
                    edgesOnEdgeNew, &
                    weightsOnEdgeNew, &
                    dvEdgeNew, &
                    dcEdgeNew, &
                    angleEdgeNew, &
                    areaCellNew, &
                    areaTriangleNew, &
                    cellsOnCellNew, &
                    verticesOnCellNew, &
                    verticesOnEdgeNew, &
                    edgesOnVertexNew, &
                    cellsOnVertexNew, &
                    kiteAreasOnVertexNew &
                   )

call write_netcdf_finalize

if (expand_from_unit_sphere) then
   xCellNew = xCellNew / sphere_radius
   yCellNew = yCellNew / sphere_radius
   zCellNew = zCellNew / sphere_radius
   xEdgeNew = xEdgeNew / sphere_radius
   yEdgeNew = yEdgeNew / sphere_radius
   zEdgeNew = zEdgeNew / sphere_radius
   xVertexNew = xVertexNew / sphere_radius
   yVertexNew = yVertexNew / sphere_radius
   zVertexNew = zVertexNew / sphere_radius
   dcEdgeNew = dcEdgeNew / sphere_radius
   dvEdgeNew = dvEdgeNew / sphere_radius
   areaCellNew = areaCellNew / (sphere_radius)**2
   areaTriangleNew = areaTriangleNew / (sphere_radius)**2
   kiteAreasOnVertexNew = kiteAreasOnVertexNew / (sphere_radius)**2
endif

end subroutine write_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Step 4: Check define_kmt routine for bottomDepth and kmt (maxLevelCell) variables
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine define_kmt
implicit none
real (kind=4), allocatable, dimension(:) :: x,y, work_kmt
real (kind=4), allocatable, dimension(:,:) :: ztopo
integer :: nx, ny, inx, iny, ix, iy, kmt_neighbor_max
real :: pi, dtr, zdata, rlon, rlat, r, ymin, ymax, xmin, xmax
real :: latmin, latmax, lonmin, lonmax, ridgeDepth, maxdc
logical :: flag, kmt_flag

!pi = 4.0*atan(1.0)
!dtr = pi / 180.0

allocate(kmt(nCells))
kmt = 0


do iCell=1,nCells
   if (keepCellMask(iCell) > 0) then
     kmt(iCell) = 1  ! any positive value will keep the cell, 0 cells will be culled
   endif
enddo


!if (cut_domain_from_sphere) then
!   latmin = -30*dtr
!   latmax = +30*dtr
!   lonmin = +10*dtr
!   lonmax = +70*dtr
!   write(6,*) ' lat min ', latmin
!   write(6,*) ' lat max ', latmax
!   where(latCell.lt.latmin) kmt = 0
!   where(latCell.gt.latmax) kmt = 0
!   where(lonCell.lt.lonmin) kmt = 0
!   where(lonCell.gt.lonmax) kmt = 0
!endif

!if (solid_boundary_in_y) then
!   ymin = minval(yCell)
!   write(6,*) ' minimum yCell ', ymin
!   ymax = maxval(yCell)
!   write(6,*) ' maximum yCell ', ymax
!   where(yCell.lt.1.001*ymin) kmt = 0
!   where(yCell.gt.0.999*ymax) kmt = 0
!endif

!if (solid_boundary_in_x) then
!   maxdc = maxval(dcEdge)
!   xmin = minval(xCell)
!   write(6,*) ' minimum xCell ', xmin
!   xmax = maxval(xCell)
!   write(6,*) ' maximum xCell ', xmax
!   where(xCell.lt.xmin+maxdc/1.5) kmt = 0
!   where(xCell.gt.xmax-maxdc/1.5) kmt = 0
!endif

    
allocate(work_kmt(nCells))
work_kmt = 0.0
where(kmt.eq.0) work_kmt=1.0
write(6,*) 'number of cells culled ',sum(work_kmt)
deallocate(work_kmt)


!! Eliminate isolated ocean cells, and make these isolated deep cells
!! flush with the deepest neighbor.
!do iCell=1,nCells
!   kmt_neighbor_max = 0
!   do j=1,nEdgesOnCell(iCell)
!      iCell1 = cellsOnCell(j,iCell)
!      kmt_neighbor_max = max(kmt_neighbor_max,kmt(iCell1))
!   enddo
!   if (kmt(iCell).gt.kmt_neighbor_max) then
!      kmt(iCell) = kmt_neighbor_max
!      bottomDepth(iCell) = refBottomDepth(kmt(iCell))
!   endif
!enddo

!if(eliminate_inland_seas) then
!call eliminateLoops(nCells,nEdges,nVertices,maxEdges,vertexDegree, &
!                    nEdgesOnCell, cellsOnCell, verticesOnEdge, cellsOnVertex, edgesOnCell, lonCell, latCell, &
!                    xCell, yCell, zCell, xEdge, yEdge, zEdge, xVertex, yVertex, zVertex, &
!                    KMT)
!endif

!! do not allow land or PBCs in top layers
!k = min(top_layers_without_land,nVertLevelsMOD)
!where(kmt.gt.0.and.kmt.le.k) 
!   bottomDepth = refBottomDepth(k)
!   kmt=k
!endwhere

end subroutine define_kmt



subroutine define_mapping
implicit none

allocate(cellMap(nCells))
allocate(edgeMap(nEdges))
allocate(vertexMap(nVertices))
cellMap = 0; edgeMap = 0; vertexMap = 0

j=1
do i=1,nCells
if(kmt(i).ne.0) then
    cellMap(i) = j
    j=j+1
endif
write(10,*) i, cellMap(i)
enddo

j=1
do i=1,nEdges
iCell1 = cellsOnEdge(1,i)
iCell2 = cellsOnEdge(2,i)
if(kmt(iCell1).ne.0.or.kmt(iCell2).ne.0) then
    edgeMap(i)=j
    j=j+1
endif
write(11,*) i,edgeMap(i)
enddo

j=1
do i=1,nVertices
iCell1 = cellsOnVertex(1,i)
iCell2 = cellsOnVertex(2,i)
iCell3 = cellsOnVertex(3,i)
if(kmt(iCell1).ne.0.or.kmt(iCell2).ne.0.or.kmt(iCell3).ne.0) then
    vertexMap(i)=j
    j=j+1
endif
write(12,*) i,vertexMap(i)
enddo

nCellsNew = 0
do i=1,nCells
if(cellMap(i).ne.0) nCellsNew = nCellsNew + 1
enddo

nEdgesNew = 0
do i=1,nEdges
if(edgeMap(i).ne.0) nEdgesNew = nEdgesNew + 1
enddo

nVerticesNew = 0
do i=1,nVertices
if(vertexMap(i).ne.0) nVerticesNew = nVerticesNew + 1
enddo

write(6,*) ' mesh mapping found '
write(6,*)  nCells, nCellsNew
write(6,*)  nEdges, nEdgesNew
write(6,*)  nVertices, nVerticesNew

allocate(indexToCellIDNew(nCellsNew))
allocate(indexToEdgeIDNew(nEdgesNew))
allocate(indexToVertexIDNew(nVerticesNew))
indextoCellIDNew = 0; indexToEdgeIDNew = 0; indexToVertexIDNew = 0

do i=1,nCellsNew
indexToCellIDNew(i)=i
enddo

do i=1,nEdgesNew
indexToEdgeIDNew(i)=i
enddo

do i=1,nVerticesNew
indexToVertexIDNew(i)=i
enddo

end subroutine define_mapping


subroutine map_vectors
implicit none

allocate(xCellNew(nCellsNew))
allocate(yCellNew(nCellsNew))
allocate(zCellNew(nCellsNew))
allocate(normalsNew(3,nEdgesNew))
allocate(latCellNew(nCellsNew))
allocate(lonCellNew(nCellsNew))
allocate(meshDensityNew(nCellsNew))
allocate(meshSpacingNew(nCellsNew))
allocate(xEdgeNew(nEdgesNew))
allocate(yEdgeNew(nEdgesNew))
allocate(zEdgeNew(nEdgesNew))
allocate(latEdgeNew(nEdgesNew))
allocate(lonEdgeNew(nEdgesNew))
allocate(xVertexNew(nVerticesNew))
allocate(yVertexNew(nVerticesNew))
allocate(zVertexNew(nVerticesNew))
allocate(latVertexNew(nVerticesNew))
allocate(lonVertexNew(nVerticesNew))
allocate(dcEdgeNew(nEdgesNew))
allocate(dvEdgeNew(nEdgesNew))
allocate(angleEdgeNew(nEdgesNew))
allocate(areaCellNew(nCellsNew))
allocate(areaTriangleNew(nVerticesNew))




xCellNew=0; yCellNew=0; zCellNew=0; latCellNew=0; lonCellNew=0; meshDensityNew=1.0; meshSpacingNew=0.0
xEdgeNew=0; yEdgeNew=0; zEdgeNew=0; latEdgeNew=0; lonEdgeNew=0
xVertexNew=0; yVertexNew=0; zVertexNew=0; latVertexNew=0; lonVertexNew=0


do i=1,nCells
jNew = cellMap(i)
if(jNew.ne.0) then
    xCellNew(jNew)=xCell(i)
    yCellNew(jNew)=yCell(i)
    zCellNew(jNew)=zCell(i)
    latCellNew(jNew)=latCell(i)
    lonCellNew(jNew)=lonCell(i)
    meshDensityNew(jNew)=meshDensity(i)
    areaCellNew(jNew)=areaCell(i)
endif
enddo

do i=1,nEdges
jNew = edgeMap(i)
if(jNew.ne.0) then
    xEdgeNew(jNew)=xEdge(i)
    yEdgeNew(jNew)=yEdge(i)
    zEdgeNew(jNew)=zEdge(i)
    latEdgeNew(jNew)=latEdge(i)
    lonEdgeNew(jNew)=lonEdge(i)
    dcEdgeNew(jNew) = dcEdge(i)
    dvEdgeNew(jNew) = dvEdge(i)
    angleEdgeNew(jNew) = angleEdge(i)
endif
enddo

do i=1,nVertices
jNew = vertexMap(i)
if(jNew.ne.0) then
    xVertexNew(jNew)=xVertex(i)
    yVertexNew(jNew)=yVertex(i)
    zVertexNew(jNew)=zVertex(i)
    latVertexNew(jNew)=latVertex(i)
    lonVertexNew(jNew)=lonVertex(i)
    areaTriangleNew(jNew)=areaTriangle(i)
endif
enddo

deallocate(xCell)
deallocate(yCell)
deallocate(zCell)
deallocate(latCell)
deallocate(lonCell)
deallocate(meshDensity)
deallocate(xEdge)
deallocate(yEdge)
deallocate(zEdge)
deallocate(latEdge)
deallocate(lonEdge)
deallocate(xVertex)
deallocate(yVertex)
deallocate(zVertex)
deallocate(latVertex)
deallocate(lonVertex)
deallocate(dcEdge)
deallocate(dvEdge)

end subroutine map_vectors



subroutine map_connectivity
implicit none

allocate(cellsOnEdgeNew(TWONew,nEdgesNew))
allocate(boundaryEdgeNew(nVertLevelsNew,nEdgesNew))
allocate(flipVerticesOnEdgeOrdering(nEdgesNew))
cellsOnEdgeNew(:,:) = 0
boundaryEdgeNew(:,:) = 0
flipVerticesOnEdgeOrdering(:) = 0
do iEdge=1,nEdges
if(edgeMap(iEdge).eq.0) cycle
iEdgeNew = edgeMap(iEdge)
iCell1 = cellsOnEdge(1,iEdge)
iCell2 = cellsOnEdge(2,iEdge)
iCell1New = cellMap(iCell1)
iCell2New = cellMap(iCell2)
cellsOnEdgeNew(1,iEdgeNew) = iCell1New
cellsOnEdgeNew(2,iEdgeNew) = iCell2New
if(iCell1New.eq.0.or.iCell2New.eq.0) boundaryEdgeNew(:,iEdgeNew) = 1
if(iCell1New.eq.0.and.iCell2New.eq.0) stop "cellsOnEdge"
if(iCell1New.eq.0) then
    cellsOnEdgeNew(1,iEdgeNew) = iCell2New
    cellsOnEdgeNew(2,iEdgeNew) = iCell1New
    flipVerticesOnEdgeOrdering(iEdgeNew) = 1
endif
enddo
deallocate(cellsOnEdge)

allocate(verticesOnEdgeNew(TWONew,nEdgesNew))
allocate(boundaryVertexNew(nVertLevelsNew,nVerticesNew))
verticesOnEdgeNew(:,:) = 0
boundaryVertexNew(:,:) = 0
do iEdge=1,nEdges
if(edgeMap(iEdge).eq.0) cycle
iEdgeNew = edgeMap(iEdge)
iVertex1 = VerticesOnEdge(1,iEdge)
iVertex2 = VerticesOnEdge(2,iEdge)
iVertex1New = vertexMap(iVertex1)
iVertex2New = vertexMap(iVertex2)
if(iVertex1New.eq.0.or.iVertex2New.eq.0) stop "verticesOnEdge"
if(flipVerticesOnEdgeOrdering(iEdgeNew).eq.0) then
  verticesOnEdgeNew(1,iEdgeNew) = iVertex1New
  verticesOnEdgeNew(2,iEdgeNew) = iVertex2New
else
  verticesOnEdgeNew(1,iEdgeNew) = iVertex2New
  verticesOnEdgeNew(2,iEdgeNew) = iVertex1New
endif
if(boundaryEdgeNew(1,iEdgeNew).eq.1) then
    boundaryVertexNew(:,iVertex1New)=1
    boundaryVertexNew(:,iVertex2New)=1
endif
enddo
deallocate(verticesOnEdge)

allocate(nEdgesOnEdgeNew(nEdgesNew))
allocate(edgesOnEdgeNew(maxEdges2,nEdgesNew))
allocate(weightsOnEdgeNew(maxEdges2,nEdgesNew))
nEdgesOnEdgeNew(:) = 0
edgesOnEdgeNew(:,:) = 0
weightsOnEdgeNew(:,:) = 0.0
do iEdge=1,nEdges
if(edgeMap(iEdge).eq.0) cycle
iEdgeNew = edgeMap(iEdge)
if(boundaryEdgeNew(1,iEdgeNew).eq.1) then
    nEdgesOnEdgeNew(iEdgeNew) = 0
    edgesOnEdgeNew(:,iEdgeNew) = 0
    weightsOnEdgeNew(:,iEdgeNew) = 0.0
else
    nEdgesOnEdgeNew(iEdgeNew) = nEdgesOnEdge(iEdge)
    do i=1,nEdgesOnEdgeNew(iEdgeNew)
    jEdge = edgesOnEdge(i,iEdge)
    jEdgeNew = edgeMap(jEdge)
    if(jEdgeNew.eq.0) stop "jEdgeNew"
    edgesOnEdgeNew(i,iEdgeNew)=jEdgeNew
    weightsOnEdgeNew(i,iEdgeNew) = weightsOnEdge(i,iEdge)
    enddo
endif
enddo
deallocate(nEdgesOnEdge)
deallocate(edgesOnEdge)
deallocate(weightsOnEdge)

allocate(cellsOnCellNew(maxEdges,nCellsNew))
allocate(nEdgesOnCellNew(nCellsNew))
cellsOnCellNew = 0
nEdgesOnCellNew = 0
do iCell=1,nCells
if(cellMap(iCell).eq.0) cycle
iCellNew = cellMap(iCell)
nEdgesOnCellNew(iCellNew)=nEdgesOnCell(iCell)
do i=1,nEdgesOnCellNew(iCellNew)
j = cellsOnCell(i,iCell)
jNew = cellMap(j)
cellsOnCellNew(i,iCellNew) = jNew
enddo
enddo
deallocate(cellsOnCell)
deallocate(nEdgesOnCell)

allocate(edgesOnCellNew(maxEdgesNew,nCellsNew))
edgesOnCellNew(:,:) = 0
meshSpacingNew(:) = 0.0
do iCell=1,nCells
if(cellMap(iCell).eq.0) cycle
iCellNew = cellMap(iCell)
do i=1,nEdgesOnCellNew(iCellNew)
j = edgesOnCell(i,iCell)
jNew = edgeMap(j)
if(jNew.eq.0) stop "edgesOnCell"
edgesOnCellNew(i,iCellNew) = jNew
meshSpacingNew(iCellNew) = meshSpacingNew(iCellNew) + dcEdgeNew(jNew)/nEdgesOnCellNew(iCellNew)
enddo
enddo
deallocate(edgesOnCell)

allocate(verticesOnCellNew(maxEdgesNew,nCellsNew))
verticesOnCellNew(:,:)=0
do iCell=1,nCells
if(cellMap(iCell).eq.0) cycle
iCellNew = cellMap(iCell)
do i=1,nEdgesOnCellNew(iCellNew)
j=verticesOnCell(i,iCell)
jNew = vertexMap(j)
if(jNew.eq.0) stop "verticesOnCell"
verticesOnCellNew(i,iCellNew) = jNew
enddo
enddo
deallocate(verticesOnCell)

allocate(cellsOnVertexNew(vertexDegreeNew,nVerticesNew))
allocate(kiteAreasOnVertexNew(vertexDegreeNew,nVerticesNew))
cellsOnVertexNew = 0
kiteAreasOnVertexNew = 0
do iVertex=1,nVertices
if(vertexMap(iVertex).eq.0) cycle
iVertexNew = vertexMap(iVertex)
do i=1,vertexDegree
j=cellsOnVertex(i,iVertex)
jNew=cellMap(j)
if(jNew.eq.0) then
    kiteAreasOnVertexNew(i,iVertexNew)=0
else
    kiteAreasOnVertexNew(i,iVertexNew)=kiteAreasOnVertex(i,iVertex)
endif
cellsOnVertexNew(i,iVertexNew)=jNew
enddo
enddo
deallocate(cellsOnVertex)
deallocate(kiteAreasOnVertex)

areaTriangleNew = 0
do iVertex=1,nVerticesNew
do i=1,vertexDegree
areaTriangleNew(iVertex) = areaTriangleNew(iVertex) + kiteAreasOnVertexNew(i,iVertex)
enddo
enddo

allocate(edgesOnVertexNew(vertexDegreeNew, nVerticesNew))
edgesOnVertexNew = 0
do iVertex=1,nVertices
if(vertexMap(iVertex).eq.0) cycle
iVertexNew = vertexMap(iVertex)
do i=1,vertexDegree
j=edgesOnVertex(i,iVertex)
jNew=edgeMap(j)
edgesOnVertexNew(i,iVertexNew)=jNew
enddo
enddo
deallocate(edgesOnVertex)

! find normals
normalsNew = 0.0
do iEdge=1,nEdgesNew
cell1 = cellsOnEdgeNew(1,iEdge)
cell2 = cellsOnEdgeNew(2,iEdge)
if(cell1.eq.0.or.cell2.eq.0) cycle
c1(1) = xCellNew(cell1); c1(2) = yCellNew(cell1); c1(3) = zCellNew(cell1)
c2(1) = xCellNew(cell2); c2(2) = yCellNew(cell2); c2(3) = zCellNew(cell2)
distance = sqrt( (c1(1)-c2(1))**2 + (c1(2)-c2(2))**2 + (c1(3)-c2(3))**2 )

if(on_a_sphere.eq.'YES') then
    normalsNew(1,iEdge) = c2(1) - c1(1)
    normalsNew(2,iEdge) = c2(2) - c1(2)
    normalsNew(3,iEdge) = c2(3) - c1(3)
    distance = sqrt( (c1(1)-c2(1))**2 + (c1(2)-c2(2))**2 + (c1(3)-c2(3))**2 )
    normalsNew(:,iEdge) = normalsNew(:,iEdge) / distance
else
!    if(distance.gt.0.5*Lx) then
!        write(6,*) ' periodic edge ', iEdge, distance
!        write(6,10) '          c1   ', c1(:)
!        write(6,10) '          c2   ', c2(:)
!        r = c2(1) - c1(1)
!        if(r.gt.0) c2(1) = c2(1) - Lx
!        if(r.lt.0) c2(1) = c2(1) + Lx
!        distance = sqrt( (c1(1)-c2(1))**2 + (c1(2)-c2(2))**2 + (c1(3)-c2(3))**2 )
!        write(6,*) ' periodic edge fix ', iEdge, r, distance
!    endif
    normalsNew(1,iEdge) = c2(1) - c1(1)
    normalsNew(2,iEdge) = c2(2) - c1(2)
    normalsNew(3,iEdge) = c2(3) - c1(3)
    distance = sqrt( (c1(1)-c2(1))**2 + (c1(2)-c2(2))**2 + (c1(3)-c2(3))**2 )
    normalsNew(:,iEdge) = normalsNew(:,iEdge) / distance
endif
enddo
10 format(a20,3e15.5)

end subroutine map_connectivity


end program crop_landice_grid
