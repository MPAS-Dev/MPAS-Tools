program map_to_basin

use read_netcdf
use read_topo
use read_TS
use read_MONTHLY
use write_netcdf
use utilities
use cullLoops

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Program: basin.F
!
! This program is meant to add land to grids, as well as initial conditions.
! 
! This program is used to take a specific mesh, and remove Cells from it
! It can be used to change a planar grid into a Channel or a basin grid, or to 
! Change a spherical grid into a Limited area spherical grid.
!
! How to use:
! Step 1: Link namelist.basin to the correct namelist file. 
! Step 2: Change parameters and flags in namelist file as needed.
! Step 3: Check get_init_conditions routine for initial T&S, thickness, etc.
! Step 4: Check define_kmt routine for bottomDepth and kmt (maxLevelCell) variables.
! Step 5: Check get_dz routine for hZLevel variable.
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

real, allocatable, dimension(:) :: fEdge, fVertex, bottomDepth, work1
real, allocatable, dimension(:,:) :: normalVelocityForcing
real, allocatable, dimension(:,:,:) :: normalVelocity, tangentialVelocity, layerThickness
real, allocatable, dimension(:,:,:) :: density

integer nlon, nlat, ndepth
real(kind=4), allocatable, dimension(:) :: t_lon, t_lat, depth_t
real(kind=4), allocatable, dimension(:,:) :: mTEMP, mSALT
real(kind=4), allocatable, dimension(:,:,:) :: TEMP, SALT
real(kind=4), allocatable, dimension(:,:) :: TAUX, TAUY

real(kind=4), allocatable, dimension(:,:,:) :: SST_MONTHLY, SSS_MONTHLY
real(kind=4), allocatable, dimension(:,:,:) :: TAUX_MONTHLY, TAUY_MONTHLY

real, dimension(40) :: dz
integer :: nMonths = 1

   real (kind=8) :: ymid, ytmp, ymax, xmid, xloc, yloc, pert, ymin, distance, r, c1(3), c2(3)
   real (kind=8) :: latmid, lattmp, latmax, latmin
   integer :: cell1, cell2
real ::    eos_linear_alpha, eos_linear_beta, eos_linear_Tref, eos_linear_Sref, eos_linear_densityref

! new grid variables
real, allocatable, dimension(:) :: hZLevel, refBottomDepth
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

real, allocatable, dimension(:) :: fEdgeNew, fVertexNew, bottomDepthNew
real, allocatable, dimension(:,:) :: normalVelocityForcingNew
real, allocatable, dimension(:,:) :: windStressMonthlyNew
real, allocatable, dimension(:,:,:) :: normalVelocityNew, layerThicknessNew
real, allocatable, dimension(:,:,:) :: densityNew, temperatureNew, salinityNew, tracer1New
real, allocatable, dimension(:) :: temperatureRestoreNew, salinityRestoreNew
real, allocatable, dimension(:,:) :: temperatureRestoreMonthlyNew, salinityRestoreMonthlyNew

! mapping variables
integer, allocatable, dimension(:) :: kmt, maxLevelCellNew
integer, allocatable, dimension(:) :: cellMap, edgeMap, vertexMap

! work variables
integer :: i,j,jNew,k,jEdge,jEdgeNew,iVertex1New,iVertex2New,iCell1New,iCell2New
integer :: iCell, iCell1, iCell2, iCell3, iEdge, iVertex, iVertex1, iVertex2
integer :: iCellNew, iEdgeNew, iVertexNew, ndata, jCell1, jCell2, jCell, iter
real :: xin, yin, zin, ulon, ulat, ux, uy, uz, rlon, rlat, temp_t, temp_s

integer :: iMonth
character(len=80) :: fileNameT, fileNameS, fileNameU

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Namelist variables
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Variables in namelist file
character (len=32) :: on_a_sphere, zLevel_thickness,bottom_topography, initial_conditions
logical :: expand_from_unit_sphere, eliminate_inland_seas, load_woce_IC, &
   write_OpenDX_flag, monthly_forcing, check_mesh, &
   cut_domain_from_sphere, solid_boundary_in_y, solid_boundary_in_x, &
   amplify_acc_wind

integer :: nVertLevelsMOD, top_layers_without_land 
real (kind=8) :: sphere_radius, layer_thickness_total_max, normalVelocityForcing_max, f0, beta, omega, Lx, &
  uniform_T, uniform_S, uniform_tracer1, &
  amp_wind_factor, amp_wind_center_lat, amp_wind_spread_lat

! specify namelist
namelist /basin/ nVertLevelsMOD, on_a_sphere, sphere_radius, &
   expand_from_unit_sphere, &
   zLevel_thickness, bottom_topography, initial_conditions, &
   eliminate_inland_seas, load_woce_IC, write_OpenDX_flag, monthly_forcing, check_mesh, &
   cut_domain_from_sphere, solid_boundary_in_y, solid_boundary_in_x, &
   top_layers_without_land, layer_thickness_total_max, normalVelocityForcing_max, f0, beta, omega, Lx, &
   uniform_T, uniform_S, uniform_tracer1, &
   amplify_acc_wind, amp_wind_factor, amp_wind_center_lat, amp_wind_spread_lat

! Default namelist values.  Default set for realistic global IC.
nVertLevelsMOD = 40
on_a_sphere = 'YES'
sphere_radius = 6.37122e6
expand_from_unit_sphere = .true.

! zLevel thickness options: 
! 'POP_40_zLevel', 'equally_spaced', 'zero'
zLevel_thickness = 'POP_40_zLevel'

! bottom topography options: 
! 'realistic_ETOPO', 'flat_bottom', 'Ilicak2_overflow', 'DOME_3D_overflow'
bottom_topography = 'realistic_ETOPO'

! initial temperature and salinity options: 
! 'realistic_WOCE', 'uniform_TS', 'Ilicak2_overflow', 'Ilicak2_overflow_sigma', 'DOME_3D_overflow', 
! 'internal_wave_Ilicak3', 'baroclinic_channel_Ilicak4'
initial_conditions = 'realistic_WOCE'
uniform_T = 10.0
uniform_S = 34.0
uniform_tracer1 = 1.0

eliminate_inland_seas=.true.
load_woce_IC = .true.
write_OpenDX_flag = .false.
monthly_forcing = .true.
check_mesh = .true.
cut_domain_from_sphere = .false.
solid_boundary_in_y = .false.
solid_boundary_in_x = .false.

! Set the number of top layers that are not allowed to have land, usually three.
top_layers_without_land = 3

layer_thickness_total_max = 2000.0 ! total layer thickness, for equally spaced case
normalVelocityForcing_max = 0.1 ! max wind stress, N/m2
f0 = -1.1e-4  ! Coriolis parameter
beta = 1.4e-11 
omega = 7.29212e-5  ! rotation rate of earth

! This needs to be changed for correct periodic boundaries
! Lx is the TOTAL domain width, and needs to be exact for correct periodic
! boundaries in x.
Lx = 3200.0e3  ! 40x80km=3200km

! amplify wind stress in acc, for study with Valis
amplify_acc_wind = .false.
amp_wind_factor = 2.0
amp_wind_center_lat = -35.0
amp_wind_spread_lat = 3.0

! Read in namelist
   open(20,file='namelist.basin',status='old')
   read(20,basin)
   close(20)

allocate (hZLevel(nVertLevelsMOD), refBottomDepth(nVertLevelsMOD))


if(monthly_forcing) then
    nMonths = 12
else
    nMonths = 1
end if

fileNameT = 'TS/annual/woce_t_ann.3600x2431x42interp.r4.nc'
fileNameS = 'TS/annual/woce_s_ann.3600x2431x42interp.r4.nc'
fileNameU = 'TS/annual/ws.old_ncep_1958-2000avg.interp3600x2431.nc'

! get to work
write(6,*) ' starting'
write(6,*)

! get depth profile for later
write(6,*) ' calling get_dz'
call get_dz

! get grid
write(6,*) ' calling read_grid'
write(6,*)
call read_grid
write(6,*) ' xCell 1: ',minval(xCell), maxval(xCell)

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

if (load_woce_IC) then
   write(6,*) ' getting woce t and s '

   call read_TS_init(nlon, nlat, ndepth,fileNameT, fileNameS, fileNameU)
   write(6,*) ' TS INIT ', nlon, nlat, ndepth
   allocate(t_lon(nlon), t_lat(nlat), depth_t(ndepth), TEMP(nlon,nlat,ndepth), SALT(nlon,nlat,ndepth))
   allocate(TAUX(nlon,nlat), TAUY(nlon,nlat))
   allocate(mTEMP(nlat,ndepth), mSALT(nlat,ndepth))
   call read_TS_fields(t_lon, t_lat, depth_t, TEMP, SALT, TAUX, TAUY)
   call read_TS_finalize()

   do k=1,ndepth
     ndata = 0; temp_t=0; temp_s=0
     do j=1,nlat
     do i=1,nlon
       if(TEMP(i,j,k).gt.-10.0) then
         ndata = ndata + 1
         temp_t = temp_t + TEMP(i,j,k)
         temp_s = temp_s + SALT(i,j,k)
       endif
     enddo
     enddo
     mTEMP(:,k) = temp_t / float(ndata)
     mSALT(:,k) = temp_s / float(ndata)
     write(6,*) ndata,mTemp(1,k),mSalt(1,k)
   enddo

endif

if(monthly_forcing) then
  allocate(SST_MONTHLY(nlon,nlat,nMonths), SSS_MONTHLY(nlon,nlat,nMonths))
  allocate(TAUX_MONTHLY(nlon,nlat,nMonths), TAUY_MONTHLY(nlon,nlat,nMonths))
  SST_MONTHLY=0; SSS_MONTHLY=0; TAUX_MONTHLY=0; TAUY_MONTHLY=0
  iMonth=1
  fileNameT='TS/SST.shf.normal_year+Hurrell.monthly.01.interp3600x2431.nc'
  fileNameS='TS/SSS.sfwf.CORE_SSS+precip.monthly01.interp3600x2431.nc'
  fileNameU='TS/TAUIJ.ws.old_ncep_1958-2000avg.mon.01.interp3600x2431.nc'
  call read_MONTHLY_init(nlon, nlat, ndepth,fileNameT, fileNameS, fileNameU)
  write(6,*) nlon,nlat,ndepth
  call read_MONTHLY_fields(SST_MONTHLY(:,:,iMonth), SSS_MONTHLY(:,:,iMonth), TAUX_MONTHLY(:,:,iMonth), TAUY_MONTHLY(:,:,iMonth))
  call read_MONTHLY_finalize()
  
  iMonth=2
  fileNameT='TS/SST.shf.normal_year+Hurrell.monthly.02.interp3600x2431.nc'
  fileNameS='TS/SSS.sfwf.CORE_SSS+precip.monthly02.interp3600x2431.nc'
  fileNameU='TS/TAUIJ.ws.old_ncep_1958-2000avg.mon.02.interp3600x2431.nc'
  call read_MONTHLY_init(nlon, nlat, ndepth,fileNameT, fileNameS, fileNameU)
  write(6,*) nlon,nlat,ndepth
  call read_MONTHLY_fields(SST_MONTHLY(:,:,iMonth), SSS_MONTHLY(:,:,iMonth), TAUX_MONTHLY(:,:,iMonth), TAUY_MONTHLY(:,:,iMonth))
  call read_MONTHLY_finalize()
  
  iMonth=3
  fileNameT='TS/SST.shf.normal_year+Hurrell.monthly.03.interp3600x2431.nc'
  fileNameS='TS/SSS.sfwf.CORE_SSS+precip.monthly03.interp3600x2431.nc'
  fileNameU='TS/TAUIJ.ws.old_ncep_1958-2000avg.mon.03.interp3600x2431.nc'
  call read_MONTHLY_init(nlon, nlat, ndepth,fileNameT, fileNameS, fileNameU)
  write(6,*) nlon,nlat,ndepth
  call read_MONTHLY_fields(SST_MONTHLY(:,:,iMonth), SSS_MONTHLY(:,:,iMonth), TAUX_MONTHLY(:,:,iMonth), TAUY_MONTHLY(:,:,iMonth))
  call read_MONTHLY_finalize()
  
  iMonth=4
  fileNameT='TS/SST.shf.normal_year+Hurrell.monthly.04.interp3600x2431.nc'
  fileNameS='TS/SSS.sfwf.CORE_SSS+precip.monthly04.interp3600x2431.nc'
  fileNameU='TS/TAUIJ.ws.old_ncep_1958-2000avg.mon.04.interp3600x2431.nc'
  call read_MONTHLY_init(nlon, nlat, ndepth,fileNameT, fileNameS, fileNameU)
  write(6,*) nlon,nlat,ndepth
  call read_MONTHLY_fields(SST_MONTHLY(:,:,iMonth), SSS_MONTHLY(:,:,iMonth), TAUX_MONTHLY(:,:,iMonth), TAUY_MONTHLY(:,:,iMonth))
  call read_MONTHLY_finalize()
  
  iMonth=5
  fileNameT='TS/SST.shf.normal_year+Hurrell.monthly.05.interp3600x2431.nc'
  fileNameS='TS/SSS.sfwf.CORE_SSS+precip.monthly05.interp3600x2431.nc'
  fileNameU='TS/TAUIJ.ws.old_ncep_1958-2000avg.mon.05.interp3600x2431.nc'
  call read_MONTHLY_init(nlon, nlat, ndepth,fileNameT, fileNameS, fileNameU)
  write(6,*) nlon,nlat,ndepth
  call read_MONTHLY_fields(SST_MONTHLY(:,:,iMonth), SSS_MONTHLY(:,:,iMonth), TAUX_MONTHLY(:,:,iMonth), TAUY_MONTHLY(:,:,iMonth))
  call read_MONTHLY_finalize()
  
  iMonth=6
  fileNameT='TS/SST.shf.normal_year+Hurrell.monthly.06.interp3600x2431.nc'
  fileNameS='TS/SSS.sfwf.CORE_SSS+precip.monthly06.interp3600x2431.nc'
  fileNameU='TS/TAUIJ.ws.old_ncep_1958-2000avg.mon.06.interp3600x2431.nc'
  call read_MONTHLY_init(nlon, nlat, ndepth,fileNameT, fileNameS, fileNameU)
  write(6,*) nlon,nlat,ndepth
  call read_MONTHLY_fields(SST_MONTHLY(:,:,iMonth), SSS_MONTHLY(:,:,iMonth), TAUX_MONTHLY(:,:,iMonth), TAUY_MONTHLY(:,:,iMonth))
  call read_MONTHLY_finalize()
  
  iMonth=7
  fileNameT='TS/SST.shf.normal_year+Hurrell.monthly.07.interp3600x2431.nc'
  fileNameS='TS/SSS.sfwf.CORE_SSS+precip.monthly07.interp3600x2431.nc'
  fileNameU='TS/TAUIJ.ws.old_ncep_1958-2000avg.mon.07.interp3600x2431.nc'
  call read_MONTHLY_init(nlon, nlat, ndepth,fileNameT, fileNameS, fileNameU)
  write(6,*) nlon,nlat,ndepth
  call read_MONTHLY_fields(SST_MONTHLY(:,:,iMonth), SSS_MONTHLY(:,:,iMonth), TAUX_MONTHLY(:,:,iMonth), TAUY_MONTHLY(:,:,iMonth))
  call read_MONTHLY_finalize()
  
  iMonth=8
  fileNameT='TS/SST.shf.normal_year+Hurrell.monthly.08.interp3600x2431.nc'
  fileNameS='TS/SSS.sfwf.CORE_SSS+precip.monthly08.interp3600x2431.nc'
  fileNameU='TS/TAUIJ.ws.old_ncep_1958-2000avg.mon.08.interp3600x2431.nc'
  call read_MONTHLY_init(nlon, nlat, ndepth,fileNameT, fileNameS, fileNameU)
  write(6,*) nlon,nlat,ndepth
  call read_MONTHLY_fields(SST_MONTHLY(:,:,iMonth), SSS_MONTHLY(:,:,iMonth), TAUX_MONTHLY(:,:,iMonth), TAUY_MONTHLY(:,:,iMonth))
  call read_MONTHLY_finalize()
  
  iMonth=9
  fileNameT='TS/SST.shf.normal_year+Hurrell.monthly.09.interp3600x2431.nc'
  fileNameS='TS/SSS.sfwf.CORE_SSS+precip.monthly09.interp3600x2431.nc'
  fileNameU='TS/TAUIJ.ws.old_ncep_1958-2000avg.mon.09.interp3600x2431.nc'
  call read_MONTHLY_init(nlon, nlat, ndepth,fileNameT, fileNameS, fileNameU)
  write(6,*) nlon,nlat,ndepth
  call read_MONTHLY_fields(SST_MONTHLY(:,:,iMonth), SSS_MONTHLY(:,:,iMonth), TAUX_MONTHLY(:,:,iMonth), TAUY_MONTHLY(:,:,iMonth))
  call read_MONTHLY_finalize()
  
  iMonth=10
  fileNameT='TS/SST.shf.normal_year+Hurrell.monthly.10.interp3600x2431.nc'
  fileNameS='TS/SSS.sfwf.CORE_SSS+precip.monthly10.interp3600x2431.nc'
  fileNameU='TS/TAUIJ.ws.old_ncep_1958-2000avg.mon.10.interp3600x2431.nc'
  call read_MONTHLY_init(nlon, nlat, ndepth,fileNameT, fileNameS, fileNameU)
  write(6,*) nlon,nlat,ndepth
  call read_MONTHLY_fields(SST_MONTHLY(:,:,iMonth), SSS_MONTHLY(:,:,iMonth), TAUX_MONTHLY(:,:,iMonth), TAUY_MONTHLY(:,:,iMonth))
  call read_MONTHLY_finalize()
  
  iMonth=11
  fileNameT='TS/SST.shf.normal_year+Hurrell.monthly.11.interp3600x2431.nc'
  fileNameS='TS/SSS.sfwf.CORE_SSS+precip.monthly11.interp3600x2431.nc'
  fileNameU='TS/TAUIJ.ws.old_ncep_1958-2000avg.mon.11.interp3600x2431.nc'
  call read_MONTHLY_init(nlon, nlat, ndepth,fileNameT, fileNameS, fileNameU)
  write(6,*) nlon,nlat,ndepth
  call read_MONTHLY_fields(SST_MONTHLY(:,:,iMonth), SSS_MONTHLY(:,:,iMonth), TAUX_MONTHLY(:,:,iMonth), TAUY_MONTHLY(:,:,iMonth))
  call read_MONTHLY_finalize()
  
  iMonth=12
  fileNameT='TS/SST.shf.normal_year+Hurrell.monthly.12.interp3600x2431.nc'
  fileNameS='TS/SSS.sfwf.CORE_SSS+precip.monthly12.interp3600x2431.nc'
  fileNameU='TS/TAUIJ.ws.old_ncep_1958-2000avg.mon.12.interp3600x2431.nc'
  call read_MONTHLY_init(nlon, nlat, ndepth,fileNameT, fileNameS, fileNameU)
  write(6,*) nlon,nlat,ndepth
  call read_MONTHLY_fields(SST_MONTHLY(:,:,iMonth), SSS_MONTHLY(:,:,iMonth), TAUX_MONTHLY(:,:,iMonth), TAUY_MONTHLY(:,:,iMonth))
  call read_MONTHLY_finalize()
end if

! generate initial conditions
call get_init_conditions

! dump new grid to netCDF
write(6,*) ' calling write_grid'
write(6,*)
call write_grid

! dump graph for partioning
write(6,*) ' call write_graph'
write(6,*)
call write_graph

! write OpenDx file
if (write_OpenDX_flag) then
   write(6,*) ' calling write_OpenDX'
   write(6,*)
   call write_OpenDX(        on_a_sphere, &
                             nCellsNew, &
                             nVerticesNew, &
                             nEdgesNew, &
                             vertexDegreeNew, &
                             maxEdgesNew, &
                             xCellNew, &
                             yCellNew, &
                             zCellNew, &
                             xVertexNew, &
                             yVertexNew, &
                             zVertexNew, &
                             xEdgeNew, &
                             yEdgeNew, &
                             zEdgeNew, &
                             nEdgesOnCellNew, &
                             verticesOnCellNew, &
                             verticesOnEdgeNew, &
                             cellsOnVertexNew, &
                             edgesOnCellNew, &
                             areaCellNew, &
                             maxLevelCellNew, &
                             meshDensityNew, &
                             bottomDepthNew, &
                             temperatureNew(1,1,:), &
                             kiteAreasOnVertexNew )
endif

!do iCell=1,nCellsNew
  !ulon = 1.0; ulat = 0.0
  !xin = xCellNew(iCell); yin = yCellNew(iCell); zin = zCellNew(iCell)
  !call transform_from_lonlat_to_xyz(xin, yin, zin, ulon, ulat, ux, uy, uz)
  !if(abs(ux).lt.1.0e-10) ux=0.0
  !if(abs(uy).lt.1.0e-10) uy=0.0
  !if(abs(uz).lt.1.0e-10) uz=0.0
  !write(20,10) ux, uy, uz
  !10 format(3e25.10)
!enddo
  
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Step 3: Check get_init_conditions routine for initial T&S, thickness, etc.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_init_conditions
implicit none
real :: halfwidth, dtr, pi, p(3), q(3), xin, yin, zin, ulon, ulat, stress, n1, n2, distance, r, temp_t, temp_s
real :: dotProd, rho_ref, rho_delta, work, factor, r1, deltay, b
real :: y0_embayment, surfaceTemperature, bottomTemperature, betaTemperature,bottomMinTemp, &
   y_a,y_0,A_0, midDepth(nVertLevelsMod), maxMidDepth
real :: x_0, x_1,x_2,x_3,width,cff1, deltaTemperature
integer :: iTracer, ix, iy, ndata, i, j, k, ixt, iyt, ncull, jcount, iNoData, kdata(nVertLevelsMod), iMonth
logical :: flag_lat

pi = 4.0*atan(1.0)
dtr = pi/180.0

! defaults
layerThicknessNew = 100.0
temperatureNew = 1.0
salinityNew = 1.0
tracer1New = 1.0
normalVelocityNew = 0
normalVelocityForcingNew = 0
densityNew = 1025.0

if (initial_conditions.eq.'uniform_TS') then

   temperatureNew = uniform_T !10.0
   salinityNew = uniform_S !34.0
   tracer1New = uniform_tracer1 !1.0
   normalVelocityNew = 0.0
   normalVelocityForcingNew = 0.0

   do iCell=1,nCellsNew
      do k=1, nVertLevelsMod
         layerThicknessNew(1,k,iCell) = hZLevel(k)
      enddo
   enddo  

elseif (initial_conditions.eq.'Ilicak2_overflow') then

   do i = 1,nCellsNew
      if(yCellNew(i) < 20000) then
         temperatureNew(1,:,i) = 10.0
      else
         temperatureNew(1,:,i) = 20.0
      endif
   enddo

   salinityNew(1,:,:) = 35.0
   Tracer1New(1,:,:) = 1.0
   normalVelocityNew = 0.0
   normalVelocityForcingNew = 0.0

   do iCell=1,nCellsNew
      do k=1,nVertLevelsMod
         layerThicknessNew(1,k,iCell) = hZLevel(k)
      enddo
   enddo  

elseif (initial_conditions.eq.'Ilicak2_overflow_sigma') then

   do i = 1,nCellsNew
      if(yCellNew(i) < 20000) then
         temperatureNew(1,:,i) = 10.0
      else
         temperatureNew(1,:,i) = 20.0
      endif
   enddo

   salinityNew(1,:,:) = 35.0
   Tracer1New(1,:,:) = 1.0
   normalVelocityNew = 0.0
   normalVelocityForcingNew = 0.0

   do iCell=1,nCellsNew
      do k=1,nVertLevelsMod
         layerThicknessNew(1,k,iCell) = bottomDepthNew(iCell) / nVertLevelsMOD
      enddo
   enddo  

elseif (initial_conditions.eq.'internal_wave_Ilicak3') then

   salinityNew(1,:,:) = 35.0
   normalVelocityNew = 0.0
   normalVelocityForcingNew = 0.0

   surfaceTemperature = 20.1
   bottomTemperature = 10.1

   midDepth(1) =  hZLevel(1)/2.0
   do k=2,nVertLevelsMod
      midDepth(k) = midDepth(k-1) + 0.5*(hZLevel(k-1) + hZLevel(k))
   enddo  

   maxMidDepth = maxval(midDepth)
   do k = 1, nVertLevelsMOD
      temperatureNew(1,k,:) = (surfaceTemperature - bottomTemperature) &
         * ((maxMidDepth - midDepth(k))/maxMidDepth) + bottomTemperature
   enddo

   A_0 = 2.0
   y_0 = 125.0e3
   y_a = 50.0e3
   do i = 1, nCellsNew
      if ( abs(yCellNew(i) - y_0) < y_a) then
         do k = 2, nVertLevelsMOD
            betaTemperature = -A_0*cos(0.5*pi*(yCellNew(i)-y_0)/y_a) &
               *sin(pi*refBottomDepth(k-1)/(maxMidDepth - midDepth(1))) 
            temperatureNew(1, k, i) = temperatureNew(1, k, i) + betaTemperature
         end do
      endif
   enddo

   do iCell=1,nCellsNew
      do k=1,nVertLevelsMod
         layerThicknessNew(1,k,iCell) = hZLevel(k)
      enddo
   enddo 
 
elseif (initial_conditions.eq.'baroclinic_channel_Ilicak4') then

   salinityNew(1,:,:) = 35.0
   normalVelocityNew = 0.0
   normalVelocityForcingNew = 0.0

   midDepth(1) =  hZLevel(1)/2.0
   do k=2,nVertLevelsMod
      midDepth(k) = midDepth(k-1) + 0.5*(hZLevel(k-1) + hZLevel(k))
   enddo  

   maxMidDepth = maxval(midDepth)

   ! Set up stratification on northern half
   surfaceTemperature = 13.1
   bottomTemperature = 10.1
   deltaTemperature = 1.2
   do k = 1, nVertLevelsMOD
      temperatureNew(1,k,:) = bottomTemperature &
         + (surfaceTemperature - bottomTemperature) * ((-midDepth(k)+maxMidDepth)/maxMidDepth)
   enddo

   y_0 = 250.0e3
   x_0 = 0.0e3
   x_1 = 160.0e3
   x_2 = 110.0e3
   x_3 = 130.0e3
   width = 40.0e3
   do i = 1, nCellsNew
      ! sine wave defines southern boundary between two temperatures.
      cff1 = width * sin (6.0 * 3.141592 * (xCellNew(i) - x_0)/(x_1 - x_0))

      ! stratification on south half
      if( yCellNew(i) < y_0 - cff1 ) then
         do k = 1, nVertLevelsMOD
            temperatureNew(1,k,i) = temperatureNew(1,k,i) - deltaTemperature
         end do

      ! linear interpolation between two halves
      else if( yCellNew(i) .ge. y_0 - cff1 .and. yCellNew(i) .le. y_0 - cff1+width) then
         do k = 1, nVertLevelsMOD
            temperatureNew(1,k,i) = temperatureNew(1,k,i) &
               - deltaTemperature*(1.0 -( yCellNew(i) - (y_0 - cff1)) / (1.0 * width))
         end do
      endif
   enddo

   ! Add an additional perturbation to southern end of third wave.
   do i = 1, nCellsNew
      cff1 = 0.5 * width * sin(1.0 * 3.141592 * (xCellNew(i) - x_2)/(x_3 - x_2))
      if( yCellNew(i) .ge. y_0 - cff1-0.5*width .and. &
          yCellNew(i) .le. y_0 - cff1+0.5*width .and. &
          xCellNew(i) .ge. x_2 .and. xCellNew(i) .le. x_3) then
         do k = 1, nVertLevelsMOD
            temperatureNew(1,k,i) = temperatureNew(1,k,i) + 0.3 * (1.0 - ( (yCellNew(i)-(y_0-cff1))/(0.5*width)))
         end do
      endif
   end do

   do iCell=1,nCellsNew
      do k=1,nVertLevelsMod
         layerThicknessNew(1,k,iCell) = hZLevel(k)
      enddo
   enddo 

elseif (initial_conditions.eq.'DOME_3D_overflow') then

   y0_embayment = 600.0e3 ! y location of beginning of embayment
   Tracer1New(1,:,:) = 0.0
   normalVelocityNew = 0.0
   normalVelocityForcingNew = 0.0
   do k=1,nVertLevelsMOD
      temperatureNew(1,k,:) = 20-0.5*k
      salinityNew(1,k,:) = 35.0
   enddo

   ! Change embayment fluid to have tracer 1 and cold water.
   ! This is only used for a plug of cold water in the initial conditions,
   ! rather than forced cold water at the inlet.
   !do iCell=1,nCellsNew
   !   if (yCellNew(iCell).gt.y0_embayment) then
   !      temperatureNew(1,:,iCell) = 0.0
   !      tracer1New(1,:,iCell) = 1.0
   !   endif
   !enddo

   do iCell=1,nCellsNew
      do k=1,nVertLevelsMod
         layerThicknessNew(1,k,iCell) = hZLevel(k)
      enddo
   enddo  

elseif (initial_conditions.eq.'SOMA_TS') then

  temperatureNew = -99.0
  salinityNew = -99.0
  rho_ref=1000.0
  rho_delta=5.0
  do iCell=1,nCellsNew
  do k = 1,maxLevelCellNew(iCell)
    layerThicknessNew(1,k,iCell) = hZLevel(k)
    if(k.eq.1) r1 = -refBottomDepth(k)/2.0
    if(k.ne.1) r1 = -(refBottomDepth(k)+refBottomDepth(k-1))/2.0
    work = rho_ref - (1.0-0.05)*rho_delta*tanh(r1/300) - 0.05*rho_delta*r1/2500
    densityNew(1,k,iCell) = work
    factor = (rho_ref-work)/2.5e-1
    temperatureNew(1,k,iCell) = 15.0 + factor
    salinityNew(1,k,iCell) = 0.0  ! salinity
  enddo
  enddo

  normalVelocityForcingNew = 0.0
  do iEdge=1,nEdgesNew
    xin = xEdgeNew(iEdge)
    yin = yEdgeNew(iEdge)
    zin = zEdgeNew(iEdge)
    rlon = lonEdgeNew(iEdge)
    rlat = latEdgeNew(iEdge)
  
    b=1.25e6
    deltay = sphere_radius * ( rlat - 35.0*dtr)
    factor = 1.0-0.5*deltay/b
    r1 = factor * 0.1 * exp( -(deltay/b)**2 ) * cos(pi*deltay/b)
    ulon = r1
    ulat = 0.0
    call transform_from_lonlat_to_xyz(xin,yin,zin,ulon,ulat,ux,uy,uz)
    if(boundaryEdgeNew(1,iEdge).eq.1) then
      normalVelocityForcingNew(1,iEdge) = 0.0
    else
      iCell1 = cellsOnEdgeNew(1,iEdge)
      iCell2 = cellsOnEdgeNew(2,iEdge)
      p(1) = xCellNew(iCell1); p(2) = yCellNew(iCell1); p(3) = zCellNew(iCell1)
      q(1) = xCellNew(iCell2); q(2) = yCellNew(iCell2); q(3) = zCellNew(iCell2)
      q = q - p
      call unit_vector_in_3space(q)
      normalVelocityForcingNew(1,iEdge) = ux*q(1) + uy*q(2) + uz*q(3)
    endif

   enddo

elseif (initial_conditions.eq.'isopycnal_3layer') then

   fEdgeNew(:) = 0.0
   fVertexNew(:) = 0.0
   bottomDepthNew(:) = 0.0
   normalVelocityNew(:,:,:) = 0.0

 ! setting for three levels - Set h values for isopycnal system
   write(6,*) ' setting three levels for isopycnal system'
       layerThicknessNew(1,1,:) = 500.0
       layerThicknessNew(1,2,:) = 1250.0
       layerThicknessNew(1,3,:) = 3250.0
       bottomDepthNew(:) = -( layerThicknessNew(1,1,:) + layerThicknessNew(1,2,:) + layerThicknessNew(1,3,:) )

   ! Noise is meant to make the flow unstable at some point 
   ! Not needed for all simulations
   write(6,*) ' adding noise to layer thickness'
   r = 0.0
   do i=1,nCellsNew
      work1(i) = float(i) / float(nCellsNew)
      call random_number(work1(i))
      r = r + work1(i)
   enddo
   r = r/float(nCells)
   work1(:) = work1(:) - r
   layerThicknessNew(1,1,:) = layerThicknessNew(1,1,:) + 1.0*work1(:)
   layerThicknessNew(1,2,:) = layerThicknessNew(1,2,:) - 1.0*work1(:)

   ! Specify Density values for isopycnal levels
   write(6,*) ' setting density - depricate soon'
   densityNew(1,:,:) = 1010.0
   densityNew(1,2,:) = 1011.0
   densityNew(1,3,:) = 1012.0

   eos_linear_alpha = 2.55e-1
   eos_linear_beta = 7.64e-1
   eos_linear_Tref = 19.0
   eos_linear_Sref = 35.0
   eos_linear_densityref = 1025.022

   ! set salinity for isopycnal levels
   salinityNew = eos_linear_Sref

   ! set temperature for isopycnal levels.  Just invert linear eos.
   write(6,*) ' setting temperature'
   do k=1,nVertLevelsMOD
     temperatureNew(1,k,:) = eos_linear_Tref + (eos_linear_densityref - densityNew(1,k,:))/eos_linear_alpha
   enddo

   ! set forcing for isopycnal levels
   write(6,*) 'setting normalVelocityForcing - wind forcing'
   normalVelocityForcingNew = 0.0
   if(on_a_sphere.eq.'YES') then
       latmin = -60*dtr
       latmax = -10*dtr
       latmid = -35*dtr
       latmin = minval(latEdgeNew)
       latmax = maxval(latEdgeNew)
       latmid = (latmin+latmax)/2.0
       r = 10.0*dtr

       write(6,*) 'normalVelocityForcing info', latmin, latmax, latmid, r
       do i = 1,nEdgesNew
           lattmp = latEdgeNew(i)
           iCell1 = cellsOnEdgeNew(1,i)
           iCell2 = cellsOnEdgeNew(2,i)
           if(iCell1>0.and.iCell2>0) then
               pert =  normalVelocityForcing_max * exp(-(lattmp-latmid)**2/(r**2))

               ulat = latEdgeNew(i)
               ulon = lonEdgeNew(i) + 0.05

               call convert_lx(xin, yin, zin, 1.0, ulat, ulon)

               xin = xin - xEdgeNew(i)
               yin = yin - yEdgeNew(i)
               zin = zin - zEdgeNew(i)

               dotProd = sqrt(xin**2 + yin**2 + zin**2)
               xin = xin/dotProd
               yin = yin/dotProd
               zin = zin/dotProd

               dotProd = normalsNew(1,i)*xin + normalsNew(2,i)*yin + normalsNew(3,i)*zin

               normalVelocityForcingNew(1,i) = pert * dotProd
               write(8,*) lattmp,pert,dotProd
           endif
       enddo
   else
       ymin = minval(yEdgeNew)
       ymax = maxval(yEdgeNew)
       r = 3.0e5
       ymid = (ymax+ymin)/2
       do i = 1,nEdgesNew
           ytmp = yEdgeNew(i)
           iCell1 = cellsOnEdgeNew(1,i)
           iCell2 = cellsOnEdgeNew(2,i)
           if(iCell1>0.and.iCell2>0) then
               pert =  normalVelocityForcing_max * exp(-(ytmp-ymid)**2/(r**2))
               write(8,*) ytmp,pert
               normalVelocityForcingNew(1,i) = pert * normalsNew(1,i)
           endif
       enddo
   endif
   write(6,*) ' normalVelocityForcingNew ', minval(normalVelocityForcingNew), maxval(normalVelocityForcingNew)

elseif (initial_conditions.eq.'realistic_WOCE') then


normalVelocityForcingNew = 0.0
windStressMonthlyNew = 0.0
do iEdge=1,nEdgesNew
  xin = xEdgeNew(iEdge)
  yin = yEdgeNew(iEdge)
  zin = zEdgeNew(iEdge)
  rlon = lonEdgeNew(iEdge)/dtr
  rlat = latEdgeNew(iEdge)/dtr
  ix = nint(rlon/0.1 - 0.05) + nlon + 1
  ix = mod(ix,nlon)+1
  iy = nlat
  do jcount=1,nlat
   if(t_lat(jcount).gt.rlat) then
    iy = jcount
    exit
   endif
  enddo
  ulon = TAUX(ix,iy)
  ulat = TAUY(ix,iy)
  !write(6,*) rlon, t_lon(ix), rlat, t_lat(iy)

  call transform_from_lonlat_to_xyz(xin,yin,zin,ulon,ulat,ux,uy,uz)
  if(boundaryEdgeNew(1,iEdge).eq.1) then
    normalVelocityForcingNew(1,iEdge) = 0.0
  else
    iCell1 = cellsOnEdgeNew(1,iEdge)
    iCell2 = cellsOnEdgeNew(2,iEdge)
    p(1) = xCellNew(iCell1); p(2) = yCellNew(iCell1); p(3) = zCellNew(iCell1)
    q(1) = xCellNew(iCell2); q(2) = yCellNew(iCell2); q(3) = zCellNew(iCell2)
    q = q - p
    call unit_vector_in_3space(q)
    normalVelocityForcingNew(1,iEdge) = ux*q(1) + uy*q(2) + uz*q(3)
  endif

  if(monthly_forcing) then
    do iMonth=1,nMonths
    ulon = TAUX_MONTHLY(ix,iy,iMonth)
    ulat = TAUY_MONTHLY(ix,iy,iMonth)
  ! if(abs(ulon).gt.1.0.or.abs(ulat).gt.1.0) then
  !    ulon=0.0
  !    ulat=0.0
  ! endif
    call transform_from_lonlat_to_xyz(xin,yin,zin,ulon,ulat,ux,uy,uz)
    if(boundaryEdgeNew(1,iEdge).eq.1) then
      windStressMonthlyNew(iMonth,iEdge) = 0.0
    else
      iCell1 = cellsOnEdgeNew(1,iEdge)
      iCell2 = cellsOnEdgeNew(2,iEdge)
      p(1) = xCellNew(iCell1); p(2) = yCellNew(iCell1); p(3) = zCellNew(iCell1)
      q(1) = xCellNew(iCell2); q(2) = yCellNew(iCell2); q(3) = zCellNew(iCell2)
      q = q - p
      call unit_vector_in_3space(q)
    ! repeat
      windStressMonthlyNew(iMonth,iEdge) = ux*q(1) + uy*q(2) + uz*q(3)
    ! windStressMonthlyNew(iMonth,iEdge) = normalVelocityForcingNew(1,iEdge)
    endif
    enddo
  else
    windStressMonthlyNew(:,:) = 0.0
  end if
enddo

! for acc runs, increase wind strength for Southern Ocean
if (amplify_acc_wind) then
  print *, 'amplify_acc_wind, amp_wind_factor, amp_wind_center_lat, amp_wind_spread_lat'
  print *, amplify_acc_wind, amp_wind_factor, amp_wind_center_lat, amp_wind_spread_lat
  do iEdge=1,nEdgesNew
     windStressMonthlyNew(:,iEdge) = windStressMonthlyNew(:,iEdge) &
       * (1.0 + (amp_wind_factor-1.0)*0.5 &
          *(1.0-tanh( (latEdgeNew(iEdge)/dtr-amp_wind_center_lat)/amp_wind_spread_lat) ) )
  enddo
endif
 
!set tracers at a first guess
temperatureNew = -99.0
salinityNew = -99.0
do iCell=1,nCellsNew
do k = 1,maxLevelCellNew(iCell)
  temperatureNew(1,k,iCell) = 20.0 - 10.0*k/nVertLevelsMod
  salinityNew(1,k,iCell) = 34.0  ! salinity
enddo
enddo

! update T and S field with WOCE data
if(load_woce_IC) then
iNoData = 0
do iCell=1,nCellsNew
  layerThicknessNew(1,:,iCell) = dz(:)
  ! if(mod(iCell,100).eq.0) write(6,*) 'load_woce_IC t and s',iCell
  rlon = lonCellNew(iCell)/dtr
  rlat = latCellNew(iCell)/dtr
  ix = nint(rlon/0.1 - 0.05) + nlon + 1
  ix = mod(ix,nlon)+1
  iy = nlat
  do j=1,nlat
   if(t_lat(j).gt.rlat) then
    iy = j
    exit
   endif
  enddo
  do k=1,maxLevelCellNew(iCell)
    ndata = 0; temp_t = 0; temp_s = 0; kdata(:) = 0
    do i=-15,15
      ixt = ix + 8*i
      if(ixt.lt.1) then
        ixt = ixt + nlon
      elseif(ixt.gt.nlon) then
        ixt = ixt - nlon
      endif
      do j=-15,15
        iyt = iy + 8*j
        flag_lat = .true.
        if(iyt.lt.1.or.iyt.gt.nlat) then
          iyt = 1
          flag_lat = .false.
        endif
        if(TEMP(ixt,iyt,k).gt.-10.0.and.flag_lat) then
          ndata = ndata + 1
          temp_t = temp_t + TEMP(ixt,iyt,k)
          temp_s = temp_s + SALT(ixt,iyt,k)
        endif
      enddo
    enddo

    if(ndata.gt.0) then
      temperatureNew(1,k,iCell) = temp_t / float(ndata)
      salinityNEW(1,k,iCell) = temp_s / float(ndata)
      kdata(k) = 1
    else
      if(k.eq.1) iNoData = iNoData + 1
      if(k.ge.3) then
        if(kdata(k-1).eq.1) maxLevelCellNew(iCell) = k-1
      endif
    endif

  enddo

enddo

! do a couple of smoothing passes
do iter=1,5
do iCell=1,nCellsNew
do k=1,maxLevelCellNew(iCell)
  ndata=1
  temp_t = temperatureNew(1,k,iCell)
  temp_s = salinityNew(1,k,iCell)
  do j=1,nEdgesOnCellNew(iCell)
    jCell = cellsOnCellNew(j,iCell)
    if(jCell.gt.0) then
    if(maxLevelCellNew(jCell).ge.k) then
      temp_t = temp_t + temperatureNew(1,k,jCell)
      temp_s = temp_s + salinityNew(1,k,jCell)
      ndata = ndata + 1
    endif
    endif
  enddo
  temperatureNew(1,k,iCell) = temp_t / ndata
  salinityNew(1,k,iCell) = temp_s / ndata
enddo
enddo
write(6,*) maxval(temperatureNew(1,1,:)),maxval(salinityNew(1,1,:))
enddo

write(6,*) iNoData, nCellsNew

temperatureRestoreNew(:) = temperatureNew(1,1,:)
salinityRestoreNew(:) = salinityNew(1,1,:)

if(monthly_forcing) then
  do iMonth=1,nMonths
  iNoData = 0
  do iCell=1,nCellsNew
    ! if(mod(iCell,100).eq.0) write(6,*) 'load_woce_IC t and s RESTORE',iCell
    rlon = lonCellNew(iCell)/dtr
    rlat = latCellNew(iCell)/dtr
    ix = nint(rlon/0.1 - 0.05) + nlon + 1
    ix = mod(ix,nlon)+1
    iy = nlat
    do j=1,nlat
     if(t_lat(j).gt.rlat) then
      iy = j
      exit
     endif
    enddo  ! j
    k=1
      ndata = 0; temp_t = 0; temp_s = 0
      do i=-15,15
        ixt = ix + 8*i
        if(ixt.lt.1) then
          ixt = ixt + nlon
        elseif(ixt.gt.nlon) then
          ixt = ixt - nlon
        endif
        do j=-15,15
          iyt = iy + 8*j
          flag_lat = .true.
          if(iyt.lt.1.or.iyt.gt.nlat) then
            iyt = 1
            flag_lat = .false.
          endif
          if(SST_MONTHLY(ixt,iyt,iMonth).gt.-10.0.and.flag_lat) then
            ndata = ndata + 1
            temp_t = temp_t + SST_MONTHLY(ixt,iyt,iMonth)
            temp_s = temp_s + SSS_MONTHLY(ixt,iyt,iMonth)
          endif
        enddo !j
       enddo  !i
  
      if(ndata.gt.0) then
        temperatureRestoreMonthlyNew(iMonth,iCell) = temp_t / float(ndata)
        salinityRestoreMonthlyNew(iMonth,iCell) = temp_s / float(ndata)
      else
        temperatureRestoreMonthlyNew(iMonth,iCell) = temperatureNew(1,1,iCell)
        salinityRestoreMonthlyNew(iMonth,iCell) = salinityNew(1,1,iCell)
      endif
  
    enddo  ! iCell
  enddo    ! iMonth
  
  ! do a couple of smoothing passes
  do iter=1,5
  do iCell=1,nCellsNew
    k=1
    ndata=1
    temp_t = temperatureRestoreMonthlyNew(iMonth,iCell)
    temp_s = salinityRestoreMonthlyNew(iMonth,iCell)
    do j=1,nEdgesOnCellNew(iCell)
      jCell = cellsOnCellNew(j,iCell)
      if(jCell.gt.0) then
      if(maxLevelCellNew(jCell).ge.k) then
        temp_t = temp_t + temperatureRestoreMonthlyNew(iMonth,iCell)
        temp_s = temp_s + salinityRestoreMonthlyNew(iMonth,iCell)
        ndata = ndata + 1
      endif
      endif
    enddo  ! j
    temperatureRestoreMonthlyNew(iMonth,iCell) = temp_t / ndata
    salinityRestoreMonthlyNew(iMonth,iCell) = temp_s / ndata
  enddo  ! iCell
  enddo  ! iter
else
  temperatureRestoreMonthlyNew(:,:) = 0.0
  salinityRestoreMonthlyNew(:,:) = 0.0
end if

endif  ! load_woce_IC

!repeat
!do iMonth=1,12
!  temperatureRestoreMonthlyNew(iMonth,:) = temperatureRestoreNew(:)
!  salinityRestoreMonthlyNew(iMonth,:) = salinityRestoreNew(:)
!enddo

else 

   print *, ' Incorrect choice of initial_conditions: ',initial_conditions
   stop

endif  ! initial_conditions


if (trim(initial_conditions).ne.'realistic_WOCE') then

      ! set coriolis parameter for grid
   write(6,*) ' setting Coriolis parameter'
   if(on_a_sphere.eq.'YES') then
       do i = 1,nVerticesNew
          fVertexNew(i) = 2.0 * omega * sin(latVertexNew(i))
       enddo

       do i = 1,nEdgesNew
          fEdgeNew(i) = 2.0 * omega * sin(latEdgeNew(i))
       enddo
   else
       do i = 1,nVerticesNew
          fVertexNew(i) = f0 + (yVertexNew(i) - ymid) * beta
       enddo

       do i = 1,nEdgesNew
          fEdgeNew(i) = f0 + (yEdgeNew(i) - ymid) * beta
       enddo
   endif

endif

write(6,*) ' done get_init_conditions'

end subroutine get_init_conditions


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

allocate(fEdge(nEdges))
allocate(fVertex(nVertices))
allocate(bottomDepth(nCells))
allocate(work1(nCells))
allocate(normalVelocityForcing(nVertLevels,nEdges))
allocate(normalVelocity(1,nVertLevels,nEdges))
allocate(tangentialVelocity(1,nVertLevels,nEdges))
allocate(layerThickness(1,nVertLevels,nCells))
allocate(density(1,nVertLevels,nCells))

xCell=0; yCell=0; zCell=0; latCell=0; lonCell=0; meshDensity=1.0
xEdge=0; yEdge=0; zEdge=0; latEdge=0; lonEdge=0
xVertex=0; yVertex=0; zVertex=0; latVertex=0; lonVertex=0

indexToCellID=0; indexToEdgeID=0; indexToVertexID=0
cellsOnEdge=0; nEdgesOnCell=0; edgesOnCell=0
edgesOnEdge=0; weightsOnEdge=0
angleEdge=0; areaCell=0; areaTriangle=0
cellsOnCell=0; verticesOnCell=0; verticesOnEdge=0
edgesOnVertex=0; cellsOnVertex=0; kiteAreasOnVertex=0

fEdge=0; fVertex=0; bottomDepth=0; normalVelocityForcing=0; work1=0
normalVelocity=0; tangentialVelocity=0; layerThickness=0; density=0


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
                    fEdge, &
                    fVertex, &
                    bottomDepth, &
                    normalVelocity, &
                    tangentialVelocity, &
                    layerThickness &
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
write(6,*) ' fEdge : ', minval(fEdge), maxval(fEdge)
write(6,*) ' fVertex : ', minval(fVertex), maxval(fVertex)
write(6,*) ' bottomDepth : ', minval(bottomDepth), maxval(bottomDepth)
write(6,*) ' normalVelocity : ', minval(normalVelocity), maxval(normalVelocity)
write(6,*) ' tangentialVelocity : ', minval(tangentialVelocity), maxval(tangentialVelocity)
write(6,*) ' layerThickness : ', minval(layerThickness), maxval(layerThickness)

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
                nMonths, &
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
                    maxLevelCellNew, &
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
                    kiteAreasOnVertexNew, &
                    fEdgeNew, &
                    fVertexNew, &
                    bottomDepthNew, &
                    boundaryEdgeNew, &
                    boundaryVertexNew, &
                    normalVelocityForcingNew, &
                    windStressMonthlyNew, &
                    normalVelocityNew, &
                    layerThicknessNew, & 
                    densityNew, &
                    temperatureNew, &
                    salinityNew, &
                    tracer1New, &
                    temperatureRestoreNew, &
                    salinityRestoreNew, &
                    temperatureRestoreMonthlyNew, &
                    salinityRestoreMonthlyNew, &
                    refBottomDepth &
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
integer :: kEmbayment
real :: pi, dtr, zdata, rlon, rlat, r, ymin, ymax, xmin, xmax
real :: latmin, latmax, lonmin, lonmax, ridgeDepth, maxdc
real :: b, H0, phi, gamma, distance, rlonC, rlatC
real :: y0_embayment, embayment_depth, embayment_slope, z
logical :: flag, kmt_flag

pi = 4.0*atan(1.0)
dtr = pi / 180.0

allocate(kmt(nCells))
kmt = 0

if (bottom_topography.eq.'realistic_ETOPO') then

    nx = 10800
    ny = 5400
    allocate(x(nx))
    allocate(y(ny))
    allocate(ztopo(nx,ny))
    x = 0.0
    y = 0.0
    ztopo = 0.0
    write(6,*) ' ztopo ', minval(ztopo), maxval(ztopo)
    call read_topo_init( inx, iny)
    if(inx.ne.nx) stop ' nx topo'
    if(iny.ne.ny) stop ' ny topo'
    call read_topo_fields(x,y,ztopo)
    call read_topo_finalize()
    write(6,*) minval(x), maxval(x), x(1)
    write(6,*) minval(y), maxval(y), y(1)
    write(6,*) minval(ztopo), maxval(ztopo)

    do iCell=1,nCells

        ! Convert from radians to degrees
        rlon = lonCell(iCell) / dtr
        rlat = latCell(iCell) / dtr

        ! Find nearest coordinate in topo file.
        ! This is 1/30th degree topo data, so multiply degrees by 30 and
        ! round to get index.
        ix = nint((rlon+180)*30) + nx
        ix = mod(ix,nx)+1
        iy = nint((rlat+90 )*30)
        ix = max(1,ix); ix = min(nx,ix)
        iy = max(1,iy); iy = min(ny,iy)

        zdata = ztopo(ix,iy)

        ! zdata is less than zero for ocean points.
        if(zdata.lt.0.0) then
            zdata = -zdata
            bottomDepth(iCell) = zdata
            r = 0
            kmt_flag=.false.
            do k=1,nVertLevelsMod
                if(.not.kmt_flag) then
                    r = r + dz(k)
                    if(r.gt.zdata) then
                        kmt(iCell) = k
                        kmt_flag = .true.
                    endif
                endif
            enddo

            ! zdata is deeper than deepest cell
            if (kmt(iCell).eq.0) then
               kmt(iCell)=nVertLevelsMod
               bottomDepth(iCell) = refBottomDepth(nVertLevelsMod)
            endif

            ! write(6,*) kmt(iCell)

        endif

        ! if(zdata.lt.0.0) kmt(iCell) = nVertLevelsMod

    enddo

    deallocate(x)
    deallocate(y)
    deallocate(ztopo)

elseif (bottom_topography.eq.'Ilicak2_overflow') then

    kmt = nVertLevelsMOD

    ridgeDepth = 500.0
    do iCell = 1,nCells
      ! From Mehmet Ilicak:
      ! depth=2000
      ! val1 = 500 is top of ridge
      ! h(i,j) = val1 + 0.5*(depth-val1) * (1.0+TANH((lon(i,j)-40000.0)/7000.0)) 
      bottomDepth(iCell) = ridgeDepth + 0.5*(layer_thickness_total_max-ridgeDepth) * (1.0+tanh((yCell(iCell)-40000.0)/7000.0)) 

      if (bottomDepth(iCell).lt.0.0.or. &
          bottomDepth(iCell).gt.refBottomDepth(nVertLevelsMOD)) then
         print *, 'error: bottomDepth cannot be less than zero or greater than refBottomDepth(nVertLevels)'
         print *, 'iCell, bottomDepth(iCell):', iCell, bottomDepth(iCell)
         exit
      end if

      do k=1,nVertLevelsMOD
         if (bottomDepth(iCell).le.refBottomDepth(k)) then
            kmt(iCell) = k
            exit
         endif
      end do

    enddo

elseif (bottom_topography.eq.'Ilicak2_overflow_sigma') then

    ridgeDepth = 500.0
    do iCell = 1,nCells
      ! From Mehmet Ilicak:
      ! depth=2000
      ! val1 = 500 is top of ridge
      ! h(i,j) = val1 + 0.5*(depth-val1) * (1.0+TANH((lon(i,j)-40000.0)/7000.0)) 
      bottomDepth(iCell) = ridgeDepth + 0.5*(layer_thickness_total_max-ridgeDepth) * (1.0+tanh((yCell(iCell)-40000.0)/7000.0)) 

      if (bottomDepth(iCell).lt.0.0.or. &
          bottomDepth(iCell).gt.refBottomDepth(nVertLevelsMOD)) then
         print *, 'error: bottomDepth cannot be less than zero or greater than refBottomDepth(nVertLevels)'
         print *, 'iCell, bottomDepth(iCell):', iCell, bottomDepth(iCell)
         exit
      end if
  
    enddo

    ! for sigma coordinates, set kmt to the max level.
    kmt = nVertLevelsMOD

elseif (bottom_topography.eq.'DOME_3D_overflow') then

    y0_embayment = 600.0e3 ! y location of beginning of embayment
    embayment_depth =600 ! depth of embayment, m
    embayment_slope = 0.01 ! slope

    ! The first level where z is deeper than zmid, set k-level of embayment
    do k=2,nVertLevelsMod
       if ((refBottomDepth(k)+refBottomDepth(k-1))/2.0.gt.embayment_depth) then
          kEmbayment = k
          exit
       endif
    enddo
    print *, 'kEmbayment',kEmbayment
    do iCell=1,nCells
      if (yCell(iCell).gt.y0_embayment) then
         ! note: I made embayment 200km wide for the 50km rez hex case, so 
         ! it is always wider than 1 cell.  For higher rez, change this
         ! to 1700 to 1800.
         if (xCell(iCell).lt.1600.0e3.or. &
             xCell(iCell).gt.1800.0e3 ) then
            kmt(iCell)=0
         else
            kmt(iCell) = kEmbayment
         endif
         bottomDepth(iCell) = embayment_depth
        
      else
         ! compute depth based on y coordinate
         bottomDepth(iCell) = min(embayment_depth + embayment_slope*(y0_embayment - yCell(iCell)),refBottomDepth(nVertLevelsMOD))
         ! default depth is 
         kmt(iCell) = nVertLevelsMod
         do k=2,nVertLevelsMod
            ! The first level where z is deeper than zmid, set kmt.
            if ((refBottomDepth(k)+refBottomDepth(k-1))/2.0.gt.bottomDepth(iCell)) then
               kmt(iCell) = k
               exit
            endif
         enddo
       endif

       ! If using no partial bottom cells, use this:
       if (kmt(iCell).gt.0) then
          bottomDepth(iCell) = refBottomDepth(kmt(iCell))
       endif
       
    enddo



elseif (bottom_topography.eq.'SOMA_Circular_Basin') then

    pi = 4.0*atan(1.0)
    dtr = pi / 180.0
    b = 1.25e6
    phi = 0.1
    H0 = 2500.0
    gamma = -0.4

    kmt = 0

    do iCell=1,nCells
        rlon = lonCell(iCell)
        rlat = latCell(iCell)
        rlonC = 0.0
        rlatC = 35.0*dtr

        r = sqrt( sin(0.5*(rlatC-rlat))**2 +  &
                cos(rlat)*cos(rlatC)*sin(0.5*(rlonC-rlon))**2 )
        distance = 2.*sphere_radius*asin(r)

        r = 1.0 - distance**2 / b**2

        if(r.gt.gamma) then
          zdata = -100.0 - (H0-100.0)/2.0 * (1.0 + tanh(r/phi))
        else
          zdata = 100.0
        endif
        bottomDepth(iCell) = -zdata


        if(zdata.lt.0.0) then
            zdata = -zdata
            r = 0
            kmt_flag=.false.
            do k=1,nVertLevelsMod
                if(.not.kmt_flag) then
                    r = r + dz(k)
                    if(r.gt.zdata) then
                        kmt(iCell) = k
                        kmt_flag = .true.
                    endif
                endif
            enddo
            if(kmt(iCell).eq.0) kmt(iCell)=nVertLevelsMod
        endif

    enddo

elseif (bottom_topography.eq.'flat_bottom') then

   kmt = nVertLevelsMOD
   bottomDepth = refBottomDepth(nVertLevelsMOD)

else 

   print *, ' Incorrect choice of bottom_topography: ',bottom_topography
   stop

endif

if (cut_domain_from_sphere) then
   latmin = -30*dtr
   latmax = +30*dtr
   lonmin = +10*dtr
   lonmax = +70*dtr
   write(6,*) ' lat min ', latmin
   write(6,*) ' lat max ', latmax
   where(latCell.lt.latmin) kmt = 0
   where(latCell.gt.latmax) kmt = 0
   where(lonCell.lt.lonmin) kmt = 0
   where(lonCell.gt.lonmax) kmt = 0
endif

if (solid_boundary_in_y) then
   ymin = minval(yCell)
   write(6,*) ' minimum yCell ', ymin
   ymax = maxval(yCell)
   write(6,*) ' maximum yCell ', ymax
   where(yCell.lt.1.001*ymin) kmt = 0
   where(yCell.gt.0.999*ymax) kmt = 0
endif

if (solid_boundary_in_x) then
   maxdc = maxval(dcEdge)
   xmin = minval(xCell)
   write(6,*) ' minimum xCell ', xmin
   xmax = maxval(xCell)
   write(6,*) ' maximum xCell ', xmax
   where(xCell.lt.xmin+maxdc/1.5) kmt = 0
   where(xCell.gt.xmax-maxdc/1.5) kmt = 0
endif

    
allocate(work_kmt(nCells))
work_kmt = 0.0
where(kmt.eq.0) work_kmt=1.0
write(6,*) 'number of cells culled ',sum(work_kmt)
deallocate(work_kmt)


! Eliminate isolated ocean cells, and make these isolated deep cells
! flush with the deepest neighbor.
do iCell=1,nCells
   kmt_neighbor_max = 0
   do j=1,nEdgesOnCell(iCell)
      iCell1 = cellsOnCell(j,iCell)
      kmt_neighbor_max = max(kmt_neighbor_max,kmt(iCell1))
   enddo
   if (kmt(iCell).gt.kmt_neighbor_max) then
      kmt(iCell) = kmt_neighbor_max
      bottomDepth(iCell) = refBottomDepth(kmt(iCell))
   endif
enddo

if(eliminate_inland_seas) then
call eliminateLoops(nCells,nEdges,nVertices,maxEdges,vertexDegree, &
                    nEdgesOnCell, cellsOnCell, verticesOnEdge, cellsOnVertex, edgesOnCell, lonCell, latCell, &
                    xCell, yCell, zCell, xEdge, yEdge, zEdge, xVertex, yVertex, zVertex, &
                    KMT)
endif

! do not allow land or PBCs in top layers
k = min(top_layers_without_land,nVertLevelsMOD)
where(kmt.gt.0.and.kmt.le.k) 
   bottomDepth = refBottomDepth(k)
   kmt=k
endwhere

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
allocate(maxLevelCellNew(nCellsNew))

allocate(fEdgeNew(nEdgesNew))
allocate(fVertexNew(nVerticesNew))
allocate(bottomDepthNew(nCellsNew))
allocate(normalVelocityForcingNew(nVertLevelsNew,nEdgesNew))
allocate(windStressMonthlyNew(nMonths,nEdgesNew))
allocate(normalVelocityNew(1,nVertLevelsNew,nEdgesNew))
allocate(layerThicknessNew(1,nVertLevelsNew,nCellsNew))
allocate(densityNew(1,nVertLevelsNew,nCellsNew))
allocate(temperatureNew(1,nVertLevelsNew,nCellsNew))
allocate(salinityNew(1,nVertLevelsNew,nCellsNew))
allocate(tracer1New(1,nVertLevelsNew,nCellsNew))

allocate(temperatureRestoreNew(nCellsNew))
allocate(salinityRestoreNew(nCellsNew))

allocate(temperatureRestoreMonthlyNew(nMonths,nCellsNew))
allocate(salinityRestoreMonthlyNew(nMonths,nCellsNew))


xCellNew=0; yCellNew=0; zCellNew=0; latCellNew=0; lonCellNew=0; meshDensityNew=1.0; meshSpacingNew=0.0
xEdgeNew=0; yEdgeNew=0; zEdgeNew=0; latEdgeNew=0; lonEdgeNew=0
xVertexNew=0; yVertexNew=0; zVertexNew=0; latVertexNew=0; lonVertexNew=0

fEdgeNew=0; fVertexNew=0; bottomDepthNew=0; normalVelocityForcingNew=0; windStressMonthlyNew = 0
normalVelocityNew=0; layerThicknessNew=0; densityNew=0
temperatureNew=0; salinityNew=0; tracer1New=0;

temperatureRestoreNew = 0.0
salinityRestoreNew = 0.0

temperatureRestoreMonthlyNew = 0.0
salinityRestoreMonthlyNew = 0.0


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
    maxLevelCellNew(jNew) = kmt(i)
    bottomDepthNew(jNew) = bottomDepth(i)
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
    fEdgeNew(jNew) = fEdge(i)
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
    fVertexNew(jNew)=fVertex(i)
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
!deallocate(bottomDepth)

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
    if(distance.gt.0.5*Lx) then
        write(6,*) ' periodic edge ', iEdge, distance
        write(6,10) '          c1   ', c1(:)
        write(6,10) '          c2   ', c2(:)
        r = c2(1) - c1(1)
        if(r.gt.0) c2(1) = c2(1) - Lx
        if(r.lt.0) c2(1) = c2(1) + Lx
        distance = sqrt( (c1(1)-c2(1))**2 + (c1(2)-c2(2))**2 + (c1(3)-c2(3))**2 )
        write(6,*) ' periodic edge fix ', iEdge, r, distance
    endif
    normalsNew(1,iEdge) = c2(1) - c1(1)
    normalsNew(2,iEdge) = c2(2) - c1(2)
    normalsNew(3,iEdge) = c2(3) - c1(3)
    distance = sqrt( (c1(1)-c2(1))**2 + (c1(2)-c2(2))**2 + (c1(3)-c2(3))**2 )
    normalsNew(:,iEdge) = normalsNew(:,iEdge) / distance
endif
enddo
10 format(a20,3e15.5)

end subroutine map_connectivity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Step 5: Check get_dz routine for hZLevel variable
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_dz
integer k

if (zLevel_thickness.eq.'equally_spaced') then

   write(6,*) ' equally spaced zLevels'
   do i = 1, nVertLevelsMOD
       hZLevel(i) =  layer_thickness_total_max / nVertLevelsMOD
   end do

elseif(zLevel_thickness.eq.'SOMA_40_Level') then

   dz( 1) = 4.6074
   dz( 2) = 5.1112
   dz( 3) = 5.6694
   dz( 4) = 6.2880
   dz( 5) = 6.9732
   dz( 6) = 7.7322
   dz( 7) = 8.5730
   dz( 8) = 9.5041
   dz( 9) = 10.5353
   dz(10) = 11.6770
   dz(11) = 12.9411
   dz(12) = 14.3404
   dz(13) = 15.8894
   dz(14) = 17.6038
   dz(15) = 19.5011
   dz(16) = 21.6005
   dz(17) = 23.9235
   dz(18) = 26.4934
   dz(19) = 29.3363
   dz(20) = 32.4807
   dz(21) = 35.9585
   dz(22) = 39.8044
   dz(23) = 44.0570
   dz(24) = 48.7588
   dz(25) = 53.9568
   dz(26) = 59.7026
   dz(27) = 66.0534
   dz(28) = 73.0722
   dz(29) = 80.8284
   dz(30) = 89.3985
   dz(31) = 98.8672
   dz(32) = 109.3274
   dz(33) = 120.8818
   dz(34) = 133.6436
   dz(35) = 147.7376
   dz(36) = 163.3011
   dz(37) = 180.4858
   dz(38) = 199.4585
   dz(39) = 220.4032
   dz(40) = 243.5224

   hZLevel = dz

elseif(zLevel_thickness.eq.'POP_40_zLevel') then

  dz( 1) = 1001.244   !   5.006218       10.01244
  dz( 2) = 1011.258   !   15.06873       20.12502
  dz( 3) = 1031.682   !   25.28342       30.44183
  dz( 4) = 1063.330   !   35.75848       41.07513
  dz( 5) = 1107.512   !   46.61269       52.15025
  dz( 6) = 1166.145   !   57.98098       63.81171
  dz( 7) = 1241.928   !   70.02135       76.23099
  dz( 8) = 1338.612   !   82.92405       89.61711
  dz( 9) = 1461.401   !   96.92412       104.2311
  dz(10) = 1617.561   !   112.3189       120.4067
  dz(11) = 1817.368   !   129.4936       138.5804
  dz(12) = 2075.558   !   148.9582       159.3360
  dz(13) = 2413.680   !   171.4044       183.4728
  dz(14) = 2863.821   !   197.7919       212.1110
  dz(15) = 3474.644   !   229.4842       246.8575
  dz(16) = 4320.857   !   268.4617       290.0660
  dz(17) = 5516.812   !   317.6501       345.2342
  dz(18) = 7230.458   !   381.3865       417.5388
  dz(19) = 9674.901   !   465.9133       514.2878
  dz(20) = 13003.92   !   579.3074       644.3270
  dz(21) = 17004.89   !   729.3514       814.3759
  dz(22) = 20799.33   !   918.3725       1022.369
  dz(23) = 23356.94   !   1139.154       1255.939
  dz(24) = 24527.19   !   1378.574       1501.210
  dz(25) = 24898.04   !   1625.701       1750.191
  dz(26) = 24983.22   !   1875.107       2000.023
  dz(27) = 24997.87   !   2125.012       2250.002
  dz(28) = 24999.79   !   2375.000       2500.000
  dz(29) = 24999.98   !   2625.000       2749.999
  dz(30) = 25000.00   !   2874.999       2999.999
  dz(31) = 25000.00   !   3124.999       3249.999
  dz(32) = 25000.00   !   3374.999       3499.999
  dz(33) = 25000.00   !   3624.999       3749.999
  dz(34) = 25000.00   !   3874.999       3999.999
  dz(35) = 25000.00   !   4124.999       4249.999
  dz(36) = 25000.00   !   4374.999       4499.999
  dz(37) = 25000.00   !   4624.999       4749.999
  dz(38) = 25000.00   !   4874.999       4999.999
  dz(39) = 25000.00   !   5124.999       5249.999
  dz(40) = 25000.00   !   5374.999       5499.999

  dz = dz / 100.0

  hZLevel = dz

elseif(zLevel_thickness.eq.'isopycnal_3layer') then

  dz( 1) =   500.0
  dz( 2) =  1250.0
  dz( 3) =  3250.0

  hZLevel = dz

else 

   print *, ' Incorrect choice of zLevel_thickness: ',zLevel_thickness
   stop

endif

  refBottomDepth(1) = hZLevel(1)
  do k = 2,nVertLevelsMod
     refBottomDepth(k) = refBottomDepth(k-1) + hZLevel(k)
  end do

  write(6,*) '    k    hZLevel       refBottomDepth'
  do k=1,nVertLevelsMod
    write(6,'(i5,2f10.2)') k,hZLevel(k), refBottomDepth(k)
  enddo
  write(6,*)

end subroutine get_dz

end program map_to_basin
