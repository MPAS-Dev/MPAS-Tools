
program smooth_topography

! on compy
! module load gcc;module load netcdf
!gfortran smooth_topo_before_init.F90 -o a.out -I /share/apps/netcdf/4.6.3/gcc/4.8.5/include/ -L /share/apps/netcdf/4.6.3/gcc/4.8.5/lib/ -lnetcdf -lnetcdff

  use netcdf

  implicit none

  integer, parameter :: &
       RKIND = selected_real_kind(13)

  ! input parameters
  character(len=1000) :: &
       filename_depth_in, &
       filename_depth_out, &
       filename_mpas_mesh, &
       filename_diagnostic

  real(kind=RKIND) :: &
       layerBottom, dz, &
       layerTop, partialCell_dz,  &
       maxBottomThicknessOrig,  &
       distanceLimit, & ! (km)
       stdDeviation     ! (km)

  integer :: &
       numSmoothingPasses, smoothingPass, minLevelForSmoothing

  namelist /smooth/ &
       filename_depth_in, &
       filename_depth_out, &
       filename_mpas_mesh, &
       distanceLimit, &
       stdDeviation, &
       numSmoothingPasses

  ! mpas grid
  integer :: &
       nCells, &
       nEdges, &
       iCellFull, &
       maxEdges

  real(kind=RKIND), dimension(:), allocatable :: &
       xCell, &
       yCell, &
       zCell, &
       lonCell, &
       latCell, &
       areaCell, &
       bedElevationOrig, &
       bedElevationNew, &
       bedElevationPrev, &
       landIceDraftOrig, &
       landIceDraftNew, &
       landIceDraftPrev, &
       landIceThkOrig, &
       landIceThkNew, &
       landIceThkPrev, &
       dcEdge

  integer, dimension(:), allocatable :: &
       nEdgesOnCell

  integer, dimension(:,:), allocatable :: &
       cellsOnCell, &
       edgesOnCell

  ! temporary variables
  integer :: &
       is, &
       ia, &
       iCell, &
       k, &
       start, &
       count, &
       ncid, &
       varid_BedOrig, &
       varid_BedNew, &
       varid_IceDraftOrig, &
       varid_IceDraftNew, &
       varid_IceThkOrig, &
       varid_IceThkNew

  real(kind=RKIND), dimension(:), allocatable :: &
       weights, &
       cellDistance, &
       distanceSpread

  integer, dimension(:), allocatable :: &
       iCellsSpread, &
       cellStatus

  real(kind=RKIND), parameter :: &
       km_to_m = 1000.0_RKIND

  ! diagnostic variables
  real(kind=RKIND), dimension(:), allocatable :: &
       weightsDiagnostic

  ! read in namelist
  write(*,*) "Reading namelist..."
  open(11,file="smooth_depth_in",status='old')
  read(11,nml=smooth)
  close(11)

  ! load the MPAS mesh file
  call load_mesh_file()

  call load_topo_file()

  ! allocate temporary arrays
  write(*,*) "Allocate temporary variables..."
  allocate(weights(nCells))
  allocate(cellDistance(nCells))
  allocate(distanceSpread(nCells))
  allocate(iCellsSpread(nCells))
  allocate(cellStatus(nCells))

  write(*,*) "Perform smoothing..."

  bedElevationPrev = bedElevationOrig
  bedElevationNew  = bedElevationOrig  !  for cells that do not get smoothed

  landIceDraftPrev = landIceDraftOrig
  landIceDraftNew  = landIceDraftOrig  !  for cells that do not get smoothed

  landIceThkPrev = landIceThkOrig
  landIceThkNew  = landIceThkOrig  !  for cells that do not get smoothed

! convert to degrees
  latCell = latCell*180./3.14159
  lonCell = lonCell*180./3.14159

  do smoothingPass = 1, numSmoothingPasses

    write(*,*)' smoothing pass ',smoothingPass

    ! perform the smoothing
    do iCell = 1, nCells

      call smooth_topo(&
        iCell, &
        distanceLimit * km_to_m, &
        nEdgesOnCell, &
        cellsOnCell, &
        edgesOnCell, &
        dcEdge, &
        cellStatus, &
        cellDistance, &
        distanceSpread)
    enddo ! iCell

    bedElevationPrev = bedElevationNew
    landIceDraftPrev = landIceDraftNew
    landIceThkPrev = landIceThkNew
  enddo ! smoothingPass

  ! output the depth file
  call create_depth_file(filename_depth_out, ncid, &
    varid_BedOrig, varid_BedNew, varid_IceDraftOrig, varid_IceDraftNew,  &
    varid_IceThkOrig, varid_IceThkNew)

  ! close the depth file
  call close_depth_file(ncid)

contains

  !----------------------------------------------------------------

  function gaussian_weight(distance, stdDeviation) result(weight)

    real(kind=RKIND), intent(in) :: &
         distance, &
         stdDeviation

    real(kind=RKIND) :: weight

    weight = exp(-1.0_RKIND * (distance**2 / (2.0_RKIND * stdDeviation**2)))

  end function gaussian_weight

  !----------------------------------------------------------------

  subroutine smooth_topo(&
       iCell_Start, &
       distanceLimit, &
       nEdgesOnCell, &
       cellsOnCell, &
       edgesOnCell, &
       dcEdge, &
       cellStatus, &
       cellDistance, &
       distanceOut)

    integer, intent(in) :: &
         iCell_Start

    real(kind=RKIND), intent(in) :: &
         distanceLimit

    integer, dimension(:), intent(in) :: &
         nEdgesOnCell

    integer, dimension(:,:), intent(in) :: &
         cellsOnCell, &
         edgesOnCell

    real(kind=RKIND), dimension(:), intent(in) :: &
         dcEdge

    integer, dimension(:), intent(out) :: &
         cellStatus

    real(kind=RKIND), dimension(:), intent(out) :: &
         cellDistance, &
         distanceOut

    integer :: &
         nCellsOut

    integer, dimension(:), allocatable :: &
         iCellsPrev, &
         iCellsOut, &
         iCellsNext

    integer :: &
         nCellsPrev, &
         nCellsNext, &
         nLayer, &
         iCellPrev, &
         iCellOnCell, &
         iCellNext, &
         iEdge

    real(kind=RKIND) :: &
         cellDistanceToAdd, &
         weight, &
         smoothedBedSum, &
         smoothedIceDraftSum, &
         smoothedIceThkSum, &
         weightSum

    allocate(iCellsPrev(nCells))
    allocate(iCellsNext(nCells))
    allocate(iCellsOut(nCells))

    ! set the internal variables
    cellStatus   = -1
    cellDistance = -1.0_RKIND

    ! set first cell
    cellStatus(iCell_Start) = 0
    nCellsPrev = 1
    iCellsPrev(1) = iCell_Start

    ! first output cell is the starting cell
    nCellsOut = 1
    iCellsOut(1) = iCell_Start
    distanceOut(1) = 0.0_RKIND

    ! initialize sums
    smoothedBedSum = 0.0_RKIND
    smoothedIceDraftSum = 0.0_RKIND
    smoothedIceThkSum = 0.0_RKIND
    weightSum = 0.0_RKIND

    ! loop over cell layers from original cell
    nLayer = 0
    do while (nCellsPrev > 0)

       ! reset number of cells found this layer
       nCellsNext = 0

       ! loop over cells defined in the previous iteration
       do iCellPrev = 1, nCellsPrev

          ! loop over neighbours of these previous cells
          do iCellOnCell = 1, nEdgesOnCell(iCellsPrev(iCellPrev))

             ! get the iCell of the next potential cell in the next cells
             iCellNext = cellsOnCell(iCellOnCell,iCellsPrev(iCellPrev))

             ! get the edge index for this crossing
             iEdge = edgesOnCell(iCellOnCell,iCellsPrev(iCellPrev))

             cellDistanceToAdd = cellDistance(iCellsPrev(iCellPrev)) + dcEdge(iEdge)

             ! check to see if we need to add it to the next array
             if (cellStatus(iCellNext) == -1 .and. cellDistanceToAdd < distanceLimit) then

                ! count how many on the next list
                nCellsNext = nCellsNext + 1

                ! add this new cell to the next list
                iCellsNext(nCellsNext) = iCellNext

                ! update the status of the cell
                cellStatus(iCellNext) = nLayer

                ! calculate the distance to this cell
                cellDistance(iCellNext) = cellDistanceToAdd

                ! output
                nCellsOut = nCellsOut + 1
                iCellsOut(nCellsOut) = iCellNext
                distanceOut(nCellsOut) = cellDistanceToAdd

                weight = gaussian_weight(cellDistanceToAdd, stdDeviation * km_to_m)
                weightSum = weightSum + weight
                smoothedBedSum = smoothedBedSum + weight*bedElevationPrev(iCellNext)
                smoothedIceDraftSum = smoothedIceDraftSum + weight*landIceDraftPrev(iCellNext)
                smoothedIceThkSum = smoothedIceThkSum + weight*landIceThkPrev(iCellNext)

             endif ! cellStatus(iCellNext) == -1

          enddo ! iCellOnCell

       enddo ! iCellPrev

       ! swap next and prev
       nCellsPrev = nCellsNext

       iCellsPrev(1:nCellsNext) = iCellsnext(1:nCellsNext)

       ! increment the layer number
       nLayer = nLayer + 1

    enddo ! nCellsNext > 0

    bedElevationNew(iCell_Start) = smoothedBedSum/weightSum
    landIceDraftNew(iCell_Start) = smoothedIceDraftSum/weightSum
    landIceThkNew(iCell_Start) = smoothedIceThkSum/weightSum

    deallocate(iCellsPrev)
    deallocate(iCellsNext)

  end subroutine smooth_topo

  !----------------------------------------------------------------

  subroutine create_depth_file(filename, ncid, &
    varid_BedOrig, varid_BedNew, varid_IceDraftOrig, varid_IceDraftNew,  &
    varid_IceThkOrig, varid_IceThkNew)

    character(len=*), intent(in) :: &
         filename

    integer, intent(out) :: &
         ncid, &
         varid_BedOrig, &
         varid_BedNew, &
         varid_IceDraftOrig, &
         varid_IceDraftNew,  &
         varid_IceThkOrig, &
         varid_IceThkNew

    integer :: &
         status

    integer :: &
         dimid_nCells, &
         dimid_ni_a, &
         dimid_nj_a, &
         dimid_nv_a, &
         dimid_src_grid_rank, &
         dimid_n_b, &
         dimid_ni_b, &
         dimid_nj_b, &
         dimid_nv_b, &
         dimid_dst_grid_rank, &
         dimid_n_s

    integer :: &
         varid_xc_a, &
         varid_yc_a, &
         varid_xv_a, &
         varid_yv_a, &
         varid_mask_a, &
         varid_area_a, &
         varid_frac_a, &
         varid_src_grid_dims, &
         varid_xc_b, &
         varid_yc_b, &
         varid_xv_b, &
         varid_yv_b, &
         varid_mask_b, &
         varid_area_b, &
         varid_frac_b, &
         varid_dst_grid_dims

    character(len=1000) :: &
         date, &
         time, &
         zone, &
         datetime

    integer, dimension(1) :: &
         start1D, &
         count1D

    integer, dimension(2) :: &
         start2D, &
         count2D

    write(*,*) "Create depth file...", trim(filename)

    ! create
!   status = nf90_create(trim(filename), NF90_CLOBBER, ncid)
    status = nf90_create(trim(filename), NF90_64BIT_OFFSET, ncid)
    call netcdf_error(status, "create_depth_file: nf90_open")

    ! define dimensions
    ! nCells
    status = nf90_def_dim(ncid, "nCells", nCells, dimid_nCells)
    call netcdf_error(status, "create_depth_file: nf90_def_dim nCells")

    ! define variables
    ! bedElevationOrig
    status = nf90_def_var(ncid, "bed_elevationOrig", NF90_DOUBLE, (/dimid_nCells/), varid_BedOrig)
    call netcdf_error(status, "create_depth_file: nf90_def_var bedElevationOrig")

    status = nf90_put_att(ncid, varid_BedOrig, "long_name", &
     "Elevation of the bottom of the ocean. Given as a negative distance from sea level.")
    call netcdf_error(status, "create_depth_file: nf90_put_att bedElevationOrig")

    status = nf90_put_att(ncid, varid_BedOrig, "units", "m")
    call netcdf_error(status, "create_depth_file: nf90_put_att bedElevationOrig")

    ! bedElevationNew
    status = nf90_def_var(ncid, "bed_elevationNew", NF90_DOUBLE, (/dimid_nCells/), varid_BedNew)
    call netcdf_error(status, "create_depth_file: nf90_def_var bedElevationNew")

    status = nf90_put_att(ncid, varid_BedNew, "long_name", &
     "Elevation of the bottom of the ocean. Given as a negative distance from sea level.")
    call netcdf_error(status, "create_depth_file: nf90_put_att bedElevationNew")

    status = nf90_put_att(ncid, varid_BedNew, "units", "m")
    call netcdf_error(status, "create_depth_file: nf90_put_att bedElevationNew")

    ! landIceDraftOrig
    status = nf90_def_var(ncid, "landIceDraftObservedOrig", NF90_DOUBLE, (/dimid_nCells/), varid_IceDraftOrig)
    call netcdf_error(status, "create_depth_file: nf90_def_var landIceDraftOrig")

    status = nf90_put_att(ncid, varid_IceDraftOrig, "long_name", &
     "The elevation of the bottom of land ice.")
    call netcdf_error(status, "create_depth_file: nf90_put_att landIceDraftOrig")

    status = nf90_put_att(ncid, varid_IceDraftOrig, "units", "m")
    call netcdf_error(status, "create_depth_file: nf90_put_att landIceDraftOrig")

    ! landIceDraftNew
    status = nf90_def_var(ncid, "landIceDraftObservedNew", NF90_DOUBLE, (/dimid_nCells/), varid_IceDraftNew)
    call netcdf_error(status, "create_depth_file: nf90_def_var landIceDraftNew")

    status = nf90_put_att(ncid, varid_IceDraftNew, "long_name", &
     "The elevation of the bottom of land ice.")
    call netcdf_error(status, "create_depth_file: nf90_put_att landIceDraftNew")

    status = nf90_put_att(ncid, varid_IceDraftNew, "units", "m")
    call netcdf_error(status, "create_depth_file: nf90_put_att landIceDraftNew")

    ! landIceThkOrig
    status = nf90_def_var(ncid, "landIceThkObservedOrig", NF90_DOUBLE, (/dimid_nCells/), varid_IceThkOrig)
    call netcdf_error(status, "create_depth_file: nf90_def_var landIceThkOrig")

    status = nf90_put_att(ncid, varid_IceThkOrig, "long_name", &
     "The thickness of the land ice.")
    call netcdf_error(status, "create_depth_file: nf90_put_att landIceThkOrig")

    status = nf90_put_att(ncid, varid_IceThkOrig, "units", "m")
    call netcdf_error(status, "create_depth_file: nf90_put_att landIceThkOrig")

    ! landIceThkNew
    status = nf90_def_var(ncid, "landIceThkObservedNew", NF90_DOUBLE, (/dimid_nCells/), varid_IceThkNew)
    call netcdf_error(status, "create_depth_file: nf90_def_var landIceThkNew")

    status = nf90_put_att(ncid, varid_IceThkNew, "long_name", &
     "The thickness of the land ice.")
    call netcdf_error(status, "create_depth_file: nf90_put_att landIceThkNew")

    status = nf90_put_att(ncid, varid_IceThkNew, "units", "m")
    call netcdf_error(status, "create_depth_file: nf90_put_att landIceThkNew")

    ! global attributes
    status = nf90_put_att(ncid, NF90_GLOBAL, "input_depth_file", trim(filename_depth_in))
    call netcdf_error(status, "create_depth_file: nf90_put_att input_depth_file")

    status = nf90_put_att(ncid, NF90_GLOBAL, "input_mesh_file", trim(filename_mpas_mesh))
    call netcdf_error(status, "create_depth_file: nf90_put_att input_mesh_file")

    status = nf90_put_att(ncid, NF90_GLOBAL, "smoothing_method", "2D Gaussian smoothing")
    call netcdf_error(status, "create_depth_file: nf90_put_att smoothing_method")

    status = nf90_put_att(ncid, NF90_GLOBAL, "created_by", "smooth_topo_before_init")
    call netcdf_error(status, "create_depth_file: nf90_put_att created_by")

    call date_and_time(date, time, zone)
    datetime = date(1:4)//"-"//date(5:6)//"-"//date(7:8)//"_"//time(1:2)//":"//time(3:4)//":"//time(5:6)//" "//trim(zone)
    status = nf90_put_att(ncid, NF90_GLOBAL, "created_at", trim(datetime))
    call netcdf_error(status, "create_depth_file: nf90_put_att created_at")

    status = nf90_put_att(ncid, NF90_GLOBAL, "distanceLimit", distanceLimit)
    call netcdf_error(status, "create_depth_file: nf90_put_att distanceLimit")

    status = nf90_put_att(ncid, NF90_GLOBAL, "stdDeviation", stdDeviation)
    call netcdf_error(status, "create_depth_file: nf90_put_att stdDeviation")

    ! end definition phase
    status = nf90_enddef(ncid)
    call netcdf_error(status, "create_depth_file: nf90_enddef")

    ! write variables
    start1D(1) = 1
    count1D(1) = nCells

    ! bedElevationOrig
    status = nf90_put_var(ncid, varid_BedOrig, bedElevationOrig, start1D, count1D)
    call netcdf_error(status, "create_depth_file: nf90_put_var bedElevationOrig")

    ! bedElevationNew
    status = nf90_put_var(ncid, varid_BedNew, bedElevationNew, start1D, count1D)
    call netcdf_error(status, "create_depth_file: nf90_put_var bedElevationNew")

    ! landIceDraftOrig
    status = nf90_put_var(ncid, varid_IceDraftOrig, landIceDraftOrig, start1D, count1D)
    call netcdf_error(status, "create_depth_file: nf90_put_var landIceDraftOrig")

    ! landIceDraftNew
    status = nf90_put_var(ncid, varid_IceDraftNew, landIceDraftNew, start1D, count1D)
    call netcdf_error(status, "create_depth_file: nf90_put_var landIceDraftNew")

    ! landIceThkOrig
    status = nf90_put_var(ncid, varid_IceThkOrig, landIceThkOrig, start1D, count1D)
    call netcdf_error(status, "create_depth_file: nf90_put_var landIceThkOrig")

    ! landIceThkNew
    status = nf90_put_var(ncid, varid_IceThkNew, landIceThkNew, start1D, count1D)
    call netcdf_error(status, "create_depth_file: nf90_put_var landIceThkNew")

  end subroutine create_depth_file

  !----------------------------------------------------------------

  subroutine close_depth_file(ncid)

    integer, intent(in) :: &
         ncid

    integer :: &
         status

    ! close
    status = nf90_close(ncid)
    call netcdf_error(status, "close_depth_file: nf90_close ")

  end subroutine close_depth_file

  !----------------------------------------------------------------

  subroutine load_mesh_file()

    integer :: &
         status, &
         ncid, &
         dimid, &
         varid

    write(*,*) "Load mesh file..."

    ! open file
    status = nf90_open(trim(filename_mpas_mesh), NF90_NOWRITE, ncid)
    call netcdf_error(status, "load_mesh_file: nf90_open")

    ! nCells
    status = nf90_inq_dimid(ncid, "nCells", dimid)
    call netcdf_error(status, "load_mesh_file: nf90_inq_dimid nCells")

    status = nf90_inquire_dimension(ncid, dimid, len=nCells)
    call netcdf_error(status, "load_mesh_file: nf90_inquire_dimension nCells")

    ! nEdges
    status = nf90_inq_dimid(ncid, "nEdges", dimid)
    call netcdf_error(status, "load_mesh_file: nf90_inq_dimid nEdges")

    status = nf90_inquire_dimension(ncid, dimid, len=nEdges)
    call netcdf_error(status, "load_mesh_file: nf90_inquire_dimension nEdges")

    ! maxEdges
    status = nf90_inq_dimid(ncid, "maxEdges", dimid)
    call netcdf_error(status, "load_mesh_file: nf90_inq_dimid maxEdges")

    status = nf90_inquire_dimension(ncid, dimid, len=maxEdges)
    call netcdf_error(status, "load_mesh_file: nf90_inquire_dimension maxEdges")

    allocate(xCell(nCells))
    allocate(yCell(nCells))
    allocate(zCell(nCells))
    allocate(latCell(nCells))
    allocate(lonCell(nCells))
    allocate(areaCell(nCells))
    allocate(dcEdge(nEdges))
    allocate(nEdgesOnCell(nCells))
    allocate(cellsOnCell(maxEdges,nCells))
    allocate(edgesOnCell(maxEdges,nCells))

    ! xCell
    status = nf90_inq_varid(ncid, "xCell", varid)
    call netcdf_error(status, "load_mesh_file: nf90_inq_varid xCell")

    status = nf90_get_var(ncid, varid, xCell)
    call netcdf_error(status, "load_mesh_file: nf90_get_var xCell")

    ! yCell
    status = nf90_inq_varid(ncid, "yCell", varid)
    call netcdf_error(status, "load_mesh_file: nf90_inq_varid yCell")

    status = nf90_get_var(ncid, varid, yCell)
    call netcdf_error(status, "load_mesh_file: nf90_get_var yCell")

    ! zCell
    status = nf90_inq_varid(ncid, "zCell", varid)
    call netcdf_error(status, "load_mesh_file: nf90_inq_varid zCell")

    status = nf90_get_var(ncid, varid, zCell)
    call netcdf_error(status, "load_mesh_file: nf90_get_var zCell")

    ! latCell
    status = nf90_inq_varid(ncid, "latCell", varid)
    call netcdf_error(status, "load_mesh_file: nf90_inq_varid latCell")

    status = nf90_get_var(ncid, varid, latCell)
    call netcdf_error(status, "load_mesh_file: nf90_get_var latCell")

    ! lonCell
    status = nf90_inq_varid(ncid, "lonCell", varid)
    call netcdf_error(status, "load_mesh_file: nf90_inq_varid lonCell")

    status = nf90_get_var(ncid, varid, lonCell)
    call netcdf_error(status, "load_mesh_file: nf90_get_var lonCell")

    ! areaCell
    status = nf90_inq_varid(ncid, "areaCell", varid)
    call netcdf_error(status, "load_mesh_file: nf90_inq_varid areaCell")

    status = nf90_get_var(ncid, varid, areaCell)
    call netcdf_error(status, "load_mesh_file: nf90_get_var areaCell")

    ! dcEdge
    status = nf90_inq_varid(ncid, "dcEdge", varid)
    call netcdf_error(status, "load_mesh_file: nf90_inq_varid dcEdge")

    status = nf90_get_var(ncid, varid, dcEdge)
    call netcdf_error(status, "load_mesh_file: nf90_get_var dcEdge")

    ! nEdgesOnCell
    status = nf90_inq_varid(ncid, "nEdgesOnCell", varid)
    call netcdf_error(status, "load_mesh_file: nf90_inq_varid nEdgesOnCell")

    status = nf90_get_var(ncid, varid, nEdgesOnCell)
    call netcdf_error(status, "load_mesh_file: nf90_get_var nEdgesOnCell")

    ! cellsOnCell
    status = nf90_inq_varid(ncid, "cellsOnCell", varid)
    call netcdf_error(status, "load_mesh_file: nf90_inq_varid cellsOnCell")

    status = nf90_get_var(ncid, varid, cellsOnCell)
    call netcdf_error(status, "load_mesh_file: nf90_get_var cellsOnCell")

    ! edgesOnCell
    status = nf90_inq_varid(ncid, "edgesOnCell", varid)
    call netcdf_error(status, "load_mesh_file: nf90_inq_varid edgesOnCell")

    status = nf90_get_var(ncid, varid, edgesOnCell)
    call netcdf_error(status, "load_mesh_file: nf90_get_var edgesOnCell")

    ! close
    status = nf90_close(ncid)
    call netcdf_error(status, "load_mesh_file: nf90_close")

  end subroutine load_mesh_file

  !----------------------------------------------------------------

  subroutine load_topo_file()

    integer :: &
         status, &
         ncid, &
         dimid, &
         varid

    write(*,*) "Load topo file..."

    ! open file
    status = nf90_open(trim(filename_depth_in), NF90_NOWRITE, ncid)
    call netcdf_error(status, "load_topo_file: nf90_open")

    allocate(bedElevationOrig(nCells))
    allocate(bedElevationNew(nCells))
    allocate(bedElevationPrev(nCells))
    allocate(landIceDraftOrig(nCells))
    allocate(landIceDraftNew(nCells))
    allocate(landIceDraftPrev(nCells))
    allocate(landIceThkOrig(nCells))
    allocate(landIceThkNew(nCells))
    allocate(landIceThkPrev(nCells))

    ! bedElevation
    status = nf90_inq_varid(ncid, "bed_elevation", varid)
    call netcdf_error(status, "load_topo_file: nf90_inq_varid bed_elevation")

    status = nf90_get_var(ncid, varid, bedElevationOrig)
    call netcdf_error(status, "load_topo_file: nf90_get_var bed_elevation")

    ! landIceDraft
    status = nf90_inq_varid(ncid, "landIceDraftObserved", varid)
    call netcdf_error(status, "load_topo_file: nf90_inq_varid landIceDraft")

    status = nf90_get_var(ncid, varid, landIceDraftOrig)
    call netcdf_error(status, "load_topo_file: nf90_get_var landIceDraft")

    ! landIceThk
    status = nf90_inq_varid(ncid, "landIceThkObserved", varid)
    call netcdf_error(status, "load_topo_file: nf90_inq_varid landIceThk")

    status = nf90_get_var(ncid, varid, landIceThkOrig)
    call netcdf_error(status, "load_topo_file: nf90_get_var landIceThk")

    ! close
    status = nf90_close(ncid)
    call netcdf_error(status, "load_topo_file: nf90_close")

  end subroutine load_topo_file

  !----------------------------------------------------------------

  subroutine netcdf_error(status, message)

    integer, intent(in) :: &
         status

    character(len=*), intent(in) :: &
         message

    if (status /= 0) then
       write(*,*) "Netcdf error: ", status, nf90_strerror(status)
       write(*,*) trim(message)
       stop
    endif

  end subroutine netcdf_error

  !----------------------------------------------------------------

end program smooth_topography
