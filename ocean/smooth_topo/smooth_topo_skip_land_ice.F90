
program smooth_topography

! on compy
! module load gcc;module load netcdf
!gfortran smooth_topo.F90 -o a.out -I /share/apps/netcdf/4.6.3/gcc/4.8.5/include/ -L /share/apps/netcdf/4.6.3/gcc/4.8.5/lib/ -lnetcdf -lnetcdff

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
       numSmoothingPasses,  &
       minLevelForSmoothing

  real(kind=RKIND), dimension(:,:), allocatable :: &
       restingThicknessOrig,  &
       restingThicknessNew

  ! mpas grid
  integer :: &
       nCells, &
       nEdges, &
       nVertLevels, &
       iCellFull, &
       maxEdges

  real(kind=RKIND), dimension(:), allocatable :: &
       xCell, &
       yCell, &
       zCell, &
       lonCell, &
       latCell, &
       areaCell, &
       bottomDepthOrig, &
       bottomDepthNew, &
       bottomDepthPrev, &
       dcEdge

  integer, dimension(:), allocatable :: &
       landIceMask,  &
       maxLevelCellOrig, &
       maxLevelCellNew, &
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
       varid_landIceMask, &
       varid_DepthOrig, &
       varid_DepthNew, &
       varid_maxLevCellOrig, &
       varid_maxLevCellNew,  &
       varid_RestThickOrig, &
       varid_RestThickNew
  
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

  ! allocate temporary arrays
  write(*,*) "Allocate temporary variables..."
  allocate(weights(nCells))
  allocate(cellDistance(nCells))
  allocate(distanceSpread(nCells))
  allocate(iCellsSpread(nCells))
  allocate(cellStatus(nCells))

  write(*,*) "Perform smoothing..."

  bottomDepthPrev = bottomDepthOrig
  bottomDepthNew  = bottomDepthOrig  !  for cells that do not get smoothed

! convert to degrees
  latCell = latCell*180./3.14159
  lonCell = lonCell*180./3.14159

  do smoothingPass = 1, numSmoothingPasses

  write(*,*)' smoothing pass ',smoothingPass

  ! perform the smoothing
    do iCell = 1, nCells

! only smooth sufficiently deep cells
      if (maxLevelCellOrig(iCell) >= minLevelForSmoothing) then
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
      endif

    enddo ! iCell

    bottomDepthPrev = bottomDepthNew
  enddo ! smoothingPass

  ! find maximum thickness in bottom layer
  maxBottomThicknessOrig = -1.0
  do iCell = 1, nCells
    if (restingThicknessOrig(nVertLevels,iCell) > maxBottomThicknessOrig)  &
         maxBottomThicknessOrig = restingThicknessOrig(nVertLevels,iCell)
  enddo

  ! find first column that extends to the bottom
  do iCell = 1, nCells
    if (maxLevelCellOrig(iCell) == nVertLevels .and.  &
         restingThicknessOrig(nVertLevels,iCell) == maxBottomThicknessOrig) then
      iCellFull = iCell
      exit
    endif
  enddo

  ! find new maxLevelCell and restingThickness
  do iCell = 1, nCells
   if (landIceMask(iCell) == 0) then
    maxLevelCellNew(iCell) = nVertLevels ! some deep levels get set to zero due to roundoff
    layerBottom = 0.0_RKIND
    layerTop    = 0.0_RKIND
    do k = 1, nVertLevels
      restingThicknessNew(k,iCell) = 0.0_RKIND
      layerBottom = layerBottom + restingThicknessOrig(k,iCellFull)
      dz = bottomDepthNew(iCell) - layerBottom
      restingThicknessNew(k,iCell) = restingThicknessOrig(k,iCellFull)

      if (dz < 0.0_RKIND) then
        partialCell_dz = bottomDepthNew(iCell) - layerTop
        if ( partialCell_dz < 0.25*restingThicknessOrig(k,iCellFull) ) then
          bottomDepthNew(iCell) = layerTop + 0.25*restingThicknessOrig(k,iCellFull)
          restingThicknessNew(k,iCell) = 0.25*restingThicknessOrig(k,iCellFull)
        else
          restingThicknessNew(k,iCell) = abs(partialCell_dz)  ! abs for roundoff
        endif 
        maxLevelCellNew(iCell) = k

        exit
      endif 
      layerTop = layerBottom
    enddo ! k

   else  !  landIceMask = 1
    bottomDepthNew(iCell) = bottomDepthOrig(iCell)
    restingThicknessNew(:,iCell) = restingThicknessOrig(:,iCell)
    maxLevelCellNew(iCell) = maxLevelCellOrig(iCell)
   endif

  enddo ! iCell

  ! output the depth file
  call create_depth_file(filename_depth_out, ncid, &
    varid_DepthOrig, varid_DepthNew, varid_maxLevCellOrig, varid_maxLevCellNew,  &
    varid_RestThickOrig, varid_RestThickNew)

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
         smoothedTopoSum, &
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
    smoothedTopoSum = 0.0_RKIND
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
             if (landIceMask(iCellNext) == 0 .and. cellStatus(iCellNext) == -1 .and. cellDistanceToAdd < distanceLimit) then

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
                smoothedTopoSum = smoothedTopoSum + weight*bottomDepthPrev(iCellNext)

             endif ! cellStatus(iCellNext) == -1

          enddo ! iCellOnCell

       enddo ! iCellPrev

       ! swap next and prev
       nCellsPrev = nCellsNext

       iCellsPrev(1:nCellsNext) = iCellsnext(1:nCellsNext)

       ! increment the layer number
       nLayer = nLayer + 1

    enddo ! nCellsNext > 0

    bottomDepthNew(iCell_Start) = smoothedTopoSum/weightSum

    deallocate(iCellsPrev)
    deallocate(iCellsNext)

  end subroutine smooth_topo

  !----------------------------------------------------------------

  subroutine create_depth_file(filename, ncid, &
    varid_DepthOrig, varid_DepthNew, varid_maxLevCellOrig, varid_maxLevCellNew,  &
    varid_RestThickOrig, varid_RestThickNew)

    character(len=*), intent(in) :: &
         filename

    integer, intent(out) :: &
         ncid, &
         varid_DepthOrig, &
         varid_DepthNew, &
         varid_maxLevCellOrig, &
         varid_maxLevCellNew,  &
         varid_RestThickOrig, &
         varid_RestThickNew
    
    integer :: &
         status

    integer :: &
         dimid_nCells, &
         dimid_nVertLevels, &
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

    ! nVertLevels
    status = nf90_def_dim(ncid, "nVertLevels", nVertLevels, dimid_nVertLevels)
    call netcdf_error(status, "create_depth_file: nf90_def_dim nVertLevels")

    ! define variables
    ! bottomDepthOrig
    status = nf90_def_var(ncid, "bottomDepthOrig", NF90_DOUBLE, (/dimid_nCells/), varid_DepthOrig)
    call netcdf_error(status, "create_depth_file: nf90_def_var bottomDepthOrig")

    status = nf90_put_att(ncid, varid_DepthOrig, "long_name", &
     "Depth of the bottom of the ocean. Given as a positive distance from sea level.")
    call netcdf_error(status, "create_depth_file: nf90_put_att DepthOrig")

    status = nf90_put_att(ncid, varid_DepthOrig, "units", "m")
    call netcdf_error(status, "create_depth_file: nf90_put_att DepthOrig")

    ! bottomDepthNew
    status = nf90_def_var(ncid, "bottomDepthNew", NF90_DOUBLE, (/dimid_nCells/), varid_DepthNew)
    call netcdf_error(status, "create_depth_file: nf90_def_var bottomDepthNew")

    status = nf90_put_att(ncid, varid_DepthNew, "long_name", &
     "Depth of the bottom of the ocean. Given as a positive distance from sea level.")
    call netcdf_error(status, "create_depth_file: nf90_put_att DepthNew")

    status = nf90_put_att(ncid, varid_DepthNew, "units", "m")
    call netcdf_error(status, "create_depth_file: nf90_put_att DepthNew")

    ! maxLevelCellOrig
    status = nf90_def_var(ncid, "maxLevelCellOrig", NF90_INT, (/dimid_nCells/), varid_maxLevCellOrig)
    call netcdf_error(status, "create_depth_file: nf90_def_var maxLevelCellOrig")

    status = nf90_put_att(ncid, varid_maxLevCellOrig, "long_name", &
     "Index to the last active ocean cell in each column.")
    call netcdf_error(status, "create_depth_file: nf90_put_att maxLevelCellOrig")

    status = nf90_put_att(ncid, varid_maxLevCellOrig, "units", "unitless")
    call netcdf_error(status, "create_depth_file: nf90_put_att maxLevelCellOrig")

    ! maxLevelCellNew
    status = nf90_def_var(ncid, "maxLevelCellNew", NF90_INT, (/dimid_nCells/), varid_maxLevCellNew)
    call netcdf_error(status, "create_depth_file: nf90_def_var maxLevelCellNew")

    status = nf90_put_att(ncid, varid_maxLevCellNew, "long_name", &
     "Index to the last active ocean cell in each column.")
    call netcdf_error(status, "create_depth_file: nf90_put_att maxLevelCellNew")

    status = nf90_put_att(ncid, varid_maxLevCellNew, "units", "unitless")
    call netcdf_error(status, "create_depth_file: nf90_put_att maxLevelCellOrig")

    ! restingThicknessOrig
    status = nf90_def_var(ncid, "restingThicknessOrig", NF90_DOUBLE, (/dimid_nVertLevels,dimid_nCells/), varid_RestThickOrig)
    call netcdf_error(status, "create_depth_file: nf90_def_var restingThicknessOrig")

    status = nf90_put_att(ncid, varid_RestThickOrig, "long_name", &
     "Layer thickness when the ocean is at rest, i.e. without SSH or internal perturbations.")
    call netcdf_error(status, "create_depth_file: nf90_put_att RestThickOrig")

    status = nf90_put_att(ncid, varid_RestThickOrig, "units", "m")
    call netcdf_error(status, "create_depth_file: nf90_put_att RestThickOrig")

    ! restingThicknessNew
    status = nf90_def_var(ncid, "restingThicknessNew", NF90_DOUBLE, (/dimid_nVertLevels,dimid_nCells/), varid_RestThickNew)
    call netcdf_error(status, "create_depth_file: nf90_def_var restingThicknessNew")

    status = nf90_put_att(ncid, varid_RestThickNew, "long_name", &
     "Layer thickness when the ocean is at rest, i.e. without SSH or internal perturbations.")
    call netcdf_error(status, "create_depth_file: nf90_put_att RestThickNew")

    status = nf90_put_att(ncid, varid_RestThickNew, "units", "m")
    call netcdf_error(status, "create_depth_file: nf90_put_att RestThickNew")

    ! global attributes
    status = nf90_put_att(ncid, NF90_GLOBAL, "input_depth_file", trim(filename_depth_in))
    call netcdf_error(status, "create_depth_file: nf90_put_att input_depth_file")

    status = nf90_put_att(ncid, NF90_GLOBAL, "input_mesh_file", trim(filename_mpas_mesh))
    call netcdf_error(status, "create_depth_file: nf90_put_att input_mesh_file")

    status = nf90_put_att(ncid, NF90_GLOBAL, "smoothing_method", "2D Gaussian smoothing")
    call netcdf_error(status, "create_depth_file: nf90_put_att smoothing_method")

    status = nf90_put_att(ncid, NF90_GLOBAL, "created_by", "smooth_runoff.exe")
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

    ! bottomDepthOrig
    status = nf90_put_var(ncid, varid_DepthOrig, bottomDepthOrig, start1D, count1D)
    call netcdf_error(status, "create_depth_file: nf90_put_var DepthOrig")
    
    ! bottomDepthNew
    status = nf90_put_var(ncid, varid_DepthNew, bottomDepthNew, start1D, count1D)
    call netcdf_error(status, "create_depth_file: nf90_put_var DepthNew")
    
    ! maxLevelCellOrig
    status = nf90_put_var(ncid, varid_maxLevCellOrig, maxLevelCellOrig, start1D, count1D)
    call netcdf_error(status, "create_depth_file: nf90_put_var maxLevelCellOrig")
    
    ! maxLevelCellNew
    status = nf90_put_var(ncid, varid_maxLevCellNew, maxLevelCellNew, start1D, count1D)
    call netcdf_error(status, "create_depth_file: nf90_put_var maxLevelCellNew")
    
    start2D(1) = 1
    start2D(2) = 1
    count2D(2) = nCells
    count2D(1) = nVertLevels

    ! restingThicknessOrig
    status = nf90_put_var(ncid, varid_RestThickOrig, restingThicknessOrig, start2D, count2D)
    call netcdf_error(status, "create_depth_file: nf90_put_var restingThicknessOrig")
    
    ! restingThicknessNew
    status = nf90_put_var(ncid, varid_RestThickNew, restingThicknessNew, start2D, count2D)
    call netcdf_error(status, "create_depth_file: nf90_put_var restingThicknessNew")
    
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

    ! nVertLevels
    status = nf90_inq_dimid(ncid, "nVertLevels", dimid)
    call netcdf_error(status, "load_mesh_file: nf90_inq_dimid nVertLevels")
    
    status = nf90_inquire_dimension(ncid, dimid, len=nVertLevels)
    call netcdf_error(status, "load_mesh_file: nf90_inquire_dimension nVertLevels")

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
    allocate(bottomDepthOrig(nCells))
    allocate(bottomDepthNew(nCells))
    allocate(bottomDepthPrev(nCells))
    allocate(maxLevelCellOrig(nCells))
    allocate(maxLevelCellNew(nCells))
    allocate(restingThicknessOrig(nVertLevels,nCells))
    allocate(restingThicknessNew(nVertLevels,nCells))
    allocate(landIceMask(nCells))

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

    ! bottomDepth
    status = nf90_inq_varid(ncid, "bottomDepth", varid)
    call netcdf_error(status, "load_mesh_file: nf90_inq_varid bottomDepth")
    
    status = nf90_get_var(ncid, varid, bottomDepthOrig)
    call netcdf_error(status, "load_mesh_file: nf90_get_var bottomDepth")

    ! maxLevelCell
    status = nf90_inq_varid(ncid, "maxLevelCell", varid)
    call netcdf_error(status, "load_mesh_file: nf90_inq_varid maxLevelCell")
    
    status = nf90_get_var(ncid, varid, maxLevelCellOrig)
    call netcdf_error(status, "load_mesh_file: nf90_get_var maxLevelCell")

    ! restingThickness
    status = nf90_inq_varid(ncid, "restingThickness", varid)
    call netcdf_error(status, "load_mesh_file: nf90_inq_varid restingThickness")
    
    status = nf90_get_var(ncid, varid, restingThicknessOrig)
    call netcdf_error(status, "load_mesh_file: nf90_get_var restingThickness")

    ! landIceMask
    status = nf90_inq_varid(ncid, "landIceMask", varid)
    call netcdf_error(status, "load_mesh_file: nf90_inq_varid landIceMask")
    
    status = nf90_get_var(ncid, varid, landIceMask)
    call netcdf_error(status, "load_mesh_file: nf90_get_var landIceMask")

    ! close
    status = nf90_close(ncid)
    call netcdf_error(status, "load_mesh_file: nf90_close")
    
  end subroutine load_mesh_file
  
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
