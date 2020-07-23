program fix_regrid_output

  use netcdf

  implicit none

  integer, parameter :: RKIND = selected_real_kind(13)

  integer :: &
       status, &
       ncid, &
       dimid, &
       varid, &
       nCells, &
       maxEdges, &
       iCell, &
       iCellOnCell, &
       iCellNeighbour, &
       it, &
       nAvg, &
       iCellNearest

  real(kind=RKIND) :: &
       newVal

  integer, parameter :: &
       nMax = 100

  real(kind=RKIND), dimension(:), allocatable :: &
       iceFraction, &
       mpasArrayOut1, &
       mpasArrayOut2, &
       icePresence, &
       latCell, &
       lonCell, &
       icePresenceModify, &
       icePresenceModifyAtlantic, &
       icePresenceModifyPacific, &
       icePresenceModifyNorth, &
       icePresenceModifySouth

  integer, dimension(:), allocatable :: &
       iceFractionMask, &
       nEdgesOnCell, &
       mpasArrayOutMask1, &
       mpasArrayOutMask2

  integer, dimension(:,:), allocatable :: &
       cellsOnCell

  character(len=256) :: &
       filenameIn, &
       filenameMesh, &
       filenameOut

  real(kind=RKIND), parameter :: &
       pi = 3.142_RKIND, &
       degreesToRadians = pi / 180.0_RKIND

  if (command_argument_count() /= 3) then
     write(*,*) "Usage: fix_regrid_output.exe filenameIn filenameMesh filenameOut"
     stop
  endif

  call get_command_argument(1, filenameIn)
  call get_command_argument(2, filenameMesh)
  call get_command_argument(3, filenameOut)

  ! load data
  write(*,*) "Read data file..."
  status = nf90_open(trim(filenameIn), NF90_NOWRITE, ncid)
  call netcdf_error(status, "nf90_open: "//trim(filenameIn))

  ! get dimensions
  status = nf90_inq_dimid(ncid, "nCells", dimid)
  call netcdf_error(status, "nf90_inq_dimid: nCells")

  status = nf90_inquire_dimension(ncid, dimid, len=nCells)
  call netcdf_error(status, "nf90_inquire_dimension: nCells")

  write(*,*) " nCells: ", nCells

  ! get variables
  allocate(iceFraction(nCells))

  status = nf90_inq_varid(ncid, "iceFraction", varid)
  call netcdf_error(status, "nf90_inq_varid: iceFraction")

  status = nf90_get_var(ncid, varid, iceFraction)
  call netcdf_error(status, "nf90_get_var: iceFraction")

  allocate(iceFractionMask(nCells))

  status = nf90_inq_varid(ncid, "iceFractionMask", varid)
  call netcdf_error(status, "nf90_inq_varid: iceFractionMask")

  status = nf90_get_var(ncid, varid, iceFractionMask)
  call netcdf_error(status, "nf90_get_var: iceFractionMask")

  ! close
  status = nf90_close(ncid)
  call netcdf_error(status, "nf90_close")

  ! load mesh
  write(*,*) "Read mesh file..."
  status = nf90_open(trim(filenameMesh), NF90_NOWRITE, ncid)
  call netcdf_error(status, "nf90_open: "//trim(filenameMesh))

  status = nf90_inq_dimid(ncid, "maxEdges", dimid)
  call netcdf_error(status, "nf90_inq_dimid: maxEdges")

  status = nf90_inquire_dimension(ncid, dimid, len=maxEdges)
  call netcdf_error(status, "nf90_inquire_dimension: maxEdges")

  allocate(nEdgesOnCell(nCells))

  status = nf90_inq_varid(ncid, "nEdgesOnCell", varid)
  call netcdf_error(status, "nf90_inq_varid: nEdgesOnCell")

  status = nf90_get_var(ncid, varid, nEdgesOnCell)
  call netcdf_error(status, "nf90_get_var: nEdgesOnCell")

  allocate(cellsOnCell(maxEdges,nCells))

  status = nf90_inq_varid(ncid, "cellsOnCell", varid)
  call netcdf_error(status, "nf90_inq_varid: cellsOnCell")

  status = nf90_get_var(ncid, varid, cellsOnCell)
  call netcdf_error(status, "nf90_get_var: cellsOnCell")

  allocate(latCell(nCells))

  status = nf90_inq_varid(ncid, "latCell", varid)
  call netcdf_error(status, "nf90_inq_varid: latCell")

  status = nf90_get_var(ncid, varid, latCell)
  call netcdf_error(status, "nf90_get_var: latCell")

  allocate(lonCell(nCells))

  status = nf90_inq_varid(ncid, "lonCell", varid)
  call netcdf_error(status, "nf90_inq_varid: lonCell")

  status = nf90_get_var(ncid, varid, lonCell)
  call netcdf_error(status, "nf90_get_var: lonCell")

  status = nf90_close(ncid)
  call netcdf_error(status, "nf90_close")

  ! fill in missing

  allocate(mpasArrayOut1(nCells))
  allocate(mpasArrayOutMask1(nCells))
  allocate(mpasArrayOut2(nCells))
  allocate(mpasArrayOutMask2(nCells))

  mpasArrayOut1 = iceFraction
  mpasArrayOutMask1 = iceFractionMask
  mpasArrayOut2 = iceFraction
  mpasArrayOutMask2 = iceFractionMask

  write(*,*) "fill missing..."
  do it = 1, nMax

     write(*,*) it, " of ", nMax

     do iCell = 1, nCells

        if (mpasArrayOutMask1(iCell) == 0) then

           newVal = 0.0
           nAvg = 0

           do iCellOnCell = 1, nEdgesOnCell(iCell)

              iCellNeighbour = cellsOnCell(iCellOnCell,iCell)

              if (iCellNeighbour > 0 .and. mpasArrayOutMask1(iCellNeighbour) == 1) then

                 newVal = newVal + mpasArrayOut1(iCellNeighbour)
                 nAvg = nAvg + 1

              endif

           enddo

           if (nAvg > 0) then
              mpasArrayOut2(iCell) = newVal / real(nAvg,RKIND)
              mpasArrayOutMask2(iCell) = 1
           endif

        endif

     enddo

     mpasArrayOutMask1 = mpasArrayOutMask2
     mpasArrayOut1     = mpasArrayOut2

  enddo

  allocate(icePresence(nCells))

  do iCell = 1, nCells
     if (mpasArrayOut1(iCell) > 0.0001) then
        icePresence(iCell) = 1
     else
        icePresence(iCell) = 0
     endif
  enddo

  allocate(icePresenceModify(nCells))
  allocate(icePresenceModifyAtlantic(nCells))
  allocate(icePresenceModifyPacific(nCells))
  allocate(icePresenceModifyNorth(nCells))
  allocate(icePresenceModifySouth(nCells))

  ! remove non-connected oceans
  iCellNearest = find_nearest_cell(0.0_RKIND, 0.0_RKIND, latCell, lonCell, nCells, "north")
  call flood_fill(icePresence, icePresenceModifyAtlantic, 2.0_RKIND, nCells, iCellNearest, nEdgesOnCell, cellsOnCell)
  where(latCell < 0.0) icePresenceModifyAtlantic = 2.0
  where(icePresenceModifyAtlantic < 1.5_RKIND) icePresenceModifyAtlantic=1.0_RKIND
  where(icePresenceModifyAtlantic > 1.5_RKIND) icePresenceModifyAtlantic=0.0_RKIND

  iCellNearest = find_nearest_cell(0.0_RKIND, 180.0_RKIND*degreesToRadians, latCell, lonCell, nCells, "north")
  call flood_fill(icePresence, icePresenceModifyPacific, 2.0_RKIND, nCells, iCellNearest, nEdgesOnCell, cellsOnCell)
  where(latCell < 0.0) icePresenceModifyPacific = 2.0
  where(icePresenceModifyPacific < 1.5_RKIND) icePresenceModifyPacific=1.0_RKIND
  where(icePresenceModifyPacific > 1.5_RKIND) icePresenceModifyPacific=0.0_RKIND

  icePresenceModifyNorth = 0.0_RKIND
  where(icePresenceModifyAtlantic == 1.0_RKIND .and. icePresenceModifyPacific == 1.0_RKIND) icePresenceModifyNorth = 1.0_RKIND

  iCellNearest = find_nearest_cell(0.0_RKIND, 0.0_RKIND, latCell, lonCell, nCells, "south")
  call flood_fill(icePresence, icePresenceModifySouth, 2.0_RKIND, nCells, iCellNearest, nEdgesOnCell, cellsOnCell)
  where(latCell > 0.0) icePresenceModifySouth = 2.0
  where(icePresenceModifySouth < 1.5_RKIND) icePresenceModifySouth=1.0_RKIND
  where(icePresenceModifySouth > 1.5_RKIND) icePresenceModifySouth=0.0_RKIND

  do iCell = 1, nCells
     if (icePresenceModifyNorth(iCell) == 1.0_RKIND .or. &
         icePresenceModifySouth(iCell)    == 1.0_RKIND) then
        icePresenceModify(iCell) = 1.0_RKIND
     else
        icePresenceModify(iCell) = 0.0_RKIND
     endif
  enddo ! iCell
  icePresence = icePresenceModify

  ! set contigous sea ice
  !iCellNearest = find_nearest_cell(90.0_RKIND*degreesToRadians, 0.0_RKIND, latCell, lonCell, nCells, "global")
  !call flood_fill(icePresence, icePresenceModify, 3.0_RKIND, nCells, iCellNearest, nEdgesOnCell, cellsOnCell)
  !icePresence = icePresenceModify

  !iCellNearest = find_nearest_cell(-90.0_RKIND*degreesToRadians, 0.0_RKIND, latCell, lonCell, nCells, "global")
  !call flood_fill(icePresence, icePresenceModify, 3.0_RKIND, nCells, iCellNearest, nEdgesOnCell, cellsOnCell)
  !icePresence = icePresenceModify

  !where(icePresence < 2.5_RKIND) icePresence=0.0_RKIND
  !where(icePresence > 0.5_RKIND) icePresence=1.0_RKIND

  ! output
  write(*,*) "Write output..."
  status = nf90_create(trim(filenameOut), NF90_NOCLOBBER, ncid)
  call netcdf_error(status, "nf90_create: "//trim(filenameOut))

  status = nf90_def_dim(ncid, "nCells", nCells, dimid)
  call netcdf_error(status, "nf90_def_dim: nCells")

  status = nf90_def_var(ncid, "icePresence", NF90_DOUBLE, (/dimid/), varid)
  call netcdf_error(status, "nf90_def_var: icePresence")

  status = nf90_enddef(ncid)
  call netcdf_error(status, "nf90_enddef")

  status = nf90_put_var(ncid, varid, icePresence)
  call netcdf_error(status, "nf90_put_var: icePresence")

  status = nf90_close(ncid)
  call netcdf_error(status, "nf90_close")

contains

  !-----------------------------------------------------------------------------

  subroutine netcdf_error(status, message)

    integer, intent(in) :: &
         status

    character(len=*), intent(in) :: &
         message

    if (status /= 0) then

       !write(*,*) "Error: ", nf90_strerror(status)
       write(*,*) trim(message)
       stop

    endif

  end subroutine netcdf_error

  !-----------------------------------------------------------------------------

  subroutine get_position_from_latlon(lat, lon, x, y, z)

    real(kind=RKIND), intent(in) :: &
         lat, &
         lon

    real(kind=RKIND), intent(out) :: &
         x, &
         y, &
         z

    x = cos(lat) * cos(lon)
    y = cos(lat) * sin(lon)
    z = sin(lat)

  end subroutine get_position_from_latlon

  !-----------------------------------------------------------------------------

  function find_nearest_cell(lat, lon, latCell, lonCell, nCells, type) result(iCellNearest)

    real(kind=RKIND), intent(in) :: &
         lat, &
         lon

    real(kind=RKIND), dimension(:), intent(in) :: &
         latCell, &
         lonCell

    integer, intent(in) :: &
         nCells

    character(len=*), intent(in) :: &
         type

    integer :: &
         iCellNearest

    integer :: &
         iCell

    real(kind=RKIND) :: &
         distance, &
         minDistance

    real(kind=RKIND) :: &
         x, y, z, &
         x0, y0, z0

    minDistance = 10.0_RKIND
    iCellNearest = -1

    call get_position_from_latlon(lat, lon, x0, y0, z0)

    if (type == "global") then

       do iCell = 1, nCells

          call get_position_from_latlon(latCell(iCell), lonCell(iCell), x, y, z)

          distance = (x0-x)**2 + (y0-y)**2 + (z0-z)**2

          if (distance < minDistance) then
             minDistance = distance
             iCellNearest = iCell
          endif

       enddo ! iCell

    else if (type == "north") then

       do iCell = 1, nCells

          if (latCell(iCell) > 0.0) then

             call get_position_from_latlon(latCell(iCell), lonCell(iCell), x, y, z)

             distance = (x0-x)**2 + (y0-y)**2 + (z0-z)**2

             if (distance < minDistance) then
                minDistance = distance
                iCellNearest = iCell
             endif

          endif ! north

       enddo ! iCell

    else if (type == "south") then

       do iCell = 1, nCells

          if (latCell(iCell) < 0.0) then

             call get_position_from_latlon(latCell(iCell), lonCell(iCell), x, y, z)

             distance = (x0-x)**2 + (y0-y)**2 + (z0-z)**2

             if (distance < minDistance) then
                minDistance = distance
                iCellNearest = iCell
             endif

          endif ! south

       enddo ! iCell

    endif

  end function find_nearest_cell

  !-----------------------------------------------------------------------------

  subroutine flood_fill(arrayIn, arrayOut, fillValue, nCells, iCellStart, nEdgesOnCell, cellsOnCell)

    real(kind=RKIND), dimension(:), intent(in) :: &
         arrayIn

    real(kind=RKIND), dimension(:), intent(out) :: &
         arrayOut

    real(kind=RKIND), intent(in) :: &
         fillValue

    integer, intent(in) :: &
         nCells, &
         iCellStart

    integer, dimension(:), intent(in) :: &
         nEdgesOnCell

    integer, dimension(:,:), intent(in) :: &
         cellsOnCell

    integer, dimension(:), allocatable :: &
         prevFilliCells, &
         nextFilliCells

    integer :: &
         nPrevFills, &
         nNextFills, &
         iCellFill, &
         iCellPrev, &
         iCellOnCell

    real(kind=RKIND) :: &
         startValue

    startValue = arrayIn(iCellStart)

    arrayOut = arrayIn

    allocate(prevFilliCells(nCells))
    allocate(nextFilliCells(nCells))

    arrayOut(iCellStart) = fillValue
    prevFilliCells(1) = iCellStart
    nPrevFills = 1
    nNextFills = 1

    do while (nNextFills > 0)

       nNextFills = 0

       do iCellFill = 1, nPrevFills

          iCellPrev = prevFilliCells(iCellFill)

          do iCellOnCell = 1, nEdgesOnCell(iCellPrev)

             iCellNeighbour = cellsOnCell(iCellOnCell,iCellPrev)

             if (arrayOut(iCellNeighbour) /= fillValue .and. &
                 arrayOut(iCellNeighbour) == startValue) then

                nNextFills = nNextFills + 1
                nextFilliCells(nNextFills) = iCellNeighbour
                arrayOut(iCellNeighbour) = fillValue

             endif

          enddo

       enddo

       nPrevFills = nNextFills
       prevFilliCells(1:nPrevFills) = nextFilliCells(1:nNextFills)

    enddo

    deallocate(prevFilliCells)

  end subroutine flood_fill

  !-----------------------------------------------------------------------------

end program fix_regrid_output
