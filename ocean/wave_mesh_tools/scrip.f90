PROGRAM scrip

     USE NETCDF
     USE read_write_gmsh
     USE WMSCRPMD
     
     CHARACTER(100) :: waves_mesh_file
     CHARACTER(100) :: waves_scrip_file
     TYPE(grid_type) :: mesh
     REAL*8, DIMENSION(:,:), ALLOCATABLE :: xyb
     INTEGER, DIMENSION(:,:), ALLOCATABLE :: trigp
     REAL*8, DIMENSION(:), ALLOCATABLE :: grid_center_lon
     REAL*8, DIMENSION(:), ALLOCATABLE :: grid_center_lat
     REAL*8, DIMENSION(:,:), ALLOCATABLE :: grid_corner_lon
     REAL*8, DIMENSION(:,:), ALLOCATABLE :: grid_corner_lat
     LOGICAL, DIMENSION(:), ALLOCATABLE :: grid_mask
     INTEGER, DIMENSION(:), ALLOCATABLE :: grid_imask
     INTEGER, DIMENSION(:), ALLOCATABLE :: grid_dims
     INTEGER :: grid_size, grid_corners, grid_rank
     INTEGER :: NCID, IERR
     INTEGER :: grid_size_dimid, grid_rank_dimid, grid_corners_dimid
     INTEGER :: grid_center_lat_varid, grid_center_lon_varid
     INTEGER :: grid_corner_lat_varid, grid_corner_lon_varid
     INTEGER :: grid_area_varid, grid_imask_varid 
     INTEGER :: grid_dims_varid
     CHARACTER (:), ALLOCATABLE :: GRID_UNITS, GRID_NAME

     NAMELIST / inputs / waves_mesh_file
     NAMELIST / outputs / waves_scrip_file

     OPEN(UNIT=15, FILE='scrip.nml')
     READ(UNIT=15, NML=inputs)
     READ(UNIT=15, NML=outputs)
     CLOSE(15)
     
     CALL read_gmsh_file(TRIM(waves_mesh_file), mesh)
     
     ALLOCATE(xyb(mesh%nn,3))
     DO i = 1,mesh%nn
       xyb(i,1) = mesh%xy(1,i)
       xyb(i,2) = mesh%xy(2,i)
       xyb(i,3) = mesh%depth(i)
     ENDDO
     
     ALLOCATE(trigp(mesh%ne,3))
     DO i = 1,mesh%ne
       DO j = 1,3
         trigp(i,j) = mesh%ect(j,i)
       ENDDO
     ENDDO
     
     
     CALL GET_SCRIP_INFO(mesh%ne, mesh%nn, xyb, trigp, &
                         grid_center_lon, grid_center_lat, &
                         grid_corner_lon, grid_corner_lat, grid_mask, &
                         grid_dims, grid_size, grid_corners, grid_rank)

     GRID_UNITS='degrees' ! the other option is radians...we don't use this
     GRID_NAME='src' ! this is not used, except for netcdf output


     IERR = NF90_CREATE(TRIM(waves_scrip_file), NF90_NETCDF4, NCID)
     IERR = NF90_DEF_DIM(NCID, 'grid_size', GRID_SIZE, grid_size_dimid)
     IERR = NF90_DEF_DIM(NCID, 'grid_corners', GRID_CORNERS, grid_corners_dimid)
     IERR = NF90_DEF_DIM(NCID, 'grid_rank', GRID_RANK, grid_rank_dimid)

     IERR = NF90_DEF_VAR(NCID, 'grid_center_lat', NF90_DOUBLE, &
                         (/grid_size_dimid/),grid_center_lat_varid)
     IERR = NF90_DEF_VAR(NCID, 'grid_center_lon', NF90_DOUBLE, &
                         (/grid_size_dimid/),grid_center_lon_varid)
     IERR = NF90_DEF_VAR(NCID, 'grid_corner_lat', NF90_DOUBLE, &
                         (/grid_corners_dimid,grid_size_dimid/), &
                         grid_corner_lat_varid)
     IERR = NF90_DEF_VAR(NCID, 'grid_corner_lon', NF90_DOUBLE, &
                         (/grid_corners_dimid,grid_size_dimid/), &
                         grid_corner_lon_varid)
     IERR = NF90_DEF_VAR(NCID, 'grid_imask', NF90_INT, &
                         (/grid_size_dimid/),grid_imask_varid)
     IERR = NF90_DEF_VAR(NCID, 'grid_dims', NF90_INT, &
                         (/grid_rank_dimid/),grid_dims_varid)
     IERR = NF90_ENDDEF(NCID)

     ALLOCATE(GRID_IMASK(GRID_DIMS(1)))
     GRID_IMASK = 0
     DO I = 1,GRID_DIMS(1)
       IF (GRID_MASK(I)) THEN
         GRID_IMASK(I) = 1
       ENDIF
     ENDDO

     IERR = NF90_PUT_ATT(NCID,grid_center_lat_varid,'units',GRID_UNITS)
     IERR = NF90_PUT_ATT(NCID,grid_center_lon_varid,'units',GRID_UNITS)
     IERR = NF90_PUT_ATT(NCID,grid_corner_lat_varid,'units',GRID_UNITS)
     IERR = NF90_PUT_ATT(NCID,grid_corner_lon_varid,'units',GRID_UNITS)
     IERR = NF90_PUT_ATT(NCID,grid_imask_varid,'units','unitless')

     IERR = NF90_PUT_VAR(NCID,grid_center_lat_varid,GRID_CENTER_LAT)
     IERR = NF90_PUT_VAR(NCID,grid_center_lon_varid,GRID_CENTER_LON)
     IERR = NF90_PUT_VAR(NCID,grid_corner_lat_varid,GRID_CORNER_LAT)
     IERR = NF90_PUT_VAR(NCID,grid_corner_lon_varid,GRID_CORNER_LON)
     IERR = NF90_PUT_VAR(NCID,grid_imask_varid,GRID_IMASK)
     IERR = NF90_PUT_VAR(NCID,grid_dims_varid,GRID_DIMS)
     IERR = NF90_CLOSE(NCID)

END PROGRAM scrip
