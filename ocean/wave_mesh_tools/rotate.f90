      PROGRAM rotate

      USE rotate_mod
      USE read_write_gmsh

      IMPLICIT NONE

      REAL(rp) :: LON_POLE,LAT_POLE 
      INTEGER :: NP
      REAL(rp), DIMENSION(:), ALLOCATABLE :: LON,LAT
      REAL(rp), DIMENSION(:), ALLOCATABLE :: LONR,LATR
      REAL(rp), DIMENSION(:), ALLOCATABLE :: UVEC,VVEC 
      REAL(rp), DIMENSION(:), ALLOCATABLE :: UVECR,VVECR 
      REAL(rp), DIMENSION(:), ALLOCATABLE :: ANGLED
      INTEGER :: ncid,lon_dimid,lat_dimid,t_dimid
      INTEGER :: lon_varid,lat_varid
      INTEGER :: u_varid,v_varid
      INTEGER :: nlon,nlat,ntime
      INTEGER :: lon_pole_varid,lat_pole_varid
      REAL(rp), DIMENSION(:), ALLOCATABLE :: lon_nc,lat_nc
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: u_nc,v_nc
      INTEGER :: i,j,k,n
      CHARACTER (:), ALLOCATABLE :: wind_file
      CHARACTER (:), ALLOCATABLE :: mesh_file
      CHARACTER (:), ALLOCATABLE :: mesh_file_out
      TYPE(grid_type) :: grid
      LOGICAL :: file_exists

      NAMELIST / inputs / LON_POLE, LAT_POLE, wind_file, mesh_file, mesh_file_out

      LON_POLE = -42.8906d0 
      LAT_POLE = 72.3200d0
      wind_file = 'wnd10mx0.5.gdas.200506.ww3.nc'
      mesh_file = 'mesh_shallow_4000_edit_mv_nd.msh'
      mesh_file_out = 'mesh_shallow_4000_edit_mv_nd_RTD.msh'
      
      OPEN(UNIT=15, FILE='rotate.nml')
      READ(UNIT=15, NML=inputs)
      CLOSE(15)

      PRINT*, 'LON_POLE = ',LON_POLE
      PRINT*, 'LAT_POLE = ',LAT_POLE
      PRINT*, 'wind_file = ',wind_file
      PRINT*, 'mesh_file = ',mesh_file
      PRINT*, 'mesh_file_out = ',mesh_file_out


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      ! Rotate mesh coordinates 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      CALL read_gmsh_file(mesh_file,grid)
      ALLOCATE(LON(grid%nn),LAT(grid%nn))
      ALLOCATE(LONR(grid%nn),LATR(grid%nn))
      ALLOCATE(ANGLED(grid%nn))

      LONR = 0d0
      LATR = 0d0
      ANGLED = 0d0
      DO i = 1,grid%nn
        LON(i) = grid%xy(1,i)
        LAT(i) = grid%xy(2,i)
      ENDDO

      CALL W3LLTOEQ ( LAT , LON, LATR, LONR,     &              
                      ANGLED, LAT_POLE, LON_POLE, grid%nn )             

      DO i = 1,grid%nn
        IF (LONR(i) > 180d0) THEN
          grid%xy(1,i) = LONR(i) - 360d0
        ELSE
          grid%xy(1,i) = LONR(i)
        ENDIF
        grid%xy(2,i) = LATR(i)
      ENDDO

      !DO i = 1,grid%nn
      !  PRINT("(I7,2(F12.6,F12.6,15x),F12.6)"), i, LON(i), LAT(i), LONR(i), LATR(i), ANGLED(i)
      !ENDDO


      CALL write_gmsh_file(mesh_file_out,grid)

      OPEN(UNIT=10,NAME='angled.d')
      DO i = 1,grid%nn
        WRITE(10,"(3(F12.6))") LON(i), LAT(i), ANGLED(i)
      ENDDO
      CLOSE(10)

      DEALLOCATE(LON,LAT)
      DEALLOCATE(LONR,LATR)
      DEALLOCATE(ANGLED)
    
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      ! Read wind data
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      INQUIRE(FILE=wind_file,EXIST=file_exists)
      IF (.not. file_exists) THEN
        PRINT*, "Wind file does not exist"
        STOP
      ENDIF

      CALL check(NF90_OPEN(wind_file,NF90_WRITE,ncid))      

      IF (NF90_INQ_VARID(ncid,'LON_POLE',lon_pole_varid) .eq. NF90_NOERR) THEN
        PRINT*, "Winds already rotated"
        STOP
      ENDIF

      !CALL check(NF90_INQ_DIMID(ncid,'lon',lon_dimid))
      !CALL check(NF90_INQ_DIMID(ncid,'lat',lat_dimid))
      CALL check(NF90_INQ_DIMID(ncid,'longitude',lon_dimid))
      CALL check(NF90_INQ_DIMID(ncid,'latitude',lat_dimid))
      CALL check(NF90_INQ_DIMID(ncid,'time',t_dimid))

      CALL check(NF90_INQUIRE_DIMENSION(ncid,lon_dimid,len=nlon))
      CALL check(NF90_INQUIRE_DIMENSION(ncid,lat_dimid,len=nlat))
      CALL check(NF90_INQUIRE_DIMENSION(ncid,t_dimid,len=ntime))
      PRINT*, nlon,nlat,ntime

      CALL check(NF90_INQ_VARID(ncid,'longitude',lon_varid))
      CALL check(NF90_INQ_VARID(ncid,'latitude',lat_varid))
      !CALL check(NF90_INQ_VARID(ncid,'U_GRD_L103',u_varid))
      !CALL check(NF90_INQ_VARID(ncid,'V_GRD_L103',v_varid))
      CALL check(NF90_INQ_VARID(ncid,'u10',u_varid))
      CALL check(NF90_INQ_VARID(ncid,'v10',v_varid))

      ALLOCATE(lon_nc(nlon),lat_nc(nlat))
      CALL check(NF90_GET_VAR(ncid,lon_varid,lon_nc))
      CALL check(NF90_GET_VAR(ncid,lat_varid,lat_nc))
      ALLOCATE(u_nc(nlon,nlat,ntime),v_nc(nlon,nlat,ntime))
      CALL check(NF90_GET_VAR(ncid,u_varid,u_nc))
      CALL check(NF90_GET_VAR(ncid,v_varid,v_nc))

      NP = nlon*nlat
      ALLOCATE(LON(NP), LAT(NP))
      ALLOCATE(LONR(NP), LATR(NP))
      ALLOCATE(ANGLED(NP))
      ALLOCATE(UVEC(NP),VVEC(NP))
      ALLOCATE(UVECR(NP),VVECR(NP))


      LONR = 0d0
      LATR = 0d0
      UVECR = 0d0
      VVECR = 0d0

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      ! Rotate wind coordinates 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      n = 1
      DO i = 1,nlat
        DO j = 1,nlon

          IF (lon_nc(j) > 180d0) THEN
            LON(n) = lon_nc(j) - 360d0
          ELSE
            LON(n) = lon_nc(j)
          ENDIF

          LAT(n) = lat_nc(i)

          n = n+1
        ENDDO
      ENDDO


      CALL W3LLTOEQ ( LAT , LON, LATR, LONR,     &              
                      ANGLED, LAT_POLE, LON_POLE, NP )             


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      ! Rotate wind vectors
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      DO k = 1,ntime
        PRINT*, 't = ',k

        n = 1
        DO i = 1,nlat
          DO j = 1,nlon

            UVEC(n) = u_nc(j,i,k)
            VVEC(n) = v_nc(j,i,k)

            n = n+1
          ENDDO
        ENDDO

        UVECR = UVEC
        VVECR = VVEC
        CALL W3XYRTN ( NP, UVECR, VVECR, -ANGLED)

        n = 1
        DO i = 1,nlat
          DO j = 1,nlon

            u_nc(j,i,k) = UVECR(n)
            v_nc(j,i,k) = VVECR(n)

            if (k == 1) then
              PRINT*, sqrt(UVEC(n)**2 + VVEC(n)**2), sqrt(UVECR(n)**2 + VVECR(n)**2)
            endif

            n = n + 1
          ENDDO
        ENDDO

      ENDDO


      CALL check(NF90_PUT_VAR(ncid,u_varid,u_nc))
      CALL check(NF90_PUT_VAR(ncid,v_varid,v_nc))

      CALL check(NF90_REDEF(ncid))
      CALL check(NF90_DEF_VAR(ncid,'LON_POLE',NF90_DOUBLE,lon_pole_varid))
      CALL check(NF90_DEF_VAR(ncid,'LAT_POLE',NF90_DOUBLE,lat_pole_varid))
      CALL check(NF90_ENDDEF(ncid))
      CALL check(NF90_PUT_VAR(ncid,lon_pole_varid,LON_POLE))
      CALL check(NF90_PUT_VAR(ncid,lat_pole_varid,LAT_POLE))

      CALL check(NF90_CLOSE(ncid))


      END PROGRAM rotate
