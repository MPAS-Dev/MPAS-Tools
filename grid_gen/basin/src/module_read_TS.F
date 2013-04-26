module read_TS
 
   integer :: rd_ncid, rd_ncids, rd_ncidu
   integer :: rdDimIDt_lon
   integer :: rdDimIDt_lat
   integer :: rdDimIDdepth_t
   integer :: rdVarIDt_lon
   integer :: rdVarIDt_lat
   integer :: rdVarIDdepth_t
   integer :: rdVarIDTEMP
   integer :: rdVarIDSALT
   integer :: rdVarIDTAUX
   integer :: rdVarIDTAUY
 
   integer :: rdLocalt_lon
   integer :: rdLocalt_lat
   integer :: rdLocaldepth_t
 
   contains
 
   subroutine read_TS_init(nx, ny, nz, fileNameT, fileNameS, fileNameU)
 
      implicit none
 
      include 'netcdf.inc'
 
      character(len=80), intent(in) :: fileNameT, fileNameS, fileNameU
      integer, intent(out) :: nx, ny, nz
 
      integer :: nferr, nferrs, nferru
 
      nferr = nf_open(fileNameT, NF_SHARE, rd_ncid)
      write(6,*) ' nferr ', nferr, rd_ncid
 
      !
      ! Get IDs for variable dimensions
      !
      nferr = nf_inq_dimid(rd_ncid, 't_lon', rdDimIDt_lon)
      write(6,*) ' nferr ', nferr, rdDimIDt_lon
      nferr = nf_inq_dimlen(rd_ncid, rdDimIDt_lon, rdLocalt_lon)
      write(6,*) ' nferr ', nferr, rdLocalt_lon
      nferr = nf_inq_dimid(rd_ncid, 't_lat', rdDimIDt_lat)
      write(6,*) ' nferr ', nferr, rdDimIDt_lat
      nferr = nf_inq_dimlen(rd_ncid, rdDimIDt_lat, rdLocalt_lat)
      write(6,*) ' nferr ', nferr, rdLocalt_lat
      nferr = nf_inq_dimid(rd_ncid, 'depth_t', rdDimIDdepth_t)
      write(6,*) ' nferr ', nferr, rdDimIDdepth_t
      nferr = nf_inq_dimlen(rd_ncid, rdDimIDdepth_t, rdLocaldepth_t)
      write(6,*) ' nferr ', nferr, rdLocaldepth_t

      nx = rdLocalt_lon
      ny = rdLocalt_lat
      nz = rdLocaldepth_t

      write(6,*) nx, ny, nz
 
      !
      ! Get IDs for variables
      !
      nferr = nf_inq_varid(rd_ncid, 't_lon', rdVarIDt_lon)
      write(6,*) ' nferr ', nferr, rdVarIDt_lon
      nferr = nf_inq_varid(rd_ncid, 't_lat', rdVarIDt_lat)
      write(6,*) ' nferr ', nferr, rdVarIDt_lat
      nferr = nf_inq_varid(rd_ncid, 'depth_t', rdVarIDdepth_t)
      write(6,*) ' nferr ', nferr, rdVarIDdepth_t
      nferr = nf_inq_varid(rd_ncid, 'TEMP', rdVarIDTEMP)
      write(6,*) ' nferr ', nferr, rdVarIDTEMP

      nferrs = nf_open(fileNameS, NF_SHARE, rd_ncids)
      nferrs = nf_inq_varid(rd_ncids, 'SALT', rdVarIDSALT)
      write(6,*) ' nferrs ', nferrs, rdVarIDSALT

      nferru = nf_open(fileNameU, NF_SHARE, rd_ncidu)
      nferru = nf_inq_varid(rd_ncidu, 'TAUX', rdVarIDTAUX)
      nferru = nf_inq_varid(rd_ncidu, 'TAUY', rdVarIDTAUY)
      write(6,*) ' nferru ', nferru, rdVarIDTAUX, rdVarIDTAUY
 
   end subroutine read_TS_init
 
   subroutine read_TS_fields(t_lon, t_lat, depth_t, TEMP, SALT, TAUX, TAUY)
 
      implicit none
 
      include 'netcdf.inc'
 
      real (kind=4), dimension(:), intent(out) :: t_lon, t_lat, depth_t
      real (kind=4), dimension(:,:,:), intent(out) :: TEMP, SALT
      real (kind=4), dimension(:,:), intent(out) :: TAUX, TAUY

      integer, dimension(1) :: start1, count1
      integer, dimension(2) :: start2, count2
      integer, dimension(3) :: start3, count3
      integer, dimension(4) :: start4, count4

      integer :: nferr, nferrs, nferru

      start1(1) = 1
      count1(1) = rdLocalt_lon
      nferr = nf_get_vara_real(rd_ncid, rdVarIDt_lon, start1, count1, t_lon)
      write(6,*) ' nferr ', nferr, rd_ncid, rdVarIDt_lon

      start1(1) = 1
      count1(1) = rdLocalt_lat
      nferr = nf_get_vara_real(rd_ncid, rdVarIDt_lat, start1, count1, t_lat)
      write(6,*) ' nferr ', nferr, rd_ncid, rdVarIDt_lat

      start1(1) = 1
      count1(1) = rdLocaldepth_t
      nferr = nf_get_vara_real(rd_ncid, rdVarIDdepth_t, start1, count1, depth_t)
      write(6,*) ' nferr ', nferr, rd_ncid, rdVarIDdepth_t

      start3(1) = 1
      start3(2) = 1
      start3(3) = 1
      count3(1) = rdLocalt_lon
      count3(2) = rdLocalt_lat
      count3(3) = rdLocaldepth_t
      nferr = nf_get_vara_real(rd_ncid, rdVarIDTEMP, start3, count3, TEMP)
      write(6,*) ' nferr ', nferr, rd_ncid, rdVarIDTEMP
      write(6,*) ' temperature' , minval(TEMP), maxval(TEMP)

      start3(1) = 1
      start3(2) = 1
      start3(3) = 1
      count3(1) = rdLocalt_lon
      count3(2) = rdLocalt_lat
      count3(3) = rdLocaldepth_t
      nferrs = nf_get_vara_real(rd_ncids, rdVarIDSALT, start3, count3, SALT)
      write(6,*) ' nferrs ', nferrs, rd_ncids, rdVarIDSALT
      write(6,*) ' salinity' , minval(SALT), maxval(SALT)

      start2(1) = 1
      start2(2) = 1
      count2(1) = rdLocalt_lon
      count2(2) = rdLocalt_lat
      nferru = nf_get_vara_real(rd_ncidu, rdVarIDTAUX, start2, count2, TAUX)
      nferru = nf_get_vara_real(rd_ncidu, rdVarIDTAUY, start2, count2, TAUY)
      write(6,*) ' nferru ', nferru, rd_ncidu, rdVarIDTAUX, rdVarIDTAUY
      write(6,*) ' TAUX' , minval(TAUX), maxval(TAUX)
      write(6,*) ' TAUY' , minval(TAUY), maxval(TAUY)


   end subroutine read_TS_fields
 
 
   subroutine read_TS_finalize()
 
      implicit none
 
      include 'netcdf.inc'
 
      integer :: nferr, nferrs, nferru
 
      nferr = nf_close(rd_ncid)
      write(6,*) ' nferr ', nferr

      nferrs = nf_close(rd_ncids)
      write(6,*) ' nferrs ', nferrs

      nferru = nf_close(rd_ncidu)
      write(6,*) ' nferru ', nferru


   end subroutine read_TS_finalize
 
end module read_TS
