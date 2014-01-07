module read_U
 
   integer :: rd_ncidu
   integer :: rdDimIDu_lon
   integer :: rdDimIDu_lat
   integer :: rdDimIDdepth_t
   integer :: rdVarIDu_lon
   integer :: rdVarIDu_lat
   integer :: rdVarIDdepth_t
   integer :: rdVarIDTAUX
   integer :: rdVarIDTAUY
 
   integer :: rdLocalu_lon
   integer :: rdLocalu_lat
   integer :: rdLocaldepth_t
 
   contains
 
   subroutine read_U_init(nx, ny, nz, fileNameU)
 
      implicit none
 
      include 'netcdf.inc'
 
      character(len=80), intent(in) :: fileNameU
      integer, intent(out) :: nx, ny, nz
 
      integer :: nferru
 
      nferru = nf_open(fileNameU, NF_SHARE, rd_ncidu)
      write(6,*) ' nferru ', nferru, rd_ncidu
 
      !
      ! Get IDs for variable dimensions
      !
      nferru = nf_inq_dimid(rd_ncidu, 'u_lon', rdDimIDu_lon)
      write(6,*) ' nferru ', nferru, rdDimIDu_lon
      nferru = nf_inq_dimlen(rd_ncidu, rdDimIDu_lon, rdLocalu_lon)
      write(6,*) ' nferru ', nferru, rdLocalu_lon
      nferru = nf_inq_dimid(rd_ncidu, 'u_lat', rdDimIDu_lat)
      write(6,*) ' nferru ', nferru, rdDimIDu_lat
      nferru = nf_inq_dimlen(rd_ncidu, rdDimIDu_lat, rdLocalu_lat)
      write(6,*) ' nferru ', nferru, rdLocalu_lat
      nferru = nf_inq_dimid(rd_ncidu, 'depth_t', rdDimIDdepth_t)
      write(6,*) ' nferru ', nferru, rdDimIDdepth_t
      nferru = nf_inq_dimlen(rd_ncidu, rdDimIDdepth_t, rdLocaldepth_t)
      write(6,*) ' nferru ', nferru, rdLocaldepth_t

      nx = rdLocalu_lon
      ny = rdLocalu_lat
      nz = rdLocaldepth_t

      write(6,*) nx, ny, nz
 
      !
      ! Get IDs for variables
      !
      nferru = nf_inq_varid(rd_ncidu, 'u_lon', rdVarIDu_lon)
      write(6,*) ' nferru ', nferru, rdVarIDu_lon
      nferru = nf_inq_varid(rd_ncidu, 'u_lat', rdVarIDu_lat)
      write(6,*) ' nferru ', nferru, rdVarIDu_lat
      nferru = nf_inq_varid(rd_ncidu, 'depth_t', rdVarIDdepth_t)
      write(6,*) ' nferru ', nferru, rdVarIDdepth_t

      nferru = nf_inq_varid(rd_ncidu, 'TAUX', rdVarIDTAUX)
      nferru = nf_inq_varid(rd_ncidu, 'TAUY', rdVarIDTAUY)
      write(6,*) ' nferru ', nferru, rdVarIDTAUX, rdVarIDTAUY
 
   end subroutine read_U_init
 
   subroutine read_U_fields(u_lon, u_lat, depth_t, TAUX, TAUY)
 
      implicit none
 
      include 'netcdf.inc'
 
      real (kind=4), dimension(:), intent(out) :: u_lon, u_lat, depth_t
      real (kind=4), dimension(:,:), intent(out) :: TAUX, TAUY

      integer, dimension(1) :: start1, count1
      integer, dimension(2) :: start2, count2
      integer, dimension(3) :: start3, count3
      integer, dimension(4) :: start4, count4

      integer :: nferru

      start1(1) = 1
      count1(1) = rdLocalu_lon
      nferru = nf_get_vara_real(rd_ncidu, rdVarIDu_lon, start1, count1, u_lon)
      write(6,*) ' nferru ', nferru, rd_ncidu, rdVarIDu_lon

      start1(1) = 1
      count1(1) = rdLocalu_lat
      nferru = nf_get_vara_real(rd_ncidu, rdVarIDu_lat, start1, count1, u_lat)
      write(6,*) ' nferru ', nferru, rd_ncidu, rdVarIDu_lat

      start1(1) = 1
      count1(1) = rdLocaldepth_t
      nferru = nf_get_vara_real(rd_ncidu, rdVarIDdepth_t, start1, count1, depth_t)
      write(6,*) ' nferru ', nferru, rd_ncidu, rdVarIDdepth_t

      start2(1) = 1
      start2(2) = 1
      count2(1) = rdLocalu_lon
      count2(2) = rdLocalu_lat
      nferru = nf_get_vara_real(rd_ncidu, rdVarIDTAUX, start2, count2, TAUX)
      nferru = nf_get_vara_real(rd_ncidu, rdVarIDTAUY, start2, count2, TAUY)
      write(6,*) ' nferru ', nferru, rd_ncidu, rdVarIDTAUX, rdVarIDTAUY
      write(6,*) ' TAUX' , minval(TAUX), maxval(TAUX)
      write(6,*) ' TAUY' , minval(TAUY), maxval(TAUY)


   end subroutine read_U_fields
 
 
   subroutine read_U_finalize()
 
      implicit none
 
      include 'netcdf.inc'
 
      integer :: nferru
 
      nferru = nf_close(rd_ncidu)
      write(6,*) ' nferru ', nferru


   end subroutine read_U_finalize
 
end module read_U
