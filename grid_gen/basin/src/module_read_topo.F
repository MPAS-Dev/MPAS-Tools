module read_topo
 
   integer :: rd_ncid
   integer :: rdDimIDnx
   integer :: rdDimIDny
   integer :: rdVarIDz
   integer :: rdVarIDx
   integer :: rdVarIDy
 
   integer :: rdLocalnx
   integer :: rdLocalny
 
   contains
 
   subroutine read_topo_init( nx, ny)
 
      implicit none
 
      include 'netcdf.inc'
 
      integer, intent(out) :: nx, ny
 
      integer :: nferr
 
 
      nferr = nf_open('topo/ETOPO2v2c_f4.nc', NF_SHARE, rd_ncid)
      write(6,*) ' nferr ', nferr, rd_ncid
 
      !
      ! Get IDs for variable dimensions
      !
      nferr = nf_inq_dimid(rd_ncid, 'x', rdDimIDnx)
      write(6,*) ' nferr ', nferr, rdDimIDnx
      nferr = nf_inq_dimlen(rd_ncid, rdDimIDnx, rdLocalnx)
      write(6,*) ' nferr ', nferr, rdLocalnx
      nferr = nf_inq_dimid(rd_ncid, 'y', rdDimIDny)
      write(6,*) ' nferr ', nferr, rdDimIDny
      nferr = nf_inq_dimlen(rd_ncid, rdDimIDny, rdLocalny)
      write(6,*) ' nferr ', nferr, rdLocalny

      nx = rdLocalnx
      ny = rdLocalny

      write(6,*) nx, ny
 
      !
      ! Get IDs for variables
      !
      nferr = nf_inq_varid(rd_ncid, 'x', rdVarIDx)
      write(6,*) ' nferr ', nferr, rdVarIDx
      nferr = nf_inq_varid(rd_ncid, 'y', rdVarIDy)
      write(6,*) ' nferr ', nferr, rdVarIDy
      nferr = nf_inq_varid(rd_ncid, 'z', rdVarIDz)
      write(6,*) ' nferr ', nferr, rdVarIDz
 
   end subroutine read_topo_init
 
 
   subroutine read_topo_fields(x,y,z)
 
      implicit none
 
      include 'netcdf.inc'
 
      real (kind=4), dimension(:), intent(out) :: x,y
      real (kind=4), dimension(:,:), intent(out) :: z

      integer, dimension(1) :: start1, count1
      integer, dimension(2) :: start2, count2
      integer, dimension(3) :: start3, count3
      integer, dimension(4) :: start4, count4

      integer :: nferr

      start1(1) = 1
      count1(1) = rdLocalnx
      nferr = nf_get_vara_real(rd_ncid, rdVarIDx, start1, count1, x)
      write(6,*) ' nferr ', nferr, rd_ncid, rdVarIDx

      start1(1) = 1
      count1(1) = rdLocalny
      nferr = nf_get_vara_real(rd_ncid, rdVarIDy, start1, count1, y)
      write(6,*) ' nferr ', nferr, rd_ncid, rdVarIDy

      start2(1) = 1
      start2(2) = 1
      count2(1) = rdLocalnx
      count2(2) = rdLocalny
      nferr = nf_get_vara_real(rd_ncid, rdVarIDz, start2, count2, z)
      write(6,*) ' nferr ', nferr, rd_ncid, rdVarIDz, rdLocalnx
 
   end subroutine read_topo_fields
 
 
   subroutine read_topo_finalize()
 
      implicit none
 
      include 'netcdf.inc'
 
      integer :: nferr
 
      nferr = nf_close(rd_ncid)
      write(6,*) ' nferr ', nferr

 
   end subroutine read_topo_finalize
 
end module read_topo
