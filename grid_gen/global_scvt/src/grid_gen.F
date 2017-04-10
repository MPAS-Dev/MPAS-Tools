program grid_gen

   use grid_params
   use grid_constants
   use data_types
   use grid_gen_utils
   use voronoi_utils
   use stripack
   use scvt
   use grid_meta

   implicit none

   real(kind=RKIND) :: dlat, dlon
   real(kind=RKIND) :: dl
   real(kind=RKIND) :: d1, d2, d3, d4
   integer :: p1, p2, p3
   integer :: if
   character (len=80) :: frame_name

   real(kind=RKIND) :: pi
   real(kind=RKIND) :: area_per_sample, nhexs, sum_nhexs, hex_area
   type (geo_point) :: p

   integer :: i, j, k, nb, ier

   real(kind=RKIND), allocatable, dimension(:)       :: rlat, rlon, vclat, vclon, x, y, z, xc, yc, zc
   integer, allocatable, dimension(:)    :: list, lptr, listc, lend
   integer, allocatable, dimension(:,:)  :: ltri

   integer :: n, nrow, ntmax, nvc
   integer :: tr1, tr2, tr3, tr4
   integer :: i1, i2, i3, k1, k2

   type (adjacency_list) :: alist, clist

   real(kind=RKIND) :: lat1, lon1, lat2, lon2, lat3, lon3, latc, lonc


   call read_namelist()


   pi = 4.0_RKIND*atan(1.0_RKIND)

   area_per_sample = 4.0_RKIND * pi * 6370000**2.0_RKIND / 6000000.0_RKIND
   sum_nhexs = 0.0_RKIND
   write(0,'(a,f10.1)') 'Computing an estimate for the required number of cells to reach dx=', min_dx
   do if = 1,5
      nhexs = 0.0_RKIND
      do i=1,6000000
         call random_point(p)
         d1 = density_for_point(p)
         dl = min_dx / (d1 ** 0.25)
         hex_area = sqrt(3.0_RKIND) / 2.0_RKIND * dl**2.0_RKIND
         nhexs = nhexs + area_per_sample / hex_area
      end do
!      write(0,'(a,i2,a,i)') 'Estimate ',if,' for required # hexs:', nint(nhexs)
      sum_nhexs = sum_nhexs + nhexs
      write(0,'(a,i3,a)',advance='no') ' ...',if*20,'%'
   end do
   write(0,*) ' '
   write(0,*) 'Estimated # hexs:', nint(sum_nhexs/5.0_RKIND)
   write(0,*) ' '


   n = np

   ntmax = 6*n
   nvc = ntmax
   nrow = 6
   allocate(rlat(n))
   allocate(rlon(n))
   allocate(x(n))
   allocate(y(n))
   allocate(z(n))
   allocate(ltri(nrow,ntmax))
   allocate(list(nvc))
   allocate(lptr(nvc))
   allocate(lend(n))
   allocate(listc(nvc))
   allocate(vclat(nvc))
   allocate(vclon(nvc))
   allocate(xc(nvc))
   allocate(yc(nvc))
   allocate(zc(nvc))


   write(0,*) ' '
   write(0,*) 'Reading generating points from a file'
   open(22,file='locs.dat',form='formatted',status='old')
   if (locs_as_xyz) then
      read(22,*)
      do i=1,n
         read(22,'(10x,f22.10,f23.10,f23.10)') x(i), y(i), z(i)
      end do
      call trans_inv(x, y, z, rlat, rlon, n)
   else
      do i=1,n
         read(22,'(f13.10,1x,f13.10)') rlat(i), rlon(i)
      end do
   end if
   close(22)


   !
   ! Compute Voronoi corners
   !
   write(0,*) ' '
   write(0,*) 'Computing Voronoi corners'
   call compute_vc(rlat, rlon, n, nrow, ntmax, list, lptr, lend, listc, vclat, vclon, nvc)


   !
   ! Form SCVT
   !
   call TRANS (n, rlat, rlon, x, y, z)
   call TRANS (nvc, vclat, vclon, xc, yc, zc)
   write(frame_name,'(a)') 'scvt_initial.ps'
   open(32,file=trim(frame_name),form='formatted',status='unknown')
   call vrplot(32, 8.0_RKIND, 0.0_RKIND, 0.0_RKIND, 90.0_RKIND ,N, X,Y,Z, 2*n-4,LISTC,LPTR,LEND,XC,YC,ZC,'(spherical centroidal voronoi tessellation)',.false.,IER)
   close(32)

   call scvt_solve(n, lend, rlat, rlon, nvc, list, lptr, if)

   call compute_vc(rlat, rlon, n, nrow, ntmax, list, lptr, lend, listc, vclat, vclon, nvc)
   call TRANS (n, rlat, rlon, x, y, z)
   call TRANS (nvc, vclat, vclon, xc, yc, zc)
   write(frame_name,'(a)') 'scvt_final.ps'
   open(32,file=trim(frame_name),form='formatted',status='unknown')
   call vrplot(32, 8.0_RKIND, 0.0_RKIND, 0.0_RKIND, 90.0_RKIND ,N, X,Y,Z, 2*n-4,LISTC,LPTR,LEND,XC,YC,ZC,'(spherical centroidal voronoi tessellation)',.false.,IER)
   close(32)


   write(0,*) ' '
   write(0,*) 'Deriving grid metadata and writing output'
   call convert_adjacency_list(n, lend, nvc, list, lptr, alist)
   call convert_corner_list(n, lend, nvc, listc, lptr, clist)
   call write_grid(rlat, rlon, n, vclat, vclon, nvc, alist, clist)


   deallocate(rlat)
   deallocate(x)
   deallocate(y)
   deallocate(z)
   deallocate(rlon)
   deallocate(ltri)
   deallocate(list)
   deallocate(lptr)
   deallocate(lend)
   deallocate(listc)
   deallocate(vclat)
   deallocate(vclon)
   deallocate(xc)
   deallocate(yc)
   deallocate(zc)

end program grid_gen
