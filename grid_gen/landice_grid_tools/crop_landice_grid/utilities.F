module utilities

contains

subroutine write_OpenDX(    on_a_sphere, &
                            nCells, &
                            nVertices, &
                            nEdges, &
                            vertexDegree, &
                            maxEdges, &
                            xCell, &
                            yCell, &
                            zCell, &
                            xVertex, &
                            yVertex, &
                            zVertex, &
                            xEdge, &
                            yEdge, &
                            zEdge, &
                            nEdgesOnCell, &
                            verticesOnCell, &
                            verticesOnEdge, &
                            cellsOnVertex, &
                            edgesOnCell, &
                            areaCell, &
                            maxLevelCell, &
                            meshSpacing, &
                            depthCell, &
                            SST, &
                            kiteAreasOnVertex )

      implicit none

      character (len=16), intent(in) :: on_a_sphere
      integer, intent(in) :: nCells, nVertices, vertexDegree, nEdges, maxEdges
      real (kind=8), dimension(nCells), intent(inout) :: xCell
      real (kind=8), dimension(nCells), intent(inout) :: yCell
      real (kind=8), dimension(nCells), intent(inout) :: zCell
      real (kind=8), dimension(nVertices), intent(inout) :: xVertex
      real (kind=8), dimension(nVertices), intent(inout) :: yVertex
      real (kind=8), dimension(nVertices), intent(inout) :: zVertex
      real (kind=8), dimension(nEdges), intent(inout) :: xEdge
      real (kind=8), dimension(nEdges), intent(inout) :: yEdge
      real (kind=8), dimension(nEdges), intent(inout) :: zEdge
      integer, dimension(nCells), intent(in) :: nEdgesOnCell
      integer, dimension(maxEdges,nCells), intent(in) :: verticesOnCell
      integer, dimension(maxEdges,nCells), intent(in) :: edgesOnCell
      integer, dimension(2,nEdges), intent(in) :: verticesOnEdge
      integer, dimension(vertexDegree, nVertices), intent(in) :: cellsOnVertex
      integer, dimension(nCells), intent(in) :: maxLevelCell
      real (kind=8), dimension(nCells), intent(in) :: areaCell
      real (kind=8), dimension(nCells), intent(in) :: depthCell, SST, meshSpacing
      real (kind=8), dimension(vertexDegree,nVertices), intent(in) :: kiteAreasOnVertex

      character(len=80) :: a, b, c, d, e, f
      integer :: i, j, k, nVerticesTotal, iEdge, iLoop, iFace, Vert(4), Edge(4), iVertex, i1, i2, jp1
      integer :: nKitesTotal, iCell, iEdge1, iEdge2, iVertex11, iVertex12, iVertex21, iVertex22, ksave
      real (kind=8) :: x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4, xscale, work(nCells), work1(nCells), work2(nCells)
      real (kind=8) :: xv, yv, zv, xc, yc, zc, dist
      logical (kind=8) :: eflag

      if(on_a_sphere.eq.'NO              ') then
         write(6,*) ' write_dx, not on a sphere '
      endif

      xscale = 1.00
      xCell = xCell*xscale
      yCell = yCell*xscale
      zCell = zCell*xscale
      xVertex = xVertex*xscale
      yVertex = yVertex*xscale
      zVertex = zVertex*xscale
      xEdge = xEdge*xscale
      yEdge = yEdge*xscale
      zEdge = zEdge*xscale

      write(6,*) 'xCell', minval(xCell), maxval(xCell)
      write(6,*) ' nCells', nCells
      write(6,*) ' nEdges', nEdges
      write(6,*) ' nVertices', nVertices
      write(6,*) ' nEdgesOnCell',minval(nEdgesOnCell), maxval(nEdgesOnCell)

      open(unit=1,file='dx/vector.dx',form='formatted',status='unknown')

      a = trim('object "positions list" class array type float rank 1 shape 3 items')
      b = trim('ascii data file vector.position.data')
      write(1,10) a, nCells
      write(1,10) b
      write(1,*)

      a = trim('object 0  class array type float rank 1 shape 3 items')
      b = trim('ascii data file vector.data')
      c = trim('attribute "dep" string "positions"')
      write(1,10) a, nCells
      write(1,10) b
      write(1,10) c
      write(1,*)

      a = trim('object "vector" class field')
      b = trim('component "positions"     "positions list"')
      c = trim('component "data"           0')
      write(1,10) a
      write(1,10) b
      write(1,10) c

      close(1)

      open(unit=14,file='dx/vector.position.data',form='formatted',status='unknown')
      do i=1,nCells
       write(14,22) xCell(i), yCell(i), zCell(i)
      enddo
      close(14)



      nVerticesTotal = 0
      do i=1,nCells
       nVerticesTotal = nVerticesTotal + nEdgesOnCell(i)
      enddo
      write(6,*) 'total number of vertices', nVerticesTotal

      open(unit=1,file='dx/ocean.dx',form='formatted',status='unknown')

      a = trim('object "positions list" class array type float rank 1 shape 3 items')
      b = trim('ascii data file ocean.position.data')
      write(1,10) a, nVerticesTotal
      write(1,10) b
      write(1,*)
      10 format(a70,i10)
      
      a = trim('object "edge list" class array type int rank 0 items')
      b = trim('ascii data file ocean.edge.data')
      c = trim('attribute "ref" string "positions"')
      write(1,10) a, nVerticesTotal
      write(1,10) b
      write(1,10) c
      write(1,*)
      
      a = trim('object "loops list" class array type int rank 0 items')
      b = trim('ascii data file ocean.loop.data')
      c = trim('attribute "ref" string "edges"')
      write(1,10) a, nCells
      write(1,10) b
      write(1,10) c
      write(1,*)
      
      a = trim('object "face list" class array type int rank 0 items')
      b = trim('ascii data file ocean.face.data')
      c = trim('attribute "ref" string "loops"')
      write(1,10) a, nCells
      write(1,10) b
      write(1,10) c
      write(1,*)
      
      a = trim('object 0  class array type float rank 0 items')
      b = trim('data file ocean.meshSpacing.data')
      c = trim('attribute "dep" string "faces"')
      write(1,10) a, nCells
      write(1,10) b
      write(1,10) c
      write(1,*)
      
      a = trim('object "area" class field')
      b = trim('component "positions"     "positions list"')
      c = trim('component "edges"         "edge list"')
      d = trim('component "loops"         "loops list"')
      e = trim('component "faces"         "face list"')
      f = trim('component "data"           0')
      write(1,10) a
      write(1,10) b
      write(1,10) c
      write(1,10) d
      write(1,10) e
      write(1,10) f

      close(1)

     
      work2 = meshSpacing
      work1 = depthCell
      work = SST

      open(unit= 8,file='dx/ocean.meshSpacing.data',form='formatted',status='unknown')
      open(unit= 9,file='dx/ocean.depth.data',form='formatted',status='unknown')
      open(unit=10,file='dx/ocean.area.data',form='formatted',status='unknown')
      open(unit=11,file='dx/ocean.face.data',form='formatted',status='unknown')
      open(unit=12,file='dx/ocean.loop.data',form='formatted',status='unknown')
      open(unit=13,file='dx/ocean.edge.data',form='formatted',status='unknown')
      open(unit=14,file='dx/ocean.position.data',form='formatted',status='unknown')

      iLoop = 0
      iEdge = 0
      do i=1,nCells
       write(8,20) work2(i)
       write(9,20) work1(i)
       write(10,20) work(i)
       write(11,21) i-1
       write(12,21) iLoop
       iLoop = iLoop + nEdgesOnCell(i)

       eflag = .false.
       do j=1,nEdgesOnCell(i)
         k = verticesOnCell(j,i)
         xv = xVertex(k); yv = yVertex(k); zv = zVertex(k)
         xc = xCell(i); yc = yCell(i); zc = zCell(i)
         dist = sqrt( (xc-xv)**2 + (yc-yv)**2 + (zc-zv)**2 )
         if(dist.gt.5.0e5.and.on_a_sphere.eq.'NO              ') then
           eflag = .true.
         endif
       enddo

       if(eflag) then

       do j=1,nEdgesOnCell(i)
         write(13,21) iEdge
         iEdge = iEdge + 1
         k = verticesOnCell(j,i)
         xv = xVertex(k); yv = yVertex(k); zv = zVertex(k)
         xc = xCell(i); yc = yCell(i); zc = zCell(i)
         dist = sqrt( (xc-xv)**2 + (yc-yv)**2 + (zc-zv)**2 )
         if(dist.gt.5.0e5) then
            write(14,22) xc, yc, zc
         else
            write(14,22) xv, yv, zv
         endif
       enddo

       else

       do j=1,nEdgesOnCell(i)
         write(13,21) iEdge
         iEdge = iEdge + 1
         k = verticesOnCell(j,i)
         if(k.le.0) write(6,*) ' vert1 ',k, verticesOnCell(:,i)
         write(14,22) xVertex(k), yVertex(k), zVertex(k)
         write(15,23) j,i,k,xVertex(k), yVertex(k), zVertex(k)
       enddo
      endif
      enddo

 20   format(e20.10)
 21   format(i20)
 22   format(3e20.10)
 23   format(3i8, 3e20.10)

      close(9)
      close(10)
      close(11)
      close(12)
      close(13)
      close(14)

  !   nVerticesTotal = 0
  !   nKitesTotal = 0
  !   do i=1,nCells
  !    nKitesTotal = nKitesTotal + nEdgesOnCell(i)
  !   enddo
  !   nVerticesTotal = nKitesTotal*4
  !   write(6,*) nKitesTotal, nVerticesTotal

  !   open(unit=1,file='dx/kite.dx',form='formatted',status='unknown')

  !   a = trim('object "positions list" class array type float rank 1 shape 3 items')
  !   b = trim('ascii data file kite.position.data')
  !   write(1,10) a, nVerticesTotal
  !   write(1,10) b
  !   write(1,*)

  !   a = trim('object "edge list" class array type int rank 0 items')
  !   b = trim('ascii data file kite.edge.data')
  !   c = trim('attribute "ref" string "positions"')
  !   write(1,10) a, nVerticesTotal
  !   write(1,10) b
  !   write(1,10) c
  !   write(1,*)

  !   a = trim('object "loops list" class array type int rank 0 items')
  !   b = trim('ascii data file kite.loop.data')
  !   c = trim('attribute "ref" string "edges"')
  !   write(1,10) a, nKitesTotal
  !   write(1,10) b
  !   write(1,10) c
  !   write(1,*)

  !   a = trim('object "face list" class array type int rank 0 items')
  !   b = trim('ascii data file kite.face.data')
  !   c = trim('attribute "ref" string "loops"')
  !   write(1,10) a, nKitesTotal
  !   write(1,10) b
  !   write(1,10) c
  !   write(1,*)

  !   a = trim('object 0  class array type float rank 0 items')
  !   b = trim('data file kite.area.data')
  !   c = trim('attribute "dep" string "faces"')
  !   write(1,10) a, nKitesTotal
  !   write(1,10) b
  !   write(1,10) c
  !   write(1,*)

  !   a = trim('object "area" class field')
  !   b = trim('component "positions"     "positions list"')
  !   c = trim('component "edges"         "edge list"')
  !   d = trim('component "loops"         "loops list"')
  !   e = trim('component "faces"         "face list"')
  !   f = trim('component "data"           0')
  !   write(1,10) a
  !   write(1,10) b
  !   write(1,10) c
  !   write(1,10) d
  !   write(1,10) e
  !   write(1,10) f

  !   close(1)

  !   open(unit=10,file='dx/kite.area.data',form='formatted',status='unknown')
  !   open(unit=11,file='dx/kite.face.data',form='formatted',status='unknown')
  !   open(unit=12,file='dx/kite.loop.data',form='formatted',status='unknown')
  !   open(unit=13,file='dx/kite.edge.data',form='formatted',status='unknown')
  !   open(unit=14,file='dx/kite.position.data',form='formatted',status='unknown')

  !   iLoop = 0
  !   iEdge = 0
  !   iFace = 0

  !   do iCell=1,nCells
  !     do j=1,nEdgesOnCell(iCell)
  !        iEdge1 = edgesOnCell(j,iCell)
  !        jp1 = j+1
  !        if(j.eq.nEdgesOnCell(iCell)) jp1=1
  !        iEdge2 = edgesOnCell(jp1,iCell)

  !        iVertex11 = verticesOnEdge(1,iEdge1)
  !        iVertex21 = verticesOnEdge(2,iEdge1)
  !        iVertex12 = verticesOnEdge(1,iEdge2)
  !        ivertex22 = verticesOnEdge(2,iEdge2)

  !        if(iVertex11.eq.iVertex12.or.iVertex11.eq.iVertex22) then
  !           iVertex = iVertex11
  !        elseif(iVertex21.eq.iVertex12.or.iVertex21.eq.iVertex22)  then
  !           iVertex = iVertex21
  !        else
  !           write(6,*) iVertex11, iVertex21, iVertex12, iVertex22
  !           stop
  !        endif

  !        ksave = 0
  !        do k=1,vertexDegree
  !          if(cellsOnVertex(k,iVertex).eq.iCell) ksave=k
  !        enddo
  !        if(ksave.eq.0) then 
  !           write(6,*) ' can not find iCell'
  !           write(6,*) cellsOnVertex(:,iVertex)
  !           write(6,*) iCell
  !           write(6,*) iEdge1, iEdge2
  !           write(6,*) iVertex11, iVertex21, iVertex21, iVertex22
  !           write(6,*) iVertex
  !           stop
  !         endif

  !        write(11,21) iFace
  !        write(12,21) iLoop
  !        iFace = iFace + 1
  !        iLoop = iLoop + 4
  !        do k=1,4
  !          write(13,21) iEdge
  !          iEdge = iEdge + 1
  !        enddo
 !
 !         x1 = xCell(iCell)    ; y1 = yCell(iCell)    ; z1 = zCell(iCell)
 !         x2 = xEdge(iEdge1)   ; y2 = yEdge(iEdge1)   ; z2 = zEdge(iEdge1)
 !         x3 = xVertex(iVertex); y3 = yVertex(iVertex); z3 = zVertex(iVertex)
 !         x4 = xEdge(iEdge2)   ; y4 = yEdge(iEdge2)   ; z4 = zEdge(iEdge2)
 !
 !         write(14,22) x1, y1, z1
 !         write(14,22) x2, y2, z2
 !         write(14,22) x3, y3, z3
 !         write(14,22) x4, y4, z4
 !         write(10,22) kiteAreasOnVertex(ksave,iVertex)

 !      enddo
 !    enddo

end subroutine write_OpenDX


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE CONVERT_LX
!
! Convert (lat,lon) to an (x, y, z) location on a sphere with specified radius.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine convert_lx(x, y, z, radius, lat, lon)

   implicit none

   real, intent(in) :: radius
   real, intent(in) :: lat, lon
   real, intent(out) :: x, y, z

   z = radius * sin(lat)
   x = radius * cos(lon) * cos(lat)
   y = radius * sin(lon) * cos(lat)

end subroutine convert_lx


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE CONVERT_XL
!
! Convert (x, y, z) to a (lat, lon) location on a sphere with
!    radius sqrt(x^2 + y^2 + z^2).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine convert_xl(x, y, z, lat,lon)

   implicit none

   real, intent(in) :: x, y, z
   real, intent(out) :: lat, lon

   real :: dl, clat, pii, rtod
   real :: eps
   parameter (eps=1.e-10)

   pii = 2.*asin(1.0)
   rtod=180./pii
   dl = sqrt(x*x + y*y + z*z)

   lat = asin(z/dl)

!  check for being close to either pole

   if (abs(x) > eps) then

      if (abs(y) > eps) then

         lon = atan(abs(y/x))

         if ((x <= 0.) .and. (y >= 0.)) then
            lon = pii-lon
         else if ((x <= 0.) .and. (y < 0.)) then
            lon = lon+pii
         else if ((x >= 0.) .and. (y <= 0.)) then
            lon = 2*pii-lon
         end if

      else ! we're either on longitude 0 or 180

         if (x > 0) then
            lon = 0.
         else
            lon = pii
         end if

      end if

   else if (abs(y) > eps) then

      if (y > 0) then
         lon = pii/2.
      else
         lon = 3.*pii/2.
      end if

   else  ! we are at a pole

      lon = 0.

   end if

end subroutine convert_xl


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine transform_from_lonlat_to_xyz(xin, yin, zin, ulon, ulat, ux, uy, uz)
!
!  transform vector measured in latitude/longitude space to a vector measured in x,y,z
!
!     INTENT(IN)
!     xin = x position
!     yin = y position
!     zin = z position
!     ulon = east component of vector
!     ulat = north component of vector
!
!     INTENT(OUT)
!     ux = x component of vector
!     uy = y component of vector
!     uz = z component of vector
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

implicit none
real, intent(in) :: xin, yin, zin, ulon, ulat
real, intent(out) :: ux, uy, uz
real :: h(3,3), p(3), q(3), g(3), X1(3,3), X2(3,3), trans_X2_to_X1(3,3), r
integer :: i,j,k
logical :: l_Pole
real, parameter :: epsvt = 1.0e-10

!-----------------------------------------------------------------------
! define the e1, e2, and e3 directions
!-----------------------------------------------------------------------
        X1(1,1) = 1.0; X1(1,2) = 0.0; X1(1,3) = 0.0
        X1(2,1) = 0.0; X1(2,2) = 1.0; X1(2,3) = 0.0
        X1(3,1) = 0.0; X1(3,2) = 0.0; X1(3,3) = 1.0

!-----------------------------------------------------------------------
! find the vectors (measured in X1) that point in the local
!   east (h(1,:)), north (h(2,:)), and vertical (h(3,:)) direction
!-----------------------------------------------------------------------
        h(3,1) = xin; h(3,2) = yin; h(3,3) = zin
        call unit_vector_in_3space(h(3,:))

!-----------------------------------------------------------------------
! g(:) is a work array and holds the vector pointing to the North Pole.
! measured in X1
!-----------------------------------------------------------------------
              g(:) = X1(3,:)

!-----------------------------------------------------------------------
! determine if the local vertical hits a pole
!-----------------------------------------------------------------------
              l_Pole = .false.
              r = g(1)*h(3,1) + g(2)*h(3,2) + g(3)*h(3,3)
              r = abs(r) + epsvt
              if(r.gt.1.0) then
                l_Pole = .true.
                h(3,:) = h(3,:) + epsvt
                call unit_vector_in_3space(h(3,:))
              endif

!-----------------------------------------------------------------------
! find the vector that is perpendicular to the local vertical vector
! and points in the direction of of the North pole, this defines the local
! north direction. measured in X1
!-----------------------------------------------------------------------
              call vector_on_tangent_plane ( h(3,:), g(:), h(2,:) )

!-----------------------------------------------------------------------
! take the cross product of the local North direction and the local vertical
! to find the local east vector. still in X1
!-----------------------------------------------------------------------
              call cross_product_in_3space ( h(2,:), h(3,:), h(1,:) )

!-----------------------------------------------------------------------
! put these 3 vectors into a matrix X2
!-----------------------------------------------------------------------
              X2(1,:) = h(1,:)              ! local east     (measured in X1)
              X2(2,:) = h(2,:)              ! local north    (measured in X1)
              X2(3,:) = h(3,:)              ! local vertical (measured in X1)

!-----------------------------------------------------------------------
! compute the transformation matrix
!-----------------------------------------------------------------------
              trans_X2_to_X1(:,:) = matmul(X1,transpose(X2))

!-----------------------------------------------------------------------
! transform (ulon, ulat) into (x,y,z)
!-----------------------------------------------------------------------
              p(1) = ulon; p(2) = ulat; p(3) = 0
              g(:) = matmul(trans_X2_to_X1(:, :), p(:))
              ux = g(1); uy = g(2); uz = g(3)

end subroutine transform_from_lonlat_to_xyz


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine transform_from_xyz_to_lonlat(xin, yin, zin, ux, uy, uz, ulon, ulat)
!
!  transform vector measured in x,y,z space to a vector measured in latitude/longitude space
!
!     INTENT(IN)
!     xin = x position
!     yin = y position
!     zin = z position
!     ux = x component of vector
!     uy = y component of vector
!     uz = z component of vector
!
!     INTENT(OUT)
!     ulon = east component of vector
!     ulat = north component of vector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

implicit none
real, intent(in) :: xin, yin, zin, ux, uy, uz
real, intent(out) :: ulon, ulat
real :: h(3,3), p(3), q(3), g(3), X1(3,3), X2(3,3), trans_X1_to_X2(3,3), r
integer :: i,j,k
logical :: l_Pole
real, parameter :: epsvt = 1.0e-10

!-----------------------------------------------------------------------
! define the e1, e2, and e3 directions
!-----------------------------------------------------------------------
        X1(1,1) = 1.0; X1(1,2) = 0.0; X1(1,3) = 0.0
        X1(2,1) = 0.0; X1(2,2) = 1.0; X1(2,3) = 0.0
        X1(3,1) = 0.0; X1(3,2) = 0.0; X1(3,3) = 1.0

!-----------------------------------------------------------------------
! find the vectors (measured in X1) that point in the local
!   east (h(1,:)), north (h(2,:)), and vertical (h(3,:)) direction
!-----------------------------------------------------------------------
        h(3,1) = xin; h(3,2) = yin; h(3,3) = zin
        call unit_vector_in_3space(h(3,:))

!-----------------------------------------------------------------------
! g(:) is a work array and holds the vector pointing to the North Pole.
! measured in X1
!-----------------------------------------------------------------------
              g(:) = X1(3,:)

!-----------------------------------------------------------------------
! determine if the local vertical hits a pole
!-----------------------------------------------------------------------
              l_Pole = .false.
              r = g(1)*h(3,1) + g(2)*h(3,2) + g(3)*h(3,3)
              r = abs(r) + epsvt
              if(r.gt.1.0) then
                l_Pole = .true.
                h(3,:) = h(3,:) + epsvt
                call unit_vector_in_3space(h(3,:))
              endif

!-----------------------------------------------------------------------
! find the vector that is perpendicular to the local vertical vector
! and points in the direction of of the North pole, this defines the local
! north direction. measured in X1
!-----------------------------------------------------------------------
              call vector_on_tangent_plane ( h(3,:), g(:), h(2,:) )

!-----------------------------------------------------------------------
! take the cross product of the local North direction and the local vertical
! to find the local east vector. still in X1
!-----------------------------------------------------------------------
              call cross_product_in_3space ( h(2,:), h(3,:), h(1,:) )

!-----------------------------------------------------------------------
! put these 3 vectors into a matrix X2
!-----------------------------------------------------------------------
              X2(1,:) = h(1,:)              ! local east     (measured in X1)
              X2(2,:) = h(2,:)              ! local north    (measured in X1)
              X2(3,:) = h(3,:)              ! local vertical (measured in X1)

!-----------------------------------------------------------------------
! compute the transformation matrix
!-----------------------------------------------------------------------
              trans_X1_to_X2(:,:) = matmul(X2,transpose(X1))

!-----------------------------------------------------------------------
! transform (ulon, ulat) into (x,y,z)
!-----------------------------------------------------------------------
              p(1) = ux; p(2) = uy; p(3) = uz
              g(:) = matmul(trans_X1_to_X2(:, :), p(:))
              ulon = g(1); ulat= g(2);

end subroutine transform_from_xyz_to_lonlat

!======================================================================
! BEGINNING OF UNIT_VECTOR_IN_3SPACE
!======================================================================
        subroutine unit_vector_in_3space (p_1)

!-----------------------------------------------------------------------
! PURPOSE : normalize p_1 to unit length and overwrite p_1
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! intent(inout)
!-----------------------------------------------------------------------
        real , intent(inout) ::                         &
                        p_1 (:)

!-----------------------------------------------------------------------
! local
!-----------------------------------------------------------------------
        real  :: length

        length = SQRT (p_1(1)**2 + p_1(2)**2 + p_1(3)**2 )
        length = 1.0/length
        p_1(1) = p_1(1)*length
        p_1(2) = p_1(2)*length
        p_1(3) = p_1(3)*length

        end subroutine unit_vector_in_3space
!======================================================================
! END OF UNIT_VECTOR_IN_3SPACE
!======================================================================

!======================================================================
! BEGINNING OF CROSS_PRODUCT_IN_3SPACE
!======================================================================
        subroutine cross_product_in_3space(p_1,p_2,p_out)

!-----------------------------------------------------------------------
! PURPOSE: compute p_1 cross p_2 and place in p_out
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! intent(in)
!-----------------------------------------------------------------------
        real , intent(in) ::                            &
                        p_1 (:),                                      &
                        p_2 (:)

!-----------------------------------------------------------------------
! intent(out)
!-----------------------------------------------------------------------
        real , intent(out) ::                           &
                        p_out (:)

        p_out(1) = p_1(2)*p_2(3)-p_1(3)*p_2(2)
        p_out(2) = p_1(3)*p_2(1)-p_1(1)*p_2(3)
        p_out(3) = p_1(1)*p_2(2)-p_1(2)*p_2(1)

        end subroutine cross_product_in_3space
!======================================================================
! END OF CROSS_PRODUCT_IN_3SPACE
!======================================================================

!======================================================================
! BEGINNING OF VECTOR_ON_TANGENT_PLANE
!======================================================================
        subroutine vector_on_tangent_plane(p_1, p_2, p_out)

!-----------------------------------------------------------------------
! PURPOSE : given two points measured in (x,y,z) and lying on
!       the unit sphere, find the vector (p_out) that lies on the plane
!       perpendicular to the p_1 vector and points in the direction of
!       the projection of p_2 onto the tangent plane.
!
! NOTE : p_1 and p_2 are assumed to be of unit length
! NOTE : p_out is normalized to unit length
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! intent(in)
!-----------------------------------------------------------------------
        real , intent(in) ::                            &
                        p_1 (:),                                      &
                        p_2 (:)

!-----------------------------------------------------------------------
! intent(out)
!-----------------------------------------------------------------------
        real , intent(out) ::                           &
                        p_out (:)

!-----------------------------------------------------------------------
! local
!-----------------------------------------------------------------------
        real  ::                                        &
                        work (3), t1(3), t2(3)

!       work (1) = - p_1(2) * ( -p_1(2) * p_2(1) + p_1(1) * p_2(2) )   &
!                  + p_1(3) * (  p_1(3) * p_2(1) - p_1(1) * p_2(3) )

!       work (2) = + p_1(1) * ( -p_1(2) * p_2(1) + p_1(1) * p_2(2) )   &
!                  - p_1(3) * ( -p_1(3) * p_2(2) + p_1(2) * p_2(3) )

!       work (3) = - p_1(1) * (  p_1(3) * p_2(1) - p_1(1) * p_2(3) )   &
!                  + p_1(2) * ( -p_1(3) * p_2(2) + p_1(2) * p_2(3) )


        t1(:) = p_2(:) - p_1(:)
        t2(:) = p_1

        call unit_vector_in_3space (t1)
        call unit_vector_in_3space (t2)

        call cross_product_in_3space(t1(:), t2(:), work(:))
        call unit_vector_in_3space (work)
        call cross_product_in_3space(t2(:),work(:),p_out(:))
        call unit_vector_in_3space (p_out)

        end subroutine vector_on_tangent_plane
!======================================================================
! END OF VECTOR_ON_TANGENT_PLANE
!======================================================================

end module utilities
