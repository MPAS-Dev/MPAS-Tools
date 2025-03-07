  module mpas_geometry_utilities

   use precision

   contains
     
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION MPAS_SPHERE_ANGLE
   !
   ! Computes the angle between arcs AB and AC, given points A, B, and C
   ! Equation numbers w.r.t. http://mathworld.wolfram.com/SphericalTrigonometry.html
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real (dp) function mpas_sphere_angle(ax, ay, az, bx, by, bz, cx, cy, cz)!{{{

      implicit none

      real (dp), intent(in) :: ax, ay, az, bx, by, bz, cx, cy, cz

      real (dp) :: a, b, c          ! Side lengths of spherical triangle ABC

      real (dp) :: ABx, ABy, ABz    ! The components of the vector AB
      real (dp) :: ACx, ACy, ACz    ! The components of the vector AC

      real (dp) :: Dx               ! The i-components of the cross product AB x AC
      real (dp) :: Dy               ! The j-components of the cross product AB x AC
      real (dp) :: Dz               ! The k-components of the cross product AB x AC

      real (dp) :: s                ! Semiperimeter of the triangle
      real (dp) :: sin_angle


      a = mpas_arc_length(bx, by, bz, cx, cy, cz)
      b = mpas_arc_length(ax, ay, az, cx, cy, cz)
      c = mpas_arc_length(ax, ay, az, bx, by, bz)

      ABx = bx - ax
      ABy = by - ay
      ABz = bz - az

      ACx = cx - ax
      ACy = cy - ay
      ACz = cz - az

      Dx =   (ABy * ACz) - (ABz * ACy)
      Dy = -((ABx * ACz) - (ABz * ACx))
      Dz =   (ABx * ACy) - (ABy * ACx)

      s = 0.5_dp*(a + b + c)
!      sin_angle = sqrt((sin(s-b)*sin(s-c))/(sin(b)*sin(c)))   ! Eqn. (28)
      sin_angle = sqrt(min(1.0_dp,max(0.0_dp,(sin(s-b)*sin(s-c))/(sin(b)*sin(c)))))   ! Eqn. (28)

      if ((Dx*ax + Dy*ay + Dz*az) >= 0.0) then
         mpas_sphere_angle =  2.0 * asin(max(min(sin_angle,1.0_dp),-1.0_dp))
      else
         mpas_sphere_angle = -2.0 * asin(max(min(sin_angle,1.0_dp),-1.0_dp))
      end if

   end function mpas_sphere_angle!}}}


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION MPAS_PLANE_ANGLE
   !
   ! Computes the angle between vectors AB and AC, given points A, B, and C, and
   !   a vector (u,v,w) normal to the plane.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real (dp) function mpas_plane_angle(ax, ay, az, bx, by, bz, cx, cy, cz, u, v, w)!{{{

      implicit none

      real (dp), intent(in) :: ax, ay, az, bx, by, bz, cx, cy, cz, u, v, w

      real (dp) :: ABx, ABy, ABz    ! The components of the vector AB
      real (dp) :: mAB              ! The magnitude of AB
      real (dp) :: ACx, ACy, ACz    ! The components of the vector AC
      real (dp) :: mAC              ! The magnitude of AC

      real (dp) :: Dx               ! The i-components of the cross product AB x AC
      real (dp) :: Dy               ! The j-components of the cross product AB x AC
      real (dp) :: Dz               ! The k-components of the cross product AB x AC

      real (dp) :: cos_angle

      ABx = bx - ax
      ABy = by - ay
      ABz = bz - az
      mAB = sqrt(ABx**2.0_dp + ABy**2.0_dp + ABz**2.0_dp)

      ACx = cx - ax
      ACy = cy - ay
      ACz = cz - az
      mAC = sqrt(ACx**2.0_dp + ACy**2.0_dp + ACz**2.0_dp)


      Dx =   (ABy * ACz) - (ABz * ACy)
      Dy = -((ABx * ACz) - (ABz * ACx))
      Dz =   (ABx * ACy) - (ABy * ACx)

      cos_angle = (ABx*ACx + ABy*ACy + ABz*ACz) / (mAB * mAC)

      if ((Dx*u + Dy*v + Dz*w) >= 0.0_dp) then
         mpas_plane_angle =  acos(max(min(cos_angle,1.0_dp),-1.0_dp))
      else
         mpas_plane_angle = -acos(max(min(cos_angle,1.0_dp),-1.0_dp))
      end if

   end function mpas_plane_angle!}}}


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION MPAS_ARC_LENGTH
   !
   ! Returns the length of the great circle arc from A=(ax, ay, az) to
   !    B=(bx, by, bz). It is assumed that both A and B lie on the surface of the
   !    same sphere centered at the origin.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real (dp) function mpas_arc_length(ax, ay, az, bx, by, bz)!{{{

      implicit none

      real (dp), intent(in) :: ax, ay, az, bx, by, bz

      real (dp) :: r, c
      real (dp) :: cx, cy, cz

      cx = bx - ax
      cy = by - ay
      cz = bz - az

!      r = ax*ax + ay*ay + az*az
!      c = cx*cx + cy*cy + cz*cz
!
!      arc_length = sqrt(r) * acos(1.0 - c/(2.0*r))

      r = sqrt(ax*ax + ay*ay + az*az)
      c = sqrt(cx*cx + cy*cy + cz*cz)
!      arc_length = sqrt(r) * 2.0 * asin(c/(2.0*r))
      mpas_arc_length = r * 2.0_dp * asin(c/(2.0_dp*r))

   end function mpas_arc_length!}}}


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTine mpas_arc_bisect
   !
   ! Returns the point C=(cx, cy, cz) that bisects the great circle arc from
   !   A=(ax, ay, az) to B=(bx, by, bz). It is assumed that A and B lie on the
   !   surface of a sphere centered at the origin.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine mpas_arc_bisect(ax, ay, az, bx, by, bz, cx, cy, cz)!{{{

      implicit none

      real (dp), intent(in) :: ax, ay, az, bx, by, bz
      real (dp), intent(out) :: cx, cy, cz

      real (dp) :: r           ! Radius of the sphere
      real (dp) :: d

      r = sqrt(ax*ax + ay*ay + az*az)

      cx = 0.5_dp*(ax + bx)
      cy = 0.5_dp*(ay + by)
      cz = 0.5_dp*(az + bz)

      if (cx == 0._dp .and. cy == 0._dp .and. cz == 0._dp) then
         write(6,*) 'arc_bisect: A and B are diametrically opposite'
         stop
      else
         d = sqrt(cx*cx + cy*cy + cz*cz)
         cx = r * cx / d
         cy = r * cy / d
         cz = r * cz / d
      end if

   end subroutine mpas_arc_bisect!}}}


!***********************************************************************
!
!  routine mpas_distance_plane
!
!> \brief   Calculates distance between two points on a plane.
!> \author  Matthew Hoffman
!> \date    5 February 2015
!> \details
!>  This routine calculates distance between two points on a plane.
!>  It does not work for periodic meshes.
!-----------------------------------------------------------------------
   real(dp) function mpas_distance_plane(a, b)!{{{
      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      real(dp), dimension(3), intent(in) :: a, b  !< Input: 3d (x,y,z) points between which to calculate distance

      mpas_distance_plane = sqrt(sum((a - b)**2))
   end function mpas_distance_plane !}}}


!***********************************************************************
!
!  routine mpas_triangle_signed_area
!
!> \brief   Calculates area of a triangle, whether on a sphere or a plane.
!> \author  Matthew Hoffman
!> \date    13 January 2015
!> \details
!>  This routine calculates the area of a triangle whether on a sphere or a plane.
!>  Note this does not handle triangles spanning planar periodic meshes because mpas_triangle_signed_area_plane does not!
!-----------------------------------------------------------------------
   real(dp) function mpas_triangle_signed_area(a, b, c, on_a_sphere, radius)
      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      real(dp), dimension(3), intent(in) :: a, b, c  !< Input: 3d (x,y,z) points forming the triangle for which to get the area
      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      logical :: on_a_sphere
      real(dp), dimension(3) :: normalvec
      real(dp) :: radius

      ! call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)
      ! call mpas_pool_get_config(meshPool, 'sphere_radius', radius)

      on_a_sphere = .true.
      if (on_a_sphere) then
         mpas_triangle_signed_area = mpas_triangle_signed_area_sphere(a, b, c, radius)
      else
         normalvec = (/ 0, 0, 1 /)
         mpas_triangle_signed_area = mpas_triangle_signed_area_plane(a, b, c, normalvec)
      endif
   end function mpas_triangle_signed_area


!***********************************************************************
!
!  routine mpas_triangle_signed_area_plane
!
!> \brief   Calculates signed area of a triangle in a plane
!> \author  Matthew Hoffman
!> \date    13 January 2015
!> \details
!>  This routine calculates the area of a triangle in a plane.
!>  Uses cross product.  Signed area will be positive if the vertices are oriented counterclockwise.
!>  Note this does not handle triangles spanning periodic meshes!
!-----------------------------------------------------------------------
   real(dp) function mpas_triangle_signed_area_plane(a, b, c, normalvec)!{{{
      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      real(dp), dimension(3), intent(in) :: a, b, c  !< Input: 3d (x,y,z) points forming the triangle for which to calculate the area
      real(dp), dimension(3), intent(in) :: normalvec  !< Input: 3d vector indicating the normal direction for the plane for assigning a sign to the area
      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      real(dp), dimension(3) :: ab, ac, crossprod, triangleNormal

      ab = b - a
      ac = c - a
      call mpas_cross_product_in_r3(ab, ac, crossprod)
      if (mpas_vec_mag_in_r3(crossprod) == 0.0_dp) then
         mpas_triangle_signed_area_plane = 0.0_dp
      else
         triangleNormal = crossprod / mpas_vec_mag_in_r3(crossprod)
         mpas_triangle_signed_area_plane = 0.5_dp * (mpas_vec_mag_in_r3(crossprod)) *  &
              sum(triangleNormal * normalvec)
      endif
    end function mpas_triangle_signed_area_plane !}}}

!***********************************************************************
!
!  routine mpas_triangle_signed_area_sphere
!
!> \brief   Calculates area of a triangle on a sphere
!> \author  Matthew Hoffman
!> \date    13 January 2015
!> \details
!>  This routine calculates the area of a triangle on the surface of a sphere.
!>  Uses the spherical analog of Heron's formula.
!>  Copied from mesh generator.  A CCW winding angle is positive.
!-----------------------------------------------------------------------
   real(dp) function mpas_triangle_signed_area_sphere(a, b, c, radius)!{{{
      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      real(dp), dimension(3), intent(in) :: a, b, c  !< Input: 3d (x,y,z) points forming the triangle in which to calculate the bary weights
      real(dp), intent(in) :: radius  !< sphere radius
      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      real(dp) :: ab, bc, ca, semiperim, tanqe
      real(dp), dimension(3) :: ablen, aclen, Dlen

      ab = mpas_arc_length(a(1), a(2), a(3), b(1), b(2), b(3))/radius
      bc = mpas_arc_length(b(1), b(2), b(3), c(1), c(2), c(3))/radius
      ca = mpas_arc_length(c(1), c(2), c(3), a(1), a(2), a(3))/radius
      semiperim = 0.5 * (ab + bc + ca)

      tanqe = sqrt(max(0.0_dp,tan(0.5_dp * semiperim) * tan(0.5_dp * (semiperim - ab)) &
                   * tan(0.5_dp * (semiperim - bc)) * tan(0.5_dp * (semiperim - ca))))

      mpas_triangle_signed_area_sphere = 4.0_dp * radius * radius * atan(tanqe)

      ! computing correct signs (in similar fashion to mpas_sphere_angle)
      ablen(1) = b(1) - a(1)
      ablen(2) = b(2) - a(2)
      ablen(3) = b(3) - a(3)

      aclen(1) = c(1) - a(1)
      aclen(2) = c(2) - a(2)
      aclen(3) = c(3) - a(3)

      dlen(1) =   (ablen(2) * aclen(3)) - (ablen(3) * aclen(2))
      dlen(2) = -((ablen(1) * aclen(3)) - (ablen(3) * aclen(1)))
      dlen(3) =   (ablen(1) * aclen(2)) - (ablen(2) * aclen(1))

      if ((Dlen(1)*a(1) + Dlen(2)*a(2) + Dlen(3)*a(3)) < 0.0_dp) then
        mpas_triangle_signed_area_sphere = -mpas_triangle_signed_area_sphere
      end if

   end function mpas_triangle_signed_area_sphere !}}}

!***********************************************************************
!
!  routine mpas_vec_mag_in_r3
!
!> \brief   MPAS 3D vector magnitude routine
!> \author  Matt Hoffman
!> \date    13 Jan. 2015
!> \details
!> This routine calculates the magnitude of a 3d vector.
!-----------------------------------------------------------------------
  real (dp) function mpas_vec_mag_in_r3(xin)!{{{
    implicit none
    real (dp), dimension(3), intent(in) :: xin !< Input: Vector
    mpas_vec_mag_in_r3 = sqrt(xin(1)**2 + xin(2)**2 + xin(3)**2)
  end function mpas_vec_mag_in_r3!}}}

!***********************************************************************
!
!  routine mpas_cross_product_in_r3
!
!> \brief   MPAS 3D cross product routine
!> \author  Xylar Asay-Davis
!> \date    03/28/13
!> \details
!> This routine computes the cross product of two input vectors.
!-----------------------------------------------------------------------
  subroutine mpas_cross_product_in_r3(p_1,p_2,p_out)!{{{
    real (dp), intent(in) :: p_1 (3) !< Input: Vector 1
    real (dp), intent(in) :: p_2 (3) !< Input: Vector 2
    real (dp), intent(out) :: p_out (3) !< Output: Cross product of vector 1 and vector 2

    p_out(1) = p_1(2)*p_2(3)-p_1(3)*p_2(2)
    p_out(2) = p_1(3)*p_2(1)-p_1(1)*p_2(3)
    p_out(3) = p_1(1)*p_2(2)-p_1(2)*p_2(1)
  end subroutine mpas_cross_product_in_r3!}}}

end module mpas_geometry_utilities
