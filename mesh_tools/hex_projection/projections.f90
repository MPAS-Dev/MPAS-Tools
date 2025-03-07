module projections

  use precision
  use projection_setup
  
  contains

subroutine xy_to_lambert_conformal_projection( xpos, ypos, zpos, latitude, longitude, mapFactor, npts )

  implicit none

  integer :: npts
  real(dp), dimension(npts), intent(inout)  :: xpos, ypos
  real(dp), dimension(npts), intent(out) :: zpos, latitude, longitude, mapFactor

  ! xpos and ypos are points on a plane that are to be projected using the lambert conformal
  ! projection to the sphere.  (x,y) = 0 is the center of the projection

  ! we determine the latitude and longitude from the xpos and ypos values on the plane, and then
  ! reset the xpos, ypos and zpos to be the coordinates of the sphere

  ! hardwire the CONUS config for HRRR and RRFS

  real(dp) :: reference_latitude
  real(dp) :: standard_parallel_1
  real(dp) :: standard_parallel_2
  real(dp) :: earth_radius
  real(dp) :: reference_longitude
  real(dp) :: standard_longitude

  real(dp) :: rho
  real(dp) :: rho_0
  real(dp) :: lambda
  real(dp) :: F_factor
  real(dp) :: n_factor
  real(dp) :: phi
  real(dp) :: pii
  real(dp) :: degrees_to_radians, radians_to_degrees
  real(dp) :: theta

  real(dp) :: F_over_rho
  real(dp) :: factor_map
  real(dp) :: xshift

  integer :: i


  pii = 2.0_dp*asin(1.0_dp)
  degrees_to_radians = pii/180.0_dp
  radians_to_degrees = 1.0_dp/degrees_to_radians

  reference_longitude = reference_longitude_degrees * degrees_to_radians
  standard_longitude = standard_longitude_degrees * degrees_to_radians
  reference_latitude = reference_latitude_degrees * degrees_to_radians
  standard_parallel_1 = standard_parallel_1_degrees * degrees_to_radians
  standard_parallel_2 = standard_parallel_2_degrees * degrees_to_radians
  earth_radius = earth_radius_km * 1000.0_dp  ! meters

  write(6,*) " projection parameters "
  write(6,*) " standard_longitude (degrees) ", standard_longitude_degrees
  write(6,*) " reference_longitude (degrees) ", reference_longitude_degrees
  write(6,*) " reference latitude (degrees)  ", reference_latitude_degrees
  write(6,*) " standard_parallel_1 (degrees) ", standard_parallel_1_degrees
  write(6,*) " standard_parallel_2 (degrees) ", standard_parallel_2_degrees
  write(6,*) " Earth radius (km)             ", earth_radius_km
  
  write(6,*) " min and max xpos ",minval(xpos),maxval(xpos)
  write(6,*) " min and max ypos ",minval(ypos),maxval(ypos)

  ! compute some constant factors given the lambert conformal parameters

  n_factor = log(cos(standard_parallel_1)/cos(standard_parallel_2))
  if(n_factor == 0.0) then
     n_factor = sin(standard_parallel_1)
  else
     n_factor = n_factor/log( tan(0.25_dp*pii+0.5_dp*standard_parallel_2)/tan(0.25_dp*pii+0.5_dp*standard_parallel_1) )
  end if

  F_factor = tan(0.25_dp*pii+0.5_dp*standard_parallel_1)**n_factor
  F_factor = cos(standard_parallel_1)*F_factor/n_factor
  rho_0 = (1./tan(0.25_dp*pii+0.5_dp*reference_latitude))**n_factor
  rho_0 = F_factor*earth_radius*rho_0

  xshift = rho_0*sin(n_factor*(reference_longitude -  standard_longitude))
  write(6,*) " xshift = ",xshift
  
  do i=1,npts
     xpos(i) = xpos(i) + xshift
     rho = sign(1.0_dp, n_factor)*sqrt(xpos(i)**2 + (rho_0 - ypos(i))**2)
     theta = atan( xpos(i)/(rho_0 - ypos(i)) )
     lambda = standard_longitude + theta/n_factor
     F_over_rho = earth_radius*F_factor/rho
     phi = 2.0_dp*atan((F_over_rho)**(1.0_dp/n_factor)) - 0.5_dp*pii

     factor_map = cos(standard_parallel_1)*(tan(0.25_dp*pii + 0.5_dp*standard_parallel_1)**n_factor)
     factor_map = factor_map/(cos(phi)*( (tan(0.25_dp*pii+0.5_dp*phi))**n_factor ))

     if(lambda < 0.0_dp) lambda = lambda + pii*2.0_dp
     
     mapFactor(i) = factor_map
     latitude(i) = phi
     longitude(i) = lambda
     xpos(i) = earth_radius*cos(phi)*cos(lambda)
     ypos(i) = earth_radius*cos(phi)*sin(lambda)
     zpos(i) = earth_radius*sin(phi)

     lambda = lambda*radians_to_degrees
     phi = phi*radians_to_degrees

!     write(6,*) " pt, longitude and latitude ",i,lambda, phi
!     write(6,*) " "
  end do
  write(6,*) " min and max mapFactor ",minval(mapFactor(:)), maxval(mapFactor(:))
     
end subroutine xy_to_lambert_conformal_projection

end module projections
