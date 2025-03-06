module projection_setup

  !  use, intrinsic :: iso_fortran_env, dp=>real64
  use precision
  
  real(dp) :: cell_spacing_km
  real(dp) :: mesh_length_x_km
  real(dp) :: mesh_length_y_km
  real(dp) :: earth_radius_km

  character(len=20) :: projection_type

  real(dp) :: reference_longitude_degrees
  real(dp) :: reference_latitude_degrees
  real(dp) :: standard_longitude_degrees
  real(dp) :: standard_parallel_1_degrees
  real(dp) :: standard_parallel_2_degrees

  namelist /mesh/ cell_spacing_km, mesh_length_x_km, mesh_length_y_km, earth_radius_km
  namelist /projection/ projection_type
  namelist /lambert_conformal/ reference_longitude_degrees, reference_latitude_degrees, &
                               standard_longitude_degrees,                              &
                               standard_parallel_1_degrees, standard_parallel_2_degrees

contains

  subroutine read_namelist()

  open(unit=10,file='namelist.projections',status='old')
  read(10, nml=mesh)
  read(10, nml=projection)
  read(10, nml=lambert_conformal)

end subroutine read_namelist

end module projection_setup
