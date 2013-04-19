module grid_params

   integer :: np
   logical :: locs_as_xyz
   logical :: l2_conv, inf_conv
   integer :: n_scvt_iterations
   real :: eps
   real :: min_dx

   contains

   subroutine read_namelist()

      implicit none

      integer :: funit
      real :: pi

      namelist /domains/ np, locs_as_xyz, n_scvt_iterations, eps, l2_conv, inf_conv, min_dx

      pi = 4.0*atan(1.0)

      funit = 21 

      np = 40962
      locs_as_xyz = .true.
      n_scvt_iterations = 0
      eps = 0.0000000001
      l2_conv = .true.
      inf_conv = .false.
      min_dx = 120000.0

      open(funit,file='namelist.input',status='old',form='formatted')
      read(funit,domains)
      close(funit)

      if(l2_conv) then
          inf_conv = .false.
      endif

   end subroutine read_namelist

end module grid_params
