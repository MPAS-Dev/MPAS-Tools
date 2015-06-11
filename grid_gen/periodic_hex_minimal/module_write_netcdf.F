module write_netcdf
 
   integer :: wr_ncid
   integer :: wrDimIDnCells
   integer :: wrDimIDnEdges
   integer :: wrDimIDnVertices
   integer :: wrDimIDmaxEdges
   integer :: wrDimIDmaxEdges2
   integer :: wrDimIDTWO
   integer :: wrDimIDvertexDegree
   integer :: wrVarIDxCell
   integer :: wrVarIDyCell
   integer :: wrVarIDzCell
   integer :: wrVarIDxVertex
   integer :: wrVarIDyVertex
   integer :: wrVarIDzVertex
   integer :: wrVarIDcellsOnVertex
 
   integer :: wrLocalnCells
   integer :: wrLocalnVertices
 
   contains
 
   subroutine write_netcdf_init( &
                               nCells, &
                               nEdges, &
                               nVertices, &
                               maxEdges, &
                               vertexDegree, &
                               dc, &
                               nx, &
                               ny )
 
      implicit none
 
      include 'netcdf.inc'
 
      integer, intent(in) :: nCells
      integer, intent(in) :: nEdges
      integer, intent(in) :: nVertices
      integer, intent(in) :: maxEdges
      integer, intent(in) :: vertexDegree
      real (kind=8), intent(in) :: dc
      integer, intent(in) :: nx
      integer, intent(in) :: ny
 
      integer :: nferr
      integer, dimension(10) :: dimlist
      character (len=16) :: on_a_sphere
      character (len=16) :: is_periodic
      real (kind=8) :: sphere_radius
      real (kind=8) :: x_offset, y_offset
 
 
      wrLocalnCells = nCells
      wrLocalnVertices = nVertices

      on_a_sphere = 'NO'
      is_periodic = 'YES'
      sphere_radius = 0.0
      x_offset = (nx) * dc
      y_offset = (ny) * (dc * sqrt(3.0)) / 2.0
 
      nferr = nf_create('grid.nc', IOR(NF_CLOBBER,NF_64BIT_OFFSET), wr_ncid)
 
      !
      ! Define dimensions
      !
      nferr = nf_def_dim(wr_ncid, 'nCells', nCells, wrDimIDnCells)
      nferr = nf_def_dim(wr_ncid, 'nEdges', nEdges, wrDimIDnEdges)
      nferr = nf_def_dim(wr_ncid, 'nVertices', nVertices, wrDimIDnVertices)
      nferr = nf_def_dim(wr_ncid, 'maxEdges', maxEdges, wrDimIDmaxEdges)
      nferr = nf_def_dim(wr_ncid, 'maxEdges2', 2*maxEdges, wrDimIDmaxEdges2)
      nferr = nf_def_dim(wr_ncid, 'TWO', 2, wrDimIDTWO)
      nferr = nf_def_dim(wr_ncid, 'vertexDegree', vertexDegree, wrDimIDvertexDegree)


      !
      ! Define attributes
      !
      nferr = nf_put_att_text(wr_ncid, NF_GLOBAL, 'on_a_sphere', 16, on_a_sphere)
      nferr = nf_put_att_text(wr_ncid, NF_GLOBAL, 'is_periodic', 16, is_periodic)
      nferr = nf_put_att_double(wr_ncid, NF_GLOBAL, 'sphere_radius', NF_DOUBLE, 1, sphere_radius)
      nferr = nf_put_att_double(wr_ncid, NF_GLOBAL, 'x_offset', NF_DOUBLE, 1, x_offset)
      nferr = nf_put_att_double(wr_ncid, NF_GLOBAL, 'y_offset', NF_DOUBLE, 1, y_offset)

 
      !
      ! Define variables
      !
      dimlist( 1) = wrDimIDnCells
      nferr = nf_def_var(wr_ncid, 'xCell', NF_DOUBLE,  1, dimlist, wrVarIDxCell)
      dimlist( 1) = wrDimIDnCells
      nferr = nf_def_var(wr_ncid, 'yCell', NF_DOUBLE,  1, dimlist, wrVarIDyCell)
      dimlist( 1) = wrDimIDnCells
      nferr = nf_def_var(wr_ncid, 'zCell', NF_DOUBLE,  1, dimlist, wrVarIDzCell)
      dimlist( 1) = wrDimIDnVertices
      nferr = nf_def_var(wr_ncid, 'xVertex', NF_DOUBLE,  1, dimlist, wrVarIDxVertex)
      dimlist( 1) = wrDimIDnVertices
      nferr = nf_def_var(wr_ncid, 'yVertex', NF_DOUBLE,  1, dimlist, wrVarIDyVertex)
      dimlist( 1) = wrDimIDnVertices
      nferr = nf_def_var(wr_ncid, 'zVertex', NF_DOUBLE,  1, dimlist, wrVarIDzVertex)
      dimlist( 1) = wrDimIDvertexDegree
      dimlist( 2) = wrDimIDnVertices
      nferr = nf_def_var(wr_ncid, 'cellsOnVertex', NF_INT,  2, dimlist, wrVarIDcellsOnVertex)
 
      nferr = nf_enddef(wr_ncid)
 
   end subroutine write_netcdf_init
 
 
   subroutine write_netcdf_fields( &
                                  xCell, &
                                  yCell, &
                                  zCell, &
                                  xVertex, &
                                  yVertex, &
                                  zVertex, &
                                  cellsOnVertex )
 
      implicit none
 
      include 'netcdf.inc'
 
      real (kind=8), dimension(:), optional, intent(in) :: xCell
      real (kind=8), dimension(:), optional, intent(in) :: yCell
      real (kind=8), dimension(:), optional, intent(in) :: zCell
      real (kind=8), dimension(:), optional, intent(in) :: xVertex
      real (kind=8), dimension(:), optional, intent(in) :: yVertex
      real (kind=8), dimension(:), optional, intent(in) :: zVertex
      integer, dimension(:,:), optional, intent(in) :: cellsOnVertex

      integer :: nferr
      integer, dimension(1) :: start1, count1
      integer, dimension(2) :: start2, count2
      integer, dimension(3) :: start3, count3
      integer, dimension(4) :: start4, count4
 
      start1(1) = 1
 
      start2(1) = 1
      start2(2) = 1
 
      start3(1) = 1
      start3(2) = 1
      start3(3) = 1
 
      start4(1) = 1
      start4(2) = 1
      start4(3) = 1
      start4(4) = 1
 
      if(present(xCell))then      
      start1(1) = 1
      count1( 1) = wrLocalnCells
      nferr = nf_put_vara_double(wr_ncid, wrVarIDxCell, start1, count1, xCell)
      endif
 
      if(present(yCell))then      
      start1(1) = 1
      count1( 1) = wrLocalnCells
      nferr = nf_put_vara_double(wr_ncid, wrVarIDyCell, start1, count1, yCell)
      endif
 
      if(present(zCell))then      
      start1(1) = 1
      count1( 1) = wrLocalnCells
      nferr = nf_put_vara_double(wr_ncid, wrVarIDzCell, start1, count1, zCell)
      endif
 
      if(present(xVertex))then  
      start1(1) = 1
      count1( 1) = wrLocalnVertices
      nferr = nf_put_vara_double(wr_ncid, wrVarIDxVertex, start1, count1, xVertex)
      endif
 
      if(present(yVertex))then      
      start1(1) = 1
      count1( 1) = wrLocalnVertices
      nferr = nf_put_vara_double(wr_ncid, wrVarIDyVertex, start1, count1, yVertex)
      endif
 
      if(present(zVertex))then      
      start1(1) = 1
      count1( 1) = wrLocalnVertices
      nferr = nf_put_vara_double(wr_ncid, wrVarIDzVertex, start1, count1, zVertex)
      endif
 
      if(present(cellsOnVertex))then
      start2(2) = 1
      count2( 1) = 3
      count2( 2) = wrLocalnVertices
      nferr = nf_put_vara_int(wr_ncid, wrVarIDcellsOnVertex, start2, count2, cellsOnVertex)
      endif
 
   end subroutine write_netcdf_fields
 
 
   subroutine write_netcdf_finalize()
 
      implicit none
 
      include 'netcdf.inc'
 
      integer :: nferr
 
      nferr = nf_close(wr_ncid)
 
   end subroutine write_netcdf_finalize
 
end module write_netcdf
