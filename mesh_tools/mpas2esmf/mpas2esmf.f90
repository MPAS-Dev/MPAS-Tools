! mpas2esmf.F90 - Create an ESMF and SCRIP file from an MPAS grid
!
!

module read_mesh

   contains
   
   subroutine read_mpas_mesh(filename, &
                             nCells, nVertices, maxEdges, &
                             latCell, lonCell, latVertex, lonVertex, nEdgesOnCell, verticesOnCell, areaCell, meshDensity, &
                             sphere_radius)
   
      use netcdf
   
      implicit none
   
      character (len=*), intent(in) :: filename
      integer, intent(inout) :: nCells, nVertices, maxEdges
      double precision, dimension(:), pointer :: latCell, lonCell, latVertex, lonVertex, areaCell, meshDensity
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: verticesOnCell
      double precision, intent(inout) :: sphere_radius
   
   
      integer :: ncid, nCellsID, nVerticesID, maxEdgesID, latCellID, lonCellID, &
                 latVertexID, lonVertexID, nEdgesOnCellID, verticesOnCellID, areaCellID, meshDensityID, status
   
      status = nf90_open(path=trim(filename), mode=nf90_nowrite, ncid=ncid)
      if (status /= nf90_noerr) then
          write(0,*) 'mpas2esmf: Error occured when opening MPAS grid: '//filename
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
   
      status = nf90_inq_dimid(ncid, 'nCells',    nCellsID)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error when getting dimid of 'nCells'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_inq_dimid(ncid, 'nVertices', nVerticesID)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error when getting dimid of 'nVertices'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_inq_dimid(ncid, 'maxEdges',  maxEdgesID)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error when getting dimid of 'maxEdges'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
   
      status = nf90_inquire_dimension(ncid, nCellsID,    len=nCells)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error on inquire dimension of 'nCellsID'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_inquire_dimension(ncid, nVerticesID, len=nVertices)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error on inquire dimension of 'nVerticesID'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_inquire_dimension(ncid, maxEdgesID,  len=maxEdges)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error on inquire dimension of 'maxEdgesID'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
   
      allocate(latCell(nCells))
      allocate(lonCell(nCells))
      allocate(latVertex(nVertices))
      allocate(lonVertex(nVertices))
      allocate(nEdgesOnCell(nCells))
      allocate(verticesOnCell(maxEdges,nCells))
      allocate(areaCell(nCells))
      allocate(meshDensity(nCells))
   
      status = nf90_inq_varid(ncid, 'latCell',   latCellID)
      status = nf90_inquire_dimension(ncid, maxEdgesID,  len=maxEdges)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error on inquire varid of 'latCell'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_inq_varid(ncid, 'lonCell',   lonCellID)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error on inquire varid of 'lonCell'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_inq_varid(ncid, 'latVertex', latVertexID)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error on inquire varid of 'latVertex'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_inq_varid(ncid, 'lonVertex', lonVertexID)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error on inquire varid of 'lonVertex'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_inq_varid(ncid, 'nEdgesOnCell', nEdgesOnCellID)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error on inquire varid of 'nEdgesOnCell'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_inq_varid(ncid, 'verticesOnCell', verticesOnCellID)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error on inquire varid of 'verticesOnCellID'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_inq_varid(ncid, 'areaCell', areaCellID)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error on inquire varid of 'areaCell'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_inq_varid(ncid, 'meshDensity', meshDensityID)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error on inquire varid of 'meshDensity'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
   
      status = nf90_get_var(ncid, latCellID, latCell)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error on get var of 'latCell'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_get_var(ncid, lonCellID, lonCell)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error on get var of 'lonCell'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_get_var(ncid, latVertexID, latVertex)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error on get var of 'latVertex'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_get_var(ncid, lonVertexID, lonVertex)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error on get var of 'lonVertex'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_get_var(ncid, nEdgesOnCellID, nEdgesOnCell)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error on get var of 'nEdgesOnCell'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_get_var(ncid, verticesOnCellID, verticesOnCell)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error on get var of 'verticesOnCell'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_get_var(ncid, areaCellID, areaCell)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error on get var of 'areaCell'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_get_var(ncid, meshDensityID, meshDensity)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error on get var of 'meshDensity'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_get_att(ncid, NF90_GLOBAL, 'sphere_radius', sphere_radius)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error on get attribute of 'sphere_radius'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
   
      status = nf90_close(ncid)
   
   end subroutine read_mpas_mesh
   
end module read_mesh


module write_desc

   contains

   subroutine write_scrip_mesh(filename, &
                              title, datestring, &
                              nCells, maxEdges, &
                              latCell, lonCell, &
                              grid_area, rrfac, latVerticesOnCell, lonVerticesOnCell, grid_imask)
   
      use netcdf
   
      implicit none
   
      character (len=*), intent(in) :: filename
      character (len=*), intent(in) :: title
      character (len=*), intent(in) :: datestring
      integer, intent(inout) :: nCells, maxEdges
      double precision, dimension(:), pointer :: latCell, lonCell, grid_area, rrfac
      double precision, dimension(:,:), pointer :: latVerticesOnCell, lonVerticesOnCell
      integer, dimension(:), pointer :: grid_imask
   
   
      integer :: ncid, grid_sizeID, grid_cornersID, grid_rankID, status
      integer :: grid_areaID, rrfacID, grid_center_latID, grid_center_lonID, &
                 grid_corner_lonID, grid_corner_latID, grid_imaskID, grid_dimsID
      integer, dimension(1) :: id1
      integer, dimension(2) :: id2


      status = nf90_create(trim(filename), NF90_64BIT_OFFSET, ncid)
      if (status /= nf90_noerr) then
          write(0,*) 'mpas2esmf: Error occured in nf90_create for mpas_scrip.nc'
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_def_dim(ncid, 'grid_size', nCells, grid_sizeID)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_def_dim for 'grid_size'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      status = nf90_def_dim(ncid, 'grid_corners', maxEdges, grid_cornersID)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_def_dim for 'grid_corners'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      status = nf90_def_dim(ncid, 'grid_rank', 1, grid_rankID)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_def_dim for 'grid_rank'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      id1(1) = grid_sizeID
      status = nf90_def_var(ncid, 'grid_area', NF90_DOUBLE, id1, grid_areaID)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_def_var for 'grid_area'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      status = nf90_put_att(ncid, grid_areaID, 'units', 'radians^2')
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_put_att for 'units' for grid_area"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      status = nf90_put_att(ncid, grid_areaID, 'long_name', 'area weights')
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_put_att for 'long_name' for grid_area"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      id1(1) = grid_sizeID
      status = nf90_def_var(ncid, 'rrfac', NF90_DOUBLE, id1, rrfacID)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_def_var for 'rrfac'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      status = nf90_put_att(ncid, rrfacID, 'units', 'dimensionless')
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_put_att for 'units' for rrfac"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      status = nf90_put_att(ncid, rrfacID, 'long_name', 'normalized dc')
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_put_att for 'long_name' for rrfac"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      id1(1) = grid_sizeID
      status = nf90_def_var(ncid, 'grid_center_lat', NF90_DOUBLE, id1, grid_center_latID)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_def_var for 'grid_center_lat'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      status = nf90_put_att(ncid, grid_center_latID, 'units', 'radians')
      if (status /= nf90_noerr) then
          write(0,*) 'mpas2esmf: Error occured in nf90_put_att for "units" for grid_center_lat'
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      id1(1) = grid_sizeID
      status = nf90_def_var(ncid, 'grid_center_lon', NF90_DOUBLE, id1, grid_center_lonID)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_def_var for 'grid_center_lon'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      status = nf90_put_att(ncid, grid_center_lonID, 'units', 'radians')
      if (status /= nf90_noerr) then
          write(0,*) 'mpas2esmf: Error occured in nf90_put_att for "units" for grid_center_lon'
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      id2(1) = grid_cornersID
      id2(2) = grid_sizeID
      status = nf90_def_var(ncid, 'grid_corner_lon', NF90_DOUBLE, id2, grid_corner_lonID)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_def_var for 'grid_corner_lon'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      status = nf90_put_att(ncid, grid_corner_lonID, 'units', 'radians')
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_put_att for 'units' for 'grid_corner_lon'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
!      status = nf90_def_var_fill(ncid, grid_corner_lonID, 0, -9999.)

      id2(1) = grid_cornersID
      id2(2) = grid_sizeID
      status = nf90_def_var(ncid, 'grid_corner_lat', NF90_DOUBLE, id2, grid_corner_latID)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_def_var for 'grid_corner_lat'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      status = nf90_put_att(ncid, grid_corner_latID, 'units', 'radians')
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_put_att for 'units' for 'grid_corner_lat'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
!      status = nf90_def_var_fill(ncid, grid_corner_latID, 0, -9999.)

      id1(1) = grid_sizeID
      status = nf90_def_var(ncid, 'grid_imask', NF90_INT, id1, grid_imaskID)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_def_var for 'grid_imask'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
!      status = nf90_def_var_fill(ncid, grid_imaskID, 0, -9999.)

      id1(1) = grid_rankID
      status = nf90_def_var(ncid, 'grid_dims', NF90_INT, id1, grid_dimsID)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_def_var for 'grid_dims'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      status = nf90_put_att(ncid, NF90_GLOBAL, 'title', title)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_put_att for 'title' for global attribute"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      status = nf90_put_att(ncid, NF90_GLOBAL, 'history', 'Created by the mpas2esmf utility, '//trim(datestring))
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_put_att for 'history' for global attribute"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_enddef(ncid)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured when calling enddef for 'mpas_scrip.nc'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_put_var(ncid, grid_areaID, grid_area)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_put_var for 'grid_area'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      status = nf90_put_var(ncid, rrfacID, rrfac)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_put_var for 'rrfac'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      status = nf90_put_var(ncid, grid_center_latID, latCell)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_put_var for 'grid_center_lat'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      status = nf90_put_var(ncid, grid_center_lonID, lonCell)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_put_var for 'grid_center_lon'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      status = nf90_put_var(ncid, grid_corner_latID, latVerticesOnCell)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_put_var for 'grid_corner_lat'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      status = nf90_put_var(ncid, grid_corner_lonID, lonVerticesOnCell)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_put_var for 'grid_corner_lon'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      status = nf90_put_var(ncid, grid_imaskID, grid_imask)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_put_var for 'grid_imask'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      status = nf90_put_var(ncid, grid_dimsID, nCells)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_put_var for 'grid_dims'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_close(ncid) 
      if (status /= nf90_noerr) then
          write(0,*) 'mpas2esmf: Error occured in nf90_close for mpas_scrip.nc'
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
   
   end subroutine write_scrip_mesh


   subroutine write_esmf_mesh(filename, &
                              input_file, title, datestring, &
                              nCells, nVertices, maxEdges, &
                              centerCoords, nodeCoords, elementConn, nEdgesOnCell, &
                              grid_area, rrfac, grid_imask)
   
      use netcdf
   
      implicit none
   
      character (len=*), intent(in) :: filename
      character (len=*), intent(in) :: input_file
      character (len=*), intent(in) :: title
      character (len=*), intent(in) :: datestring
      integer, intent(inout) :: nCells, nVertices, maxEdges
      double precision, dimension(:), pointer :: grid_area, rrfac
      double precision, dimension(:,:), pointer :: centerCoords, nodeCoords
      integer, dimension(:,:), pointer :: elementConn
      integer, dimension(:), pointer :: grid_imask, nEdgesOnCell
   
   
      integer :: ncid, nVerticesID, nCellsID, maxNodePElementID, coordDimID, status
      integer :: nodeCoordsID, elementConnID, numElementConnID, centerCoordsID, &
                 elementAreaID, elementMaskID, elementRefinementRatioID
      integer, dimension(1) :: id1
      integer, dimension(2) :: id2


      status = nf90_create(trim(filename), NF90_64BIT_OFFSET, ncid)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_create for esmf file"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_def_dim(ncid, 'nodeCount', nVertices, nVerticesID)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_def_dim for 'nodeCount'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_def_dim(ncid, 'elementCount', nCells, nCellsID)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_def_dim for 'elementCount'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_def_dim(ncid, 'maxNodePElement', maxEdges, maxNodePElementID)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_def_dim for 'maxNodePElementID'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_def_dim(ncid, 'coordDim', 2, coordDimID)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_def_dim for 'coordDim'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      id2(1) = coordDimID
      id2(2) = nVerticesID
      status = nf90_def_var(ncid, 'nodeCoords', NF90_DOUBLE, id2, nodeCoordsID)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_def_var for 'nodeCoords'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_put_att(ncid, nodeCoordsID, 'units', 'degrees')
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_put_att for 'units' for 'nodeCoordsID'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      id2(1) = maxNodePElementID
      id2(2) = nCellsID
      status = nf90_def_var(ncid, 'elementConn', NF90_INT, id2, elementConnID)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_def_var for 'elementConn'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_put_att(ncid, elementConnID, 'long_name', &
                            'Node Indices that define the element connectivity')
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_put_att for 'long_name' for 'elementConnID'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      status = nf90_put_att(ncid, elementConnID, '_FillValue', -1)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_put_att for '_FillValue' for 'elementConnID'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      id1(1) = nCellsID
      status = nf90_def_var(ncid, 'numElementConn', NF90_BYTE, id1, numElementConnID)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_def_var for 'numElementConn'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      status = nf90_put_att(ncid, numElementConnID, 'long_name', 'Number of nodes per element')
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_put_att for 'long_name' for 'numElementConnID'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      id2(1) = coordDimID
      id2(2) = nCellsID
      status = nf90_def_var(ncid, 'centerCoords', NF90_DOUBLE, id2, centerCoordsID)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_def_var for 'centerCoords'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_put_att(ncid, centerCoordsID, 'units', 'degrees')
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_put_att for 'units' for 'centerCoordsID'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      id1(1) = nCellsID
      status = nf90_def_var(ncid, 'elementArea', NF90_DOUBLE, id1, elementAreaID)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_def_var for 'elementArea'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_put_att(ncid, elementAreaID, 'units', 'radians^2')
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_put_att for 'units' for 'elementAreaID'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_put_att(ncid, elementAreaID, 'long_name', 'area weights')
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_put_att for 'long_name' for 'elementAreaID'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      id1(1) = nCellsID
      status = nf90_def_var(ncid, 'elementRefinementRatio', NF90_DOUBLE, id1, elementRefinementRatioID)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_def_var for 'elementRefinementRatio'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_put_att(ncid, elementRefinementRatioID, 'units', 'dimensionless')
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_put_att for 'units' for 'elementRefinementRatioID'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_put_att(ncid, elementRefinementRatioID, 'long_name', 'normalized dc')
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_put_att for 'long_name' for 'elementRefinementRatioID'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      id1(1) = nCellsID
      status = nf90_def_var(ncid, 'elementMask', NF90_INT, id1, elementMaskID)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_def_var for 'elementMask'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_put_att(ncid, NF90_GLOBAL, 'gridType', 'unstructured')
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_put_att for global attribute 'gridType'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      status = nf90_put_att(ncid, NF90_GLOBAL, 'version', '0.9')
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_put_att for global attribute 'version'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      status = nf90_put_att(ncid, NF90_GLOBAL, 'inputFile', trim(input_file))
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_put_att for global attribute 'inputFile'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      status = nf90_put_att(ncid, NF90_GLOBAL, 'timeGenerated', trim(datestring))
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_put_att for global attribute 'timeGenerated" 
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      status = nf90_put_att(ncid, NF90_GLOBAL, 'history', 'Created by the mpas2esmf utility')
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_put_att for global attribute 'history'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_enddef(ncid)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_enddef esmf file"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_put_var(ncid, nodeCoordsID, nodeCoords)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_put_var for 'nodeCoords'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      status = nf90_put_var(ncid, elementConnID, elementConn)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_put_var for 'elementConn'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      status = nf90_put_var(ncid, numElementConnID, nEdgesOnCell)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_put_var for 'numElementConn'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      status = nf90_put_var(ncid, centerCoordsID, centerCoords)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_put_var for 'centerCoords'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      status = nf90_put_var(ncid, elementAreaID, grid_area)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_put_var for 'elementArea'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      status = nf90_put_var(ncid, elementRefinementRatioID, rrfac)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_put_var for 'elementRefinementRatio'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
      status = nf90_put_var(ncid, elementMaskID, grid_imask)
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_put_var for 'elementMask'"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if

      status = nf90_close(ncid) 
      if (status /= nf90_noerr) then
          write(0,*) "mpas2esmf: Error occured in nf90_close"
          write(0,*) trim(nf90_strerror(status))
          stop 
      end if
   
   end subroutine write_esmf_mesh
   
end module write_desc


program mpas2esmf

   use read_mesh
   use write_desc

   implicit none

   integer :: nCells, nVertices, maxEdges
   integer :: iCell, iVtx
   double precision, dimension(:), pointer :: latCell, lonCell, latVertex, lonVertex, grid_area, meshDensity
   integer, dimension(:), pointer :: nEdgesOnCell
   integer, dimension(:,:), pointer :: verticesOnCell, elementConn
   double precision, dimension(:,:), pointer :: latVerticesOnCell, lonVerticesOnCell, &
                                                centerCoords, nodeCoords 
   double precision, dimension(:), pointer :: rrfac
   double precision :: dcMin, dcMax
   double precision :: sphere_radius
   integer, dimension(:), pointer :: grid_imask
   character (len=1024) :: input_file_name
   character (len=1024) :: title
   character (len=1024) :: datestring

   if (command_argument_count() /= 3) then
      write(0,*) 'Usage: mpas2esmf <mpas grid file> <title> <date>'
      stop
   end if

   call getarg(1,input_file_name)
   call getarg(2,title)
   call getarg(3,datestring)

   write(0,'(A)',advance='no') "mpas2esmf: Reading MPAS mesh ... "

   call read_mpas_mesh(trim(input_file_name), &
                       nCells, nVertices, maxEdges, &
                       latCell, lonCell, latVertex, lonVertex, nEdgesOnCell, verticesOnCell, grid_area, meshDensity, &
                       sphere_radius)

   if (sphere_radius /= 1.0) then
       write(0,*) "ERROR: Please only run mpas2esmf on MPAS files with sphere_radius == 1.0"
       write(0,*) "ERROR: Any grid.nc MPAS are on an unit sphere"
       stop
   end if

   write(0,*) "DONE!"

   maxEdges = maxval(nEdgesOnCell)

   write(0,'(A)',advance='no') "mpas2esmf: Allocating and creating fields for SCRIP file ... "

   allocate(latVerticesOnCell(maxEdges,nCells))
   allocate(lonVerticesOnCell(maxEdges,nCells))
   allocate(centerCoords(2,nCells))
   allocate(nodeCoords(2,nVertices))
   allocate(elementConn(maxEdges,nCells))
   allocate(grid_imask(nCells))
   allocate(rrfac(nCells))

   ! convert meshDensity to normalized dc (cell separation distance, aka dx)
   do iCell=1,nCells
      meshDensity(iCell) = 1./meshDensity(iCell)**0.25
   end do

   ! find min and max dc
   dcMin = meshDensity(1)
   dcMax = meshDensity(1)
   do iCell=2,nCells
     dcMin = min(dcMin, meshDensity(iCell))
     dcMax = max(dcMax, meshDensity(iCell))
   end do

   ! compute CAM refinement factor
   do iCell=1,nCells
     rrfac(iCell) = dcMax/meshDensity(iCell)
   end do

   do iCell=1,nCells
      do iVtx=1,nEdgesOnCell(iCell)
         latVerticesOnCell(iVtx,iCell) = latVertex(verticesOnCell(iVtx,iCell))
         lonVerticesOnCell(iVtx,iCell) = lonVertex(verticesOnCell(iVtx,iCell))
         elementConn(iVtx,iCell) = verticesOnCell(iVtx,iCell)
      end do
      do iVtx=nEdgesOnCell(iCell)+1,maxEdges
         latVerticesOnCell(iVtx,iCell) = latVerticesOnCell(iVtx-1,iCell)
         lonVerticesOnCell(iVtx,iCell) = lonVerticesOnCell(iVtx-1,iCell)
         elementConn(iVtx,iCell)       = -1
      end do
      centerCoords(1,iCell) = lonCell(iCell)
      centerCoords(2,iCell) = latCell(iCell)
   end do

   do iVtx=1,nVertices
      nodeCoords(1,iVtx) = lonVertex(iVtx)
      nodeCoords(2,iVtx) = latVertex(iVtx)
   end do

   grid_imask(:) = 1

   write(0,*) "DONE!"

   write(0,'(A)',advance='no') "mpas2esmf: Writing SCRIP file ... "

   call write_scrip_mesh('mpas_scrip.nc', &
                        title, datestring, &
                        nCells, maxEdges, &
                        latCell, lonCell, &
                        grid_area, rrfac, latVerticesOnCell, lonVerticesOnCell, grid_imask)

   write(0,*) "DONE!"
   write(0,'(A)',advance='no') "mpas2esmf: Creating fields for ESMF files ... "

   do iVtx=1,nVertices
      nodeCoords(1,iVtx) = nodeCoords(1,iVtx) * 90.0 / asin(1.0)
      nodeCoords(2,iVtx) = nodeCoords(2,iVtx) * 90.0 / asin(1.0)
   end do

   do iCell=1,nCells
      centerCoords(1,iCell) = centerCoords(1,iCell) * 90.0 / asin(1.0)
      centerCoords(2,iCell) = centerCoords(2,iCell) * 90.0 / asin(1.0)
   end do

   write(0,*) "DONE!"
   write(0,'(A)',advance='no') "mpas2esmf: Writing ESMF files ... "

   call write_esmf_mesh('mpas_esmf.nc', &
                        input_file_name, title, datestring, &
                        nCells, nVertices, maxEdges, &
                        centerCoords, nodeCoords, elementConn, nEdgesOnCell, &
                        grid_area, rrfac, grid_imask)

   write(0,*) "DONE!"

   deallocate(latCell)
   deallocate(lonCell)
   deallocate(latVertex)
   deallocate(lonVertex)
   deallocate(nEdgesOnCell)
   deallocate(verticesOnCell)

   deallocate(latVerticesOnCell)
   deallocate(lonVerticesOnCell)
   deallocate(centerCoords)
   deallocate(nodeCoords)
   deallocate(elementConn)
   deallocate(grid_imask)
   deallocate(grid_area)

   write(0,*) "mpas2emsf: FINISHED!"

end program mpas2esmf
