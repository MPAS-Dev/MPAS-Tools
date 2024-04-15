      PROGRAM cull_waves_mesh

      USE netcdf
      USE in_cell_mod, ONLY: in_cell_init, pt_in_cell, check, pi, lonlat2xyz, R
      USE globals, ONLY: rp, nverts
      USE write_vtk, ONLY: write_vtk_file
      USE read_write_gmsh, ONLY: write_header, write_nodes, write_elements
      USE edge_connectivity_mod, ONLY: elements_per_node, find_edge_pairs, find_interior_edges, &
                                       find_element_edges, find_neighbor_elements, find_boundary_nodes
      USE fix_elements, ONLY: fix_single_node_connections_across_islands, &
                              fix_element_node_orientation, &
                              create_new_ect, &
                              remove_unconnected_nodes, &
                              add_elements_to_ect

      IMPLICIT NONE

      INTEGER :: waves_ncid
      INTEGER :: ocean_ncid
      INTEGER :: nCells_dimid, nVertices_dimid
      INTEGER :: cellsOnVertex_varid, lonCell_varid, latCell_varid
      INTEGER :: bottomDepth_varid

      CHARACTER(100) :: waves_mesh_file
      CHARACTER(100) :: waves_mesh_culled_vtk
      CHARACTER(100) :: waves_mesh_culled_gmsh
      CHARACTER(100) :: ocean_mesh_file

      INTEGER :: ne_waves, nn_waves
      INTEGER :: ne_new, nn_new
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ect_waves
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ect_new
      REAL(rp), DIMENSION(:), ALLOCATABLE :: x_waves, y_waves
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: xy_waves
      REAL(rp), DIMENSION(:), ALLOCATABLE :: x_new, y_new
      REAL(rp), DIMENSION(:), ALLOCATABLE :: depth_waves
      REAL(rp), DIMENSION(:), ALLOCATABLE :: depth_new
      REAL(rp) :: xy(2)
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: xyz
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: xy_new
      REAL(rp) :: x1, x2, x3
      REAL(rp) :: y1, y2, y3
      REAL(rp) :: x2x1, x1x3, area

      INTEGER :: ncells_ocean
      REAL(rp), DIMENSION(:), ALLOCATABLE :: depth_ocean

      INTEGER :: el, nd, ged
      INTEGER :: i,j
      INTEGER :: in_cell
      INTEGER :: tmp
      INTEGER :: mnepn
      INTEGER :: ned
      INTEGER, DIMENSION(:), ALLOCATABLE :: nepn
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: epn
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ged2el
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ged2nn
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ged2led      
      INTEGER, DIMENSION(:), ALLOCATABLE :: el_type
      INTEGER :: nied
      INTEGER, DIMENSION(:), ALLOCATABLE :: iedn
      INTEGER, DIMENSION(:), ALLOCATABLE :: ed_type
      INTEGER, DIMENSION(:), ALLOCATABLE :: recv_edge
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: el2ged
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: el2el
      INTEGER :: nbed 
      INTEGER, DIMENSION(:), ALLOCATABLE :: bedn
      INTEGER :: nbnd
      INTEGER, DIMENSION(:), ALLOCATABLE :: bndn
      INTEGER, DIMENSION(:), ALLOCATABLE :: new_node_numbers
      INTEGER, DIMENSION(:), ALLOCATABLE :: nd_flag
      INTEGER, DIMENSION(:), ALLOCATABLE :: elem_flag
      INTEGER :: nfill_elements
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: fill_elements
      INTEGER, DIMENSION(:), ALLOCATABLE  :: ed_flag
      INTEGER :: nd_count
      INTEGER :: nd_start, nd_next
      INTEGER, DIMENSION(:), ALLOCATABLE :: nd_found

      INTEGER, DIMENSION(:), ALLOCATABLE :: wave_node_in_ocean
      INTEGER, DIMENSION(:), ALLOCATABLE :: wave_elements_keep

      NAMELIST / inputs / waves_mesh_file, ocean_mesh_file
      NAMELIST / output / waves_mesh_culled_vtk, waves_mesh_culled_gmsh
 
      OPEN(UNIT=15, FILE='cull_waves_mesh.nml')
      READ(UNIT=15, NML=inputs)
      READ(UNIT=15, NML=output)
      CLOSE(15)
 
      PRINT*, 'waves_mesh_file = ', TRIM(ADJUSTL(waves_mesh_file))
      PRINT*, 'ocean_mesh_file = ', TRIM(ADJUSTL(ocean_mesh_file))
      PRINT*, 'waves_mesh_file_vtk = ', TRIM(ADJUSTL(waves_mesh_culled_vtk))
      PRINT*, 'waves_mesh_file_gmsh = ', TRIM(ADJUSTL(waves_mesh_culled_gmsh))
      CALL check(NF90_OPEN(TRIM(ADJUSTL(waves_mesh_file)), NF90_NOWRITE, waves_ncid)) 
      CALL check(NF90_INQ_DIMID(waves_ncid, 'nCells', nCells_dimid))
      CALL check(NF90_INQ_DIMID(waves_ncid, 'nVertices', nVertices_dimid))

      CALL check(NF90_INQUIRE_DIMENSION(waves_ncid, nCells_dimid, len=nn_waves))
      CALL check(NF90_INQUIRE_DIMENSION(waves_ncid, nVertices_dimid, len=ne_waves))

      CALL check(NF90_INQ_VARID(waves_ncid, 'cellsOnVertex', cellsOnVertex_varid))
      CALL check(NF90_INQ_VARID(waves_ncid, 'lonCell', lonCell_varid))
      CALL check(NF90_INQ_VARID(waves_ncid, 'latCell', latCell_varid))

      ALLOCATE(x_waves(nn_waves), y_waves(nn_waves))
      CALL check(NF90_GET_VAR(waves_ncid, lonCell_varid, x_waves))
      CALL check(NF90_GET_VAR(waves_ncid, latCell_varid, y_waves))
      ALLOCATE(ect_waves(3, ne_waves))
      CALL check(NF90_GET_VAR(waves_ncid, cellsOnVertex_varid, ect_waves))
      CALL check(NF90_CLOSE(waves_ncid))

      CALL check(NF90_OPEN(ocean_mesh_file, NF90_NOWRITE, ocean_ncid))
      CALL check(NF90_INQ_DIMID(ocean_ncid, 'nCells', nCells_dimid))
      CALL check(NF90_INQUIRE_DIMENSION(ocean_ncid, nCells_dimid, len=ncells_ocean))
      CALL check(NF90_INQ_VARID(ocean_ncid, 'bottomDepth', bottomDepth_varid))
      ALLOCATE(depth_ocean(ncells_ocean))
      CALL check(NF90_GET_VAR(ocean_ncid, bottomDepth_varid, depth_ocean))
      CALL check(NF90_CLOSE(ocean_ncid))


      
      CALL in_cell_init(TRIM(ADJUSTL(ocean_mesh_file)))
  
      ALLOCATE(xy_waves(2,nn_waves))
      DO i = 1,nn_waves
        xy_waves(1,i) = x_waves(i)
        xy_waves(2,i) = y_waves(i)
      ENDDO


      ALLOCATE(wave_node_in_ocean(nn_waves))
      wave_node_in_ocean = 0
      ALLOCATE(wave_elements_keep(ne_waves))
      wave_elements_keep = 0
      ALLOCATE(depth_waves(nn_waves))
      depth_waves = -2d0

      ! Determine which elements to keep based on nodes found in ocean mesh
      ne_new = 0
elems:DO el = 1,ne_waves

         IF (mod(el,1000) == 0) THEN
           PRINT*, el
         ENDIF

   nds:  DO i = 1,3
         
          nd = ect_waves(i,el)

          IF (nd == 0) THEN
            wave_elements_keep(el) = 0
            EXIT nds
          ENDIF

          IF (wave_node_in_ocean(nd) == 1) THEN
            wave_elements_keep(el) = 1
            ne_new = ne_new + 1
            CYCLE nds
          ENDIF

          xy(1) = x_waves(nd)
          xy(2) = y_waves(nd)

          CALL pt_in_cell(xy,in_cell)

          IF (in_cell > 0) THEN
            wave_node_in_ocean(nd) = 1
            depth_waves(nd) = depth_ocean(in_cell)
            IF (wave_elements_keep(el) == 0) THEN
              wave_elements_keep(el) = 1
              ne_new = ne_new + 1
            ENDIF
          ENDIF    

        ENDDO nds
      ENDDO elems


      ! Create new element connectivity table
      CALL create_new_ect(ne_waves,ect_waves,wave_elements_keep,ne_new,ect_new)
      PRINT*, ""
      PRINT*, "Elements in original waves mesh: ", ne_waves
      PRINT*, "Elements in culled waves mesh:", ne_new

      ! Find number of elements per node
      ALLOCATE(el_type(ne_new))
      el_type = 1
      CALL elements_per_node(ne_new,nn_waves,nverts,el_type,ect_new,nepn,mnepn,epn)

      ! Build new coordinate list
      DEALLOCATE(ect_waves)
      ALLOCATE(ect_waves(3,ne_new))
      ect_waves = ect_new
      CALL remove_unconnected_nodes(nn_waves,nepn,xy_waves,ne_new,ect_waves,nn_new,xy_new,ect_new,depth_waves,depth_new)
      PRINT*, "Nodes in original waves mesh: ", nn_waves
      PRINT*, "Nodes in culled waves mesh:", nn_new


      ! Convert to Cartesian coordinates for vtk output
      ALLOCATE(xyz(3,nn_new))
      DO i = 1,nn_new
        CALL lonlat2xyz(xy_new(1,i),xy_new(2,i),R,xyz(1,i),xyz(2,i),xyz(3,i))
      ENDDO

      ! Convert lon,lat to degrees and adjust lon to -180,180 range
      DO i = 1,nn_new
        xy_new(1,i) = xy_new(1,i)*180d0/pi
        xy_new(2,i) = xy_new(2,i)*180d0/pi
        IF (xy_new(1,i) > 180d0) THEN
          xy_new(1,i) = xy_new(1,i) - 360d0
        ENDIF
      ENDDO

      ! Compute edge connectivity and indentify single node connections
      DEALLOCATE(nepn,epn)
      CALL elements_per_node(ne_new,nn_new,nverts,el_type,ect_new,nepn,mnepn,epn)
      CALL find_edge_pairs(ne_new,nverts,el_type,ect_new,nepn,epn,ned,ged2el,ged2nn,ged2led)
      CALL find_interior_edges(ned,ged2el,nied,iedn,ed_type,recv_edge,nbed,bedn)
      CALL find_element_edges(ne_new,ned,ged2el,ged2led,el2ged)
      CALL find_neighbor_elements(ne_new,ned,ged2el,ged2led,el2el)
      CALL find_boundary_nodes(nn_new,nbed,bedn,ged2nn,nbnd,bndn)
      ALLOCATE(elem_flag(ne_new))
      CALL fix_single_node_connections_across_islands(nn_new,nbnd,bndn,nbed,bedn,nepn,epn,ged2nn,nd_flag,elem_flag)

      ! Find nodes in order to add elements to correct single node connections
      PRINT*, "Fixing single node connections"
      ALLOCATE(fill_elements(3,ne_new))
      ALLOCATE(nd_found(nbed))
      ALLOCATE(ed_flag(ned))
      nfill_elements = 0
      DO nd = 1,nn_new
        IF (nd_flag(nd) > 2) THEN
          nd_start = nd
          nd_next = nd
          nd_count = 1
          nd_found(nd_count) = nd_next
          ed_flag = 0
    link: DO 
    search: DO i = 1,nbed

              ged = bedn(i)
              IF (ed_flag(ged) == 1) THEN
                CYCLE search
              ENDIF

              IF (ged2nn(1,ged) == nd_next) THEN
                nd_next = ged2nn(2,ged)
                ed_flag(ged) = 1
                nd_count = nd_count + 1
                nd_found(nd_count) = nd_next
              ELSE IF (ged2nn(2,ged) == nd_next) THEN
                nd_next = ged2nn(1,ged)
                ed_flag(ged) = 1
                nd_count = nd_count + 1
                nd_found(nd_count) = nd_next
              ENDIF

              IF (ed_flag(ged) == 1) THEN
                IF (nd_next == nd_start) THEN
                  EXIT link
                ENDIF
                EXIT search
              ENDIF

            ENDDO search
          ENDDO link

          nfill_elements = nfill_elements + 1
          fill_elements(1,nfill_elements) = nd_found(1)
          fill_elements(2,nfill_elements) = nd_found(2)    
          fill_elements(3,nfill_elements) = nd_found(nd_count-1)    
          PRINT*, "Adding element: ", (fill_elements(j,nfill_elements), j = 1,3)

        ENDIF 
      ENDDO

      CALL add_elements_to_ect(ne_new,ect_new,nfill_elements,fill_elements)
      CALL fix_element_node_orientation(ne_new,xy_new,ect_new)

     

      ! Write mesh
      CALL write_vtk_file(TRIM(ADJUSTL(waves_mesh_culled_vtk)),nn_new,xyz,ne_new,ect_new,depth_new)

      CALL write_header(TRIM(ADJUSTL(waves_mesh_culled_gmsh)))
      CALL write_nodes(nn_new,xy_new,depth_new)
      CALL write_elements(ne_new,ect_new)

      END PROGRAM cull_waves_mesh
