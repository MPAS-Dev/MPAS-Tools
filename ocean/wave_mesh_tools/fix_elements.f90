      MODULE fix_elements 

      USE globals, ONLY: rp

      IMPLICIT NONE

      CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      
      SUBROUTINE flag_problem_elements(ne,ect,nepn,el2el,keep_element) 

      IMPLICIT NONE
 
      INTEGER, INTENT(IN) :: ne
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ect
      INTEGER, DIMENSION(:), INTENT(IN) :: nepn
      INTEGER, DIMENSION(:,:), INTENT(IN) :: el2el
      INTEGER, DIMENSION(:), INTENT(INOUT) :: keep_element

      INTEGER :: i,el,nd,led
      INTEGER :: nneigh
      INTEGER :: nconnds

      DO el = 1,ne

        ! Determine how many edge neighbors an element has
        nneigh = 0
        DO led = 1,3
          IF (el2el(el,led) > 0) THEN
            nneigh = nneigh + 1
          ENDIF
        ENDDO

        ! Determine how many vertices are shared by more than one element
        nconnds = 0
        DO i = 1,3
          nd = ect(i,el)
          IF (nepn(nd) > 1) THEN
            nconnds = nconnds + 1
          ENDIF
        ENDDO

        ! Flag elements with no edge neighbors
        IF (nneigh == 0) THEN
          keep_element(el) = 0          
        ENDIF

        ! Flag elements with two boundary edges and three shared verticies
        IF (nneigh == 1 .and. nconnds == 3) THEN
          keep_element(el) = 0
        ENDIF

      ENDDO

      RETURN 
      END SUBROUTINE flag_problem_elements

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE flag_single_element_passages(ne,nn,ect,nbnd,bndn,keep_element)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ne
      INTEGER, INTENT(IN) :: nn
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ect
      INTEGER :: nbnd
      INTEGER, DIMENSION(:), INTENT(IN) :: bndn
      INTEGER, DIMENSION(:), INTENT(INOUT) :: keep_element

      INTEGER :: el,nd
      INTEGER :: bnd_count
      INTEGER, DIMENSION(:), ALLOCATABLE :: bnd_flag

      ! Compute list of boundary node flags
      ALLOCATE(bnd_flag(nn))
      bnd_flag = 0
      DO nd = 1,nbnd
        bnd_flag(bndn(nd)) = 1
      ENDDO

      ! Flag elements with 3 boundary nodes
      DO el = 1,ne
        bnd_count = 0
        DO nd = 1,3
          bnd_count = bnd_count + bnd_flag(ect(nd,el))
        ENDDO
        IF (bnd_count ==3) THEN
          keep_element(el) = 0
        ENDIF
      ENDDO

      RETURN
      END SUBROUTINE flag_single_element_passages

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      
      SUBROUTINE flag_isolated_element_patches(ne,el2el,keep_element)
     
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ne
      INTEGER, DIMENSION(:,:), INTENT(IN) :: el2el
      INTEGER, DIMENSION(:), INTENT(INOUT) :: keep_element

      INTEGER :: i,j,k,el
      INTEGER :: el_next,el_check
      INTEGER :: nel_found,found
      INTEGER :: pos
      INTEGER, DIMENSION(:), ALLOCATABLE :: el_found

      ALLOCATE(el_found(ne))
      DO el = 1,ne
      
        el_next = el
        nel_found = 1
        pos = 1
        el_found(nel_found) = el_next
 search:DO i = 1,200
  elems:  DO j = 1,3

            ! find elements connected to el
            el_check = el2el(el_next,j)             

            IF (el_check == 0 ) THEN
              CYCLE elems
            ENDIF

            IF (keep_element(el_check) == 0) THEN
              CYCLE elems
            ENDIF

            ! check if the element has been added to the queue 
            found = 0                            
      quen: DO k = 1,nel_found
              IF (el_found(k) == el_check) THEN
                found= 1
                EXIT quen
              ENDIF
            ENDDO quen

            ! add element to the queue if it hasn't been already
            IF (found == 0 ) THEN         
              nel_found = nel_found + 1               
              el_found(nel_found) = el_check
            ENDIF
  
          ENDDO elems

          ! exit when all elements have been added
          IF (pos == nel_found) THEN            
            EXIT search
          ENDIF      

          ! get ready to check the next element
          pos = pos + 1                 
          el_next = el_found(pos)
 
        ENDDO search

        ! Flag elements
        IF (nel_found < 200) THEN
          DO i = 1,nel_found
            keep_element(el_found(i)) = 0
          ENDDO
        ENDIF

      ENDDO

      END SUBROUTINE flag_isolated_element_patches

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE fix_single_node_connections_across_islands(nn,nbnd,bndn,nbed,bedn,nepn,epn,ged2nn,nd_flag,keep_element)
 
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nn
      INTEGER, INTENT(IN) :: nbnd
      INTEGER, DIMENSION(:), INTENT(IN) :: bndn
      INTEGER, INTENT(IN) :: nbed
      INTEGER, DIMENSION(:), INTENT(IN) :: bedn
      INTEGER, DIMENSION(:), INTENT(IN) :: nepn
      INTEGER, DIMENSION(:,:), INTENT(IN) :: epn
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ged2nn
      INTEGER, DIMENSION(:), ALLOCATABLE,  INTENT(OUT) :: nd_flag
      INTEGER, DIMENSION(:), INTENT(INOUT) :: keep_element

      INTEGER :: i,j
      INTEGER :: nd,ged

      ! Find nodes that are shared by more than two boundary edges
      ALLOCATE(nd_flag(nn))
      nd_flag = 0
      DO i = 1,nbnd
        nd = bndn(i)
        DO j = 1,nbed
          ged = bedn(j)
          IF (ged2nn(1,ged) == nd .OR. ged2nn(2,ged) == nd) THEN
            nd_flag(nd) = nd_flag(nd) + 1
          ENDIF
        ENDDO
      ENDDO

      ! Eliminate elements associated with these nodes
      DO nd = 1,nn
        IF (nd_flag(nd) > 2) THEN
          DO i = 1,nepn(nd)
            keep_element(epn(i,nd)) = 0
          ENDDO
        ENDIF
      ENDDO

      RETURN
      END SUBROUTINE fix_single_node_connections_across_islands

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE flag_single_element_islands(ned,nbed,bedn,ged2nn,xy,nfill_elements,fill_elements)     
   
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ned
      INTEGER, INTENT(IN) :: nbed
      INTEGER, DIMENSION(:), INTENT(IN) :: bedn
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ged2nn
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: xy
      INTEGER, INTENT(OUT) :: nfill_elements
      INTEGER, DIMENSION(:,:), INTENT(OUT) :: fill_elements
      
      INTEGER :: i,j,ed,ged
      INTEGER :: nd_start,nd_end,nd_next
      INTEGER :: ed_count
      INTEGER :: tmp
      INTEGER, DIMENSION(:), ALLOCATABLE :: ed_flag 
      INTEGER, DIMENSION(:), ALLOCATABLE :: ed_skip 
      INTEGER, DIMENSION(5) :: ed_found,nd_found          
      REAL(rp) :: x1,x2,x3,y1,y2,y3
      REAL(rp) :: area

      ALLOCATE(ed_flag(ned))
      ALLOCATE(ed_skip(ned))
      ed_skip = 0
      nfill_elements = 0
 edge:DO ed = 1,nbed
        ed_flag = 0
        ged = bedn(ed)
        IF (ed_skip(ged) == 1) THEN
          CYCLE edge
        ENDIF

        ed_flag(ged) = 1
        nd_start = ged2nn(1,ged)
        nd_end = ged2nn(2,ged)
        nd_next = nd_start
        nd_found = 0
        ed_found = 0
        nd_found(1) = nd_start
        ed_found(1) = ged
        ed_count = 1 

        ! Try to find 4 consecutive edges before wrapping around
   link:DO j = 1,4
   search:DO i = 1,nbed
            ged = bedn(i)
            IF (ed_flag(ged) == 1) THEN
              CYCLE search 
            ENDIF

            IF (ged2nn(1,ged) == nd_next) THEN
              nd_next = ged2nn(2,ged)
              ed_flag(ged) = 1
              ed_count = ed_count + 1
              nd_found(ed_count) = nd_next
              ed_found(ed_count) = ged
            ELSE IF (ged2nn(2,ged) == nd_next) THEN
              nd_next = ged2nn(1,ged)
              ed_flag(ged) = 1
              ed_count = ed_count + 1
              nd_found(ed_count) = nd_next
              ed_found(ed_count) = ged
            ENDIF

            IF (ed_flag(ged) == 1) THEN
              EXIT search
            ENDIF
            IF (nd_next == nd_end) THEN
              EXIT link
            ENDIF

          ENDDO search
        ENDDO link
    
        ! Fill element if only three edges/nodes were found
        IF (ed_count == 3) THEN
          DO i = 1,3
            ed_skip(ed_found(i)) = 1
          ENDDO
          
          ! Check to make sure node numbering is CCW
          x1 = xy(1,nd_found(1)) 
          x2 = xy(1,nd_found(2))
          x3 = xy(1,nd_found(3))
          y1 = xy(2,nd_found(1))
          y2 = xy(2,nd_found(2))
          y3 = xy(2,nd_found(3))
          area = (x1-x3)*(y2-y3)+(x3-x2)*(y1-y3)
          IF (area < 0d0) THEN
            tmp = nd_found(2)            
            nd_found(2) = nd_found(3)
            nd_found(3) = tmp
          ENDIF

          ! Add element 
          nfill_elements = nfill_elements + 1
          DO i = 1,3
            fill_elements(i,nfill_elements) = nd_found(i)
          ENDDO
        ENDIF
 
      ENDDO edge

      RETURN
      END SUBROUTINE flag_single_element_islands
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE create_new_ect(ne,ect,keep_element,ne_new,ect_new)
  
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ne
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ect
      INTEGER, DIMENSION(:), INTENT(IN) :: keep_element
      INTEGER, INTENT(OUT) :: ne_new
      INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: ect_new

      INTEGER :: el,i

      ne_new = 0
      DO el = 1,ne
        IF (keep_element(el) == 1) THEN
          ne_new = ne_new + 1
        ENDIF
      ENDDO 

      ALLOCATE(ect_new(3,ne_new))
      ne_new = 0
      DO el = 1,ne
        IF (keep_element(el) == 1) THEN
          ne_new = ne_new + 1
          DO i = 1,3
            ect_new(i,ne_new) = ect(i,el)
          ENDDO
        ENDIF
      ENDDO 
      PRINT*, ne,ne_new

      RETURN
      END SUBROUTINE create_new_ect

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE add_elements_to_ect(ne,ect,nfill_elements,fill_elements)

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: ne
      INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: ect
      INTEGER, INTENT(IN) :: nfill_elements
      INTEGER, DIMENSION(:,:), INTENT(IN) :: fill_elements

      INTEGER :: i,j
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ect_tmp

      ALLOCATE(ect_tmp(3,2*ne))
     
      DO i = 1,ne
        DO j = 1,3
          ect_tmp(j,i) = ect(j,i)
        ENDDO
      ENDDO

      DEALLOCATE(ect)
      ALLOCATE(ect(3,ne+nfill_elements))

      DO i = 1,ne
        DO j = 1,3
          ect(j,i) = ect_tmp(j,i)
        ENDDO
      ENDDO

     
      DO i = 1,nfill_elements
        DO j = 1,3
          ect(j,ne+i) = fill_elements(j,i)
        ENDDO
      ENDDO
      PRINT*, ne,ne+nfill_elements

      ne = ne + nfill_elements

      RETURN
      END SUBROUTINE add_elements_to_ect

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
     
      SUBROUTINE remove_unconnected_nodes(nn,nepn,xy,ne,ect,nn_new,xy_new,ect_new,depth,depth_new)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nn
      INTEGER, DIMENSION(:), INTENT(IN) :: nepn
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: xy
      INTEGER, INTENT(IN) :: ne
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ect
      INTEGER, INTENT(OUT) :: nn_new
      REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: xy_new 
      INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: ect_new 
      REAL(rp), DIMENSION(:), INTENT(IN), OPTIONAL :: depth
      REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: depth_new
      

      INTEGER :: nd,el
      INTEGER, DIMENSION(:), ALLOCATABLE :: keep_nd
      INTEGER :: new_node_number

      ! Find nodes that aren't connencted to any elements
      ALLOCATE(keep_nd(nn))
      nn_new = 0
      DO nd = 1,nn
        IF (nepn(nd) == 0) THEN
          keep_nd(nd) = 0
        ELSE
          nn_new = nn_new + 1
          keep_nd(nd) = nn_new
        ENDIF
      ENDDO

      ! Remove nodes from node coordinate array
      ALLOCATE(xy_new(2,nn_new))
      IF (PRESENT(depth)) THEN
        ALLOCATE(depth_new(nn_new))
      ENDIF
      DO nd = 1,nn
        new_node_number = keep_nd(nd)
        IF (new_node_number /= 0) THEN
          xy_new(1,new_node_number) = xy(1,nd)
          xy_new(2,new_node_number) = xy(2,nd)
          IF (PRESENT(depth)) THEN
            depth_new(new_node_number) = depth(nd)
          ENDIF
        ENDIF
      ENDDO 

      ! Renumber element connectivity table   
      ALLOCATE(ect_new(3,ne))
      DO el = 1,ne
        DO nd = 1,3
          new_node_number = keep_nd(ect(nd,el))
          IF (new_node_number == 0) THEN
            PRINT*, "Error: Node connected to element was removed"
            STOP
          ENDIF
          ect_new(nd,el) = new_node_number
        ENDDO
      ENDDO
     
      PRINT*, nn,nn_new

      RETURN
      END SUBROUTINE remove_unconnected_nodes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE fix_element_node_orientation(ne,xy,ect)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ne
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: xy
      INTEGER, DIMENSION(:,:), INTENT(INOUT) :: ect

      INTEGER :: el
      INTEGER :: tmp
      REAL(rp) :: x1,x2,x3
      REAL(rp) :: y1,y2,y3
      REAL(rp) :: x1x3,x2x1
      REAL(rp) :: area

      DO el = 1,ne
        x1 = xy(1,ect(1,el))
        y1 = xy(2,ect(1,el))
        x2 = xy(1,ect(2,el))
        y2 = xy(2,ect(2,el))
        x3 = xy(1,ect(3,el))
        y3 = xy(2,ect(3,el))

        x1x3 = x1-x3
        IF (abs(x1x3) > 180d0) THEN
          x1x3 = x1x3 - SIGN(360d0,x1x3)
        ENDIF

        x2x1 = x2-x1
        IF (abs(x2-x1) > 180d0) THEN
          x2x1 = x2x1 - SIGN(360d0,x2x1)
        ENDIF

        area = (y2-y1)*(x1x3)+(y3-y1)*(x2x1)

        IF (area < 0d0) THEN
          tmp = ect(2,el)
          ect(2,el) = ect(3,el)
          ect(3,el) = tmp 
        ENDIF
        
      ENDDO

      RETURN
      END SUBROUTINE fix_element_node_orientation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      END MODULE fix_elements
