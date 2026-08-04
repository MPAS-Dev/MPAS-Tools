      MODULE edge_connectivity_mod
 
      USE globals, ONLY: rp
      
      IMPLICIT NONE

      INTEGER, SAVE :: nbed
      INTEGER, SAVE, DIMENSION(:), ALLOCATABLE :: bedn ! allocated in find_interior_edges, deallocated in find_flow_edges

      CONTAINS
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       

      
      SUBROUTINE elements_per_node(ne,nn,nverts,el_type,vct,nepn,mnepn,epn)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ne
      INTEGER, INTENT(IN) :: nn
      INTEGER, DIMENSION(:), INTENT(IN) :: nverts
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      INTEGER, DIMENSION(:,:), INTENT(IN) :: vct
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: nepn
      INTEGER, INTENT(OUT) :: mnepn
      INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: epn
      
      INTEGER :: el
      INTEGER :: nd
      INTEGER :: alloc_status      
      INTEGER :: nvert  
      INTEGER :: n           
      
      

      ! count the number of elements per node
      
      ALLOCATE(nepn(nn), STAT=alloc_status)
      IF (alloc_status /= 0) THEN
        PRINT*, "Allocation error"
        STOP
      ENDIF         

      nepn(:) = 0
      DO el = 1,ne
        nvert = nverts(el_type(el))
        DO nd = 1,nvert
          n = vct(nd,el)
          nepn(n) = nepn(n) + 1
        ENDDO        
      ENDDO

      
      ! find the maximum number of elements per node

      mnepn = 0
      DO nd = 1,nn
        IF(nepn(nd) > mnepn) THEN
          mnepn = nepn(nd)
        ENDIF
      ENDDO
      
      

      ! find the elements associated with each node
      
      ALLOCATE(epn(mnepn,nn), STAT=alloc_status)
      IF (alloc_status /= 0) THEN
        PRINT*, "Allocation error"
        STOP
      ENDIF         

      nepn(:) = 0
      DO el = 1,ne
        nvert = nverts(el_type(el))
        DO nd = 1,nvert
          n = vct(nd,el)
          nepn(n) = nepn(n) + 1
          epn(nepn(n),n) = el
        ENDDO
      ENDDO      
      
      RETURN
      END SUBROUTINE elements_per_node
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 


      SUBROUTINE find_edge_pairs(ne,nverts,el_type,vct,nepn,epn,ned,ged2el,ged2nn,ged2led)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ne
      INTEGER, DIMENSION(:), INTENT(IN) :: nverts
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      INTEGER, DIMENSION(:,:), INTENT(IN) :: vct
      INTEGER, DIMENSION(:), INTENT(IN) :: nepn
      INTEGER, DIMENSION(:,:), INTENT(IN) :: epn
      INTEGER, INTENT(OUT) :: ned      
      INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: ged2el
      INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: ged2nn
      INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: ged2led
      
      INTEGER :: el
      INTEGER :: nd
      INTEGER :: alloc_status      
      INTEGER :: el1,el2
      INTEGER :: led1,led2
      INTEGER :: nvert1,nvert2
      INTEGER :: n1ed1,n2ed1
      INTEGER :: n1ed2,n2ed2
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ged2nn_temp, ged2el_temp, ged2led_temp
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: edflag
      
      
      ALLOCATE(ged2nn_temp(2,4*ne),ged2el_temp(2,4*ne),ged2led_temp(2,4*ne), STAT = alloc_status)      
      IF (alloc_status /= 0) THEN
        PRINT*, "Allocation error"
        STOP
      ENDIF             
      
      ALLOCATE(edflag(4,ne), STAT = alloc_status)
      IF (alloc_status /= 0) THEN
        PRINT*, "Allocation error"
        STOP
      ENDIF       
      
      
      ned = 0
      edflag = 0
      ged2nn_temp = 0
      ged2led_temp = 0
      ged2el_temp = 0
      
   elem1: DO el1 = 1,ne ! loop through trial elements
   
            nvert1 = nverts(el_type(el1))

 local_ed1: DO led1 = 1,nvert1 ! loop through trial edges

              IF(edflag(led1,el1) == 1) THEN ! skip if edge has already been flagged
                CYCLE local_ed1
              ENDIF

              ned = ned + 1 ! increment edge number

              n1ed1 = vct(mod(led1+0,nvert1)+1,el1) ! find nodes on trial edge
              n2ed1 = vct(mod(led1+1,nvert1)+1,el1)
              
              ged2nn_temp(1,ned) = n1ed1 ! set nodes on global edge # ned
              ged2nn_temp(2,ned) = n2ed1

              ged2led_temp(1,ned) = led1 ! set local edge number of first element sharing the edge
              ged2el_temp(1,ned) = el1   ! set the first element that shares the edge

              edflag(led1,el1) = 1  ! flag the edge so it is not repeated

       elem2: DO el = 1,nepn(n1ed1) ! loop through test elements (only those that contain node n1ed1)

                el2 = epn(el,n1ed1) ! choose a test element that contains node n1ed1 
                
                IF(el2 == el1) THEN ! skip if the test element is the same as the trial element
                  CYCLE elem2
                ENDIF
                
                nvert2 = nverts(el_type(el2))

     local_ed2: DO led2 = 1,nvert2 ! loop through local test edge numbers
                  
                  n1ed2 = vct(MOD(led2+0,nvert2)+1,el2) ! find nodes on test edge
                  n2ed2 = vct(MOD(led2+1,nvert2)+1,el2)

                  IF(((n1ed1 == n1ed2) .AND. (n2ed1 == n2ed2)) .OR. & ! check if nodes on trial edge matches test edge
                     ((n1ed1 == n2ed2) .AND. (n2ed1 == n1ed2))) THEN

                     ged2led_temp(2,ned) = led2 ! set local edge number of second element sharing the edge
                     ged2el_temp(2,ned) = el2   ! set the second element that shares the edge
                     
                     edflag(led2,el2) = 1       ! flag the edge so it is not repeated

                     EXIT elem2          
                  ENDIF

                ENDDO local_ed2

              ENDDO elem2

            ENDDO local_ed1
      
          ENDDO elem1

      ALLOCATE(ged2nn(2,ned),ged2el(2,ned),ged2led(2,ned), STAT=alloc_status)
      IF (alloc_status /= 0) THEN
        PRINT*, "Allocation error"
        STOP
      ENDIF          

      ged2nn(:,1:ned) = ged2nn_temp(:,1:ned)
      ged2el(:,1:ned) = ged2el_temp(:,1:ned)
      ged2led(:,1:ned) = ged2led_temp(:,1:ned)      
      
      RETURN
      END SUBROUTINE find_edge_pairs

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 


      SUBROUTINE find_interior_edges(ned,ged2el,nied,iedn,ed_type,recv_edge,nbou_edges,bou_edges)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ned
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ged2el
      INTEGER, INTENT(OUT) :: nied
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: iedn
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: ed_type
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: recv_edge
      INTEGER, INTENT(OUT), OPTIONAL :: nbou_edges
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: bou_edges

      INTEGER :: el
      INTEGER :: nd
      INTEGER :: alloc_status      
      INTEGER :: ged,ed
      INTEGER :: el1,el2
      INTEGER, DIMENSION(:), ALLOCATABLE :: iedn_temp
      
      ALLOCATE(iedn_temp(ned), STAT = alloc_status)
      IF (alloc_status /= 0) THEN
        PRINT*, "Allocation error"
        STOP
      ENDIF       
      
      ALLOCATE(recv_edge(ned),ed_type(ned), STAT = alloc_status)
      IF (alloc_status /= 0) THEN
        PRINT*, "Allocation error"
        STOP
      ENDIF      

      IF (.NOT. ALLOCATED(bedn)) THEN
        ALLOCATE(bedn(ned), STAT = alloc_status)
        IF (alloc_status /= 0) THEN
          PRINT*, "Allocation error"
          STOP
        ENDIF      
      ENDIF
      
      recv_edge = 1 
      ed_type = -999

      nied = 0
      nbed = 0
      iedn_temp(:) = 0
      DO ged = 1,ned
        el1 = ged2el(1,ged)
        el2 = ged2el(2,ged)
        IF ((el1 /= 0) .AND. (el2 /= 0)) THEN
          nied = nied + 1
          iedn_temp(nied) = ged
          recv_edge(ged) = 0   
          ed_type(ged) = 0
        ELSE 
          nbed = nbed + 1
          bedn(nbed) = ged
        ENDIF        
      ENDDO      
      
      ALLOCATE(iedn(nied), STAT=alloc_status)
      IF (alloc_status /= 0) THEN
        PRINT*, "Allocation error"
        STOP
      ENDIF        
      
      iedn(1:nied) = iedn_temp(1:nied)   
      
      ! ed_type : interior edges       = 0
      !           open edges           = 1
      !           no normal flow edges = 10
      !           specified flow edges = 12
      !           recieve edges        = -1
      
      IF (PRESENT(nbou_edges).AND.PRESENT(bou_edges)) THEN      
        ALLOCATE(bou_edges(nbed))
        DO ed = 1,nbed
          bou_edges(ed) = bedn(ed)
        ENDDO      
        nbou_edges = nbed
      ENDIF
      
      RETURN
      END SUBROUTINE find_interior_edges
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE find_open_edges(nope,obseg,obnds,ged2nn,nobed,obedn,ed_type,recv_edge)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: nope
      INTEGER, DIMENSION(:), INTENT(IN) :: obseg
      INTEGER, DIMENSION(:,:), INTENT(IN) :: obnds
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ged2nn
      INTEGER, INTENT(OUT) :: nobed
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: obedn
      INTEGER, DIMENSION(:), INTENT(INOUT) :: ed_type
      INTEGER, DIMENSION(:), INTENT(INOUT) :: recv_edge
      
      
      INTEGER :: el
      INTEGER :: alloc_status      
      INTEGER :: seg,nd,ed,ged
      INTEGER :: n1bed,n2bed
      INTEGER :: n1ed2,n2ed2
      INTEGER, DIMENSION(:), ALLOCATABLE :: obedn_temp
      
      ALLOCATE(obedn_temp(nbed), STAT = alloc_status)
      IF (alloc_status /= 0) THEN
        PRINT*, "Allocation error"
        STOP
      ENDIF        
      
      nobed = 0 
      obedn_temp = 0
      
      DO seg = 1,nope
        DO nd = 1,obseg(seg)-1
        
          n1bed = obnds(nd,seg)
          n2bed = obnds(nd+1,seg)
          
  edges1: DO ed = 1,nbed
  
            ged = bedn(ed) 
            n1ed2 = ged2nn(1,ged)
            n2ed2 = ged2nn(2,ged)
            
            IF(((n1ed2 == n1bed).AND.(n2ed2 == n2bed)).OR. &
               ((n1ed2 == n2bed).AND.(n2ed2 == n1bed))) THEN
               
              nobed = nobed + 1
              obedn_temp(nobed) = ged
              recv_edge(ged) = 0
              ed_type(ged) = 1
              
              EXIT edges1
            ENDIF
            
          ENDDO edges1
        ENDDO
      ENDDO      
      
      ALLOCATE(obedn(nobed), STAT=alloc_status)
      IF (alloc_status /= 0) THEN
        PRINT*, "Allocation error"
        STOP
      ENDIF          
      
      obedn(1:nobed) = obedn_temp(1:nobed)
      
      RETURN
      END SUBROUTINE find_open_edges

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE find_flow_edges(nbou,fbseg,fbnds,ged2nn,nnfbed,nfbedn,nfbednn,nfbed,fbedn,recv_edge,ed_type)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: nbou
      INTEGER, DIMENSION(:,:), INTENT(IN) :: fbseg
      INTEGER, DIMENSION(:,:), INTENT(IN) :: fbnds
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ged2nn
      INTEGER, INTENT(OUT) :: nnfbed
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: nfbedn
      INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: nfbednn
      INTEGER, INTENT(OUT) :: nfbed
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: fbedn
      INTEGER, DIMENSION(:), INTENT(INOUT) :: recv_edge
      INTEGER, DIMENSION(:), INTENT(INOUT) :: ed_type      
      
      INTEGER :: el
      INTEGER :: nd
      INTEGER :: alloc_status      
      INTEGER :: seg,ed
      INTEGER :: ged,segtype
      INTEGER :: n1bed,n2bed
      INTEGER :: n1ed2,n2ed2
      INTEGER :: found
      INTEGER, DIMENSION(:), ALLOCATABLE :: nfbedn_temp
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: nfbednn_temp
      INTEGER, DIMENSION(:), ALLOCATABLE :: fbedn_temp
      
      
      ALLOCATE(nfbedn_temp(nbed),nfbednn_temp(nbed,3),fbedn_temp(nbed), STAT=alloc_status)
      IF (alloc_status /= 0) THEN
        PRINT*, "Allocation error"
        STOP
      ENDIF             
      
      
      nnfbed = 0 ! # no normal flow edge
      nfbedn_temp = 0
      nfbednn_temp = 0
      
      nfbed = 0 ! # flow specifed edges
      fbedn_temp = 0

      found = 0

      DO seg = 1,nbou
      
        segtype = fbseg(2,seg)
              
        DO nd = 1,fbseg(1,seg)-1
        
          n1bed = fbnds(nd,seg)
          n2bed = fbnds(nd+1,seg)
          found = 0 
          
  edges2: DO ed = 1,nbed
  
            ged = bedn(ed)
            n1ed2 = ged2nn(1,ged)
            n2ed2 = ged2nn(2,ged)
            
            IF(((n1ed2 == n1bed).AND.(n2ed2 == n2bed)).OR. &
               ((n1ed2 == n2bed).AND.(n2ed2 == n1bed))) THEN

              ! no normal flow edges
              IF( segtype == 0 .OR. segtype == 10 .OR. segtype == 20  .OR. &   ! land boundaries
                  segtype == 1 .OR. segtype == 11 .OR. segtype == 21 ) THEN    ! island boundaries
                  
                nnfbed = nnfbed + 1
                nfbedn_temp(nnfbed) = ged
                nfbednn_temp(nnfbed,1) = seg
                nfbednn_temp(nnfbed,2) = n1bed
                nfbednn_temp(nnfbed,3) = nd                
                recv_edge(ged) = 0
                ed_type(ged) = 10
                found = 1
                
                EXIT edges2               
              ENDIF

              ! specified normal flow edges
              IF ( segtype == 2 .OR. segtype == 12 .OR. segtype == 22 ) THEN
              
                nfbed = nfbed + 1
                fbedn_temp(nfbed) = ged
                recv_edge(ged) = 0
                ed_type(ged) = 12
                found = 1
                
                EXIT edges2
              ENDIF

            ENDIF
          ENDDO edges2
          
          IF (found == 0) THEN
!             PRINT*, seg, nd, fbseg(1,seg), segtype,myrank    
            PRINT "(4(A,I9))", "seg = ",seg, " nd = ",nd, " #nds in seg = ",fbseg(1,seg), " type = ",segtype            
            PRINT "(A)", "  edge not found"
          ENDIF
          
        ENDDO
      ENDDO
      
      ALLOCATE(nfbedn(nnfbed),nfbednn(nnfbed,3),fbedn(nfbed), STAT=alloc_status)
      DEALLOCATE(bedn, STAT=alloc_status)
      
      
      nfbedn(1:nnfbed) = nfbedn_temp(1:nnfbed)
      nfbednn(1:nnfbed,1:3) = nfbednn_temp(1:nnfbed,1:3)
      fbedn(1:nfbed) = fbedn_temp(1:nfbed)

      END SUBROUTINE find_flow_edges
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE find_recieve_edges(ned,recv_edge,nred,redn,ed_type)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ned
      INTEGER, DIMENSION(:), INTENT(IN) :: recv_edge
      INTEGER, INTENT(OUT) :: nred
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: redn
      INTEGER, DIMENSION(:), INTENT(INOUT) :: ed_type
      
      INTEGER :: el
      INTEGER :: nd
      INTEGER :: alloc_status      
      INTEGER :: ged                  

      nred = 0
      
#ifdef CMPI      
      DO ged = 1,ned
        IF(recv_edge(ged) == 1) THEN
          nred = nred + 1
        ENDIF
      ENDDO
      
      ALLOCATE(redn(nred), STAT = alloc_status)
      IF (alloc_status /= 0) THEN
        PRINT*, "Allocation error"
        STOP
      ENDIF          
      
      nred = 0
      DO ged = 1,ned
        IF(recv_edge(ged) == 1) THEN
          nred = nred + 1
          redn(nred) = ged
          ed_type(ged) = -1
        ENDIF
      ENDDO
#endif         
      
      RETURN
      END SUBROUTINE find_recieve_edges

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 


      SUBROUTINE find_element_edges(ne,ned,ged2el,ged2led,el2ged)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ne
      INTEGER, INTENT(IN) :: ned
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ged2el
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ged2led
      INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: el2ged
      
      INTEGER :: ged
      INTEGER :: el1,el2
      INTEGER :: led
      INTEGER :: alloc_status
      
      ALLOCATE(el2ged(ne,4), STAT = alloc_status)
      IF (alloc_status /= 0) THEN
        PRINT*, "Allocation error"
        STOP
      ENDIF          
      
      el2ged = 0
      
      DO ged = 1,ned
        
        el1 = ged2el(1,ged)
        el2 = ged2el(2,ged)
        
        IF (el1 /= 0) THEN
          led = ged2led(1,ged)
          el2ged(el1,led) = ged
        ENDIF
        
        IF (el2 /= 0) THEN
          led = ged2led(2,ged)
          el2ged(el2,led) = ged        
        ENDIF       
      
      ENDDO
      
      RETURN 
      END SUBROUTINE find_element_edges



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE find_neighbor_elements(ne,ned,ged2el,ged2led,el2el)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ne
      INTEGER, INTENT(IN) :: ned
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ged2el
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ged2led
      INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: el2el
      
      INTEGER :: ged
      INTEGER :: el1,el2
      INTEGER :: led1,led2
      INTEGER :: alloc_status      
      
      ALLOCATE(el2el(ne,4), STAT=alloc_status)
      IF (alloc_status /= 0) THEN
        PRINT*, "Allocation error"
        STOP
      ENDIF          
      
      el2el = 0
      
      DO ged = 1,ned
      
        el1 = ged2el(1,ged)
        el2 = ged2el(2,ged)
        
        IF (el1 /= 0 .and. el2 /= 0 ) THEN
          led1 = ged2led(1,ged)
          led2 = ged2led(2,ged)
          
          el2el(el1,led1) = el2
          el2el(el2,led2) = el1
        ENDIF
      
      ENDDO
      
      
      END SUBROUTINE
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE find_adjacent_nodes(nn,mnepn,ned,ged2nn,nadjnds,adjnds)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: nn
      INTEGER, INTENT(IN) :: mnepn
      INTEGER, INTENT(IN) :: ned
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ged2nn
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: nadjnds
      INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: adjnds
      
      INTEGER :: ed
      INTEGER :: n1,n2
      
      ALLOCATE(nadjnds(nn))
      ALLOCATE(adjnds(2*mnepn,nn))
      nadjnds = 0
      
      DO ed = 1,ned
        n1 = ged2nn(1,ed)  ! find the node numbers on each edge
        n2 = ged2nn(2,ed)
        
        nadjnds(n1) = nadjnds(n1) + 1 ! count the nodes adjacent to node n1
        nadjnds(n2) = nadjnds(n2) + 1 ! count the nodes adjacent to node n2
        
        adjnds(nadjnds(n1),n1) = n2 ! node n2 is adjacent to node n1
        adjnds(nadjnds(n2),n2) = n1 ! node n1 is adjacent to node n2                       
      ENDDO      
      
      RETURN
      END SUBROUTINE find_adjacent_nodes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE find_boundary_nodes(nn,nbed,bedn,ged2nn,nbnd,bndn)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nn
      INTEGER, INTENT(IN) :: nbed
      INTEGER, DIMENSION(:), INTENT(IN) :: bedn
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ged2nn
      INTEGER, INTENT(OUT) :: nbnd
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: bndn

      INTEGER :: ed,i
      INTEGER :: ged,nd
      
      INTEGER, DIMENSION(:), ALLOCATABLE :: nd_flag

      ALLOCATE(bndn(2*nbed),nd_flag(nn))
      nd_flag = 0
      nbnd = 0
      DO ed = 1,nbed
        ged = bedn(ed)

        DO i = 1,2
          nd = ged2nn(i,ged)
          IF (nd_flag(nd) == 0) THEN
            nbnd = nbnd + 1
            bndn(nbnd) = nd
            nd_flag(nd) = 1
          ENDIF
        ENDDO

      ENDDO

      RETURN
      END SUBROUTINE find_boundary_nodes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE min_edge_length(ned,ne,xy,ged2el,ged2nn,edlen,edlen_min)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ned
      INTEGER, INTENT(IN) :: ne
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: xy
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ged2el
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ged2nn
      REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: edlen
      REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: edlen_min
      
      INTEGER :: ed
      INTEGER :: el_in,el_ex
      INTEGER :: n1,n2
      REAL(rp) :: x1,x2,y1,y2
      INTEGER :: alloc_status 
      

      ALLOCATE(edlen(ned),edlen_min(ne), STAT=alloc_status)
      IF (alloc_status /= 0) THEN
        PRINT*, "Allocation error"
        STOP
      ENDIF     
      
      
      
      edlen_min = 1d10
      DO ed = 1,ned
        n1 = ged2nn(1,ed)
        n2 = ged2nn(2,ed)                
        
        x1 = xy(1,n1)
        x2 = xy(1,n2)
        
        y1 = xy(2,n1)
        y2 = xy(2,n2)
        
        edlen(ed) = sqrt((x2-x1)**2 + (y2-y1)**2)        
        
        el_in = ged2el(1,ed)
        el_ex = ged2el(2,ed)
        
        IF (edlen(ed) < edlen_min(el_in)) THEN
          edlen_min(el_in) = edlen(ed)
        ENDIF
        
        IF (el_ex > 0 ) THEN
          IF (edlen(ed) < edlen_min(el_ex)) THEN
            edlen_min(el_ex) = edlen(ed)
          ENDIF          
        ENDIF
        
        
      ENDDO
      
      END SUBROUTINE min_edge_length
         

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 


      SUBROUTINE print_connect_info(mnepn,ned,nied,nobed,nfbed,nnfbed,nred)
      
      IMPLICIT NONE      
      
      INTEGER, INTENT(IN) :: mnepn
      INTEGER, INTENT(IN) :: ned
      INTEGER, INTENT(IN) :: nied
      INTEGER, INTENT(IN) :: nobed
      INTEGER, INTENT(IN) :: nfbed
      INTEGER, INTENT(IN) :: nnfbed
      INTEGER, INTENT(IN) :: nred
      
      
      
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "       Edge Connectivity Information         "
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", " "        
      PRINT "(A,I7)", '   maximum elements per node:', mnepn
      PRINT "(A)", ' '            
      PRINT "(A,I7)", '   number of total edges:', ned
      PRINT "(A)", ' '  
      PRINT "(A,I7)", '   number of interior edges:', nied
      PRINT "(A)", ' '        
      PRINT "(A,I7)", '   number of open boundary edges:', nobed
      PRINT "(A)", ' '                
      PRINT "(A,I7)", '   number of specified normal boundary edges:', nfbed
      PRINT "(A,I7)", '   number of no normal flow boundary edges:', nnfbed
      PRINT "(A)", ' '        
      PRINT "(A,I7)", '   number of recieve edges:', nred
      PRINT "(A)", ' ' 
      PRINT "(A,I7)", '   number of missing edges:',ned-(nied+nobed+nfbed+nnfbed+nred)
      PRINT "(A)", ' '         
      
      
      RETURN
      END SUBROUTINE print_connect_info

      END MODULE edge_connectivity_mod
