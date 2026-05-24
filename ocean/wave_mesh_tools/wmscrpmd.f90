!/ ------------------------------------------------------------------- /
      MODULE WMSCRPMD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III                     |
!/                  |   E. Rogers, M. Dutour,           |
!/                  |   A. Roland, F. Ardhuin           |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         10-Dec-2014 |
!/                  +-----------------------------------+
!/
!/    06-Sep-2012 : Origination, transfer from WMGRIDMD ( version 4.08 )
!/    10-Dec-2014 : Add checks for allocate status      ( version 5.04 )
!/
!/    Copyright 2012 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS. 
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     Routines to determine and process grid dependencies in the 
!     multi-grid wave model.
!
!  2. Variables and types :
!
!  3. Subroutines and functions :
!
!      Name                          Type  Scope      Description
!     ----------------------------------------------------------------
!      scrip_wrapper                  Subr. Public   as the name says ... 
!      get_scrip_info_structured      Subr. Public   as the name says ... 
!      get_scrip_info_unstructured    Subr. Public   as the name says ... 
!      get_scrip_info                 Subr. Public   as the name says ... 
!      scrip_info_renormalization     Subr. Public   as the name says ... 
!      TRIANG_INDEXES                 Subr. Public   as the name says ... 
!      get_unstructured_vertex_degree Subr. Public   as the name says ... 
!      GET_BOUNDARY                   Subr. Public   as the name says ... 
!     ----------------------------------------------------------------
!
!  4. Subroutines and functions used :
!
!      Name                             Type  Module   Description
!     ----------------------------------------------------------------
!     get_unstructured_vertex_degree    Subr. W3TRIAMD Manage data structures
!     ----------------------------------------------------------------
!
!  5. Remarks :
!     
!     - The subroutines TRIANG_INDEXES, get_unstructured_vertex_degree, and 
!       GET_BOUNDARY should probably be renamed and moved to the module w3triamd
!
!  6. Switches :
!
!
!     !/S    Enable subroutine tracing.
!     !/Tn   Enable test output.
!
!  7. Source code :
!
!/ ------------------------------------------------------------------- /
!/
!/ Specify default accessibility
!/
      PUBLIC
      REAL*8 , ALLOCATABLE :: GRID_AREA(:)
!/
!/ Module private variable for checking error returns
!/
      INTEGER, PRIVATE        :: ISTAT
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE GET_SCRIP_INFO_UNSTRUCTURED (MNE, MNP, XYB, TRIGP,    &
       GRID_CENTER_LON, GRID_CENTER_LAT,                             &
       GRID_CORNER_LON, GRID_CORNER_LAT, GRID_MASK,                  &
       GRID_DIMS, GRID_SIZE, GRID_CORNERS, GRID_RANK)
!/                  +-----------------------------------+
!/                  | WAVEWATCH III                     |
!/                  | M. Dutour, A. Roland              |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         10-Dec-2014 !
!/                  +-----------------------------------+
!/
!  1. Original author :
!
!     Mathieu Dutour Sikiric, IRB & Aron Roland, Z&P
!
!  2. Last update :
!
!     See revisions.
!
!  3. Revisions :
!
!     20-Feb-2012 : Origination
!/    10-Dec-2014 : Add checks for allocate status      ( version 5.04 )
!
!  4. Copyright :
!
!  5. Purpose : 
!
!     Compute grid arrays for scrip for a specific unstructured grid
!     For interior vertices, we select for every triangle the barycenter
!     of the triangle. So to every node contained in N triangles we associate
!     a domain with N corners. Every one of those corners is contained
!     in 3 different domains.
!     For nodes that are on the boundary, we have to proceed differently.
!     For every such node, we have NEIGHBOR_PREV and NEIGHBOR_NEXT which
!     are the neighbor on each side of the boundary.
!     We put a corner on the middle of the edge. We also put a corner
!     on the node itself.
!     Note that instead of taking barycenters, we could have taken
!     Voronoi vertices, but this carries danger since Voronoi vertices
!     can be outside of the triangle. And it leaves open how to treat
!     the boundary.
!
!  6. Method :
!
!  7. Parameters, Variables and types :
!  
!  8. Called by : 
!
!     Subroutine get_scrip_info
!
!  9. Subroutines and functions used :
!
! 10. Error messages: 
!
! 11. Remarks :
!
! 12. Structure :
!
! 13. Switches :
!
! 14. Source code :
      IMPLICIT NONE
      INTEGER, INTENT(IN)      :: MNE
      INTEGER, INTENT(IN)      :: MNP
      REAL*8, INTENT(IN), DIMENSION(:,:) :: XYB
      INTEGER, INTENT(IN), DIMENSION(:,:) :: TRIGP
      REAL*8, INTENT(OUT), ALLOCATABLE :: GRID_CENTER_LON(:)
      REAL*8, INTENT(OUT), ALLOCATABLE :: GRID_CENTER_LAT(:)
      LOGICAL, INTENT(OUT), ALLOCATABLE :: GRID_MASK(:)
      REAL*8, INTENT(OUT), ALLOCATABLE :: GRID_CORNER_LON(:,:)
      REAL*8, INTENT(OUT), ALLOCATABLE :: GRID_CORNER_LAT(:,:)
      INTEGER, INTENT(OUT), ALLOCATABLE :: GRID_DIMS(:)
      INTEGER, INTENT(OUT) :: GRID_SIZE, GRID_CORNERS, GRID_RANK

      INTEGER DIRAPPROACH, DUALAPPROACH, THEAPPROACH
      INTEGER IE, IP, I, J
      INTEGER NBPLUS, NBMINUS
      INTEGER I1, I2, I3
      REAL*8 :: ELON1, ELON2, ELON3, ELON, ELONC
      REAL*8 :: ELAT1, ELAT2, ELAT3, ELAT, ELATC
      REAL *8 :: DELTALON12, DELTALON13, DELTALAT12, DELTALAT13
      REAL *8 :: THEDET
      REAL*8 :: PT(3,2)
      INTEGER, POINTER :: IOBP(:), TRIGINCD(:)
      INTEGER, POINTER :: NEIGHBOR_PREV(:), NEIGHBOR_NEXT(:)
      INTEGER, POINTER :: NBASSIGNEDCORNER(:), LISTNBCORNER(:)
      INTEGER, POINTER :: STATUS(:), NEXTVERT(:), PREVVERT(:), FINALVERT(:)
      INTEGER :: MAXCORNER, NBCORNER
      INTEGER :: IDX, IPNEXT, IPPREV, NB, INEXT, IPREV
      REAL*8, POINTER :: LON_CENT_TRIG(:), LAT_CENT_TRIG(:)
      REAL*8 :: ELONIP, ELONNEXT, ELONPREV, ELONN, ELONP
      REAL*8 :: ELATIP, ELATNEXT, ELATPREV, ELATN, ELATP
      INTEGER :: ISFINISHED, ZPREV
      INTEGER :: DODEBUG
      GRID_RANK=2
      DIRAPPROACH=1
      DUALAPPROACH=2
      THEAPPROACH=DUALAPPROACH
      IF (THEAPPROACH .EQ. DIRAPPROACH) THEN
        ALLOCATE(GRID_CENTER_LON(MNE), STAT=ISTAT)
        CALL CHECK_ALLOC_STATUS ( ISTAT )
        ALLOCATE(GRID_CENTER_LAT(MNE), STAT=ISTAT)
        CALL CHECK_ALLOC_STATUS ( ISTAT )
        ALLOCATE(GRID_CORNER_LON(3,MNE), STAT=ISTAT)
        CALL CHECK_ALLOC_STATUS ( ISTAT )
        ALLOCATE(GRID_CORNER_LAT(3,MNE), STAT=ISTAT)
        CALL CHECK_ALLOC_STATUS ( ISTAT )
        ALLOCATE(GRID_MASK(MNE), STAT=ISTAT)
        CALL CHECK_ALLOC_STATUS ( ISTAT )
        DO IE=1,MNE
          I1=TRIGP(IE,1)
          I2=TRIGP(IE,2)
          I3=TRIGP(IE,3)
          ELON1=XYB(I1,1)
          ELON2=XYB(I2,1)
          ELON3=XYB(I3,1)
          ELAT1=XYB(I1,2)
          ELAT2=XYB(I2,2)
          ELAT3=XYB(I3,2)
          ELON=(ELON1 + ELON2 + ELON3)/3
          ELAT=(ELAT1 + ELAT2 + ELAT3)/3
          GRID_CENTER_LON(IE)=ELON
          GRID_CENTER_LAT(IE)=ELAT
          GRID_CORNER_LON(1,IE)=ELON1
          GRID_CORNER_LON(2,IE)=ELON2
          GRID_CORNER_LON(3,IE)=ELON3
          GRID_CORNER_LAT(1,IE)=ELAT1
          GRID_CORNER_LAT(2,IE)=ELAT2
          GRID_CORNER_LAT(3,IE)=ELAT3
          GRID_MASK(IE)=.TRUE.
        END DO
        GRID_CORNERS=3
      END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (THEAPPROACH .EQ. DUALAPPROACH) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ALLOCATE(TRIGINCD(MNP), STAT=ISTAT)
        CALL CHECK_ALLOC_STATUS ( ISTAT )
        ALLOCATE(IOBP(MNP), STAT=ISTAT)
        CALL CHECK_ALLOC_STATUS ( ISTAT )
        ALLOCATE(NEIGHBOR_NEXT(MNP), STAT=ISTAT)
        CALL CHECK_ALLOC_STATUS ( ISTAT )
        ALLOCATE(NEIGHBOR_PREV(MNP), STAT=ISTAT)
        CALL CHECK_ALLOC_STATUS ( ISTAT )
        ALLOCATE(NBASSIGNEDCORNER(MNP), STAT=ISTAT)
        CALL CHECK_ALLOC_STATUS ( ISTAT )
        ALLOCATE(LISTNBCORNER(MNP), STAT=ISTAT)
        CALL CHECK_ALLOC_STATUS ( ISTAT )

        ALLOCATE(STATUS(MNP), STAT=ISTAT)
        CALL CHECK_ALLOC_STATUS ( ISTAT )
        ALLOCATE(NEXTVERT(MNP), STAT=ISTAT)
        CALL CHECK_ALLOC_STATUS ( ISTAT )
        ALLOCATE(PREVVERT(MNP), STAT=ISTAT)
        CALL CHECK_ALLOC_STATUS ( ISTAT )
        ALLOCATE(FINALVERT(MNP), STAT=ISTAT)
        CALL CHECK_ALLOC_STATUS ( ISTAT )
        ALLOCATE(LON_CENT_TRIG(MNE), STAT=ISTAT)
        CALL CHECK_ALLOC_STATUS ( ISTAT )
        ALLOCATE(LAT_CENT_TRIG(MNE), STAT=ISTAT)
        CALL CHECK_ALLOC_STATUS ( ISTAT )

        CALL GET_UNSTRUCTURED_VERTEX_DEGREE (MNP, MNE,             &
               TRIGP, TRIGINCD)
        CALL GET_BOUNDARY(MNP, MNE, TRIGP, IOBP,     &
               NEIGHBOR_PREV, NEIGHBOR_NEXT)

        ! Find max number of corners
        MAXCORNER=0
        DO IP=1,MNP
          !IF (IP == 53591) THEN 
          !   PRINT*, "NEIGHBOR_NEXT", NEIGHBOR_NEXT(IP)
          !   PRINT*, XYB(IP,1), XYB(IP,2), XYB(IP,3)
          !   PRINT*, TRIGINCD(IP)
          !ENDIF
          IF (NEIGHBOR_NEXT(IP) .EQ. 0) THEN
            NBCORNER=TRIGINCD(IP)
          ELSE
            NBCORNER=TRIGINCD(IP) + 3
            PRINT*, IP, IOBP(IP)
          END IF
          LISTNBCORNER(IP)=NBCORNER
          IF (NBCORNER .GT. MAXCORNER) THEN
            MAXCORNER=NBCORNER
          END IF
        END DO
        GRID_CORNERS=MAXCORNER

        ALLOCATE(GRID_CENTER_LON(MNP), STAT=ISTAT)
        CALL CHECK_ALLOC_STATUS ( ISTAT )
        ALLOCATE(GRID_CENTER_LAT(MNP), STAT=ISTAT)
        CALL CHECK_ALLOC_STATUS ( ISTAT )
        ALLOCATE(GRID_CORNER_LON(MAXCORNER,MNP), STAT=ISTAT)
        CALL CHECK_ALLOC_STATUS ( ISTAT )
        ALLOCATE(GRID_CORNER_LAT(MAXCORNER,MNP), STAT=ISTAT)
        CALL CHECK_ALLOC_STATUS ( ISTAT )
        ALLOCATE(GRID_MASK(MNP), STAT=ISTAT)
        CALL CHECK_ALLOC_STATUS ( ISTAT )

        ! Add first three corners for boundaries
        NBASSIGNEDCORNER(:)=0
        DO IP=1,MNP
          GRID_MASK(IP)=.TRUE.
          IF (NEIGHBOR_NEXT(IP) .GT. 0) THEN
            IPNEXT=NEIGHBOR_NEXT(IP)
            IPPREV=NEIGHBOR_PREV(IP)
            ELONIP=DBLE(XYB(IP,1))
            ELATIP=DBLE(XYB(IP,2))
            ELONNEXT=DBLE(XYB(IPNEXT,1))
            ELATNEXT=DBLE(XYB(IPNEXT,2))
            ELONPREV=DBLE(XYB(IPPREV,1))
            ELATPREV=DBLE(XYB(IPPREV,2))
 
            ! Periodicity fix for corner node
            IF ( ABS(ELONIP - ELONNEXT) .GT. 180.0 ) THEN
                     ELONNEXT = ELONNEXT -SIGN(360.0d0,(ELONIP - ELONNEXT))
            ENDIF
            IF ( ABS(ELONIP - ELONPREV) .GT. 180.0 ) THEN
                     ELONPREV = ELONPREV -SIGN(360.0d0,(ELONIP - ELONPREV))
            ENDIF

            ELONN=(ELONIP+ELONNEXT)/2.0
            ELATN=(ELATIP+ELATNEXT)/2.0
            ELONP=(ELONIP+ELONPREV)/2.0
            ELATP=(ELATIP+ELATPREV)/2.0

            
            GRID_CORNER_LON(1,IP)=ELONN
            GRID_CORNER_LAT(1,IP)=ELATN
            GRID_CORNER_LON(2,IP)=ELONIP
            GRID_CORNER_LAT(2,IP)=ELATIP
            GRID_CORNER_LON(3,IP)=ELONP
            GRID_CORNER_LAT(3,IP)=ELATP
            NBASSIGNEDCORNER(IP)=3
          END IF
        END DO

        ! Compute centers
        DO IP=1,MNP
          GRID_CENTER_LON(IP)=DBLE(XYB(IP,1))
          GRID_CENTER_LAT(IP)=DBLE(XYB(IP,2))
        END DO

        ! Check triangle node orientation
        ! Compute triangle centers
        NBPLUS=0
        NBMINUS=0
        DO IE=1,MNE
          I1=TRIGP(IE,1)
          I2=TRIGP(IE,2)
          I3=TRIGP(IE,3)

          PT(1,1)=DBLE(XYB(I1,1))
          PT(2,1)=DBLE(XYB(I2,1))
          PT(3,1)=DBLE(XYB(I3,1))
          PT(1,2)=DBLE(XYB(I1,2))
          PT(2,2)=DBLE(XYB(I2,2))
          PT(3,2)=DBLE(XYB(I3,2))

          CALL FIX_PERIODCITY(PT)

          ELON1 = PT(1,1)
          ELON2 = PT(2,1)
          ELON3 = PT(3,1)
          ELAT1 = PT(1,2)
          ELAT2 = PT(2,2)
          ELAT3 = PT(3,2)

          DELTALON12=ELON2 - ELON1
          DELTALON13=ELON3 - ELON1
          DELTALAT12=ELAT2 - ELAT1
          DELTALAT13=ELAT3 - ELAT1
          THEDET=DELTALON12*DELTALAT13 - DELTALON13*DELTALAT12
          IF (THEDET.GT.0) THEN
            NBPLUS=NBPLUS+1
          END IF
          IF (THEDET.LT.0) THEN
            NBMINUS=NBMINUS+1
          END IF
          ELON=(ELON1 + ELON2 + ELON3)/3.0
          ELAT=(ELAT1 + ELAT2 + ELAT3)/3.0
       
          
          LON_CENT_TRIG(IE)=ELON
          LAT_CENT_TRIG(IE)=ELAT

        END DO
        DODEBUG=0
        IF (DODEBUG.EQ.1) THEN
          print *, 'nbplus=', nbplus, ' nbminus=', nbminus
        END IF
        
        STATUS(:) = 0
        NEXTVERT(:) = 0
        PREVVERT(:) = 0
        DO IE=1,MNE
          DO I=1,3
            CALL TRIANG_INDEXES(I, INEXT, IPREV)
            IP=TRIGP(IE,I)
            IPNEXT=TRIGP(IE,INEXT)
            IPPREV=TRIGP(IE,IPREV)
            IF (STATUS(IP).EQ.0) THEN
              IF (NEIGHBOR_NEXT(IP).EQ.0) THEN
                STATUS(IP)=1
                FINALVERT(IP)=IPPREV
                PREVVERT(IP)=IPPREV
                NEXTVERT(IP)=IPNEXT
              ELSE
                IF (NEIGHBOR_PREV(IP).EQ.IPPREV) THEN
                  STATUS(IP)=1
                  PREVVERT(IP)=IPPREV
                  NEXTVERT(IP)=IPNEXT
                  FINALVERT(IP)=NEIGHBOR_NEXT(IP)
                END IF
              END IF
            END IF
          END DO
        END DO
        STATUS(:)=0
        DO
          ISFINISHED=1
          DO IE=1,MNE
            ELON=LON_CENT_TRIG(IE)
            ELAT=LAT_CENT_TRIG(IE)
            DO I=1,3
              CALL TRIANG_INDEXES(I, INEXT, IPREV)
              IP=TRIGP(IE,I)
              IPNEXT=TRIGP(IE,INEXT)
              IPPREV=TRIGP(IE,IPREV)
              IF (STATUS(IP).EQ.0) THEN
                ISFINISHED=0
                ZPREV=PREVVERT(IP)
                IF (ZPREV.EQ.IPPREV) THEN
                  IDX=NBASSIGNEDCORNER(IP)
                  IDX=IDX+1
                  GRID_CORNER_LON(IDX,IP)=ELON
                  GRID_CORNER_LAT(IDX,IP)=ELAT
                  NBASSIGNEDCORNER(IP)=IDX
                  PREVVERT(IP)=IPNEXT
                  IF (IPNEXT.EQ.FINALVERT(IP)) THEN
                    STATUS(IP)=1
                  END IF
                END IF
              END IF
            END DO
          END DO
          IF (ISFINISHED.EQ.1) THEN
            EXIT
          END IF
        END DO
        DO IP=1,MNP
          IF (NBASSIGNEDCORNER(IP).NE.LISTNBCORNER(IP)) THEN
            WRITE(*,*) 'Incoherent number at IP=', IP
            WRITE(*,*) '  NbAssignedCorner(IP)=', NbAssignedCorner(IP)
            WRITE(*,*) '  ListNbCorner(IP)=', ListNbCorner(IP)
            WRITE(*,*) '  N_N=', NEIGHBOR_NEXT(IP), 'N_P=', NEIGHBOR_PREV(IP)
            WRITE(*,*) '  TrigIncd=', TrigIncd(IP)
            !STOP 'wmscrpmd, case 2'
          END IF
        END DO

        ! if the number of corner is below threshold, we have to
        ! add some more.
        DO IP=1,MNP
          NB=NBASSIGNEDCORNER(IP)
          IF (NB .LT. MAXCORNER) THEN
            ELON=GRID_CORNER_LON(NB,IP)
            ELAT=GRID_CORNER_LAT(NB,IP)
            DO IDX=NB+1,MAXCORNER
              GRID_CORNER_LON(IDX,IP)=ELON
              GRID_CORNER_LAT(IDX,IP)=ELAT
            END DO
          END IF
        END DO

        ! Calculate area of cell
        ALLOCATE(GRID_AREA(MNP))
        DO I1=1,MNP
          GRID_AREA(I1) = 0.0
          NB=NBASSIGNEDCORNER(I1)
          DO I2=1,NB
            I3 = MOD(I2,NB)+1
            PT(1,1)=GRID_CENTER_LON(I1)
            PT(2,1)=GRID_CORNER_LON(I2,I1)
            PT(3,1)=GRID_CORNER_LON(I3,I1)
            PT(1,2)=GRID_CENTER_LAT(I1)
            PT(2,2)=GRID_CORNER_LAT(I2,I1)
            PT(3,2)=GRID_CORNER_LON(I3,I1)

            CALL FIX_PERIODCITY(PT)

            ELON1 = PT(1,1)
            ELON2 = PT(2,1)
            ELON3 = PT(3,1)
            ELAT1 = PT(1,2)
            ELAT2 = PT(2,2)
            ELAT3 = PT(3,2)

            DELTALON12=ELON2 - ELON1
            DELTALON13=ELON3 - ELON1
            DELTALAT12=ELAT2 - ELAT1
            DELTALAT13=ELAT3 - ELAT1
            THEDET=DELTALON12*DELTALAT13 - DELTALON13*DELTALAT12

            GRID_AREA(I1) = GRID_AREA(I1) + ABS(THEDET)
          ENDDO
        ENDDO

        DEALLOCATE(NBASSIGNEDCORNER, STAT=ISTAT)
        CALL CHECK_DEALLOC_STATUS ( ISTAT )
        DEALLOCATE(LISTNBCORNER, STAT=ISTAT)
        CALL CHECK_DEALLOC_STATUS ( ISTAT )
        DEALLOCATE(TRIGINCD, STAT=ISTAT)
        CALL CHECK_DEALLOC_STATUS ( ISTAT )
        DEALLOCATE(IOBP, STAT=ISTAT)
        CALL CHECK_DEALLOC_STATUS ( ISTAT )
        DEALLOCATE(NEIGHBOR_PREV, STAT=ISTAT)
        CALL CHECK_DEALLOC_STATUS ( ISTAT )
        DEALLOCATE(NEIGHBOR_NEXT, STAT=ISTAT)
        CALL CHECK_DEALLOC_STATUS ( ISTAT )
        DEALLOCATE(STATUS, STAT=ISTAT)
        CALL CHECK_DEALLOC_STATUS ( ISTAT )
        DEALLOCATE(NEXTVERT, STAT=ISTAT)
        CALL CHECK_DEALLOC_STATUS ( ISTAT )
        DEALLOCATE(PREVVERT, STAT=ISTAT)
        CALL CHECK_DEALLOC_STATUS ( ISTAT )
        DEALLOCATE(FINALVERT, STAT=ISTAT)
        CALL CHECK_DEALLOC_STATUS ( ISTAT )
        DEALLOCATE(LON_CENT_TRIG, STAT=ISTAT)
        CALL CHECK_DEALLOC_STATUS ( ISTAT )
        DEALLOCATE(LAT_CENT_TRIG, STAT=ISTAT)
        CALL CHECK_DEALLOC_STATUS ( ISTAT )

        GRID_SIZE=MNP
        GRID_RANK=1
        ALLOCATE(GRID_DIMS(GRID_RANK))
        GRID_DIMS(1) = GRID_SIZE

        DO I = 1,GRID_SIZE
          IF (GRID_CENTER_LON(I) < 0.0) THEN
            GRID_CENTER_LON(I) = GRID_CENTER_LON(I)+360.0
          ENDIF
          DO J = 1,GRID_CORNERS
            IF (GRID_CORNER_LON(J,I) < 0.0) THEN
              GRID_CORNER_LON(J,I) = GRID_CORNER_LON(J,I)+360.0 
            ENDIF
          ENDDO
        ENDDO

      END IF
      END SUBROUTINE GET_SCRIP_INFO_UNSTRUCTURED

!/ ------------------------------------------------------------------- /
      SUBROUTINE GET_SCRIP_INFO(MNE, MNP, XYB, TRIGP,                   &
     &   GRID_CENTER_LON, GRID_CENTER_LAT,                              &
     &   GRID_CORNER_LON, GRID_CORNER_LAT, GRID_MASK,                   &
     &   GRID_DIMS, GRID_SIZE, GRID_CORNERS, GRID_RANK)
!  1. Original author :
!
!     Mathieu Dutour Sikiric, IRB & Aron Roland, Z&P
!
!  2. Last update :
!
!     See revisions.
!
!  3. Revisions :
!
!     20-Feb-2012 : Origination, this is adapted from Erick Rogers
!                   code by splitting the code into sections.
!
!  4. Copyright :
!
!  5. Purpose : 
!
!     Compute grid for scrip for a specific structured grid.
!     This is adapted from Erick Rogers code by making it cleaner.
!
!  6. Method :
!
!  7. Parameters, Variables and types :
!  
!  8. Called by : 
!
!     Subroutine scrip_wrapper
!
!  9. Subroutines and functions used :
!
! 10. Error messages: 
!
! 11. Remarks :
!
! 12. Structure :
!
! 13. Switches :
!
! 14. Source code :
      IMPLICIT NONE
      INTEGER, INTENT(IN)      :: MNE 
      INTEGER, INTENT(IN)      :: MNP
      REAL*8, INTENT(IN), DIMENSION(:,:) :: XYB
      INTEGER, INTENT(IN), DIMENSION(:,:) :: TRIGP
      REAL*8, INTENT(OUT), ALLOCATABLE :: GRID_CENTER_LON(:)
      REAL*8, INTENT(OUT), ALLOCATABLE :: GRID_CENTER_LAT(:)
      LOGICAL, INTENT(OUT), ALLOCATABLE :: GRID_MASK(:)
      REAL*8, INTENT(OUT), ALLOCATABLE :: GRID_CORNER_LON(:,:)
      REAL*8, INTENT(OUT), ALLOCATABLE :: GRID_CORNER_LAT(:,:)
      INTEGER, INTENT(OUT), ALLOCATABLE :: GRID_DIMS(:)
      INTEGER, INTENT(OUT) :: GRID_SIZE, GRID_CORNERS, GRID_RANK
      REAL*8 :: DLON1, DLAT1, DLON2, DLAT2, THEDET
      INTEGER :: I, J
      INTEGER :: IC, JC, IP, CHECKSIGNS, NBPLUS, NBMINUS, NBZERO
      INTEGER :: PRINTDATA, PRINTMINMAX
      REAL*8 :: MINLON, MINLAT, MAXLON, MAXLAT
      REAL*8 :: MINLONCORNER, MAXLONCORNER, MINLATCORNER, MAXLATCORNER
      REAL*8 :: PT(3,2)
      CALL GET_SCRIP_INFO_UNSTRUCTURED (MNE, MNP, XYB, TRIGP,           &
     &   GRID_CENTER_LON, GRID_CENTER_LAT,                              &
     &   GRID_CORNER_LON, GRID_CORNER_LAT, GRID_MASK,                   &
     &   GRID_DIMS, GRID_SIZE, GRID_CORNERS, GRID_RANK)
      CHECKSIGNS=1
      IF (CHECKSIGNS.EQ.1) THEN
        NBPLUS=0
        NBMINUS=0
        NBZERO=0
        DO IP=1,GRID_SIZE
          DO IC=1,GRID_CORNERS
            IF (IC.EQ.GRID_CORNERS) THEN
              JC=1
            ELSE
              JC=IC+1
            END IF

            PT(1,1) = GRID_CENTER_LON(IP) 
            PT(1,2) = GRID_CENTER_LAT(IP)
            PT(2,1) = GRID_CORNER_LON(IC,IP)
            PT(2,2) = GRID_CORNER_LAT(IC,IP)
            PT(3,1) = GRID_CORNER_LON(JC,IP)
            PT(3,2) = GRID_CORNER_LAT(JC,IP)

            CALL FIX_PERIODCITY(PT)

            DLON1=PT(2,1)-PT(1,1)
            DLON2=PT(3,1)-PT(1,1)
            DLAT1=PT(2,2)-PT(1,2)
            DLAT2=PT(3,2)-PT(1,2)

            THEDET=DLON1*DLAT2 - DLON2*DLAT1
            IF (THEDET.GT.1d-8) THEN
              NBPLUS=NBPLUS+1
            ELSE IF (THEDET.LT.-1d-8) THEN
              NBMINUS=NBMINUS+1
            ELSE
              NBZERO=NBZERO+1
            END IF
          END DO
        END DO

        WRITE(*,*) 'SI nbPlus=', nbPlus, ' nbMinus=', nbMinus, ' nbZero=', nbZero

      END IF
      END SUBROUTINE GET_SCRIP_INFO

!/ ------------------------------------------------------------------- /
      SUBROUTINE SCRIP_INFO_RENORMALIZATION(                            &
     &   GRID_CENTER_LON, GRID_CENTER_LAT,                              &
     &   GRID_CORNER_LON, GRID_CORNER_LAT, GRID_MASK,                   &
     &   GRID_DIMS, GRID_SIZE, GRID_CORNERS, GRID_RANK,                 &
     &   CONV_DX, CONV_DY, OFFSET, GRIDSHIFT)
!  1. Original author :
!
!     Mathieu Dutour Sikiric, IRB & Aron Roland, Z&P
!     Adapted from Erick Rogers scrip_wrapper
!
!  2. Last update :
!
!     See revisions.
!
!  3. Revisions :
!
!     20-Feb-2012 : Origination
!
!  4. Copyright :
!
!  5. Purpose : 
!
!     This is adapted from Erick Rogers scrip_wrapper
!     Purpose is to rescale according to whether the grid is spherical
!     or not and to adjust by some small shift to keep SCRIP happy
!     in situations where nodes of different grids overlay
!
!  6. Method :
!
!    We apply various transformations to the longitude latitude
!    following here the transformations that were done only in
!    finite difference meshes.
!
!  7. Parameters, Variables and types :
!  
!  8. Called by : 
!
!     Subroutine WMGHGH
!
!  9. Subroutines and functions used :
!
! 10. Error messages: 
!
! 11. Remarks :
!
! 12. Structure :
!
! 13. Switches :
!
! 14. Source code :
      IMPLICIT NONE
      REAL*8, INTENT(INOUT) :: GRID_CENTER_LON(:)
      REAL*8, INTENT(INOUT) :: GRID_CENTER_LAT(:)
      LOGICAL, INTENT(IN) :: GRID_MASK(:)
      REAL*8, INTENT(INOUT) :: GRID_CORNER_LON(:,:)
      REAL*8, INTENT(INOUT) :: GRID_CORNER_LAT(:,:)
      INTEGER, INTENT(IN) :: GRID_DIMS(:)
      INTEGER, INTENT(IN) :: GRID_SIZE, GRID_CORNERS, GRID_RANK
      REAL*8 :: CONV_DX, CONV_DY, OFFSET, GRIDSHIFT
      REAL*8 DEG2RAD
      !
      INTEGER :: I, J, IP
      REAL*8 :: MINLON, MINLAT, MAXLON, MAXLAT, HLON, HLAT
      REAL*8 :: MINLONCORNER, MAXLONCORNER, MINLATCORNER, MAXLATCORNER

      DO I=1,GRID_SIZE
        GRID_CENTER_LON(I)=(GRID_CENTER_LON(I)+OFFSET)/CONV_DX +        &
     &       GRIDSHIFT
        GRID_CENTER_LAT(I)=GRID_CENTER_LAT(I)/CONV_DY +                 &
     &       GRIDSHIFT
        IF(GRID_CENTER_LON(I)>360.0) THEN
          GRID_CENTER_LON(I)=GRID_CENTER_LON(I)-360.0
        END IF
        IF(GRID_CENTER_LON(I)<000.0) THEN
          GRID_CENTER_LON(I)=GRID_CENTER_LON(I)+360.0
        END IF
        DO J=1,GRID_CORNERS
          GRID_CORNER_LON(J, I)=(GRID_CORNER_LON(J, I)+OFFSET)/CONV_DX+ &
     &       GRIDSHIFT
          GRID_CORNER_LAT(J, I)=GRID_CORNER_LAT(J, I)/CONV_DY +         &
     &       GRIDSHIFT
          IF(GRID_CORNER_LON(J,I)>360.0) THEN
            GRID_CORNER_LON(J,I)=GRID_CORNER_LON(J,I)-360.0
          END IF
          IF(GRID_CORNER_LON(J,I)<000.0) THEN
            GRID_CORNER_LON(J,I)=GRID_CORNER_LON(J,I)+360.0
          END IF
        END DO
      END DO

      END SUBROUTINE SCRIP_INFO_RENORMALIZATION

!/ ------------------------------------------------------------------- /
      SUBROUTINE TRIANG_INDEXES(I, INEXT, IPREV)
!  1. Original author :
!
!     Mathieu Dutour Sikiric, IRB & Aron Roland, Z&P
!
        INTEGER, INTENT(IN)  :: I
        INTEGER, INTENT(OUT) :: INEXT, IPREV
        IF (I.EQ.1) THEN
          INEXT=3
        ELSE
          INEXT=I-1
        END IF
        IF (I.EQ.3) THEN
          IPREV=1
        ELSE
          IPREV=I+1
        END IF
      END SUBROUTINE

!/ ------------------------------------------------------------------- /
      SUBROUTINE GET_UNSTRUCTURED_VERTEX_DEGREE (MNP, MNE, TRIGP,         &
     &          TRIGINCD)
!  Written:
!
!    20-Feb-2012
!
!  Author:
!
!    Mathieu Dutour Sikiric, IRB & Aron Roland, Z&P
!
!  Parameters:
!   Input:
!    MNP: number of nodes
!    INE: list of nodes
!    MNE: number of triangles
!   Output:
!    TrigIncd (number of triangles contained by vertices
!
!  Description:
!    this function returns the list of incidences
!    
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: MNP, MNE
      INTEGER, INTENT(IN) :: TRIGP(:,:)
      INTEGER, INTENT(OUT) :: TRIGINCD(:)
      INTEGER :: IP, IE, I
      TRIGINCD=0
      DO IE=1,MNE
        DO I=1,3
          IP=TRIGP(IE,I)
          TRIGINCD(IP)=TRIGINCD(IP) + 1
        END DO
      END DO
      END SUBROUTINE GET_UNSTRUCTURED_VERTEX_DEGREE

!/ ------------------------------------------------------------------- /
      SUBROUTINE GET_BOUNDARY(MNP, MNE, TRIGP, IOBP, NEIGHBOR_PREV,       &
     &   NEIGHBOR_NEXT)
!/                  +-----------------------------------+
!/                  | WAVEWATCH III                     |
!/                  | M. Dutour, A. Roland              |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         10-Dec-2014 !
!/                  +-----------------------------------+
!/
!  Written:
!
!    20-Feb-2012
!/    10-Dec-2014 : Add checks for allocate status      ( version 5.04 )
!
!  Author:
!
!    Mathieu Dutour Sikiric, IRB & Aron Roland, Z&P
!
!  Parameters:
!   Input:
!    MNP: number of nodes
!    TRIGP: list of nodes
!    MNE: number of triangles
!   Output:
!    NEIGHBOR
!
!  Description:
!    if a node belong to a boundary, the function
!    returns the neighbor of this point on one side.
!    if the point is interior then the value 0 is set.
!    
        IMPLICIT NONE
               
        INTEGER, INTENT(IN)             :: MNP, MNE, TRIGP(MNE,3)
        INTEGER, INTENT(INOUT)          :: IOBP(MNP)
        INTEGER, INTENT(INOUT)          :: NEIGHBOR_PREV(MNP)
        INTEGER, INTENT(INOUT)          :: NEIGHBOR_NEXT(MNP)
        
        INTEGER, POINTER :: STATUS(:)
        INTEGER, POINTER :: COLLECTED(:)
        INTEGER, POINTER :: NEXTVERT(:)
        INTEGER, POINTER :: PREVVERT(:)
        
        INTEGER :: IE, I, IP, IP2, IP3
        INTEGER :: ISFINISHED, INEXT, IPREV
        INTEGER :: IPNEXT, IPPREV, ZNEXT, ZPREV
        
        ALLOCATE(STATUS(MNP), STAT=ISTAT)
        CALL CHECK_ALLOC_STATUS ( ISTAT )
        ALLOCATE(COLLECTED(MNP), STAT=ISTAT)
        CALL CHECK_ALLOC_STATUS ( ISTAT )
        ALLOCATE(PREVVERT(MNP), STAT=ISTAT)
        CALL CHECK_ALLOC_STATUS ( ISTAT )
        ALLOCATE(NEXTVERT(MNP), STAT=ISTAT)
        CALL CHECK_ALLOC_STATUS ( ISTAT )
        IOBP = 0
        NEIGHBOR_NEXT = 0
        NEIGHBOR_PREV = 0
        
!  Now computing the next items
        STATUS = 0
        NEXTVERT = 0
        PREVVERT = 0
        DO IE=1,MNE
          DO I=1,3
            CALL TRIANG_INDEXES(I, INEXT, IPREV)
            IP=TRIGP(IE,I)
            IPNEXT=TRIGP(IE,INEXT)
            IPPREV=TRIGP(IE,IPREV)
            IF (STATUS(IP).EQ.0) THEN
              STATUS(IP)=1
              PREVVERT(IP)=IPPREV
              NEXTVERT(IP)=IPNEXT
            END IF
          END DO
        END DO
        STATUS(:)=0
        DO
          COLLECTED(:)=0
          DO IE=1,MNE
            DO I=1,3
              CALL TRIANG_INDEXES(I, INEXT, IPREV)
              IP=TRIGP(IE,I)
              IPNEXT=TRIGP(IE,INEXT)
              IPPREV=TRIGP(IE,IPREV)
              IF (STATUS(IP).EQ.0) THEN
                ZNEXT=NEXTVERT(IP)
                IF (ZNEXT.EQ.IPPREV) THEN
                  COLLECTED(IP)=1
                  NEXTVERT(IP)=IPNEXT
                  IF (NEXTVERT(IP).EQ.PREVVERT(IP)) THEN
                    STATUS(IP)=1
                  END IF
                END IF
              END IF
            END DO
          END DO

          ISFINISHED=1
          DO IP=1,MNP
            IF ((COLLECTED(IP).EQ.0).AND.(STATUS(IP).EQ.0)) THEN
              STATUS(IP)=-1
              NEIGHBOR_NEXT(IP)=NEXTVERT(IP)
            END IF
            IF (STATUS(IP).EQ.0) THEN
              ISFINISHED=0
            END IF
          END DO
          IF (ISFINISHED.EQ.1) THEN
            EXIT
          END IF
        END DO

!  Now computing the prev items
        STATUS = 0
        NEXTVERT = 0
        PREVVERT = 0
        DO IE=1,MNE
          DO I=1,3
            CALL TRIANG_INDEXES(I, INEXT, IPREV)
            IP=TRIGP(IE,I)
            IPNEXT=TRIGP(IE,INEXT)
            IPPREV=TRIGP(IE,IPREV)
            IF (STATUS(IP).EQ.0) THEN
              STATUS(IP)=1
              PREVVERT(IP)=IPPREV
              NEXTVERT(IP)=IPNEXT
            END IF
          END DO
        END DO
        STATUS(:)=0
        DO
          COLLECTED(:)=0
          DO IE=1,MNE
            DO I=1,3
              CALL TRIANG_INDEXES(I, INEXT, IPREV)
              IP=TRIGP(IE,I)
              IPNEXT=TRIGP(IE,INEXT)
              IPPREV=TRIGP(IE,IPREV)
              IF (STATUS(IP).EQ.0) THEN
                ZPREV=PREVVERT(IP)
                IF (ZPREV.EQ.IPNEXT) THEN
                  COLLECTED(IP)=1
                  PREVVERT(IP)=IPPREV
                  IF (PREVVERT(IP).EQ.NEXTVERT(IP)) THEN
                    STATUS(IP)=1
                  END IF
                END IF
              END IF
            END DO
          END DO

          ISFINISHED=1
          DO IP=1,MNP
            IF ((COLLECTED(IP).EQ.0).AND.(STATUS(IP).EQ.0)) THEN
              STATUS(IP)=-1
              NEIGHBOR_PREV(IP)=PREVVERT(IP)     ! new code
            END IF
            IF (STATUS(IP).EQ.0) THEN
              ISFINISHED=0
            END IF
          END DO
          IF (ISFINISHED.EQ.1) THEN
            EXIT
          END IF
        END DO
!  Now making checks
        DO IP=1,MNP
          IP2=NEIGHBOR_NEXT(IP)
          IF (IP2.GT.0) THEN
            IP3=NEIGHBOR_PREV(IP2)
            IF (ABS(IP3 - IP).GT.0) THEN
              WRITE(*,*) 'IP=', IP, ' IP2=', IP2, ' IP3=', IP3
              WRITE(*,*) 'We have a dramatic inconsistency'
              STOP 'wmscrpmd, case 3'
            END IF
          END IF
        END DO
!   Now assigning the boundary IOBP array
        DO IP=1,MNP
          IF (STATUS(IP).EQ.-1 .AND. IOBP(IP) .EQ. 0) THEN
            IOBP(IP)=1
          END IF
        END DO

        DEALLOCATE(STATUS, STAT=ISTAT)
        CALL CHECK_DEALLOC_STATUS ( ISTAT )
        DEALLOCATE(COLLECTED, STAT=ISTAT)
        CALL CHECK_DEALLOC_STATUS ( ISTAT )
        DEALLOCATE(NEXTVERT, STAT=ISTAT)
        CALL CHECK_DEALLOC_STATUS ( ISTAT )
        DEALLOCATE(PREVVERT, STAT=ISTAT)
        CALL CHECK_DEALLOC_STATUS ( ISTAT )
        
      END SUBROUTINE
!/ ------------------------------------------------------------------- /
      SUBROUTINE FIX_PERIODCITY(PT)

!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  | Steven Brus                       |
!/                  | Ali Abdolali                      |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         21-May-2020 |
!/                  +-----------------------------------+
!/
!/    21-May-2020 : Origination.                        ( version 6.07 )
!/    
!/
!  1. Purpose :
!
!     Adjust element longitude coordinates for elements straddling the
!     dateline with distance of ~360 degrees
!
!  2. Method :
!
!     Detect if element has nodes on both sides of dateline and adjust
!     coordinates so that all nodes have the same sign
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
      IMPLICIT NONE
      REAL*8, INTENT(INOUT) :: PT(3,2)
!     ----------------------------------------------------------------
!
!     Local variables.
!     ----------------------------------------------------------------
      INTEGER :: I
      INTEGER :: R1GT180, R2GT180, R3GT180
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!

!  5. Called by :
!
!      Name         Type  Module   Description
!     ----------------------------------------------------------------
!      GET_SCRIP_INFO_UNSTRUCTURED Subr. WMSCRPMD  Element center calculation 
!      GET_SCRIP_INFO              Subr. WMSCRPMD  Check signs 
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!       None.
!
!  7. Remarks :
!
!  8. Structure :
!
!  9. Switches :
!
! 10. Source code :     
!/ ------------------------------------------------------------------- /

       R1GT180 = MERGE(1, 0, ABS(PT(3,1)-PT(2,1)).GT.180)
       R2GT180 = MERGE(1, 0, ABS(PT(1,1)-PT(3,1)).GT.180)
       R3GT180 = MERGE(1, 0, ABS(PT(2,1)-PT(1,1)).GT.180)
! if R1GT180+R2GT180+R3GT180 .eq. 0 the element does not cross the
! dateline
! if R1GT180+R2GT180+R3GT180 .eq. 1 the element contains the pole
! if R1GT180+R2GT180+R3GT180 .eq. 2 the element crosses the dateline

      IF ( R1GT180 + R2GT180 == 2 ) THEN
             PT(3,1)=PT(3,1)-SIGN(360.0d0,(PT(3,1)-PT(2,1)))
      ELSE IF ( R2GT180 + R3GT180 == 2 ) THEN
             PT(1,1)=PT(1,1)-SIGN(360.0d0,(PT(1,1)-PT(2,1)))
      ELSE IF ( R1GT180 + R3GT180 == 2 ) THEN
             PT(2,1)=PT(2,1)-SIGN(360.0d0,(PT(2,1)-PT(3,1)))
      ENDIF

      RETURN
      END SUBROUTINE FIX_PERIODCITY

!/ ------------------------------------------------------------------- /
      SUBROUTINE EXTCDE ( IEXIT, UNIT, MSG, FILE, LINE, COMM )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         06-Jun-2018 |
!/                  +-----------------------------------+
!/
!/    06-Jan-1998 : Final FORTRAN 77                    ( version 1.18 )
!/    23-Nov-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!/    10-Dec-2014 : Add checks for allocate status      ( version 5.04 )
!/    11-Mar-2015 : Allow non-error exit (iexit=0)      ( version 5.04 )
!/    20-Jan-2017 : Add optional MPI communicator arg   ( version 6.02 )
!/    06-Jun-2018 : Add optional MPI                    ( version 6.04 )
!/
!  1. Purpose :
!
!     Perform a program stop with an exit code.
!
!     If exit code IEXIT=0, then it is not an error, but 
!     a stop has been requested by the calling routine:
!     wait for other processes in communicator to catch up.
!     
!     If exit code IEXIT.ne.0, then abort program w/out
!     waiting for other processes to catch up (important for example
!     when not all processes are used by WW3).
!
!  2. Method :
!
!     Machine dependent.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IEXIT   Int.   I   Exit code to be used.
!       UNIT    Int.   I   (optional) file unit to write error message
!       MSG     Str.   I   (optional) error message
!       FILE    Str.   I   (optional) name of source code file
!       LINE    Int.   I   (optional) line number in source code file
!       COMM    Int.   I   (optional) MPI communicator
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!  5. Called by :
!
!     Any.
!
!  9. Switches :
!
!     !/MPI  MPI finalize interface if active
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      IMPLICIT NONE
!
!/MPI      INCLUDE "mpif.h"
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN) :: IEXIT
      INTEGER,      INTENT(IN), OPTIONAL :: UNIT
      CHARACTER(*), INTENT(IN), OPTIONAL :: MSG
      CHARACTER(*), INTENT(IN), OPTIONAL :: FILE
      INTEGER,      INTENT(IN), OPTIONAL :: LINE
      INTEGER,      INTENT(IN), OPTIONAL :: COMM
!/
!/ ------------------------------------------------------------------- /
!/
!/MPI      INTEGER                 :: IERR_MPI
!/MPI      LOGICAL                 :: RUN
      INTEGER                 :: IUN
      CHARACTER(256)          :: LMSG = ""
      CHARACTER(6)            :: LSTR
      CHARACTER(10)           :: PREFIX = "WW3 ERROR:"
!/
!/ Set file unit for error output
!/
      IUN = 0
      IF (PRESENT(UNIT)) IUN = UNIT
!/
!/ Report error message
!/
      IF (PRESENT(MSG)) THEN
          WRITE (IUN,"(A)") PREFIX//" "//TRIM(MSG)
        END IF
!/
!/ Report context
!/
      IF ( PRESENT(FILE) ) THEN
          LMSG = TRIM(LMSG)//" FILE="//TRIM(FILE)
        END IF
      IF ( PRESENT(LINE) ) THEN
          WRITE (LSTR,'(I0)') LINE
          LMSG = TRIM(LMSG)//" LINE="//TRIM(LSTR)
        END IF
      IF ( LEN_TRIM(LMSG).GT.0 ) THEN
          WRITE (IUN,"(A)") PREFIX//TRIM(LMSG)
        END IF
!/F90      CALL EXIT ( IEXIT )
!/DUM      STOP
!/
!/ End of EXTCDE ----------------------------------------------------- /
!/
      END SUBROUTINE EXTCDE

   SUBROUTINE CHECK_ALLOC_STATUS( STAT ) 
   INTEGER, INTENT(IN) :: STAT
   IF ( STAT .NE. 0 ) &
   CALL EXTCDE ( 99, MSG="ALLOCATE FAILED")
   END SUBROUTINE CHECK_ALLOC_STATUS

   SUBROUTINE CHECK_DEALLOC_STATUS( STAT ) 
   INTEGER, INTENT(IN) :: STAT
   IF ( STAT .NE. 0 ) &
   CALL EXTCDE ( 99, MSG="DEALLOCATE FAILED")
   END SUBROUTINE CHECK_DEALLOC_STATUS

!/
!/ End of module WMSCRPMD -------------------------------------------- /
!/
      END MODULE WMSCRPMD
