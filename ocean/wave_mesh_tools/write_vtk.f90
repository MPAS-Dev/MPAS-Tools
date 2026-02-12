      MODULE write_vtk

      USE globals, ONLY: rp

      IMPLICIT NONE

      CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE write_vtk_file(filename,nn,xy,ne,ect,depth,nope,obseg,obnds,nbou,fbseg,fbnds,nsnap,solution)

      IMPLICIT NONE
      
      CHARACTER(*), INTENT(IN) :: filename
      INTEGER, INTENT(IN) :: nn
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: xy
      INTEGER, INTENT(IN) :: ne
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ect
      REAL(rp), DIMENSION(:), INTENT(IN), OPTIONAL :: depth
      INTEGER, INTENT(IN),OPTIONAL :: nope
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: obseg 
      INTEGER, DIMENSION(:,:), INTENT(IN), OPTIONAL :: obnds
      INTEGER, INTENT(IN), OPTIONAL :: nbou
      INTEGER, INTENT(IN), DIMENSION(:,:), OPTIONAL :: fbseg
      INTEGER, INTENT(IN), DIMENSION(:,:), OPTIONAL :: fbnds
      INTEGER, INTENT(IN), OPTIONAL :: nsnap
      REAL(rp), INTENT(IN), DIMENSION(:,:), OPTIONAL :: solution 

      INTEGER :: open_bou_flag
      INTEGER :: flow_bou_flag
      INTEGER :: depth_flag
      INTEGER :: solution_flag

      CHARACTER(100) :: tmp1,tmp2,tmp3
      CHARACTER(3) :: tmp4
      INTEGER :: i,j
      INTEGER :: tbnds
      INTEGER :: snap
      INTEGER :: R3

      open_bou_flag = 0
      IF (PRESENT(nope) .and. PRESENT(obseg) .and. PRESENT(obnds)) THEN
        open_bou_flag = 1
      ENDIF

      flow_bou_flag = 0
      IF (PRESENT(nbou) .and. PRESENT(fbseg) .and. PRESENT(fbnds)) THEN
        flow_bou_flag = 1
      ENDIF

      depth_flag = 0
      IF (PRESENT(depth)) THEN
        depth_flag = 1
      ENDIF

      solution_flag = 0
      IF (PRESENT(solution)) THEN
        solution_flag = 1
      ENDIF

      OPEN(UNIT=11,FILE=TRIM(filename))

      ! Write header
      WRITE(11,"(A)") "# vtk DataFile Version 3.0"
      WRITE(11,"(A)") filename 
      WRITE(11,"(A)") "ASCII"
      WRITE(11,"(A)") "DATASET UNSTRUCTURED_GRID"
      WRITE(tmp1,"(I10)") nn
      WRITE(11,"(A,1x,A,1x,A)") "POINTS", TRIM(ADJUSTL(tmp1)), "double"

      IF (SIZE(xy,1) == 3) THEN
        R3 = 1
      ELSE IF (SIZE(xy,1) == 2) THEN
        R3 = 0
      ELSE
        PRINT*, "Incorrect xy size"
        STOP
      ENDIF
   
      ! Write point coordinates
      WRITE(tmp3,"(I4)") 0
      DO i = 1,nn        
        WRITE(tmp1,"(F20.14)") xy(1,i)
        WRITE(tmp2,"(F20.14)") xy(2,i)
        IF (R3 == 1) THEN
          WRITE(tmp3,"(F20.14)") xy(3,i)
        ENDIF
        WRITE(11,"(A,1x,A,1x,A)") TRIM(ADJUSTL(tmp1)),TRIM(ADJUSTL(tmp2)),TRIM(ADJUSTL(tmp3))
      ENDDO


      ! Coutn total boundary nodes
      tbnds = 0
      IF (open_bou_flag == 1) THEN
        DO i = 1,nope
          tbnds = tbnds + obseg(i)-1
        ENDDO
      ENDIF
      IF (flow_bou_flag == 1) THEN
        DO i = 1,nbou
          tbnds = tbnds + fbseg(1,i)-1
        ENDDO
      ENDIF

      ! Write boundary and element connectivity
      WRITE(tmp1,"(I10)") ne+tbnds
      WRITE(tmp2,"(I10)") 4*(ne+tbnds)
      WRITE(11,"(A,1x,A,1x,A)") "CELLS",TRIM(ADJUSTL(tmp1)), TRIM(ADJUSTL(tmp2))
      IF (open_bou_flag == 1) THEN
        DO i = 1,nope
          DO j = 1,obseg(i)-1
            WRITE(tmp1,"(I10)") obnds(j,i)-1
            WRITE(tmp2,"(I10)") obnds(j+1,i)-1
            WRITE(11,"(I1,1x,A,1x,A)") 2,TRIM(ADJUSTL(tmp1)),TRIM(ADJUSTL(tmp2))
          ENDDO
        ENDDO
      ENDIF
      IF (flow_bou_flag == 1) THEN
        DO i = 1,nbou
          DO j = 1,fbseg(1,i)-1
            WRITE(tmp1,"(I10)") fbnds(j,i)-1
            WRITE(tmp2,"(I10)") fbnds(j+1,i)-1
            WRITE(11,"(I1,1x,A,1x,A)") 2,TRIM(ADJUSTL(tmp1)),TRIM(ADJUSTL(tmp2))
          ENDDO
        ENDDO
      ENDIF
      DO i = 1,ne
        WRITE(tmp1,"(I10)") ect(1,i)-1
        WRITE(tmp2,"(I10)") ect(2,i)-1
        WRITE(tmp3,"(I10)") ect(3,i)-1
        WRITE(11,"(I1,1x,A,1x,A,1x,A)") 3,TRIM(ADJUSTL(tmp1)), TRIM(ADJUSTL(tmp2)), TRIM(ADJUSTL(tmp3))
      ENDDO
    
      ! Write cell types 
      WRITE(tmp1,"(I10)") ne+tbnds
      WRITE(11,"(A,1x,A)") "CELL_TYPES", TRIM(ADJUSTL(tmp1))
      IF (open_bou_flag == 1) THEN
        DO i = 1,nope
          DO j = 1,obseg(i)-1
            WRITE(11,"(I1)") 3
          ENDDO
        ENDDO
      ENDIF
      IF (flow_bou_flag == 1) THEN
        DO i = 1,nbou
          DO j = 1,fbseg(1,i)-1
            WRITE(11,"(I1)") 3
          ENDDO
        ENDDO
      ENDIF
      DO i = 1,ne
        WRITE(11,"(I1)") 5
      ENDDO

      ! Write depth
      IF (depth_flag == 1) THEN
        WRITE(tmp1,"(I10)") nn
        WRITE(11,"(A,1x,A)") "POINT_DATA", TRIM(ADJUSTL(tmp1))
        WRITE(11,"(A)") "SCALARS bathymetry float 1"
        WRITE(11,"(A)") "LOOKUP_TABLE default"
        DO i = 1,nn
          WRITE(tmp1,"(F24.7)") depth(i)
          WRITE(11,"(A)") TRIM(ADJUSTL(tmp1))
        ENDDO
      ENDIF
      
      IF (solution_flag == 1) THEN
        WRITE(tmp1,"(I10)") nn
        WRITE(11,"(A,1x,A)") "POINT_DATA", TRIM(ADJUSTL(tmp1))
        DO snap = 1,nsnap
          WRITE(tmp1,"(I10)") snap 
          tmp4 = '000'
          WRITE(tmp4,"(I3.3)") snap
          WRITE(11,"(A)") "SCALARS t="//TRIM(ADJUSTL(tmp4))//" float 1"
          WRITE(11,"(A)") "LOOKUP_TABLE default"
          DO i = 1,nn
            WRITE(tmp1,"(F24.7)") solution(i,snap)
            WRITE(11,"(A)") TRIM(ADJUSTL(tmp1))
          ENDDO
        ENDDO
        !WRITE(tmp1,"(I10)") nsnap 
        !WRITE(11,"(A,1x,A)") "FIELD solution", TRIM(ADJUSTL(tmp1))
        !DO snap = 1,nsnap
        !  WRITE(tmp1,"(I10)") snap 
        !  WRITE(tmp2,"(I10)") nn
        !  WRITE(11,"(A)") "t="//TRIM(ADJUSTL(tmp1))//" 1 "//TRIM(ADJUSTL(tmp2))//" float"
        !  DO i = 1,nn
        !    WRITE(tmp1,"(F24.7)") solution(i,snap)
        !    WRITE(11,"(A)") TRIM(ADJUSTL(tmp1))
        !  ENDDO
        !ENDDO
      ENDIF

      CLOSE(11)

      RETURN
      END SUBROUTINE write_vtk_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      END MODULE write_vtk

