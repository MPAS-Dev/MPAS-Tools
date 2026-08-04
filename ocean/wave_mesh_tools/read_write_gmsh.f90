     MODULE read_write_gmsh

     USE grid_file_mod, ONLY: grid_type
     USE globals, ONLY: rp


     CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     SUBROUTINE read_gmsh_file(filename,mesh)

     IMPLICIT NONE

     CHARACTER(*), INTENT(IN) :: filename
     TYPE(grid_type), INTENT(OUT) :: mesh

     CALL read_header(filename)
     CALL read_nodes(mesh%nn,mesh%xy,mesh%depth)
     CALL read_elements(mesh%ne,mesh%ect,mesh%el_type)

     RETURN
     END SUBROUTINE read_gmsh_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     SUBROUTINE read_header(filename)

     IMPLICIT NONE

     CHARACTER(*), INTENT(IN) :: filename
     
     CHARACTER(20) :: comment

     OPEN(UNIT=14,FILE=filename)
     READ(14,*) comment
     READ(14,*) comment
     READ(14,*) comment

     RETURN
     END SUBROUTINE read_header

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
     SUBROUTINE read_nodes(nn,xy,depth)

     IMPLICIT NONE

     INTEGER, INTENT(OUT) :: nn
     REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: xy
     REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: depth
     
     CHARACTER(20) :: comment
     INTEGER :: i
     INTEGER :: j
     
     READ(14,*) comment
     READ(14,*) nn

     ALLOCATE(xy(2,nn))
     ALLOCATE(depth(nn))

     DO i = 1,nn
       READ(14,*), j, xy(1,j), xy(2,j), depth(j)
     ENDDO

     READ(14,*) comment

     RETURN
     END SUBROUTINE read_nodes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     SUBROUTINE read_elements(ne,ect,el_type)

     IMPLICIT NONE

     INTEGER, INTENT(OUT) :: ne
     INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: ect
     INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: el_type

     CHARACTER(20) :: comment
     INTEGER :: i
     INTEGER :: j
     INTEGER :: k
     INTEGER :: tag,ntags
     INTEGER, DIMENSION(15) :: lookup

     lookup(15) = 1
     lookup(2) = 3

     READ(14,*) comment
     READ(14,*) ne

     ALLOCATE(ect(3,ne))
     ALLOCATE(el_type(ne))

     DO i = 1,ne
       READ(14,*) j,el_type(j),ntags,(tag,k=1,ntags),(ect(k,j),k=1,lookup(el_type(j)))
     ENDDO

     READ(14,*) comment

     RETURN
     END SUBROUTINE read_elements

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     SUBROUTINE write_gmsh_file(filename,mesh)

     IMPLICIT NONE
     CHARACTER(*), INTENT(IN) :: filename
     TYPE(grid_type), INTENT(IN) :: mesh

     CALL write_header(filename)
     CALL write_nodes(mesh%nn,mesh%xy,mesh%depth)
     CALL write_elements(mesh%ne,mesh%ect,mesh%el_type)

     RETURN
     END SUBROUTINE write_gmsh_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     SUBROUTINE write_header(filename)

     IMPLICIT NONE
     CHARACTER(*), INTENT(IN) :: filename
     
     OPEN(UNIT=14,FILE=filename)
     WRITE(14,"(A)") "$MeshFormat"
     WRITE(14,"(A)") "2 0 8"
     WRITE(14,"(A)") "$EndMeshFormat"

     RETURN
     END SUBROUTINE write_header

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     SUBROUTINE write_nodes(nn,xy,depth)

     IMPLICIT NONE

     INTEGER, INTENT(IN) :: nn
     REAL(rp), DIMENSION(:,:), INTENT(IN) ::  xy
     REAL(rp), DIMENSION(:), INTENT(IN) ::  depth

     INTEGER :: i

     WRITE(14,"(A)") "$Nodes"
     WRITE(14,"(I0)") nn
     DO i = 1,nn
       WRITE(14,"(I0,3(F16.8))") i, xy(1,i), xy(2,i), depth(i)
     ENDDO
     WRITE(14,"(A)") "$EndNodes"

     RETURN
     END SUBROUTINE write_nodes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     SUBROUTINE write_elements(ne,ect,el_type)

     IMPLICIT NONE

     INTEGER, INTENT(IN) :: ne
     INTEGER, DIMENSION(:,:), INTENT(IN) :: ect
     INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: el_type
   
     INTEGER :: i,j

     WRITE(14,"(A)") "$Elements"
     WRITE(14,"(I0)") ne
     DO i = 1,ne
       WRITE(14,"(9(I0,2x))") i, 2, 3, 0, i, 0, ect(1,i), ect(2,i), ect(3,i)
     ENDDO
     WRITE(14,"(A)") "$EndElements"

     CLOSE(14)

     RETURN
     END SUBROUTINE write_elements

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     END MODULE read_write_gmsh
