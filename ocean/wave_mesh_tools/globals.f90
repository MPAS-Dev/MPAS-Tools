      MODULE globals

      IMPLICIT NONE

      INTEGER, PARAMETER :: rp = kind(1d0)
      INTEGER, DIMENSION(4), PARAMETER :: nverts = (/3,4,3,4/)

      TYPE :: boundary_type
        INTEGER :: nnds
        INTEGER :: bou_type
        INTEGER :: flag
        REAL(rp), DIMENSION(:,:), ALLOCATABLE :: xy
        INTEGER, DIMENSION(:), ALLOCATABLE :: nds
      END TYPE

      END MODULE globals
