      MODULE rotate_mod

      USE NETCDF
      USE globals, ONLY: rp

      REAL (rp), PARAMETER :: PI = 4d0*atan(1d0)
      REAL (rp), PARAMETER :: DERA = PI/180d0
      REAL (rp), PARAMETER :: RECIP_PI_OVER_180 = 180d0 / PI
      REAL (rp), PARAMETER :: PI_OVER_180   = PI / 180d0
      REAL (rp), PARAMETER :: DEG2RAD = PI / 180d0
      REAL (rp), PARAMETER :: RAD2DEG = 180d0 / PI
      REAL (rp) :: UNDEF = -999.9d0

      CONTAINS 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Given a radius vector (\varphi,\theta) (in rad)
      ! of z', the new z after rotation, return
      ! a transformation matrix:
      !        RTOT: (x,y,z) -> (x',y',z')  
      SUBROUTINE GET_ROTMAT_ZNVEC( RTOT, VARPHIN, THETAN )
        IMPLICIT NONE

        REAL (rp):: RTOT(:,:)
        REAL (rp):: VARPHIN, THETAN, PHIN, FPI
        REAL (rp) :: XROT(3,3),YROT(3,3),ZROT(3,3)

        YROT(1,1) = sin( THETAN ) 
        YROT(2,1) = 0d0
        YROT(3,1) = -cos( THETAN )

        YROT(1,2) = 0d0
        YROT(2,2) = 1d0
        YROT(3,2) = 0d0
       
        YROT(1,3) = cos( THETAN )
        YROT(2,3) = 0d0
        YROT(3,3) = sin( THETAN )
 
        !!!!!!!!!!!!!!!!!!!!!!!!!!

        PHIN = PI - VARPHIN 

        ZROT(1,1) = cos( PHIN )
        ZROT(2,1) = sin( PHIN )
        ZROT(3,1) = 0d0

        ZROT(1,2) = -sin( PHIN )
        ZROT(2,2) = cos( PHIN )
        ZROT(3,2) = 0d0
       
        ZROT(1,3) = 0d0
        ZROT(2,3) = 0d0
        ZROT(3,3) = 1d0

        RTOT = MATMUL(YROT,ZROT)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !RTOT = 0.0_rp ;

        !RTOT(1,1) = sin( VARPHIN ) ;
        !RTOT(2,1) = sin( THETAN )*cos( VARPHIN ) ;
        !RTOT(3,1) = cos( THETAN )*cos( VARPHIN ) ;

        !RTOT(1,2) = -cos( VARPHIN ) ;
        !RTOT(2,2) = sin( THETAN )*sin( VARPHIN ) ;
        !RTOT(3,2) = cos( THETAN )*sin( VARPHIN ) ;

        !RTOT(1,3) = 0.0_rp ;
        !RTOT(2,3) = -cos( THETAN ) ;
        !RTOT(3,3) =  sin( THETAN ) ;

        RETURN ;
      END SUBROUTINE GET_ROTMAT_ZNVEC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !
      ! Transformation matrix mapping a vector from the Cartesian to spherical
      ! coordinates 
      !
      ! [e_{r},e_{varphi},e_{theta}]^{T} = TM [ i, j, k]^{T}
      FUNCTION CT2SP_VECTRANSMAT( VARPHI, THETA ) RESULT(TM)
        IMPLICIT NONE

        REAL (rp) :: TM(3,3)
        REAL (rp), intent(in):: VARPHI, THETA


        TM(1,1) = cos( THETA )*cos( VARPHI ) ;
        TM(2,1) = -sin( VARPHI ) ;
        TM(3,1) = -sin( THETA )*cos( VARPHI ) ;

        TM(1,2) = cos( THETA )*sin( VARPHI ) ;
        TM(2,2) = cos( VARPHI ) ;
        TM(3,2) = -sin( THETA )*sin( VARPHI ) ;

        TM(1,3) = sin( THETA ) ;
        TM(2,3) = 0.0_rp ;
        TM(3,3) = cos( THETA ) ;


        RETURN ;
      END FUNCTION CT2SP_VECTRANSMAT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE SPCOORSROTS1( RTOTS, LONR, LATR, LONO, LATO, NN )
        IMPLICIT NONE

        !c dummy c!
        REAL (rp), INTENT(IN) :: RTOTS(3,3)
        REAL (rp), dimension(:), INTENT(OUT) :: LONR, LATR
        REAL (rp), dimension(:), INTENT(IN) :: LONO, LATO
        INTEGER, INTENT(IN) :: NN

        !c local c!
        INTEGER:: IP
        REAL (rp):: XP(3), XPR(3), LLO, LTO


        DO IP = 1, NN
           LLO = LONO( IP ) ;
           LTO = LATO( IP ) ;

           XP(1) = cos( LTO )*cos( LLO ) ;
           XP(2) = cos( LTO )*sin( LLO ) ;
           XP(3) = sin( LTO ) ;

           XPR = MATMUL( RTOTS, XP) ;

           LATR(IP) = atan2( XPR(3), SQRT(XPR(1)*XPR(1) &
                                       + XPR(2)*XPR(2)) )
           LONR(IP) = atan2( XPR(2), XPR(1) ) ;
        END DO

        RETURN ;
      END SUBROUTINE SPCOORSROTS1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! (e_{r},e_{\varphi},e_{\theta})^{T} = T (i,j,k)^{T}
      !  x' = R x
      ! (e_{r'},e_{\varphi'},e_{\theta'})^{T} = T' (i',j',k')^{T}
      !
      ! V'(r',\varphi',\theta') = T' R T^{T} V(r,\varphi,\theta)
      !
      ! RVELF = T' R T^{T}
      ! RVELI = (T')^{R} R^{T} T
      !
      ! Evaluated the forward and inverse transformation matrices
      ! between two spherical coordinate
      !     RVELF: VEC(LONO,LATO) -> VEC(LONR,LATR)
      !     RVELI: VEC(LONR,LATR) -> VEC(LONO,LATO)
      !
      !  lonr,latr in rad
      !
      SUBROUTINE SPVECROTSMAT( RVELF, RVELI, RTOTS, &
                               LONR, LATR, LONO, LATO, NN )
         IMPLICIT NONE

         ! Dummys !
         REAL (rp), DIMENSION(:,:,:), INTENT(OUT) :: RVELF, RVELI

         INTEGER:: NN
         REAL (rp), INTENT(IN) :: RTOTS(3,3)
         REAL (rp), DIMENSION(:), INTENT(IN) :: LONR, LATR, LONO, LATO

         ! local !
         INTEGER:: IP
         REAL (rp), DIMENSION(3,3):: TP, TM, TRT


         DO IP = 1, NN
           TRT = 0.0_rp ;

           ! T' = T(varphi',theta')
           TP = CT2SP_VECTRANSMAT( LONR(IP), LATR(IP) ) ;

           ! T = T(varphi,theta)
           TM = CT2SP_VECTRANSMAT( LONO(IP), LATO(IP) ) ;
           TM = TRANSPOSE( TM ) ; ! T^{T}

           ! T' R T^{T}
           TRT = MATMUL( RTOTS, TM ) ;
           TRT = MATMUL( TP, TRT ) ;

           RVELF(1:2,1:2,IP) = TRT(2:3,2:3) ; ! forward transform
           RVELI(:,:,IP) = TRANSPOSE( RVELF(:,:,IP) ) ; ! inverse transform
         END DO

         RETURN ;
      END SUBROUTINE SPVECROTSMAT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !    
      ! Map two vectors in spherical coordinates 
      !    
      !    VR(:,II) = SPRV(:,:,II)*VO(:,II)
      SUBROUTINE MAP2DSPVECS( VPHIR, VTHETAR, VPHIO, VTHETAO, SPRV, NN ) 
        IMPLICIT NONE 

        ! dummy !
        INTEGER, INTENT(IN) :: NN
        REAL (rp), INTENT(IN) :: SPRV(:,:,:)
        REAL (rp), DIMENSION(:), INTENT(OUT) :: VPHIR, VTHETAR
        REAL (rp), DIMENSION(:), INTENT(IN) :: VPHIO, VTHETAO

        ! local !
        INTEGER:: II
        REAL (rp):: VO(2), VR(2)

        DO  II = 1, NN
           VO = (/ VPHIO(II), VTHETAO(II) /)  ;
           VR = MATMUL( SPRV(1:2,1:2,II), VO ) ;  
     
           VPHIR(II) = VR(1) ;
           VTHETAR(II) = VR(2) ;
        END DO
     
        RETURN ;
      END SUBROUTINE MAP2DSPVECS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE CHECK_RTOTMAT( RTOTS, IERR )
        IMPLICIT NONE 

        INTEGER:: IERR 
        REAL (rp), intent(in):: RTOTS(3,3)

        !c local c!
        REAL (rp):: RTT(3,3), TR, SOFFD

        RTT = TRANSPOSE( RTOTS ) ;  

        RTT = MATMUL(RTT, RTOTS ) ;  

        ! trace of a matrix RTT 
        TR = RTT(1,1) + RTT(2,2) + RTT(3,3) ; 
     
        RTT(1,1) = 0.0 ;
        RTT(2,2) = 0.0 ;
        RTT(3,3) = 0.0 ;
     
        SOFFD = SUM( MATMUL(ABS(RTT), (/ 1.0_rp, 1.0_rp, 1.0_rp /)) )  ; 

        IERR = 0 ;
        IF ( SOFFD > 1.0e-10_rp ) THEN 
           IERR = 1 ;  
     
           PRINT*, "Error: a given rotation matrix is not orthgonal and thus invalid" ;
           STOP
        END IF

        RETURN ;
      END SUBROUTINE CHECK_RTOTMAT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Li  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!Li
!Li  Merged UM source code for rotated grid, consisting the following
!Li  original subroutines in UM 6.1
!Li    LLTOEQ1A  WCOEFF1A  and  LBCROTWINDS1
!Li  The last subroutine is modified to process only one level winds
!Li  cpp directives are removed and required header C_Pi.h inserted.
!Li         Jian-Guo Li     26 May 2005
!Li
!Li  The WCOEFF1A subroutine is merged into LLTOEQ to reduce repetition
!Li  of the same calculations. Subroutine interface changed to 
!Li  LLTOEQANGLE
!Li         Jian-GUo Li     23 Aug 2005
!Li
!Li  Subroutine W3LLTOEQ   --------------------------------------------    
!Li                                                                        
!Li  Purpose:  Calculates latitude and longitude on equatorial             
!Li            latitude-longitude (eq) grid used in regional               
!Li            models from input arrays of latitude and                    
!Li            longitude on standard grid. Both input and output           
!Li            latitudes and longitudes are in degrees.                    
!Li            Also calculate rotation angle in degree to tranform
!Li            standard wind velocity into equatorial wind.
!Li            Valid for 0<PHI_POLE<90 or new pole in N. hemisphere.
!Li                                                                        
!* Arguments:--------------------------------------------------------    
      SUBROUTINE W3LLTOEQ ( PHI, LAMBDA, PHI_EQ, LAMBDA_EQ,     &              
                      ANGLED, PHI_POLE, LAMBDA_POLE, POINTS )             
                                                                           
      IMPLICIT NONE                                                        
                                                                           
      INTEGER:: POINTS    !IN  Number of points to be processed             

      REAL (rp) :: PHI_POLE,  & !IN  Latitude of equatorial lat-lon pole
                  LAMBDA_POLE  !INOUT  Longitude of equatorial lat-lon pole
                                                                           
      REAL (rp) , DIMENSION(POINTS) ::         &
                  PHI,       & !IN  Latitude
                  LAMBDA,    & !IN  Longitude
                  ANGLED,    & !OUT turning angle in deg for standard wind
                  LAMBDA_EQ, & !OUT Longitude in equatorial lat-lon coords
                  PHI_EQ       !OUT Latitude in equatorial lat-lon coords

! Define local varables:-----------------------------------------------
      REAL(rp) :: A_LAMBDA, A_PHI, E_LAMBDA, E_PHI,                 &
                      SIN_PHI_POLE, COS_PHI_POLE,                       &
                      TERM1, TERM2, ARG, LAMBDA_ZERO, LAMBDA_POLE_KEEP
      INTEGER      :: I 

      REAL(rp), PARAMETER :: SMALL=1.0D-6


!*----------------------------------------------------------------------   

! 1. Initialise local constants
! Scale lambda pole to range -180 to 180 degs
      LAMBDA_POLE_KEEP=LAMBDA_POLE
      IF (LAMBDA_POLE.LE.-180.D0) LAMBDA_POLE=LAMBDA_POLE+360.D0
      IF (LAMBDA_POLE.GT. 180.D0) LAMBDA_POLE=LAMBDA_POLE-360.D0

! Latitude of zeroth meridian
      LAMBDA_ZERO=LAMBDA_POLE+180.D0
! Sine and cosine of latitude of eq pole
      IF (PHI_POLE >= 0.0) THEN
        SIN_PHI_POLE =  SIN(PI_OVER_180*PHI_POLE)
        COS_PHI_POLE =  COS(PI_OVER_180*PHI_POLE)
      ELSE
        SIN_PHI_POLE = -SIN(PI_OVER_180*PHI_POLE)
        COS_PHI_POLE = -COS(PI_OVER_180*PHI_POLE)
      ENDIF

! 2. Transform from standard to equatorial latitude-longitude

      DO I= 1, POINTS

! Scale longitude to range -180 to +180 degs

        A_LAMBDA=LAMBDA(I)-LAMBDA_ZERO
        IF(A_LAMBDA.GT. 180.D0) A_LAMBDA=A_LAMBDA-360.D0
        IF(A_LAMBDA.LE.-180.D0) A_LAMBDA=A_LAMBDA+360.D0

! Convert latitude & longitude to radians

        A_LAMBDA=PI_OVER_180*A_LAMBDA
        A_PHI=PI_OVER_180*PHI(I)

! Compute eq latitude using equation (4.4)

        ARG=-COS_PHI_POLE*COS(A_PHI)*COS(A_LAMBDA)   &
            +SIN_PHI_POLE*SIN(A_PHI)
        ARG=MIN(ARG, 1.D0)
        ARG=MAX(ARG,-1.D0)
        E_PHI=ASIN(ARG)
        PHI_EQ(I)=RECIP_PI_OVER_180*E_PHI

! Compute eq longitude using equation (4.6)

        TERM1 = SIN_PHI_POLE*COS(A_PHI)*COS(A_LAMBDA)   &
               +COS_PHI_POLE*SIN(A_PHI)
        TERM2 = COS(E_PHI)
        IF(TERM2 .LT. SMALL) THEN
          E_LAMBDA=0.D0
        ELSE
          ARG=TERM1/TERM2
          ARG=MIN(ARG, 1.D0)
          ARG=MAX(ARG,-1.D0)
          E_LAMBDA=RECIP_PI_OVER_180*ACOS(ARG)
          E_LAMBDA=SIGN(E_LAMBDA,A_LAMBDA)
        ENDIF

! Scale longitude to range 0 to 360 degs

        IF(E_LAMBDA.GE.360.D0) E_LAMBDA=E_LAMBDA-360.D0
        IF(E_LAMBDA.LT.  0.D0) E_LAMBDA=E_LAMBDA+360.D0
        LAMBDA_EQ(I)=E_LAMBDA

!Li  Calculate turning angle for standard wind velocity

        E_LAMBDA=PI_OVER_180*LAMBDA_EQ(I)

! Formulae used are from eqs (4.19) and (4.21)

        TERM2=SIN(E_LAMBDA)
        ARG= SIN(A_LAMBDA)*TERM2*SIN_PHI_POLE      &
            +COS(A_LAMBDA)*COS(E_LAMBDA)
        ARG=MIN(ARG, 1.D0)
        ARG=MAX(ARG,-1.D0)
        TERM1=RECIP_PI_OVER_180*ACOS(ARG)
        ANGLED(I)=SIGN(TERM1,TERM2)
!Li

      ENDDO

! Reset Lambda pole to the setting on entry to subroutine
      LAMBDA_POLE=LAMBDA_POLE_KEEP

      RETURN
      END SUBROUTINE W3LLTOEQ 
!Li  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!Li
!Li  Merged UM source code for rotated grid, consiting the following
!Li  original subroutines in UM 6.1
!Li    EQTOLL1A  WCOEFF1A  and  LBCROTWINDS1
!Li  The last subroutine is modified to process only one level winds
!Li  cpp directives are removed and required header C_Pi.h inserted.
!Li         Jian-Guo Li     26 May 2005
!Li
!Li  The WCOEFF1A subroutine is merged into EQTOLL to reduce repetition
!Li  of the same calculations. Subroutine interface changed to 
!Li  EQTOLLANGLE
!Li  First created:   Jian-GUo Li     23 Aug 2005
!Li  Last modified:   Jian-GUo Li     25 Feb 2008
!Li
!Li  Subroutine W3EQTOLL  --------------------------------------------
!Li
!Li  Purpose:  Calculates latitude and longitude on standard grid
!Li            from input arrays of latitude and longitude on
!Li            equatorial latitude-longitude (eq) grid used
!Li            in regional models. Both input and output latitudes
!Li            and longitudes are in degrees.
!Li            Also calculate rotation angle in degree to tranform
!Li            standard wind velocity into equatorial wind.
!Li            Valid for 0<PHI_POLE<90 or new pole in N. hemisphere.
!Li
!Li  Arguments:--------------------------------------------------------

      SUBROUTINE W3EQTOLL( PHI_EQ, LAMBDA_EQ, PHI, LAMBDA,   &
     &                 ANGLED, PHI_POLE, LAMBDA_POLE, POINTS )

      IMPLICIT NONE

      INTEGER:: POINTS      !IN  Number of points to be processed

      REAL(rp) :: PHI_POLE,   & !IN  Latitude of equatorial lat-lon pole
     &        LAMBDA_POLE   !IN  Longitude of equatorial lat-lon pole

      REAL(rp), DIMENSION(POINTS) ::         &
     &        PHI,       & !OUT Latitude
     &        LAMBDA,    & !OUT Longitude (0 =< LON < 360)
     &        ANGLED,    & !OUT turning angle in deg for standard wind
     &        LAMBDA_EQ, & !IN  Longitude in equatorial lat-lon coords
     &        PHI_EQ       !IN  Latitude in equatorial lat-lon coords

! Local varables:------------------------------------------------------
      REAL(rp) :: E_LAMBDA, E_PHI, A_LAMBDA, A_PHI,                 &
                      SIN_PHI_POLE, COS_PHI_POLE,                       &
                      TERM1, TERM2, ARG, LAMBDA_ZERO
      INTEGER :: I

      REAL(rp), PARAMETER :: SMALL=1.0E-6

! ----------------------------------------------------------------------

! 1. Initialise local constants

! Latitude of zeroth meridian
      LAMBDA_ZERO=LAMBDA_POLE+180.D0
! Sine and cosine of latitude of eq pole
      IF (PHI_POLE >= 0.0) THEN
        SIN_PHI_POLE =  SIN(PI_OVER_180*PHI_POLE)
        COS_PHI_POLE =  COS(PI_OVER_180*PHI_POLE)
      ELSE
        SIN_PHI_POLE = -SIN(PI_OVER_180*PHI_POLE)
        COS_PHI_POLE = -COS(PI_OVER_180*PHI_POLE)
      ENDIF

! 2. Transform from equatorial to standard latitude-longitude

      DO I= 1, POINTS

! Scale eq longitude to range -180 to +180 degs

        E_LAMBDA=LAMBDA_EQ(I)
        IF(E_LAMBDA.GT. 180.0) E_LAMBDA=E_LAMBDA-360.D0
        IF(E_LAMBDA.LT.-180.0) E_LAMBDA=E_LAMBDA+360.D0

! Convert eq latitude & longitude to radians

        E_LAMBDA=PI_OVER_180*E_LAMBDA
        E_PHI=PI_OVER_180*PHI_EQ(I)

! Compute latitude using equation (4.7)

        ARG=COS_PHI_POLE*COS(E_PHI)*COS(E_LAMBDA)    &
       &   +SIN_PHI_POLE*SIN(E_PHI)
        ARG=MIN(ARG, 1.D0)
        ARG=MAX(ARG,-1.D0)
        A_PHI=ASIN(ARG)
        PHI(I)=RECIP_PI_OVER_180*A_PHI

! Compute longitude using equation (4.8)

        TERM1 = COS(E_PHI)*SIN_PHI_POLE*COS(E_LAMBDA)   &
       &       -SIN(E_PHI)*COS_PHI_POLE
        TERM2 = COS(A_PHI)
        IF(TERM2.LT.SMALL) THEN
          A_LAMBDA=0.D0
        ELSE
          ARG=TERM1/TERM2
          ARG=MIN(ARG, 1.D0)
          ARG=MAX(ARG,-1.D0)
          A_LAMBDA=RECIP_PI_OVER_180*ACOS(ARG)
          A_LAMBDA=SIGN(A_LAMBDA,E_LAMBDA)
          A_LAMBDA=A_LAMBDA+LAMBDA_ZERO
        END IF

! Scale longitude to range 0 to 360 degs

        IF(A_LAMBDA.GE.360.0) A_LAMBDA=A_LAMBDA-360.D0
        IF(A_LAMBDA.LT.  0.0) A_LAMBDA=A_LAMBDA+360.D0
        LAMBDA(I)=A_LAMBDA

!Li  Calculate turning angle for standard wind velocity

        A_LAMBDA=PI_OVER_180*(LAMBDA(I)-LAMBDA_ZERO)

! Formulae used are from eqs (4.19) and (4.21)

        TERM2=SIN(E_LAMBDA)
        ARG=SIN(A_LAMBDA)*TERM2*SIN_PHI_POLE     &
       &           +COS(A_LAMBDA)*COS(E_LAMBDA)
        ARG=MIN(ARG, 1.D0)
        ARG=MAX(ARG,-1.D0)
        TERM1=RECIP_PI_OVER_180*ACOS(ARG)
        ANGLED(I)=SIGN(TERM1,TERM2)
!Li

      ENDDO

      RETURN
      END SUBROUTINE W3EQTOLL  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!/ ------------------------------------------------------------------- /
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3XYRTN ( NSEA, XVEC, YVEC, AnglD )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III            NOAA/NMC |
!/                  |             A. Saulter            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         01-Mar-2018 |
!/                  +-----------------------------------+
!/
!/    01-Mar-2018 : Added subroutine                   ( version 6.02 )
!
!  1. Purpose :
!     Subroutine to de-rotate x,y vectors from rotated to standard pole
!     reference system
!
!  2. Method:
!   Rotates x,y vectors anticlockwise by angle alpha in radians
!
!/ ------------------------------------------------------------------- /
      IMPLICIT NONE 
!
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN) :: NSEA        ! Number of sea points
      REAL(rp),    INTENT(IN) :: AnglD(NSEA) ! Turning angle (degrees)
      REAL(rp), INTENT(INOUT) :: XVEC(NSEA), YVEC(NSEA)
!
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER :: ISEA 
      REAL(rp)    :: XVTMP, YVTMP
!
!/ ------------------------------------------------------------------- /
! Apply the rotation
!
      DO ISEA=1, NSEA 
        IF (( XVEC(ISEA) .NE. UNDEF ) .AND. &
            ( YVEC(ISEA) .NE. UNDEF )) THEN 
          XVTMP = XVEC(ISEA)*COS(AnglD(ISEA)*DERA) + &
                   YVEC(ISEA)*SIN(AnglD(ISEA)*DERA)
          YVTMP = YVEC(ISEA)*COS(AnglD(ISEA)*DERA) - &
                   XVEC(ISEA)*SIN(AnglD(ISEA)*DERA)
          XVEC(ISEA) = XVTMP
          YVEC(ISEA) = YVTMP
        END IF
      END DO

      RETURN 
      END SUBROUTINE W3XYRTN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE check(status)

      IMPLICIT NONE
      INTEGER :: status

      IF(status /= NF90_NOERR) THEN
        PRINT("(A,A)"), "fatal error from ", TRIM(NF90_STRERROR(status))  
      ENDIF
    
      END SUBROUTINE check


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      END MODULE rotate_mod
