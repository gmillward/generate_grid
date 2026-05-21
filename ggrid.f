      SUBROUTINE GGRID (NVERT,GLAMN,GLAMX,GLOMN,GLOMX,ALTMN,ALTMX,
     +                 GPLAT,GPLON,GPALT,NLAT,NLON,NALT)
C          Written Aug 95 R. Barnes, NCAR
C          Part of the Apex magnetic coordinate interpolation routines, this
C          computes grid values using the method originally installed in
C          apxntrp.  (lat,lon,alt) limits and a grid resolution parameter
C          (NVERT) are input and coordinates are returned which bounds these
C          limits.
C          The grid values are some subset (or possibly all) of the
C          following grid, which is defined for all of space above Earth's
C          surface.  In the vertical, the quantity (Re/(Re+height)) is
C          divided into NVERT segments, where Re is the mean Earth radius.
C          In latitude, 180 degrees is divided into 3*NVERT segments.  In
C          longitude, 360 degrees is divided into 5*NVERT segments.  Only
C          the minimum portion of this larger grid that encompasses the
C          specified ranges of latitudes, longitudes, and heights is
C          retained.  Note that only the non-dipole component of the gridded
C          fields is interpolated in APXALL, APXMALL, and APXQ2G, while the
C          dipole component is handled analytically in order to improve
C          accuracy.  Since the non-dipole component becomes spatially
C          smoother and relatively less important than the dipole component
C          with increasing altitude, the grid spacing can become larger with
C          increasing altitude.  GGRID takes advantage of this fact by
C          making the vertical spacing of the grid points have equal
C          increments of (Re/(Re+height)).  The latitudinal grid spacing is
C          chosen so that at Earth's surface it is approximately the same as
C          the vertical distance between the two lowest grid points in
C          altitude.  The longitude spacing, in degrees, is chosen to be 6/5
C          of the latitude spacing, i.e., roughly equivalent in distance at
C          low and middle latitudes, though with decreasing longitude
C          spacing (in km) as the poles are approached.
C
C          INPUTS:
C            NVERT   = Spatial resolution parameter (see above for an
C                      explanation).  Increasing NVERT increases the model
C                      resolution.  NVERT must be at least 2, but in order to
C                      obtain useful results a value of at least 10 should be
C                      used; 30 or more is desirable unless computer
C                      limitations make the computation and storage of arrays
C                      of this size inconvenient (for NVERT = 30 and a global
C                      grid 0 to at least 1000 km, it took 63 minutes on a Sun
C                      4m to produce arrays for seven epochs and the resulting
C                      file is 16 MBytes).  However, one can increase the grid
C                      resolution (NVERT > 100 on a 32-bit single precision
C                      cmoputer) to the point of degraded accuracy close to
C                      the poles of quantities involving east-west gradients
C                      (B(1), G, H, W, Bhat(1), D1(1), D2(1),D3, E1, E2, E3,
C                      F1(2), and F2(2)).  Accuracy of all quantities involving
C                      gradients can be degraded at all latitudes for very
C                      large values of NVERT ( > 1000 on a 32-bit machine).
C            GLAMN   = minimum latitude for which apex coordinates will be
C                      needed (between -90. and +90., inclusive)
C            GLAMX   = maximum latitude for which apex coordinates will be
C                      needed (greater than or equal to GLAMN)
C            GLOMN   = minimum longitude for which apex coordinates will be
C                      needed (between -180. and +180., inclusive)
C            GLOMX   = maximum longitude for which apex coordinates will be
C                      needed (between -180. and +540., greater than or equal
C                      to GLOMN).  GLOMX will be lowered to GLOMN+360 if
C                      it initially exceeds that value.
C            ALTMN   = minimum altitude, in km, for which apex coordinates
C                      will be needed (greater than or equal to 0)
C            ALTMX   = maximum altitude, in km, for which apex coordinates
C                      will be needed (greater than or equal to ALTMN)
C          RETURNS:
C            GPLAT,GPLON,GPALT = Grid point latitudes, longitudes, and
C                     altitudes.  Latitudes are degrees North, longitudes
C                     degrees East, and altitudes are km msl.  These arrays
C                     must be adequately dimensioned.
C            NLAT,NLON,NALT = Number of assigned values in respective arrays
C                     GPLAT,GPLON,GPALT
C          Formal argument declarations
      DIMENSION GPLAT(*), GPLON(*), GPALT(*)

C          Local declarations:
      DOUBLE PRECISION DLON , DLAT , DIHT , DNV , DRE
      PARAMETER (RE = 6371.2 ,  DRE = 6371.2D0)
C     PARAMETER (RE = 6371.2)
C                RE = Mean Earth radius

C          Code from subroutine SETLIM in the previous version of apxntrp
C          which establishes indices defining the geographic grid (plus
C          a few modifications, such as some double precision variables).
      IF (GLAMX     .LT. GLAMN) GO TO 9100
      IF (GLOMX     .LT. GLOMN) GO TO 9200
      IF (ALTMX     .LT. ALTMN) GO TO 9300
      IF (ABS(GLAMN).GT.   90.) GO TO 9400
      IF (ABS(GLOMN).GT.  180.) GO TO 9500
      IF (ALTMN     .LT.    0.) GO TO 9600

      DNV  = DBLE (NVERT)
      DLON = 360.D0 / (5.D0*DNV)
      DLAT = 180.D0 / (3.D0*DNV)
      DIHT = 1.D0 / DNV

      LAMN =  MAX0 (INT((DBLE(GLAMN)+90.D0 )/DLAT)     , 0)
      LAMX =  MIN0 (INT((DBLE(GLAMX)+90.D0 )/DLAT+1.D0), 3*NVERT)

      LOMN  = MAX0 (INT((DBLE(GLOMN)+180.D0)/DLON)  , 0)
      GLOMX = AMIN1 (GLOMX,GLOMN+360.)
      LOMX  = MIN0 (INT((DBLE(GLOMX)+180.D0)/DLON+1.D0),10*NVERT)

      X = RE/(RE+ALTMX)/DIHT - 1.E-5
      IHTMN = AMAX1(X,1.)
      IHTMN = MIN0(IHTMN,NVERT-1)
      X = RE/(RE+ALTMN)/DIHT + 1.E-5
      I = X + 1.
      IHTMX = MIN0(I, NVERT)

      NLAT = LAMX - LAMN + 1
      NLON = LOMX - LOMN + 1
      NLON = MIN0 (NLON,5*NVERT+1)
      NALT = IHTMX - IHTMN + 1

C          Code from old MAKEXYZV which converts from indices to lat,lon,alt
      DO 110 I=1,NLAT
  110 GPLAT(I) = DLAT*DBLE(LAMN+I-1) - 90.D0
      DO 120 J=1,NLON
  120 GPLON(J) = DLON*DBLE(LOMN+J-1) - 180.D0
      DO 130 K=1,NALT
      IHT = IHTMX - K + 1
  130 GPALT(K) = DRE*(DBLE(NVERT-IHT) - 1.D-05) / (DBLE(IHT)+1.D-5)

C        Adjustments required to match test case results
      IF (GPLON(NLON-1) .GE. GLOMX) NLON = NLON - 1
      GPALT(1) = AMAX1 (GPALT(1),0.)

      RETURN

 9100 WRITE(6,'(''GGRID:  GLAMX < GLAMN'')')
      STOP
 9200 WRITE(6,'(''GGRID:  GLOMX < GLOMN'')')
      STOP
 9300 WRITE(6,'(''GGRID:  HTMAX < HTMIN'')')
      STOP
 9400 WRITE(6,'(''GGRID:  |GLATMIN| > 90.'')')
      STOP
 9500 WRITE(6,'(''GGRID:  |GLONMIN| > 180.'')')
      STOP
 9600 WRITE(6,'(''GGRID:  HTMIN < 0.'')')
      STOP
      END
