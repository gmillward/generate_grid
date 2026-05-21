      SUBROUTINE APXMKA (MSGUN, EPOCH, GPLAT,GPLON,GPALT,NLAT,NLON,NALT,
     +                  WK,LWK, IST)

C          Lunux specific settings: see IRLF

C          This contains ENTRYs that determine quantities related to Apex
C          coordinates.  They are designed to speed up coordinate determination
C          by substituting interpolation from previously calculated tables
C          for direct calculation.  They also permit calculating quantities
C	   involving gradients of Apex coordinates, such as base vectors in
C	   the Modified-Apex and Quasi-Dipole systems.  Options allow table 
C	   creation, storage, and read-back in addition to interpolation.
C
C          Tables (arrays) must be prepared prior to coordinate determination.
C          Tables may be created for one date and held only in memory by calling
C
C            APXMKA - make magnetic arrays.
C
C          Or, they may be created and written, then read back later by
C
C            APXWRA - make and write arrays
C            APXRDA - read stored arrays.
C
C          Tables may be created for multiple times (e.g., the DGRF dates)
C          and written to a file using APXWRA.  However, tables are held in
C          memory for only one time due to the potential to exceed memory
C          capacity.  Thus, when reading back multiples times, APXRDA
C          interpolates/extrapolates to a single time.
C
C          Alternatively, If APXWRA is called to make and write tables for one
C          time only, then it is not necessary to call APXRDA because tables
C          have already been initialized for that time.
C
C          Creation of global lookup tables can be time consuming, thus, it
C          is expected that this is done rarely; it might be considered an
C          installation step.  Once the tables are established, programs may
C          reference them by calling APXRDA.  In this case, the grid coordinates
C          (originally specified during the table creation) may be ascertained
C          by calling
C
C            APXGGC - get grid coordinates.
C
C          After the tables have been loaded in memory for a given date, it
C          is possible to interpolate in space to determinine various magnetic
C          parameters using the following calls:
C
C            APXALL  - get Apex coordinates (radius, lat, lon).
C            APXMALL - get Modified Apex coordinates, Quasi-Dipole
C                      coordinates and associated base vectors.
C            APXQ2G  - convert quasi-dipole to geodetic coordinates.
C
C          Details for each ENTRY introduced above are provided below.
C          Following the ENTRY summaries are comments describing EXTERNALS,
C          INSTALLATION SPECIFICS for different computers, and a HISTORY.
C
C          REFERENCE:
C          Richmond, A. D., Ionospheric Electrodynamics Using Magnetic Apex
C          Coordinates, J. Geomag. Geoelectr., 47, 191-212, 1995.
C
C          ALGORITHM:
C          When arrays are created, APXMKA calls subroutine MAKEXYZV, which
C          in turn calls subroutine APEX at each grid point to get the apex
C          radius A, the apex longitude PHIA, and the magnetic potential
C          VMP.  The cosine (CLP) and sine (SLP) of the quasi-dipole
C          latitude are computed from A.  From these can be computed
C          preliminary quantites defined as
C
C            x = cos(quasi-dipole latitude)*cos(apex longitude)
C            y = cos(quasi-dipole latitude)*sin(apex longitude)
C            z = sin(quasi-dipole latitude)
C            v = (VMP/VP)*((RE+ALT)/RE)**2
C
C          where VP is the magnitude of the magnetic potential of the
C          geomagnetic dipole at a radius of RE; ALT is altitude; and RE is
C          the mean Earth radius.  Note that all of these quantities vary
C          smoothly across the apex poles, unlike apex latitude or
C          quasi-dipole latitude, so that they can be linearly interpolated
C          near these poles.  Corresponding values of x,y,z,v for a dipole
C          field on a spherical Earth are computed analytically and
C          subtracted from the above quantities, and the results are put
C          into the 3D arrays X,Y,Z,V.  When APXALL or APXMALL is called,
C          trilinear interpolations (in latitude, longitude, and inverse
C          altitude) are carried out between the grid points.  Gradients are
C          calculated for APXMALL in order to determine the base vectors.
C          Analytic formulas appropriate for a dipole field on a spherical
C          Earth are used to determine the previously removed dipole
C          components of x,y,z,v, and their gradients and these are added
C          back to the interpolated values obtained for the non-dipole
C          components.  Finally, the apex-based coordinates and their base
C          vectors are calculated from x,y,z,v and their gradients.
C
C------------------------------------------------------------------------------
C
C       ENTRY APXMKA and APXWRA:
C          These create gridded arrays of quantities that can later be
C          linearly interpolated to any desired geographic point within
C          their range by APXALL, APXMALL, and APXQ2G.  The spatial extent
C          and resolution of the arrays in latitude, longitude, and altitude
C          can be tailored for the application at hand.  In one extreme,
C          very high interpolation accuracy can be achieved by creating a
C          2x2x2 array with very small spatial increments surrounding the
C          point of interest. (However, since vector quantities are obtained
C          by taking differences of the quantities at adjacent points,
C          computer roundoff errors will limit the benefits of making the
C          increments too small.) In another extreme, a global array from
C          the ground to an altitude of several Earth radii can be created.
C          (In this case, the spatial resolution is limited by the need to
C          maintain reasonable array sizes.) For most purposes, we have
C          found it adequate to use a global grid with dimensions 91,151,6
C          (latitude,longitude,altitude) and altitude from ground to 1274 km
C          as created by GGRID with NVERT = 30 (see file test.f).
C
C          WARNING: execution time to create these look-up tables can be
C          substantial.  It is dependent on the grid dimensions and the
C          number of Epochs to be computed.  For eight epochs and dimensions
C          (91,151,6), it took 76 minutes on a Sun SPARCstation IPX and the
C          resulting file is 15.8 Mbytes.
C
C          Use APXMKA to create the interpolation tables for a single time
C          without writing them.  Otherwise, use APXWRA create the tables
C          and save them in a file.
C
C              CALL APXMKA (MSGUN, EPOCH, GPLAT,GPLON,GPALT,NLAT,NLON,NALT,
C             +            WK,LWK, IST)
C
C              CALL APXWRA (MSGUN, FILNAM,IUN, EPOCH,NEPOCH,
C             +            GPLAT,GPLON,GPALT,NLAT,NLON,NALT, WK,LWK, IST)
C
C          INPUTS:
C            MSGUN  = Fortran unit number to write diagnostic messages.
C            EPOCH  = Time formatted as UT year and fraction; e.g., 1990.0
C                     is 0 UT 1 Jan 1990.  For APXMKA EPOCH is single valued;
C                     for APXWRA EPOCH contains NEPOCH values.
C            GPLAT,GPLON,GPALT = Grid point latitudes, longitudes, and
C                     altitudes.  NLAT values in GPLAT must be in the range
C                     -90. to +90. degrees North.  NLON values in GPLON
C                     must be in the range -270. to 270. degrees East.  NALT
C                     values in GPALT are km MSL.  Each array must be in
C                     ascending numerical order but they may have arbitrary
C                     spacing (as desired for model grid resolution).
C                     Subroutine GGRID can be used to set up these arrays.
C            NLAT,NLON,NALT = Number of latitudes, longitudes and altitudes
C                     are the single dimensions of the arrays GPLAT,GPLON,
C                     and GPALT.  Subroutine GGRID can be used to select
C                     appropriate values.  In general, increasing the number
C                     of grid points, increases the model resolution.  While
C                     array values at the grid points are exact, values at
C                     locations between grid points are estimated by linear
C                     interpolation.  One can increase the grid resolution to
C                     the point of degraded accuracy very close to the poles
C                     of quantities involving east-west gradients, viz.,
C                     B(1), G, H, W, Bhat(1), D1(1), D2(1), D3, E1, E2, E3,
C                     F1(2), and F2(2).  This is evident with dimensions
C                     301x501x15 for a global grid, 0-1000 km.
C            WK    =  Work array is used internally, i.e., there is no need
C                     to access the contents.  WK should not be altered
C                     between initialization (by APXMKA, APXWRA or APXRDA)
C                     and use (by APXGGC, APXALL, APXMALL, or APXQ2G).
C            LWK   =  Dimension of WK >=  NLAT*NLON*NALT*5 + NLAT+NLON+NALT
C
C          Additional INPUTS for APXWRA only:
C            FILNAM = file name where arrays are stored.
C            IUN    = Fortran unit number to be associated with FILNAM.
C            NEPOCH = Number of times in EPOCH.
C
C          RETURNS:
C            IST = Return status: = 0  okay
C                                 > 0  failed
C
C          Declarations for APXMKA formal arguments:
C            DIMENSION GPLAT(NLAT), GPLON(NLON), GPALT(NALT), WK(LWK)
C
C          Additional declarations for APXWRA formal arguments:
C            DIMENSION EPOCH(NEPOCH)
C            CHARACTER*(*) FILNAM
C
C------------------------------------------------------------------------------
C
C       ENTRY APXRDA:
C          Read back tables (previously created by calling APXWRA) in
C          preparation for magnetic coordinate determination (using APXALL,
C          APXMALL, APXQ2G).  If multiple dates were stored, interpolate
C          in time.
C
C              CALL APXRDA (MSGUN, FILNAM,IUN, DATE, WK,LWK, IST)
C
C          INPUTS:
C            MSGUN  = Fortran unit number to write diagnostic messages.
C            FILNAM = file name where arrays are stored.
C            IUN    = Fortran unit number to be associated with FILNAM.
C            DATE   = Time formatted as UT year and fraction; e.g., 1990.0
C                     is 0 UT 1 Jan 1990.
C            WK,LWK = Same as APXMKA
C
C          RETURNS:
C            IST = Return status: = 0  okay
C                                 = -1 non-fatal date problem
C                                 > 0  failed
C
C          Declarations for APXRDA formal arguments:
C            CHARACTER*(*) FILNAM
C            DIMENSION WK(LWK)
C
C------------------------------------------------------------------------------
C
C       ENTRY APXGGC:
C          Obtain grid coordinates from tables read back using APXRDA.
C
C              CALL APXGGC (MSGUN, WK,LWK, GPLAT,GPLON,GPALT,NLAT,NLON,NALT,
C             +            IST)
C
C          INPUTS:
C            MSGUN  = Fortran unit number to write diagnostic messages.
C            WK,LWK = Same as APXMKA
C
C          RETURNS:
C            GPLAT,GPLON,GPALT = Grid point latitudes, longitudes, and
C                                altitudes.  Units are degrees and kilometers.
C            NLAT,NLON,NALT    = Number of latitudes, longitudes and altitudes.
C            IST               = Return status, where IST=0 implies okay.
C
C          Declarations for APXGGC formal arguments:
C            DIMENSION GPLAT(NLAT), GPLON(NLON), GPALT(NALT), WK(LWK)
C
C
C------------------------------------------------------------------------------
C
C
C       ENTRY APXALL:
C          Determine Apex coordinates by interpolation from precalculated
C          arrays.  First call APXMKA, APXWRA, or APXRDA to load look-up
C          tables in memory.
C
C              CALL APXALL (GLAT,GLON,ALT, WK, A,ALAT,ALON, IST)
C
C          INPUTS:
C            GLAT = Geographic (geodetic) latitude, degrees, must be within
C                   the grid domain (GPLAT(1) <= GLAT <= GPLAT(NLAT)).
C            GLON = Geographic (geodetic) longitude, degrees, must be within
C                   one revolution of the grid domain:
C                     GPLON(1) <= GLON-360.,GLON, or GLON+360. <= GPLON(NLON))
C            ALT  = Altitude, km
C            WK   = same as entry APXMKA
C          RETURNS:
C            A    = Apex radius, normalized by Req
C            ALAT = Apex latitude, degrees
C            ALON = Apex longitude, degrees
C            IST  = Return status:  okay (0); or failure (1).
C
C          Dimensions of non-scalar arguments to APXALL:
C            WK(LWK)
C
C
C------------------------------------------------------------------------------
C
C       ENTRY APXMALL:
C          Compute Modified Apex coordinates, quasi-dipole coordinates,
C          base vectors and other parameters by interpolation from
C          precalculated arrays.  First call APXMKA, APXWRA, or APXRDA
C          to load look-up tables in memory.
C
C              CALL APXMALL (GLAT,GLON,ALT,HR, WK,
C             +             B,BHAT,BMAG,SI,
C             +             ALON,
C             +             XLATM,VMP,W,D,BE3,SIM,D1,D2,D3,E1,E2,E3,
C             +             XLATQD,F,F1,F2 , IST)
C
C          Reference:  Richmond, A. D., Ionospheric Electrodynamics Using
C          Magnetic Apex Coordinates, J. Geomag. Geoelectr., 47, 191-212, 1995.
C
C          INPUTS:
C            GLAT = Geographic (geodetic) latitude, degrees, must be within
C                   the grid domain (GPLAT(1) <= GLAT <= GPLAT(NLAT)).
C            GLON = Geographic (geodetic) longitude, degrees, must be within
C                   one revolution of the grid domain:
C                     GPLON(1) <= GLON-360.,GLON, or GLON+360. <= GPLON(NLON))
C            ALT  = Altitude, km
C            HR   = Reference altitude, km (used only for Modified Apex coords)
C            WK   = same as entry APXMKA
C          RETURNS:
C            B      = magnetic field components (east, north, up), in nT
C            BHAT   = components (east, north, up) of unit vector along
C                     geomagnetic field direction
C            BMAG   = magnitude of magnetic field, in nT
C            SI     = sin(I)
C            ALON   = apex longitude = modified apex longitude = quasi-dipole
C                     longitude, degrees
C            XLATM  = Modified Apex latitude, degrees
C            VMP    = magnetic potential, in T.m
C            W      = W of reference above, in km**2 /nT  (i.e., 10**15 m**2 /T)
C            D      = D of reference above
C            BE3    = B_e3 of reference above (= Bmag/D), in nT
C            SIM    = sin(I_m) described in reference above
C            D1,D2,D3,E1,E2,E3 = components (east, north, up) of base vectors
C                     described in reference above
C            XLATQD = quasi-dipole latitude, degrees
C            F      = F described in reference above for quasi-dipole
C                     coordinates
C            F1,F2  = components (east, north) of base vectors described in
C                     reference above for quasi-dipole coordinates
C            IST    = Return status:  okay (0); or failure (1).
C
C          Dimensions of non-scalar arguments to APXMALL:
C            GPLAT(NLAT),GPLON(NLON),GPALT(NALT),WK(LWK),
C            B(3),BHAT(3),D1(3),D2(3),D3(3), E1(3),E2(3),E3(3), F1(2),F2(2)
C
C
C------------------------------------------------------------------------------
C
C       ENTRY APXQ2G:
C          Convert from quasi-dipole to geodetic coordinates, APXQ2G
C          (input magnetic, output geodetic) is the functional inverse
C          of APXALL or APXMALL (input geodetic, output magnetic).  First
C          call APXMKA, APXWRA, or APXRDA to load look-up tables in memory.
C
C              CALL APXQ2G (QDLAT,QDLON,ALT, WK, GDLAT,GDLON, IST)
C
C          INPUTS:
C            QDLAT = quasi-dipole latitude in degrees
C            QDLON = quasi-dipole longitude in degrees
C            ALT   = altitude in km
C            WK    = same as entry APXMKA
C          RETURNS:
C            GDLAT = geodetic latitude in degrees
C            GDLON = geodetic longitude in degrees
C            IST   = Return status: =  0 okay
C                                   = -1 non-fatal, results are not as
C                                        close as PRECISE
C                                   >  0 failure
C
C          Dimensions of non-scalar arguments to APXQ2G:
C            WK(LWK)
C
C
C------------------------------------------------------------------------------
C
C          EXTERNALS:
C            apex.f   - APEX, etc. source code
C            magfld.f - COFRM, DYPOL, FELDG are the DGRF/IGRF.
C            magloctm.f, subsol.f and cossza.f - are not externals, but they
C                       are related, computing magnetic local time, the sub-
C                       solar point, and the cosine of the solar zenith angle.
C            test.f  -  also not an external, rather it is an example driver
C                       program demonstrating various call sequences.  It
C                       includes subroutine GGRID.
C
C          INSTALLATION SPECIFICS:
C            The only (known) machine dependency is parameter IRLF described
C            below.
C            Compiler default is usually to initialize memory to 0.  If this
C            is not true, an error trap involving KGMA may not work.
C
C          HISTORY:
C            Originally coded 940824 by A. D. Richmond, NCAR.  Components of
C            routines in the package came from previous work by Vincent
C            Wickwar (COFRM, FELDG in magfld.f).
C
C            Modified Sep 95 (Roy Barnes):  Changes were made with the
C            objective to allow the user to completely control definition of
C            the interpolation grid.  While doing this the order of the
C            ENTRYs in the first subroutine and arguments lists were
C            changed:  APXMKA, APXWRA, APXRDA (formerly GETARR) and the other
C            ENTRYs now include a work array (WK) which holds arrays X,Y,Z,V,
C            GPLAT,GPLON and GPALT.  Subroutine SETLIM was removed and the
C            grid creation algorithm based on NVERT originally integral to
C            GETARR and SETLIM has been extracted, but still available as
C            an example (see GGRID in file test.f).  Subroutine TSTDIM has a
C            different role, so it is now CKGP (check grid points).
C            MAKEXYZV was also changed to accomodate explicit grid point
C            arrays instead of computed values from an index, minimum and
C            delta.  Only one format is written now, so that is is possible
C            to concatinate files for different epochs.  This required
C            changing delta time logic from fixed 5 yr epochs to computed
C            differences which may vary.  Also, the ties to DGRF/IGRF dates
C            have been removed.
C
C            Modified Sep 96 (Roy Barnes):  Corrected bug in APXQ2G longitude
C            iteration.  The latitude iteration is now constrained to not
C            step beyond the (partial global) interpolation grid.  NITER was
C            increased from 9 to 14 and code to determine cos(lambda') (CLP)
C            was revised by Art Richmond to reduce truncation errors near
C            the geographic poles.
C
C            Modified Sep 97:  Corrected comments, definition of COLAT.
C
C            Modified Dec 98:  Change GLON input to try +/-360 values
C            before rejecting when out of defined grid (GDLON) range;
C            affects INTRP
C
C            Modified Feb-Mar 99:  (1) Corrected a typo (bad variable name)
C            in diagnostic print in INTRP: GLO -> GLON.  This error was
C            probably introduced Dec 98, so no-one had access to the
C            bad copy and, therefore, no announcement is needed.  (2) Also
C            modified APXMALL and APXQ2G:  When very close to a geographic
C            pole, gradients are recalculated after stepping away; in this
C            situation, the latitude input to INTRP was changed from GLAT
C            to GLATX.  This is affects gradients when within 0.1 degrees
C            of a pole (defined as the larger of GLATLIM, 89.9 deg, or the
C            second largest grid point).  (3) Changed MAKEXYZV to make X,Y,Z,
C            V constant for all longitudes (at each alt) when poleward of
C            POLA; POLA was added to /APXCON/.  (4)  Replaced definition of
C            DYLON in APXQ2G.  (5) Reduced NITER to 10 from 14.  Changes 3-5
C            fix a problem where APXQ2G calculations were failing to satisify
C            PRECISE within NITER iterations at the pole.  (6) Replace
C            XYZ2APX with revised trigonometry to avoid a truncation problem
C            at the magnetic equator.  Most changes were devised by Art
C            Richmond and installed by Roy B.
C
C            May 2000:  Relaxed acceptable input longitude check in CKGP
C            to be +/- 270 degrees.  Also updated IGRF coefficients in
C            magfld.f; now DGRF epochs are 1900, 1905, ... 1995 and IGRF
C            is for 2000-2005.
C
C            Questions about this version should be directed to Roy Barnes
C            NCAR (email: bozo@ucar.edu phone: 303/497-1230).
C
C------------------------------------------------------------------------------

      PARAMETER (XMISS=-32767. , GLATLIM=89.9 , PRECISE=7.6E-11,
     +          DATDMX=1. , DATIMX=2.5 , IRLF=4 )
C          XMISS   = value used to fill returned variables when requested
C                    point is outside array boundaries
C          GLATLIM = Extreme polar latitude allowed before changing east-west
C                    gradient calculation to avoid possible underflow at
C                    poles.  GLATLIM is superseded when the second to last
C                    grid point value is closer to a pole.
C          PRECISE = (angular distance)**2 (radians**2) of precision of
C                    transform from quasi-dipole to geodetic coordinates.
C                    7.6E-11 corresponds to an accuracy of 0.0005 degree.
C          IRLF    = Record length factor required to convert the computed
C                    length of direct access records from words to those
C                    units appropriate to this computer.  IRLF is 1 when
C                    RECL units are words, otherwise it is a function of
C                    the computer word length; e.g.,
C                      IRLF = 1 on a DEC  (32-bit words, RECL units are words)
C                      IRLF = 4 on a PC   (32-bit words, RECL units are bytes)
C                      IRLF = 4 on a Sun  (32-bit words, RECL units are bytes)
C                      IRLF = 8 on a Cray (64-bit words, RECL units are bytes)
C          DATDMX  = maximum time difference (years) allowed between the
C                    requested date and the single epoch in arrays.
C          DATIMX  = maximum time difference (years) allowed between the
C                    requested date and the closest epoch in the stored
C                    arrays (apropos multiple stored dates).

C          Formal argument declarations
      DIMENSION GPLAT(*),GPLON(*),GPALT(*), EPOCH(*), WK(*),
     +          GRADX(3), GRADY(3), GRADZ(3), GRADV(3),
     +          GRCLM(3), CLMGRP(3), RGRLP(3),
     +          B(3),BHAT(3), D1(3),D2(3),D3(3), E1(3),E2(3),E3(3),
     +          F1(2),F2(2)
      CHARACTER*(*) FILNAM

C          Local declarations
      CHARACTER CALNM*7, CMP*10, EDG*5
      SAVE KGMA, GLALMN,GLALMX, NLA,NLO,NAL, LBX,LBY,LBZ,LBV,LLA,LLO,LAL
      DATA KGMA /0/

C          Common APXDIPL is assigned in MAKEXYZV but computed in DYPOL
C            COLAT = geocentric colatitude (degrees) of north geomagnetic pole
C            ELON  = geocentric east longitude (degrees) of north geomagnetic
C                    pole
C            VP    = magnetic potential at 1 RE, south geomagnetic pole
C            CTP   = cos(colat*dtor)
C            STP   = sin(colat*dtor)
      COMMON /APXDIPL/  COLAT,ELON,VP,CTP,STP

C          Common APXCON is assigned here
C            RTOD, DTOR = radians to degrees (180/pi) and inverse
C            RE, REQ    = 6371.2, 6378.160 (Mean and equatorial Earth radius)
C            MSGU       = MSGUN to be passed to subroutines
C            POLA       = Pole angle (deg); when the geographic latitude is
C                         poleward of POLA, X,Y,Z,V are forced to be constant.
C                         for all longitudes at each altitude
      COMMON /APXCON/ RTOD,DTOR,RE,REQ,MSGU,POLA

      CALNM = 'APXMKA'
      LCN   = 6
      KGMA  = 1
      MSGU  = MSGUN
      NEPOK = 1
      IF (NLAT .LT. 2 .OR. NLON .LT. 2 .OR. NALT .LT. 2) GO TO 9100
      NLA = NLAT
      NLO = NLON
      NAL = NALT
      GO TO 40

      ENTRY APXGGC (MSGUN,WK,LWK, GPLAT,GPLON,GPALT,NLAR,NLOR,NALR,IST)
C          Sep 95 R. Barnes
      CALNM = 'APXGGC'
      LCN   = 6
      MSGU  = MSGUN
      IF (KGMA .LT. 1) GO TO 9300
      J = LLA
      DO 10 I=1,NLA
      GPLAT(I) = WK(J)
   10 J = J + 1
      DO 20 I=1,NLO
      GPLON(I) = WK(J)
   20 J = J + 1
      DO 30 I=1,NAL
      GPALT(I) = WK(J)
   30 J = J + 1
      NLAR = NLA
      NLOR = NLO
      NALR = NAL
      IST = 0
      RETURN

      ENTRY APXWRA (MSGUN, FILNAM,IUN, EPOCH,NEPOCH,
     +             GPLAT,GPLON,GPALT,NLAT,NLON,NALT, WK,LWK, IST)
C          Sep 95 R. Barnes
      CALNM = 'APXWRA'
      LCN   = 6
      KGMA  = 2
      MSGU  = MSGUN
      NEPOK = NEPOCH
      IF (NLAT .LT. 2 .OR. NLON .LT. 2 .OR. NALT .LT. 2) GO TO 9100
      NLA = NLAT
      NLO = NLON
      NAL = NALT
      GO TO 40

      ENTRY APXRDA (MSGUN, FILNAM,IUN, DATE, WK,LWK, IST)
C          Sep 95 R. Barnes
      CALNM = 'APXRDA'
      LCN   = 6
      KGMA  = 3

C          Open the read-back file with a temporary record length, get
C          the grid dimensions from the first values, then close it, (so
C          it can be reopened later with the proper LDR):
      LDR = 7*IRLF
      OPEN (IUN,FILE=FILNAM,ACCESS='direct',RECL=LDR,STATUS='old',
     +                                                       IOSTAT=IST)
      MSGU = MSGUN
      IF (IST .NE. 0) GO TO 9110
      READ  (IUN,REC=1,IOSTAT=IST) YEAR,COLAT,ELON,VP,NLA,NLO,NAL
      IF (IST .NE. 0) GO TO 9120
      CLOSE (IUN)

   40 RE   = 6371.2
      REQ  = 6378.160
      RTOD = 45./ATAN(1.)
      DTOR = 1./RTOD
      POLA = 90. - SQRT (PRECISE) * RTOD
      LFN = 0
      IF (KGMA .EQ. 1) GO TO 51
      DO 50 I=1,LEN(FILNAM)
      IF (FILNAM(I:I) .EQ. ' ') GO TO 51
   50 LFN = LFN + 1
   51 CONTINUE

C          Save grid dimensions, establish direct access rec length, and
C          determine indices into the work array.  WK contains arrays
C          X,Y,Z,V,temp,GPLAT,GPLON,GPALT where X thru tmp are dimensioned
C          (NLAT,NLON,NALT); tmp is scratch space used during read back.
      NGP = NLA*NLO*NAL
      NGM1= NGP - 1
      LDR = NGP * IRLF
      LBX = 1
      LBY = LBX + NGP
      LBZ = LBY + NGP
      LBV = LBZ + NGP
      LBT = LBV + NGP
      LLA = LBT + NGP
      LLO = LLA + NLA
      LAL = LLO + NLO
      LEG = LAL + NAL-1
      IF (LWK .LT. LEG) GO TO 9130

      IF (KGMA .EQ. 3) GO TO 200

C          Make and optionally write interpolation arrays for NEPOK times
      IF (KGMA .EQ. 2) THEN
	OPEN (IUN,FILE=FILNAM,ACCESS='direct',RECL=LDR,STATUS='new',
     +       IOSTAT=IST)
	IF (IST .NE. 0) GO TO 9115
      ENDIF

      CALL CKGP (CALNM(:LCN),MSGUN,NLAT,NLON,NALT,GPLAT,GPLON,GPALT,IST)
      IF (IST .NE. 0) RETURN
      I = LLA - 1
      DO 60 J=1,NLAT
      I = I + 1
   60 WK(I) = GPLAT(J)
      DO 70 J=1,NLON
      I = I + 1
   70 WK(I) = GPLON(J)
      DO 80 J=1,NALT
      I = I + 1
   80 WK(I) = GPALT(J)

      IF (NEPOK .LT. 1) GO TO 9140
      J = 1
      DO 100 I=1,NEPOK
      CALL MAKEXYZV (EPOCH(I),NLAT,NLON,NALT,GPLAT,GPLON,GPALT,
     +               WK(LBX),WK(LBY),WK(LBZ),WK(LBV))
      IF (KGMA .EQ. 1) GO TO 100
      WRITE (IUN,REC=J) EPOCH(I),COLAT,ELON,VP,NLAT,NLON,NALT
      WRITE (IUN,REC=J+1) (WK(K),K=LLA,LEG)
      WRITE (IUN,REC=J+2) (WK(K),K=LBX,LBX+NGM1)
      WRITE (IUN,REC=J+3) (WK(K),K=LBY,LBY+NGM1)
      WRITE (IUN,REC=J+4) (WK(K),K=LBZ,LBZ+NGM1)
      WRITE (IUN,REC=J+5) (WK(K),K=LBV,LBV+NGM1)
  100 J = J + 6
      IF (KGMA .EQ. 2) CLOSE (IUN)
      IST = 0
      GO TO 300

C          Read back interpolation arrays.  When arrays for multiple times
C          are available, interpolate using the pair bounding the desired
C          time (DATE).  Make an initial pass only to identify closest
C          available times and issue any diagnostics, then the second pass
C          to read the stored arrays (GPLAT,GPLON,GPALT,X,Y,Z,V) and do
C          the linear interpolation/extrapolation.
  200 OPEN (IUN,FILE=FILNAM,ACCESS='direct',RECL=LDR,STATUS='old',
     +      IOSTAT=IST)
      IF (IST .NE. 0) GO TO 9110

      READ (IUN,REC=1,IOSTAT=IST) TB
      IF (IST .NE. 0) GO TO 9120
      I2 = 1
      TL = TB
      IL = 1
      I = 1
  210 I = I + 6
      READ (IUN,REC=I,IOSTAT=JST) T
C         JST .NE. 0 is assumed to mean read beyond last record
      IF (JST .NE. 0) GO TO 220
      TO = TL
      TL = T
      IL = I
      IF (DATE .GT. TL) GO TO 210

  220 I1 = IL - 6

      IST = 0
      IF (TL .EQ. TB) THEN
	DATD = ABS (DATE-TB)
	IF (DATD .GT. DATDMX) THEN
	  WRITE (MSGU,9150) CALNM(1:LCN),DATE,DATD,TB
	  IF (TB .EQ. 0.) WRITE (MSGU,9155) FILNAM(1:LFN)
	  IST = -1
	ENDIF
	I1 = 1
	I2 = 0
      ELSE IF (DATE .LT. TB) THEN
	WRITE (MSGU,9160) CALNM(1:LCN),DATE,TB,FILNAM(1:LFN)
	IST = -1
      ELSE IF (DATE .GT. TL) THEN
	WRITE (MSGU,9170) CALNM(1:LCN),DATE,TL,FILNAM(1:LFN)
	IST = -1
      ELSE
	DATD = AMIN1 (DATE-TO,TL-DATE)
	IF (DATD .GT. DATIMX) THEN
	  WRITE (MSGU,9180) CALNM(1:LCN),DATE,TB,TL,FILNAM(1:LFN),DATD
	  IST = -1
	ENDIF
      ENDIF

      READ (IUN,REC=I1) YEAR1,COLAT,ELON,VP,NLAI,NLOI,NALI
      TI = YEAR1
      IF (NLAI.NE.NLA .OR. NLOI.NE.NLO .OR. NALI.NE.NAL) GO TO 9190
      READ (IUN,REC=I1+1) (WK(I),I=LLA,LEG)
      READ (IUN,REC=I1+2) (WK(I),I=LBX,LBX+NGM1)
      READ (IUN,REC=I1+3) (WK(I),I=LBY,LBY+NGM1)
      READ (IUN,REC=I1+4) (WK(I),I=LBZ,LBZ+NGM1)
      READ (IUN,REC=I1+5) (WK(I),I=LBV,LBV+NGM1)
      IF (I2 .EQ. 1) THEN
	READ (IUN,REC=I1+6) YEAR2,COLA2,ELON2,VP2,NLAI,NLOI,NALI
	TI = YEAR2
	IF (NLAI.NE.NLA .OR. NLOI.NE.NLO .OR. NALI.NE.NAL) GO TO 9190
	LE = LBT + NLA+NLO+NAL - 1
	READ (IUN,REC=I1+7) (WK(I),I=LBT,LE)
	J = LLA
	DO 230 I=LBT,LE
	IF (WK(J) .NE. WK(I)) GO TO 9200
  230   J = J + 1
	FRAC = (DATE-YEAR1) / (YEAR2-YEAR1)
	OMF  = 1. - FRAC
	LE = LBT + NGM1
	READ (IUN,REC=I1+8) (WK(I),I=LBT,LE)
	J = LBX
	DO 240 I=LBT,LE
	WK(J) =  OMF*WK(J) + FRAC*WK(I)
  240   J = J + 1
	READ (IUN,REC=I1+9) (WK(I),I=LBT,LE)
	DO 250 I=LBT,LE
	WK(J) =  OMF*WK(J) + FRAC*WK(I)
  250   J = J + 1
	READ (IUN,REC=I1+10) (WK(I),I=LBT,LE)
	DO 260 I=LBT,LE
	WK(J) =  OMF*WK(J) + FRAC*WK(I)
  260   J = J + 1
	READ (IUN,REC=I1+11) (WK(I),I=LBT,LE)
	DO 270 I=LBT,LE
	WK(J) =  OMF*WK(J) + FRAC*WK(I)
  270   J = J + 1
	COLAT = OMF*COLAT + FRAC*COLA2
	ELON  = OMF*ELON  + FRAC*ELON2
	VP    = OMF*VP    + FRAC*VP2
      ENDIF

      CTP  = COS (COLAT*DTOR)
      STP  = SIN (COLAT*DTOR)
      YEAR = DATE
      CLOSE (IUN)

C          Establish for this grid polar latitude limits beyond which east-west
C          gradients are computed differently to avoid potential underflow
  300 GLALMX = AMAX1 ( GLATLIM,WK(LLA+NLA-2))
      GLALMN = AMIN1 (-GLATLIM,WK(LLA+1))

      RETURN

C*******************************************************************************

      ENTRY APXMALL (GLAT,GLO ,ALT,HR, WK,                              !Inputs
     +              B,BHAT,BMAG,SI,                                     !Mag Fld
     +              ALON,                                               !Apx Lon
     +              XLATM,VMP,W,D,BE3,SIM,D1,D2,D3,E1,E2,E3,            !Mod Apx
     +              XLATQD,F,F1,F2 , IST)                               !Qsi-Dpl
C          940822 A. D. Richmond, Sep 95 R. Barnes

C          Test to see if WK has been initialized
      CALNM = 'APXMALL'
      LCN   = 7
      IF (KGMA .LT. 1) GO TO 9300

C          Alias geodetic longitude input in case INTRP adjusts it by +/-360
      GLON = GLO

      CALL INTRP (GLAT,GLON,ALT, WK(LBX),WK(LBY),WK(LBZ),WK(LBV),
     +           NLA,NLO,NAL,    WK(LLA),WK(LLO),WK(LAL),
     +           FX,FY,FZ,FV,
     +           DFXDTH,DFYDTH,DFZDTH,DFVDTH,DFXDLN,DFYDLN,DFZDLN,
     +           DFVDLN,DFXDH,DFYDH,DFZDH,DFVDH, CALNM(1:LCN),IST)

      IF (IST .NE. 0) THEN
	CALL SETMISS (XMISS, XLATM,ALON,VMP,B,BMAG,BE3,SIM,SI,F,D,W,
     +               BHAT,D1,D2,D3,E1,E2,E3,F1,F2)
	RETURN
      ENDIF

      CALL ADPL (GLAT,GLON,CTH,STH,FX,FY,FZ,FV,
     +          DFXDTH,DFYDTH,DFZDTH,DFVDTH,DFXDLN,DFYDLN,DFZDLN,DFVDLN)
      CALL GRADXYZV (ALT,CTH,STH,
     +          DFXDTH,DFYDTH,DFZDTH,DFVDTH,DFXDLN,DFYDLN,DFZDLN,DFVDLN,
     +          DFXDH,DFYDH,DFZDH,DFVDH,GRADX,GRADY,GRADZ,GRADV)

      IF (GLAT .GT. GLALMX .OR. GLAT .LT. GLALMN) THEN
C          If the point is very close to either the North or South
C          geographic pole, recompute the east-west gradients after
C          stepping a small distance from the pole.
	GLATX = GLALMX
	IF (GLAT .LT. 0.) GLATX = GLALMN
C990225 CALL INTRP (GLAT ,GLON,ALT, WK(LBX),WK(LBY),WK(LBZ),WK(LBV),
	CALL INTRP (GLATX,GLON,ALT, WK(LBX),WK(LBY),WK(LBZ),WK(LBV),
     +             NLA,NLO,NAL,    WK(LLA),WK(LLO),WK(LAL),
     +             FXDUM,FYDUM,FZDUM,FVDUM,
     +             DMXDTH,DMYDTH,DMZDTH,DMVDTH,DFXDLN,DFYDLN,DFZDLN,
     +             DFVDLN,DMXDH,DMYDH,DMZDH,DMVDH, CALNM(1:LCN),IST)
	CALL ADPL (GLATX,GLON,CTH,STH,FXDUM,FYDUM,FZDUM,FVDUM, DMXDTH,
     +            DMYDTH,DMZDTH,DMVDTH,DFXDLN,DFYDLN,DFZDLN,DFVDLN)
	CALL GRAPXYZV (ALT,CTH,STH, DFXDLN,
     +                DFYDLN,DFZDLN,DFVDLN,GRADX,GRADY,GRADZ,GRADV)
      ENDIF

      CALL GRADLPV (HR,ALT,FX,FY,FZ,FV,GRADX,GRADY,GRADZ,GRADV,
     +             XLATM,ALON,VMP,GRCLM,CLMGRP,XLATQD,RGRLP,B,CLM,R3_2)
      CALL BASVEC (HR,XLATM,GRCLM,CLMGRP,RGRLP,B,CLM,R3_2,
     +            BMAG,SIM,SI,F,D,W,BHAT,D1,D2,D3,E1,E2,E3,F1,F2)
      BE3 = BMAG/D

      IST = 0
      RETURN


C*******************************************************************************

      ENTRY APXALL (GLAT,GLO ,ALT, WK, A,ALAT,ALON, IST)
C          940802 A. D. Richmond, Sep 95 R. Barnes

C          Test to see if WK has been initialized
      CALNM = 'APXALL'
      LCN   = 6
      IF (KGMA .LT. 1) GO TO 9300

C          Alias geodetic longitude input in case INTRPSC adjusts it by +/-360
      GLON = GLO

      CALL INTRPSC (GLAT,GLON,ALT, WK(LBX),WK(LBY),WK(LBZ),
     +              NLA,NLO,NAL,   WK(LLA),WK(LLO),WK(LAL),
     +              FX,FY,FZ, CALNM(1:LCN), IST)
      IF (IST .NE. 0) GO TO 600

      CALL ADPLSC (GLAT,GLON,FX,FY,FZ)

      CALL XYZ2APX (ALT,FX,FY,FZ,A,ALAT,ALON,IST)
      IF (IST .EQ. 0) GO TO 601

  600 A    = XMISS
      ALAT = XMISS
      ALON = XMISS
  601 CONTINUE

      RETURN

C*******************************************************************************

      ENTRY APXQ2G (QDLAT,QDLON,ALT, WK, GDLAT,GDLON, IST)
C          940819 A. D. Richmond, Sep 95 R. Barnes, Sep 96 mod A. D. Richmond
C          Input guessed geodetic coordinates (YLAT,YLON) to INTRP and
C          compare the returned magnetic coordinates to those desired.
C          If the guess is not sufficiently close (PRECISE), make another
C          guess by moving in the direction of the gradient of a quantity
C          (DIST2) that approximates the squared angular distance between
C          the returned and desired magnetic coordinates.

C          Test to see if WK has been initialized
      CALNM = 'APXQ2G'
      LCN   = 6
      IF (KGMA .LT. 1) GO TO 9300

C          Determine quasi-cartesian coordinates on a unit sphere of the
C          desired magnetic lat,lon in quasi-dipole coordinates.
      X0 = COS (QDLAT*DTOR) * COS (QDLON*DTOR)
      Y0 = COS (QDLAT*DTOR) * SIN (QDLON*DTOR)
      Z0 = SIN (QDLAT*DTOR)

C          Initial guess:  use centered dipole, convert to geocentric coords
      CALL GM2GC (QDLAT,QDLON,YLAT,YLON)

C          Iterate until (angular distance)**2 (units: radians) is within
C          PRECISE of location (QDLAT,QDLON) on a unit sphere.
      NITER = 10
      DO 1000 ITER=1,NITER
      CALL INTRP (YLAT,YLON,ALT, WK(LBX),WK(LBY),WK(LBZ),WK(LBV),
     +           NLA,NLO,NAL,    WK(LLA),WK(LLO),WK(LAL),
     +           FX,FY,FZ,FV,
     +           DFXDTH,DFYDTH,DFZDTH,DFVDTH,DFXDLN,DFYDLN,DFZDLN,
     +           DFVDLN,DFXDH,DFYDH,DFZDH,DFVDH, CALNM(1:LCN),IST)
      IF (IST .NE. 0) GO TO 9400
      CALL ADPL (YLAT,YLON,CTH,STH,FX,FY,FZ,FV,
     +          DFXDTH,DFYDTH,DFZDTH,DFVDTH,DFXDLN,DFYDLN,DFZDLN,DFVDLN)

      DISTLON = COS(YLAT*DTOR)
      IF (YLAT .GT. GLALMX .OR. YLAT .LT. GLALMN) THEN
	GLATX = GLALMX
	IF (YLAT.LT.0.) GLATX = GLALMN
	DISTLON = COS (GLATX*DTOR)
C990225 CALL INTRP (YLAT ,YLON,ALT, WK(LBX),WK(LBY),WK(LBZ),WK(LBV),
	CALL INTRP (GLATX,YLON,ALT, WK(LBX),WK(LBY),WK(LBZ),WK(LBV),
     +             NLA,NLO,NAL,     WK(LLA),WK(LLO),WK(LAL),
     +             FXDUM,FYDUM,FZDUM,FVDUM,
     +             DMXDTH,DMYDTH,DMZDTH,DMVDTH,DFXDLN,DFYDLN,DFZDLN,
     +             DFVDLN,DMXDH,DMYDH,DMZDH,DMVDH, CALNM(1:LCN),IST)
	CALL ADPL (GLATX,YLON,CTH,STH,FXDUM,FYDUM,FZDUM,FVDUM,
     +          DMXDTH,DMYDTH,DMZDTH,DMVDTH,DFXDLN,DFYDLN,DFZDLN,DFVDLN)
      ENDIF

C          At this point, FX,FY,FZ are approximate quasi-cartesian
C          coordinates on a unit sphere for the quasi-dipole coordinates
C          corresponding to the geodetic coordinates YLAT, YLON.
C          Normalize the vector length of (FX,FY,FZ) to unity using XNORM
C          so that the resultant vector can be directly compared with the
C          target vector (X0,Y0,Z0).
      XNORM = SQRT(FX*FX + FY*FY + FZ*FZ)
      XDIF = FX/XNORM - X0
      YDIF = FY/XNORM - Y0
      ZDIF = FZ/XNORM - Z0
C          DIST2 = square of distance between normalized (FX,FY,FZ) and
C          X0,Y0,Z0.
      DIST2 = XDIF*XDIF + YDIF*YDIF + ZDIF*ZDIF

      IF (DIST2 .LE. PRECISE) GO TO 1001
C          HGRD2* = one-half of east or north gradient of DIST2 on unit sphere.
      HGRD2E =  (XDIF*DFXDLN + YDIF*DFYDLN + ZDIF*DFZDLN)/DISTLON
      HGRD2N = -(XDIF*DFXDTH + YDIF*DFYDTH + ZDIF*DFZDTH)
      HGRD2  = SQRT(HGRD2E*HGRD2E + HGRD2N*HGRD2N)
C          ANGDIST = magnitude of angular distance to be moved for new guess
C          of YLAT, YLON.
      ANGDIST = DIST2/HGRD2

C          Following spherical trigonometry moves YLAT, YLON to new location,
C          in direction of grad(DIST2), by amount ANGDIST.
      CAL = -HGRD2N/HGRD2
      SAL = -HGRD2E/HGRD2
      COSLM = COS(YLAT*DTOR)
      SLM = SIN(YLAT*DTOR)
      CAD = COS(ANGDIST)
      SAD = SIN(ANGDIST)
      SLP = SLM*CAD + COSLM*SAD*CAL

C          Old code (below) introduced truncation errors near geographic poles
C     SLP = AMIN1 (SLP,SGLAN) ; sglan = sin(wk(lla+nla-1))*dtor)
C     SLP = AMAX1 (SLP,SGLAS) ; sglas = sin(wk(lla))      *dtor)
C     CLP = SQRT(1. - SLP*SLP)
C     YLAT = ASIN(SLP)*RTOD

      CLM2 = COSLM*COSLM
      SLM2 = SLM*SLM
      SAD2 = SAD*SAD
      CAL2 = CAL*CAL
      CLP2 = CLM2 + SLM2*SAD2 - 2.*SLM*CAD*COSLM*SAD*CAL -CLM2*SAD2*CAL2
      CLP = SQRT (AMAX1(0.,CLP2))
      YLAT = ATAN2(SLP,CLP)*RTOD

C          Restrict latitude iterations to stay within the interpolation grid
C          limits, but let INTRP find any longitude exceedence.  This is only
C          an issue when the interpolation grid does not cover the entire
C          magnetic pole region.
      YLAT = AMIN1(YLAT,WK(LLA+NLA-1))
      YLAT = AMAX1(YLAT,WK(LLA))

      DYLON = ATAN2 (SAD*SAL,CAD*COSLM-SAD*SLM*CAL)*RTOD
C          Old formula (replaced Mar 99) had truncation problem near poles
C     DYLON = ATAN2 (CLP*SAD*SAL,CAD-SLM*SLP)*RTOD

      YLON  = YLON + DYLON
      IF (YLON .GT.  WK(LLO+NLO-1)) YLON = YLON - 360.
      IF (YLON .LT.  WK(LLO)      ) YLON = YLON + 360.
 1000 CONTINUE

      WRITE (MSGU,'(''APXQ2G: Warning'',I3,'' iterations only reduced th
     +e angular difference to'',/,8X,F8.5,'' degrees ('',F8.5,'' degrees
     + is the test criterion)'')') NITER, SQRT(DIST2  )*RTOD ,
     +                                    SQRT(PRECISE)*RTOD
      EDG = ' '
      IF (YLAT .EQ. WK(LLA+NLA-1)) EDG = 'north'
      IF (YLAT .EQ. WK(LLA))       EDG = 'south'
      IF (EDG .NE. ' ') WRITE (MSGU,'(''        Coordinates are on the '
     +',A,'' edge of the interpolation grid and'',/,''        latitude i
     +s constrained to stay within grid limits when iterating.'')') EDG

      IST = -1
      GO TO 1010

 1001 IST = 0

 1010 GDLAT = YLAT
      GDLON = YLON

      RETURN

C*******************************************************************************

C          Error Trap diagnostics
 9100 WRITE (MSGU,'(A,'':  NLAT,NLON or NALT < 2 '',3I8)')
     +             CALNM(1:LCN),  NLAT,NLON,NALT
      IST = 1
      RETURN
 9110 WRITE (MSGU,'(A,'': Trouble opening old file "'',A,''"'')')
     + CALNM(1:LCN), FILNAM(1:LFN)
      RETURN
 9115 WRITE (MSGU,'(A,'': Trouble opening new file "'',A,''"'')')
     + CALNM(1:LCN), FILNAM(1:LFN)
      RETURN
 9120 WRITE (MSGU,'(A,'': Trouble reading first record of '',A)')
     + CALNM(1:LCN), FILNAM(1:LFN)
      RETURN
 9130 WRITE (MSGU,'(A,'': LWK is too small; LWK must be at least'',I5,''
     + but LWK ='',I5)')  CALNM(1:LCN), LEG, LWK
      IST = 1
      RETURN
 9140 WRITE (MSGU,'(A,'':  NEPOCH must be positive; NEPOCH ='',I8)')
     +       CALNM(1:LCN), NEPOK
      IST = 1
      RETURN

 9150 FORMAT (A,': DATE (',F7.2,') differs by',F5.2,' years from the sto
     +red EPOCH (',F7.2,')')
 9155 FORMAT ('        A stored date = 0. implies "',A,'" is incorrectly
     + formatted')
 9160 FORMAT (A,': DATE (',F7.2,') is extrapolated before first EPOCH ('
     +,F7.2,') in "',A,'"')
 9170 FORMAT (A,': DATE (',F7.2,') is extrapolated after last EPOCH (',F
     +7.2,') in "',A,'"')
 9180 FORMAT (A,': DATE (',F7.2,') minimum difference from the nearest s
     +tored ',/,'        EPOCHs (',F7.2,', ',F7.2,') in "',A,'"',/,'
     +   is',F6.2,' years')

 9190 WRITE (MSGU,'(A,'': Dimensions (lat,lon,alt) read from "'',A,''" f
     +or EPOCH '',F7.2,/,''        (''I5,'','',I5,'','',I5,'') do not ma
     +tch ('',I5,'','',I5,'','',I5,'') from the'',/,''        first EPOC
     +H read ('',F7.2,'')'')') CALNM(1:LCN), FILNAM(1:LFN), TI ,
     + NLAI,NLOI,NALI,  NLA,NLO,NAL, TB
      IST = 1
      RETURN
 9200 CMP = 'latitudes'
      LC  = 9
      I1 = LLA - 1
      IT = LBT - 1
      N  = NLA
      IF (I .LT. LLO) GO TO 9201
      CMP = 'longitudes'
      LC   = 10
      I1 = I1 + NLA
      IT = IT + NLA
      N  = NLO
      IF (I .LT. LAL) GO TO 9201
      CMP = 'altitudes'
      LC  = 9
      I1 = I1 + NLO
      IT = IT + NLO
      N  = NAL
 9201 WRITE (MSGU,'(A,'': Grid '',A,'' read from "'',A,''" for EPOCH'',F
     +8.2,'' do not match the'',/,''        first EPOCH ('',F7.2,'')'',/
     +,''        First    Current'',/,(4X,2F10.3))')  CALNM(1:LCN),
     +  CMP(1:LC), FILNAM(1:LFN), TI,TB, (WK(I1+I),WK(IT+I),I=1,N)
      IST = 1
      RETURN

 9300 WRITE(MSGU,'(A,'': Must first load lookup tables by calling APXMKA
     +, APXWRA or APXRDA'')') CALNM(1:LCN)
      IST = 1
      RETURN
 9400 WRITE (MSGU,'(''APXQ2G: INTRP failed (maybe coordinates are not wi
     +thin interpolation grid)'')')
      IST = 1
      RETURN

      END

C*******************************************************************************
      SUBROUTINE INTRP (GLAT,GLON,ALT, X,Y,Z,V, NLAT,NLON,NALT,
     +                 GPLAT,GPLON,GPALT, FX,FY,FZ,FV,
     +                 DFXDTH,DFYDTH,DFZDTH,DFVDTH,DFXDLN,DFYDLN,DFZDLN,
     +                 DFVDLN,DFXDH,DFYDH,DFZDH,DFVDH, CALNM,IST)
C          Interpolation of x,y,z,v and their derivatives
C          940806 A. D. Richmond
C          INPUTS:
C            GLAT    = latitude, degrees
C            GLON    = longitude, degrees
C            ALT     = altitude, km
C            X,Y,Z,V = gridded arrays
C            NLAT,NLON,NALT = three dimensions of x,y,z,v and respective
C                      dimensions of GP___ arrays
C            GPLAT,GPLON,GPALT = grid point geographic locations
C            CALNM   = Name of calling routine (for error diagnostics)
C          OUTPUT:
C            FX = interpolated value of x
C            FY = interpolated value of y
C            FZ = interpolated value of z
C            FV = interpolated value of v
C            DFXDTH,DFYDTH,DFZDTH,DFVDTH = derivatives of x,y,z,v with
C                  respect to colatitude, in radians-1
C            DFXDLN,DFYDLN,DFZDLN,DFVDLN = derivatives of x,y,z,v with
C                  respect to longitude, in radians-1
C            DFXDH,DFYDH,DFZDH,DFVDH = derivatives of x,y,z,v with
C                  respect to altitude, in km-1
C            IST = Status =  0 = okay.
      DIMENSION X(NLAT,NLON,NALT), Y(NLAT,NLON,NALT), Z(NLAT,NLON,NALT),
     +          V(NLAT,NLON,NALT) , GPLAT(NLAT),GPLON(NLON),GPALT(NALT)
      CHARACTER*(*) CALNM

      COMMON /APXCON/ RTOD,DTOR,RE,REQ,MSGU,POLA

C          Local declarations
      DATA IO, JO, KO / 1, 1, 1 /
      SAVE IO, JO, KO

      IENTRY = 0
      GO TO 5

C*******************************************************************************
      ENTRY INTRPSC (GLAT,GLON,ALT, X,Y,Z, NLAT,NLON,NALT,
     +               GPLAT,GPLON,GPALT, FX,FY,FZ , CALNM,IST)
C          Interpolation of x,y,z
C          940803 A. D. Richmond.
C          Inputs and outputs:  for definitions, see above.

      IENTRY = 1
    5 IST = 0

      IF (GLAT .LT. GPLAT(1) .OR. GLAT .GT. GPLAT(NLAT)) GO TO 9100
      IF (ALT  .LT. GPALT(1) .OR. ALT  .GT. GPALT(NALT)) GO TO 9300
C          Accept input longitude range +/- one revolution (Dec 98)
      IF (GLON .LT. GPLON(1)   ) GLON = GLON + 360.
      IF (GLON .GT. GPLON(NLON)) GLON = GLON - 360.
      IF (GLON .LT. GPLON(1) .OR. GLON .GT. GPLON(NLON)) GO TO 9200

      I = IO
      IF (GLAT .LE. GPLAT(I)) GO TO 15
   12 I = I + 1
      IF (GPLAT(I) .LT. GLAT) GO TO 12
      I = I - 1
      GO TO 16
   14 I = I - 1
   15 IF (GPLAT(I) .GT. GLAT) GO TO 14
   16 IO = I
      DLAT = GPLAT(I+1) - GPLAT(I)
      XI   = (GLAT - GPLAT(I)) / DLAT

      J = JO
      IF (GLON .LE. GPLON(J)) GO TO 25
   22 J = J + 1
      IF (GPLON(J) .LT. GLON) GO TO 22
      J = J - 1
      GO TO 26
   24 J = J - 1
   25 IF (GPLON(J) .GT. GLON) GO TO 24
   26 JO = J
      DLON = GPLON(J+1) - GPLON(J)
      YJ   = (GLON - GPLON(J)) / DLON

      K = KO
      IF (ALT .LE. GPALT(K)) GO TO 35
   32 K = K + 1
      IF (GPALT(K) .LT. ALT) GO TO 32
      K = K - 1
      GO TO 36
   34 K = K - 1
   35 IF (GPALT(K) .GT. ALT) GO TO 34
   36 KO = K
      HTI  = RE/(RE+ALT)
      DIHT = RE/(RE+GPALT(K+1)) - RE/(RE+GPALT(K))
      ZK  = (HTI - RE/(RE+GPALT(K))) / DIHT

C          For intrpsc:
      IF (IENTRY.EQ.1) THEN
	CALL TRILINS (X(I,J,K),NLAT,NLON,XI,YJ,ZK,FX)
	CALL TRILINS (Y(I,J,K),NLAT,NLON,XI,YJ,ZK,FY)
	CALL TRILINS (Z(I,J,K),NLAT,NLON,XI,YJ,ZK,FZ)
        RETURN
      ENDIF
	
C          For intrp:
      CALL TRILIN (X(I,J,K),NLAT,NLON,XI,YJ,ZK,FX,DFXDN,DFXDE,DFXDD)
      DFXDTH = -DFXDN*RTOD/DLAT
      DFXDLN =  DFXDE*RTOD/DLON
      DFXDH  = -HTI*HTI*DFXDD/(RE*DIHT)
      CALL TRILIN (Y(I,J,K),NLAT,NLON,XI,YJ,ZK,FY,DFYDN,DFYDE,DFYDD)
      DFYDTH = -DFYDN*RTOD/DLAT
      DFYDLN =  DFYDE*RTOD/DLON
      DFYDH  = -HTI*HTI*DFYDD/(RE*DIHT)
      CALL TRILIN (Z(I,J,K),NLAT,NLON,XI,YJ,ZK,FZ,DFZDN,DFZDE,DFZDD)
      DFZDTH = -DFZDN*RTOD/DLAT
      DFZDLN =  DFZDE*RTOD/DLON
      DFZDH  = -HTI*HTI*DFZDD/(RE*DIHT)
      CALL TRILIN (V(I,J,K),NLAT,NLON,XI,YJ,ZK,FV,DFVDN,DFVDE,DFVDD)
      DFVDTH = -DFVDN*RTOD/DLAT
      DFVDLN =  DFVDE*RTOD/DLON
      DFVDH  = -HTI*HTI*DFVDD/(RE*DIHT)

      IF (NLAT .LT. 3) RETURN

C          Improve calculation of longitudinal derivatives near poles
      IF (GLAT .GE. DLAT-90.) GO TO 40
      FAC = .5*XI
      OMFAC = 1. - FAC
      XI = XI - 1.
      I = I + 1
      CALL TRILIN (X(I,J,K),NLAT,NLON,XI,YJ,ZK,DMF,DMDFDN,DMDFDE,DMDFDD)
      DFXDLN = DFXDLN*OMFAC + FAC*DMDFDE*RTOD/DLON
      CALL TRILIN (Y(I,J,K),NLAT,NLON,XI,YJ,ZK,DMF,DMDFDN,DMDFDE,DMDFDD)
      DFYDLN = DFYDLN*OMFAC + FAC*DMDFDE*RTOD/DLON
      CALL TRILIN (V(I,J,K),NLAT,NLON,XI,YJ,ZK,DMF,DMDFDN,DMDFDE,DMDFDD)
      DFVDLN = DFVDLN*OMFAC + FAC*DMDFDE*RTOD/DLON

   40 IF (GLAT .LE. 90.-DLAT) GO TO 50
      FAC = .5*(1.- XI)
      OMFAC = 1. - FAC
      XI = XI + 1.
      I = I - 1
      CALL TRILIN (X(I,J,K),NLAT,NLON,XI,YJ,ZK,DMF,DMDFDN,DMDFDE,DMDFDD)
      DFXDLN = DFXDLN*OMFAC + FAC*DMDFDE*RTOD/DLON
      CALL TRILIN (Y(I,J,K),NLAT,NLON,XI,YJ,ZK,DMF,DMDFDN,DMDFDE,DMDFDD)
      DFYDLN = DFYDLN*OMFAC + FAC*DMDFDE*RTOD/DLON
      CALL TRILIN (V(I,J,K),NLAT,NLON,XI,YJ,ZK,DMF,DMDFDN,DMDFDE,DMDFDD)
      DFVDLN = DFVDLN*OMFAC + FAC*DMDFDE*RTOD/DLON
   50 RETURN

C          Error trap diagnostics
 9100 WRITE (MSGU,'(A,'':  Latitude out of range; GPLAT(1),GLAT,GPLAT(NL
     +AT)='',3F10.3)')  CALNM,GPLAT(1),GLAT,GPLAT(NLAT)
      IST = 1
      RETURN
 9200 WRITE (MSGU,'(A,'':  Longitude out of range; GPLON(1),GLON,GPLON(N
     +LON)='',3F10.3)') CALNM,GPLON(1),GLON,GPLON(NLON)
      IST = 1
      RETURN
 9300 WRITE (MSGU,'(A,'':  Altitude out of range; GPALT(1),ALT,GPALT(NAL
     +T)='',3F10.3)')   CALNM,GPALT(1),ALT,GPALT(NALT)
      IST = 1
      RETURN

      END

C*******************************************************************************
      SUBROUTINE TRILIN (U,NLAT,NLON,XI,YJ,ZK,FU,DFUDX,DFUDY,DFUDZ)
C  Trilinear interpolation of u and its derivatives
C 940803 A. D. Richmond
C Inputs:
C   u(1,1,1) = address of lower corner of interpolation box 
C   nlat = first dimension of u from calling routine
C   nlon = second dimension of u from calling routine
C   xi = fractional distance across box in x direction 
C   yj = fractional distance across box in y direction 
C   zk = fractional distance across box in z direction 
C Outputs:
C   fu = interpolated value of u
C   dfudx = interpolated derivative of u with respect to i (x direction)
C   dfudy = interpolated derivative of u with respect to j (y direction)
C   dfudz = interpolated derivative of u with respect to k (z direction)
      DIMENSION U(NLAT,NLON,2)

      IENTRY = 0
      GOTO 5
C*******************************************************************************
      ENTRY TRILINS (U,NLAT,NLON,XI,YJ,ZK,FU)
C  Trilinear interpolation of u only
C 940803 A. D. Richmond
C Inputs and outputs:  see above for definitions
      IENTRY = 1

    5 CONTINUE
      OMXI = 1. - XI
      OMYJ = 1. - YJ
      OMZK = 1. - ZK

      FU = U(1,1,1)*OMXI*OMYJ*OMZK
     2   + U(2,1,1)*XI*OMYJ*OMZK
     3   + U(1,2,1)*OMXI*YJ*OMZK
     4   + U(1,1,2)*OMXI*OMYJ*ZK
     5   + U(2,2,1)*XI*YJ*OMZK
     6   + U(2,1,2)*XI*OMYJ*ZK
     7   + U(1,2,2)*OMXI*YJ*ZK
     8   + U(2,2,2)*XI*YJ*ZK

      IF (IENTRY.NE.0) RETURN

      DFUDX = (U(2,1,1)-U(1,1,1))*OMYJ*OMZK
     2      + (U(2,2,1)-U(1,2,1))*YJ*OMZK
     3      + (U(2,1,2)-U(1,1,2))*OMYJ*ZK
     4      + (U(2,2,2)-U(1,2,2))*YJ*ZK
      DFUDY = (U(1,2,1)-U(1,1,1))*OMXI*OMZK
     2      + (U(2,2,1)-U(2,1,1))*XI*OMZK
     3      + (U(1,2,2)-U(1,1,2))*OMXI*ZK
     4      + (U(2,2,2)-U(2,1,2))*XI*ZK
      DFUDZ = (U(1,1,2)-U(1,1,1))*OMXI*OMYJ
     2      + (U(2,1,2)-U(2,1,1))*XI*OMYJ
     3      + (U(1,2,2)-U(1,2,1))*OMXI*YJ
     4      + (U(2,2,2)-U(2,2,1))*XI*YJ
      RETURN
      END

C*******************************************************************************
      SUBROUTINE ADPL(GLAT,GLON,CTH,STH,FX,FY,FZ,FV
     1  ,DFXDTH,DFYDTH,DFZDTH,DFVDTH,DFXDLN,DFYDLN,DFZDLN,DFVDLN)
C  v is used for vr2n
C  Add-back of pseudodipole component to x,y,z,v and their derivatives.
C  940715 A. D. Richmond
C Inputs:
C   glat = latitude, degrees
C   glon = longitude, degrees
C   fx = interpolated value of x
C   fy = interpolated value of y
C   fz = interpolated value of z
C   fv = interpolated value of v
C   dfxdth,dfydth,dfzdth,dfvdth = derivatives of x,y,z,v with respect to 
C	colatitude, in radians-1
C   dfxdln,dfydln,dfzdln,dfvdln = derivatives of x,y,z,v with respect to 
C	longitude, in radians-1
C Output:
C   cth,sth = cos(colatitude), sin(colatitude)
C   fx = interpolated value of x
C   fy = interpolated value of y
C   fz = interpolated value of z
C   fv = interpolated value of v
C   dfxdth,dfydth,dfzdth,dfvdth = derivatives of x,y,z,v with respect to 
C	colatitude, in radians-1
C   dfxdln,dfydln,dfzdln,dfvdln = derivatives of x,y,z,v with respect to 
C	longitude, in radians-1
      COMMON /APXCON/ RTOD,DTOR,RE,REQ,MSGU,POLA
      COMMON /APXDIPL/ COLAT,ELON,VP,CTP,STP
      CPH = COS((GLON-ELON)*DTOR)
      SPH = SIN((GLON-ELON)*DTOR)
      CTH = SIN(GLAT*DTOR)
      STH = COS(GLAT*DTOR)
      CTM = CTP*CTH + STP*STH*CPH
      FX = FX + STH*CTP*CPH - CTH*STP
      FY = FY + STH*SPH
      FZ = FZ + CTM
      FV = FV - CTM
      DFXDTH = DFXDTH + CTP*CTH*CPH + STP*STH
      DFYDTH = DFYDTH + CTH*SPH
      DFZDTH = DFZDTH - CTP*STH + STP*CTH*CPH
      DFVDTH = DFVDTH + CTP*STH - STP*CTH*CPH
      DFXDLN = DFXDLN - CTP*STH*SPH
      DFYDLN = DFYDLN + STH*CPH
      DFZDLN = DFZDLN - STP*STH*SPH
      DFVDLN = DFVDLN + STP*STH*SPH
      RETURN
      END

C*******************************************************************************
      SUBROUTINE ADPLSC (GLAT,GLON,FX,FY,FZ)
C  Add-back of pseudodipole component to x,y,z 
C  940801 A. D. Richmond
C Inputs:
C   glat = latitude, degrees
C   glon = longitude, degrees
C   fx = interpolated value of x
C   fy = interpolated value of y
C   fz = interpolated value of z
C Output:
C   fx = interpolated value of x
C   fy = interpolated value of y
C   fz = interpolated value of z
      COMMON /APXCON/ RTOD,DTOR,RE,REQ,MSGU,POLA
      COMMON /APXDIPL/ COLAT,ELON,VP,CTP,STP
      CPH = COS((GLON-ELON)*DTOR)
      SPH = SIN((GLON-ELON)*DTOR)
      CTH = SIN(GLAT*DTOR)
      STH = COS(GLAT*DTOR)
      CTM = CTP*CTH + STP*STH*CPH
      FX = FX + STH*CTP*CPH - CTH*STP
      FY = FY + STH*SPH
      FZ = FZ + CTM
      RETURN
      END

C*******************************************************************************
      SUBROUTINE GRADXYZV (ALT,CTH,STH,
     +          DFXDTH,DFYDTH,DFZDTH,DFVDTH,DFXDLN,DFYDLN,DFZDLN,DFVDLN,
     +          DFXDH,DFYDH,DFZDH,DFVDH,GRADX,GRADY,GRADZ,GRADV)
C          Calculates east,north,up components of gradients of x,y,z,v in
C          geodetic coordinates.  All gradients are in inverse km.  Assumes
C          flatness of 1/298.25 and equatorial radius (REQ) of 6378.16 km.
C          940803 A. D. Richmond
      COMMON /APXGEOD/ RHO,DDISDTH
      DIMENSION GRADX(3),GRADY(3),GRADZ(3),GRADV(3)
      IENTRY = 0
      GOTO 5
C*******************************************************************************
      ENTRY GRAPXYZV (ALT,CTH,STH,
     +              DFXDLN,DFYDLN,DFZDLN,DFVDLN,GRADX,GRADY,GRADZ,GRADV)
C          Calculates east component of gradient near pole.
C          940803 A. D. Richmond
C          Inputs and outputs:  see above for definitions
      IENTRY = 1

    5 CONTINUE
      D2 = 40680925.E0 - 272340.E0*CTH*CTH
C          40680925. = req**2 (rounded off)
C          272340.   = req**2 * E2, where E2 = (2. - 1./298.25)/298.25
C                      is square of eccentricity of ellipsoid.
      D = SQRT(D2)
      RHO = STH*(ALT + 40680925.E0/D)
      DDDTHOD = 272340.E0*CTH*STH/D2
      DRHODTH = ALT*CTH + (40680925.E0/D)*(CTH-STH*DDDTHOD)
      DZETDTH =-ALT*STH - (40408585.E0/D)*(STH+CTH*DDDTHOD)
      DDISDTH = SQRT(DRHODTH*DRHODTH + DZETDTH*DZETDTH)
      GRADX(1) = DFXDLN/RHO
      GRADY(1) = DFYDLN/RHO
      GRADZ(1) = DFZDLN/RHO
      GRADV(1) = DFVDLN/RHO

      IF (IENTRY .NE. 0) RETURN

      GRADX(2) = -DFXDTH/DDISDTH
      GRADY(2) = -DFYDTH/DDISDTH
      GRADZ(2) = -DFZDTH/DDISDTH
      GRADV(2) = -DFVDTH/DDISDTH
      GRADX(3) = DFXDH
      GRADY(3) = DFYDH
      GRADZ(3) = DFZDH
      GRADV(3) = DFVDH

      RETURN
      END

C*******************************************************************************
      SUBROUTINE GRADLPV (HR,ALT,FX,FY,FZ,FV,GRADX,GRADY,GRADZ,GRADV,
     +              XLATM,XLONM,VMP,GRCLM,CLMGRP,QDLAT,RGRLP,B,CLM,R3_2)
C          Uses gradients of x,y,z,v to compute geomagnetic field and
C          gradients of apex latitude, longitude.
C          940819 A. D. Richmond
C          INPUT:
C            HR     = reference altitude
C            ALT    = altitude
C            FX,FY,FZ,FV = interpolated values of x,y,z,v, plus pseudodipole
C                     component
C            GRADX,GRADY,GRADZ,GRADV = interpolated gradients of x,y,z,v,
C                     including pseudodipole components (east,north,up)
C          OUTPUT:
C            XLATM  = modified apex latitude (lambda_m), degrees
C            XLONM  = apex longitude (phi_a), degrees
C            VMP    = magnetic potential, in T.m.
C            GRCLM  = grad(cos(lambda_m)), in km-1
C            CLMGRP = cos(lambda_m)*grad(phi_a), in km-1
C            QDLAT  = quasi-dipole latitude, degrees
C            RGRLP  = (re + alt)*grad(lambda')
C            B      = magnetic field, in nT
C            CLM    = cos(lambda_m)
C            R3_2   = ((re + alt)/(re + hr))**(3/2)

      DIMENSION GRADX(3),GRADY(3),GRADZ(3),GRADV(3),
     +          GRCLM(3),CLMGRP(3),RGRLP(3),B(3)

      COMMON /APXDIPL/ COLAT,ELON,VP,CTP,STP
      COMMON /APXCON/ RTOD,DTOR,RE,REQ,MSGU,POLA

      RR = RE + HR
      R  = RE + ALT
      RN = R/RE
      SQRROR = SQRT(RR/R)
      R3_2 = 1./SQRROR/SQRROR/SQRROR
      XLONM = ATAN2(FY,FX)
      CPM = COS(XLONM)
      SPM = SIN(XLONM)
      XLONM = RTOD*XLONM
      BO = VP*1.E6
C             1.E6 converts T to nT and km-1 to m-1.
      RN2 = RN*RN
      VMP = VP*FV/RN2
      B(1) = -BO*GRADV(1)/RN2
      B(2) = -BO*GRADV(2)/RN2
      B(3) = -BO*(GRADV(3)-2.*FV/R)/RN2

      X2PY2 = FX*FX + FY*FY
      XNORM = SQRT(X2PY2 + FZ*FZ)
      XLP = ATAN2(FZ,SQRT(X2PY2))
      SLP = SIN(XLP)
      CLP = COS(XLP)
      QDLAT = XLP*RTOD
      CLM = SQRROR*CLP
      IF (CLM.LE.1.) GOTO 5
      WRITE (6,*) 'Stopped in gradlpv because point lies below field lin      
     1e that peaks at reference height.'
      STOP
    5 XLATM = RTOD*ACOS(CLM)
C  If southern magnetic hemisphere, reverse sign of xlatm
      IF (SLP.LT.0.) XLATM = - XLATM
      DO 10 I=1,3
	GRCLP = CPM*GRADX(I) + SPM*GRADY(I)
	RGRLP(I) = R*(CLP*GRADZ(I) - SLP*GRCLP)
	GRCLM(I) = SQRROR*GRCLP
   10   CLMGRP(I) = SQRROR*(CPM*GRADY(I)-SPM*GRADX(I))
      GRCLM(3) = GRCLM(3) - SQRROR*CLP/(2.*R)
      RETURN
C*******************************************************************************
      ENTRY XYZ2APX (ALT,FX,FY,FZ,A,ALAT,ALON,IERR)
C          Computes apex latitude, longitude.
C          990309 A. D. Richmond
C          INPUT:
C            ALT      = altitude
C            FX,FY,FZ = interpolated values of x,y,z, plus pseudodipole
C                       component
C          OUTPUT:
C            A    = apex radius, normalized by req
C            ALAT = apex latitude, degrees
C            ALON = apex longitude, degrees
C
C          Mod (Mar 99):  Lines 19-30 are changed from the original in order
C          to avoid a truncation error near the magnetic equator.  What was
C          done is to make use of the identity
C
C                  SIN(ALAT/RTOD)**2 + COS(ALAT/RTOD)**2 = 1,
C
C          so that since
C
C                  COS(ALAT/RTOD)**2 = 1/A (Eq. 1),
C
C          then
C
C                  SIN(ALAT/RTOD)**2 = (A-1)/A (Eq. 2)
C
C          Below AM1 = A-1.  If AM1 is less than 1, use Eq. 1;
C          otherwise use Eq. 2.  Mathematically, both equations should
C          give identical results, but numerically it is better to use
C          that function ASIN or ACOS that has the smaller argument.
C          The jump from one equation to the other occurs at ALAT = 45.

      IERR  = 0
      ALON  = ATAN2(FY,FX)*RTOD
      SLP2  = FZ*FZ
      X2PY2 = FX*FX + FY*FY
      XNORM = SLP2 + X2PY2
      SLP2  = SLP2/XNORM
      CLP2  = X2PY2/XNORM
      AM1   = (RE*SLP2 + ALT)/(REQ*CLP2)
      A = 1. + AM1

      IF (AM1.LT.0.) THEN
        IERR = 1
        WRITE (6,*) 'Missing alat returned because point lies below fiel
     1d line that peaks at Earth surface.'
        RETURN
      ELSEIF (AM1.LT.1.) THEN
        ALAT = RTOD*ASIN(SQRT(AM1/A))
      ELSE
        ALAT = RTOD*ACOS(1./SQRT(A))
      ENDIF
C  If southern magnetic hemisphere, reverse sign of alat
      IF (FZ.LT.0.) ALAT = - ALAT
      RETURN
      END

C*******************************************************************************
      SUBROUTINE BASVEC (HR,XLATM,GRCLM,CLMGRP,RGRLP,B,CLM,R3_2,
     +                  BMAG,SIM,SI,F,D,W,BHAT,D1,D2,D3,E1,E2,E3,F1,F2)
C          Computes base vectors and other parameters for apex coordinates.
C          Vector components:  east, north, up
C          940801 A. D. Richmond
C          Reference:
C            Richmond, A. D., Ionospheric Electrodynamics Using Magnetic Apex
C            Coordinates, J. Geomag. Geoelectr., 47, 191-212, 1995.
C          INPUTS:
C            HR     = reference altitude
C            XLATM  = modified apex latitude, degrees
C            GRCLM  = grad(cos(lambda_m)), in km-1
C            CLMGRP = cos(lambda_m)*grad(phi_a), in km-1
C            RGRLP  = (re + altitude)*grad(lambda')
C            B      = magnetic field, in nT
C            CLM    = cos(lambda_m)
C            R3_2   = ((re + altitude)/(re + hr))**(3/2)
C          RETURNS:
C            BMAG    = magnitude of magnetic field, in nT
C            SIM     = sin(I_m) of article
C            SI      = sin(I)
C            F       = F of article
C            D       = D of article
C            W       = W of article
C            BHAT    = unit vector along geomagnetic field direction
C            D1...F2 = base vectors of article
      COMMON /APXCON/ RTOD,DTOR,RE,REQ,MSGU,POLA
      DIMENSION GRCLM(3),CLMGRP(3),RGRLP(3),B(3)
      DIMENSION BHAT(3),D1(3),D2(3),D3(3),E1(3),E2(3),E3(3),F1(2),F2(2)        
      RR = RE + HR
      SIMOSLM = 2./SQRT(4. - 3.*CLM*CLM)
      SIM = SIMOSLM*SIN(XLATM*DTOR)
      BMAG = SQRT(B(1)*B(1) + B(2)*B(2) + B(3)*B(3))
      D1DB = 0.
      D2DB = 0.
      DO 10 I=1,3
        BHAT(I) = B(I)/BMAG
        D1(I) = RR*CLMGRP(I)
        D1DB = D1DB + D1(I)*BHAT(I)
        D2(I) = RR*SIMOSLM*GRCLM(I)
   10   D2DB = D2DB + D2(I)*BHAT(I)
C Ensure that d1,d2 are exactly perpendicular to B:
      DO 15 I=1,3
        D1(I) = D1(I) - D1DB*BHAT(I)
   15   D2(I) = D2(I) - D2DB*BHAT(I)
      E3(1) = D1(2)*D2(3) - D1(3)*D2(2)
      E3(2) = D1(3)*D2(1) - D1(1)*D2(3)
      E3(3) = D1(1)*D2(2) - D1(2)*D2(1)
      D = BHAT(1)*E3(1) + BHAT(2)*E3(2) + BHAT(3)*E3(3)
      DO 20 I=1,3
        D3(I) = BHAT(I)/D
C Following step may be unnecessary, but it ensures that e3 lies along bhat.
        E3(I) = BHAT(I)*D
   20   CONTINUE
      E1(1) = D2(2)*D3(3) - D2(3)*D3(2)
      E1(2) = D2(3)*D3(1) - D2(1)*D3(3)
      E1(3) = D2(1)*D3(2) - D2(2)*D3(1)
      E2(1) = D3(2)*D1(3) - D3(3)*D1(2)
      E2(2) = D3(3)*D1(1) - D3(1)*D1(3)
      E2(3) = D3(1)*D1(2) - D3(2)*D1(1)
      W = RR*RR*CLM*ABS(SIM)/(BMAG*D)
      SI = -BHAT(3)
      F1(1) =  RGRLP(2) 
      F1(2) = -RGRLP(1)
      F2(1) = -D1(2)*R3_2
      F2(2) =  D1(1)*R3_2
      F = F1(1)*F2(2) - F1(2)*F2(1)
      RETURN
      END

C*******************************************************************************
      SUBROUTINE CKGP (CALNM,MSGUN,NLAT,NLON,NALT,GPLAT,GPLON,GPALT,IST)
C          Check grid point values tests extremes and order of the grid
C          point arrays, producing diagnostics to MSGUN and IST=1 when
C          rules have been broken.
      DIMENSION GPLAT(NLAT), GPLON(NLON), GPALT(NALT)
      CHARACTER*(*) CALNM

      IST = 1
      OLAT = -90.
      DO 10 I=1,NLAT
      IF (ABS (GPLAT(I)) .GT.  90.) GO TO 9100
      IF (     GPLAT(I)  .LT. OLAT) GO TO 9200
   10 OLAT = GPLAT(I)

      OLON = -270.
      DO 20 I=1,NLON
      IF (ABS (GPLON(I)) .GT. 270.) GO TO 9300
      IF (     GPLON(I)  .LT. OLON) GO TO 9400
   20 OLON = GPLON(I)

      OALT = 0.
      DO 30 I=1,NALT
      IF (GPALT(I) .LT.   0.) GO TO 9500
      IF (GPALT(I) .LT. OALT) GO TO 9600
   30 OALT = GPALT(I)

      IST = 0

  100 RETURN

 9100 WRITE (MSGUN,'(A,'':  |GPLAT(I)| > 90; I,GPLAT(I)'',I5,F10.3)')
     +            CALNM,                     I,GPLAT(I)
      GO TO 100
 9200 WRITE (MSGUN,'(A,'':  GPLAT(I) < GPLAT(I-1); I,GPLAT(I),GPLAT(I-1)
     +='',I5,2F10.3)')                     CALNM,  I,GPLAT(I),OLAT
      GO TO 100
 9300 WRITE (MSGUN,'(A,'':  |GPLON(I)| > 180; I,GPLON(I)'',I5,F10.3)')
     +                                 CALNM, I,GPLON(I)
      GO TO 100
 9400 WRITE (MSGUN,'(A,'':  GPLON(I) < GPLON(I-1); I,GPLON(I),GPLON(I-1)
     +='',I5,2F10.3)')                    CALNM, I,GPLON(I),OLON
      GO TO 100
 9500 WRITE (MSGUN,'(A,'':  GPALT(I) <  0; I,GPALT(I)'',I5,F10.3)')
     +                                CALNM, I,GPALT(I)
      GO TO 100
 9600 WRITE (MSGUN,'(A,'':  GPALT(I) < GPALT(I-1); I,GPALT(I),GPALT(I-1)
     +='',I5,2F10.3)')                      CALNM, I,GPALT(I),OALT
      GO TO 100
      END

C*******************************************************************************
      SUBROUTINE MAKEXYZV (EPOCH,NLAT,NLON,NALT,GPLAT,GPLON,GPALT,
     +                                                          X,Y,Z,V)
C          Sets up grid arrays for later interpolation
C          940822 A. D. Richmond, NCAR
C          INPUT:
C            EPOCH = year and fraction (e.g., 1994.50 for 1994 July 2)
C            NLAT,NLON,NALT = triple dimensions of X,Y,Z,V and respective
C                    single dimensions of GP___ arrays
C            GPLAT,GPLON,GPALT = grid point latitudes, longitudes and altitudes
C          OUTPUT:
C            X = array containing cos(lambda')cos(phi_a) less pseudodipole
C                component
C            Y = array containing cos(lambda')sin(phi_a) less pseudodipole
C                component
C            Z = array containing sin(lambda') less pseudodipole component
C            V = array containing ((magnetic potential)/vp)*((re+height)/re)**2,
C                less pseudodipole component
C
C          Modification (99 Mar):  Make X,Y,Z,V constant near the poles
C          for all GPLON(j) at each height.  Add POLA to APXCON

      DIMENSION X(NLAT,NLON,NALT), Y(NLAT,NLON,NALT), Z(NLAT,NLON,NALT),
     +          V(NLAT,NLON,NALT), GPLAT(NLAT), GPLON(NLON), GPALT(NALT)


C            POLA       = Pole angle (deg); when the geographic latitude is
C                         poleward of POLA, X,Y,Z,V are forced to be constant.
C                         for all longitudes at each altitude.  POLA is defined
C                         in APXMKA (POLA = 90. - SQRT (PRECISE) * RTOD),
C                         which currently makes POLA = 89.995
      COMMON /APXCON/ RTOD,DTOR,RE,REQ,MSGU,POLA
      COMMON /APXDIPL/ COLAT,ELON,VP,CTP,STP

      CALL COFRM (EPOCH)
      CALL DYPOL (COLAT,ELON,VP)
      CTP = COS (COLAT*DTOR)
      STP = SIN (COLAT*DTOR)
      REQORE = REQ/RE
      RQORM1 = REQORE-1.

      DO 100 I=1,NLAT
      CT = SIN (GPLAT(I)*DTOR)
      ST = COS (GPLAT(I)*DTOR)
      KPOL = 0
      IF (ABS (GPLAT(I)) .GT. POLA) KPOL = 1

      DO 100 J=1,NLON
      IF (KPOL .EQ. 0) GO TO 20
      IF (J    .EQ. 1) GO TO 20
C          KPOL = 1 (poleward of POLA) and first lon's XYZV are defined
      DO 10 K=1,NALT
      V(I,J,K) = V(I,1,K)
      X(I,J,K) = X(I,1,K)
      Y(I,J,K) = Y(I,1,K)
   10 Z(I,J,K) = Z(I,1,K)
      GO TO 100

   20 CP  = COS ((GPLON(J)-ELON)*DTOR)
      SP  = SIN ((GPLON(J)-ELON)*DTOR)
C           ctm   is pseudodipole component of z
C          -ctm   is pseudodipole component of v
C          stmcpm is pseudodipole component of x
C          stmspm is pseudodipole component of y
      CTM    = CTP*CT + STP*ST*CP
      STMCPM = ST*CTP*CP - CT*STP
      STMSPM = ST*SP

      DO 30 K=1,NALT
      CALL APEX (EPOCH,GPLAT(I),GPLON(J),GPALT(K),
     +           A,ALAT,PHIA,BMAG,XMAG,YMAG,ZDOWN,VMP)
      VNOR = VMP/VP
      RP = 1. + GPALT(K)/RE
      V(I,J,K) = VNOR*RP*RP + CTM
      REQAM1 = REQ*(A-1.)
      SLP = SQRT(AMAX1(REQAM1-GPALT(K),0.)/(REQAM1+RE))
C          Reverse sign of slp in southern magnetic hemisphere
      IF (ZDOWN.LT.0.) SLP = -SLP
      CLP = SQRT (RP/(REQORE*A-RQORM1))
      PHIAR = PHIA*DTOR
      X(I,J,K) = CLP*COS (PHIAR) - STMCPM
      Y(I,J,K) = CLP*SIN (PHIAR) - STMSPM
      Z(I,J,K) = SLP - CTM
   30 CONTINUE

  100 CONTINUE

      RETURN
      END

C*******************************************************************************

      SUBROUTINE SETMISS (XMISS
     2 ,XLATM,XLONM,VMP,B,BMAG,BE3,SIM,SI,F,D,W
     3 ,BHAT,D1,D2,D3,E1,E2,E3,F1,F2)
      DIMENSION BHAT(3),D1(3),D2(3),D3(3),E1(3),E2(3),E3(3),F1(2),F2(2)        
     2 ,B(3)
      XLATM = XMISS
      XLONM = XMISS
      VMP = XMISS
      BMAG = XMISS
      BE3 = XMISS
      SIM = XMISS
      SI = XMISS
      F = XMISS
      D = XMISS
      W = XMISS
      DO 5 I=1,3
	B(I) = XMISS
	BHAT(I) = XMISS
	D1(I) = XMISS
	D2(I) = XMISS
	D3(I) = XMISS
	E1(I) = XMISS
	E2(I) = XMISS
    5   E3(I) = XMISS
      DO 6 I=1,2
	F1(I) = XMISS
    6   F2(I) = XMISS
      RETURN
      END

C*******************************************************************************
      SUBROUTINE GM2GC (GMLAT,GMLON,GCLAT,GCLON)
C  Converts geomagnetic to geocentric coordinates.
C  940819 A. D. Richmond
C
C  Inputs:
C	gmlat = geomagnetic latitude in degrees
C	gmlon = geomagnetic longitude in degrees
C  Outputs:
C	gclat = geocentric latitude in degrees
C	gclon = geocentric longitude in degrees
C
C  Common/consts/
C	rtod, dtor = 180/pi, pi/180
C	re, req    = 6371.2, 6378.160
      COMMON /APXCON/ RTOD,DTOR,RE,REQ,MSGU,POLA
C  Common/dipol/
C       colat = geocentric colatitude of north geomagnetic pole, in degrees
C	elon  = geocentric east longitude of north geomagnetic pole, in degrees
C	vp    = magnetic potential at 1 RE, south geomagnetic pole
C	ctp   = cos(colat*dtor)
C	stp   = sin(colat*dtor)
      COMMON /APXDIPL/ COLAT,ELON,VP,CTP,STP
C
      STM = COS(GMLAT*DTOR)
      CTM = SIN(GMLAT*DTOR)
      CTC = CTP*CTM - STP*STM*COS(GMLON*DTOR)
      CTC = AMIN1(CTC,1.)
      CTC = AMAX1(CTC,-1.)
      GCLAT = ASIN(CTC)*RTOD
      GCLON = ATAN2(STP*STM*SIN(GMLON*DTOR),CTM-CTP*CTC)
      GCLON = GCLON*RTOD + ELON
      IF (GCLON.LT.-180.) GCLON = GCLON + 360.
      RETURN
      END
