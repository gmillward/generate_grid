      PROGRAM APEX2000
C          Program to set up global Apex arrays to L=4
C          Need to compile with apxntrpb4lf.f, ggrid.f, apex.f, and 
C            magfld.f, which are in directory /home/bozo/apex/src

C          FILNAM will contain gridded arrays for the epoch given by "DATE."
      CHARACTER*80 FILNAM

C          Declarations needed for APXMKA, APXWRA, or APXRDA
C            MSGUN = Fortran unit number for diagnostics
C            IUN   = Fortran unit number for I/O
C            MLAT,MLON,MALT = Maximum number of grid latitudes, longitudes
C                    and altitudes.
C            NGRF  = Number of epochs in the current DGRF/IGRF; see COFRM in
C                    file magfld.f
C MALT never needs to be larger than NVERT.  The lower the maximum
C   altitude ALTMX, the lower MALT can be.  The value of NALT returned
C   upon calling GGRID is the minimum value of MALT. 
C MALT=24 is large enough to hold altitudes up to L=4 when NVERT=30.
      PARAMETER (MSGUN=6, IUN=12, MLAT=91,MLON=151,MALT=24, NGRF=8,
     +           LWK= MLAT*MLON*MALT*5 + MLAT+MLON+MALT)
      DIMENSION GPLAT(MLAT),GPLON(MLON),GPALT(MALT),WK(LWK),TGRF(NGRF)
      DATA TGRF /1965.,1970.,1975.,1980.,1985.,1990.,1995.,2000./

C          Sample date (corresponding to 2000 July 2)
c     DATE  = 1980.00

      read(5,*) date

      NDATE = 1
C          Set up IGRF coefficients for the date desired
      CALL COFRM (DATE)

C          Name of sample write files to be created 
      FILNAM = 'Apex_grid_data'

C Create global 3-D grid out to L=4 at resolution corresponding to 
C   NVERT=30 for one time (DATE), and store in FILNAM: 
      NVERT = 30
      GLAMN = -90.
      GLAMX = 90.
      GLOMN = -180.
      GLOMX = 180.
      ALTMN = 0.
      ALTMX = 3.*6371.2
      CALL GGRID (NVERT,GLAMN,GLAMX,GLOMN,GLOMX,ALTMN,ALTMX,
     +            GPLAT,GPLON,GPALT,NLAT,NLON,NALT)

      WRITE (6,'(1X,/,''GGRID produced for NVERT='',I3,'' and limits: la
     +titude '',F8.3,'' - '',F8.3,/,40X,''longitude '',F8.3,'' - '',F8.3
     +,/,40X,'' altitude '',F8.3,'' - '',F8.3)') NVERT,GLAMN,GLAMX,
     +                                     GLOMN,GLOMX,ALTMN,ALTMX
      WRITE (6,'(9X,''NLAT ='',I4,''  Latitudes:'',5F8.3,/,(31X,5F8.3))'
     +      )    NLAT,(GPLAT(I),I=1,NLAT)
      WRITE (6,'(9X,''NLON ='',I4,'' Longitudes:'',5F8.3,/,(31X,5F8.3))'
     +      )    NLON,(GPLON(I),I=1,NLON)
      WRITE (6,'(9X,''NALT ='',I4,''  Altitudes:'',5F8.1,/,(31X,5F8.1))'
     +      )    NALT,(GPALT(I),I=1,NALT)

      CALL APXWRA (MSGUN, FILNAM,IUN, DATE,NDATE,
     +            GPLAT,GPLON,GPALT,NLAT,NLON,NALT, WK,LWK, IST)
      IF (IST .NE. 0) write(6,*) 'Error: IST= ',IST
      write(6,*) 'Finished OK'
      END

      SUBROUTINE APXMKA (MSGUN, EPOCH, GPLAT,GPLON,GPALT,NLAT,NLON,NALT,
     +                  WK,LWK, IST)

C          This translates from geographic coordinates to geomagnetic and back,
C          providing Apex coordinates (VanZandt et. al., 1972), Modified Apex
C          coordinates and Quasi-Dipole coordinates (Richmond, 1995).  This is
C          intended as a faster alternative to direct Apex calculation and it is
C          the officially supported code for Modified Apex and Quasi-Dipole
C          coordinates.  A variation of Apex coordinates, Modified Apex
C          coordinates are designed to make the magnetic latitude coordinate
C          continuous across the magnetic equator at a given reference altitude
C          and constant along magnetic field lines at other heights; hence,
C          Modified Apex coordinates are suited for organizing phenomona along
C          magnetic field lines.  In contrast, Quasi-Dipole latitude and
C          longitude are approximately constant in geographic altitude so they
C          are suited for organizing horizontally stratified phenomona such as
C          ionospheric currents.
C
C          Initially, external routines (apex.f and magfld.f) are called to
C          derive Apex coordinates from which tables are prepared here,
C          optionally saved or read back, and ultimately used for coordinate
C          conversion by interpolation.  Whereas apex.f provides Apex
C          coordinates, it does so by tracing magnetic field lines to their apex
C          which is slower than the interpolation performed here.  Moreover,
C          this routine provides two additional coordinate systems, their base
C          vectors, other quantities related to gradients, and the reverse
C          translation to geographic coordinates.
C
C          These advantages are only possible after first creating the look-up
C          tables and associated parameters.  They are created for one time and
C          held in memory by calling
C
C            APXMKA - make magnetic arrays
C
C          or they may be created and written, then read back later by calling
C
C            APXWRA - make and write arrays
C            APXRDA - read stored arrays.
C
C          APXWRA can create tables for multiple times (e.g., the IGRF dates)
C          but only one batch is held in memory because they can be large.
C          Consequently, when reading back a file with multiple times, APXRDA
C          interpolates or extrapolates to a single specified time.  Instead, if
C          APXWRA writes tables for a single time, then it is not necessary to
C          call APXRDA because tables for that time are already in memory.
C
C          Creation of tables with global extent can be time consuming so it is
C          expedient to use APXWRA in a program which is seldom executed; i.e.,
C          only when making new tables.  Afterward, other programs use APXRDA to
C          read those tables before proceeding with magnetic coordinate
C          translations.  Example programs (mkglob.f, xapxntrp.f, xglob.f) are
C          given in separate files.
C
C          After loading into memory tables for a specified time, entries here
C          provide the following translations of geographic (geodetic) latitude,
C          longitude, altitude:
C
C            APXALL  - geographic to Apex coordinates
C            APXMALL - geographic to Modified Apex and Quasi-Dipole coordinates
C            APXQ2G  - Quasi-Dipole  to geographic coordinates.
C            APXA2G  - Apex          to geographic coordinates.
C            APXM2G  - modified apex to geographic coordinates.
C
C          Geographic coordinates of points defining each table are stored with
C          the tables and may be retrieved by calling
C
C            APXGGC - get grid coordinates.
C
C          Details for each of these entries are given below after REFERENCES
C          and ALGORITHM; last are EXTERNALS, INSTALLATION SPECIFICS and
C          HISTORY.
C------------------------------------------------------------------------------
C          REFERENCES:
C
C          Richmond, A. D., Ionospheric Electrodynamics Using Magnetic Apex
C               Coordinates, J. Geomag. Geoelectr., 47, 191-212, 1995.
C          VanZandt, T.E., W.L.Clark, and J.M.Warnock, Magnetic Apex
C               Coordinates: A magnetic coordinate system for the ionospheric
C               F2 layer, J Geophys. Res., 77, 2406-2411, 1972.
C------------------------------------------------------------------------------
C          ALGORITHM:
C
C          When arrays are created, APXMKA calls subroutine MAKEXYZV, which in
C          turn calls external subroutine APEX at each grid point to get the
C          apex radius (A), the apex longitude (PHIA), and the magnetic
C          potential (VMP).  The cosine (CLP) and sine (SLP) of the Quasi-Dipole
C          latitude are computed from A.  From these preliminary quantites can
C          be computed defined as
C
C            x = cos(Quasi-Dipole latitude)*cos(apex longitude)
C            y = cos(Quasi-Dipole latitude)*sin(apex longitude)
C            z = sin(Quasi-Dipole latitude)
C            v = (VMP/VP)*((RE+ALT)/RE)**2
C
C          where VP is the magnitude of the magnetic potential of the
C          geomagnetic dipole at a radius of RE; ALT is altitude; and RE is the
C          mean Earth radius.  Note that x, y, z and v vary smoothly across the
C          apex and poles, unlike apex latitude or Quasi-Dipole latitude, so
C          that they can be linearly interpolated near these poles.
C          Corresponding values of x,y,z,v for a dipole field on a spherical
C          Earth are computed analytically and subtracted from the above
C          quantities, and the results are put into the 3D arrays X,Y,Z,V.  When
C          APXALL or APXMALL is called, trilinear interpolations (in latitude,
C          longitude, and inverse altitude) are carried out between the grid
C          points.  Gradients are calculated for APXMALL in order to determine
C          the base vectors.  Analytic formulas appropriate for a dipole field
C          on a spherical Earth are used to determine the previously removed
C          dipole components of x,y,z,v, and their gradients and these are added
C          back to the interpolated values obtained for the non-dipole
C          components.  Finally, the apex-based coordinates and their base
C          vectors are calculated from x,y,z,v and their gradients.
C------------------------------------------------------------------------------
C          APXMKA and APXWRA:
C
C          These create tables which are used later by APXALL, APXMALL, APXQ2G,
C          APXA2G or APXM2G.  The tables are computed for a grid of geographic
C          locations (latitude, longitude and altitude) which are inputs; i.e.
C          one must define a suitable grid, a compromise between global extent
C          and resolution.  The smallest possible grid (2,2,2) with small
C          spatial increments bounding the point of interest will insure high
C          interpolation accuracy at the risk of round-off errors in vector
C          quantities which are determined from differences between adjacent
C          grid point values.  Such degraded accuracy has been observed in
C          east-west gradients very close to the poles using a 301x501x15 grid
C          which has increments 0.6x0.72 (deg. lat. by lon.).  This introduces
C          the other extreme, a global array extending to several Earth radii
C          which can produce huge tables.  Examples are given separately in
C          program xapxntrp.f and the global table generator program mkglob.f;
C          for global coverage up to at least 1000 km (ALTMX=1000) we have used
C          two grids, whose resolution is specified by NVERT:
C
C                     Dimension    file size MB   Number of
C             NVERT (latxlonxalt) (32-bit words) Times (epochs)
C               30     91,121,7      15.8            8
C               40    121,201,7      32.7            8
C               40    121,201,7      36.8            9
C
C          Use APXMKA to create interpolation tables for a single time or use
C          APXWRA create the tables for one or more times and the tables in a
C          file.
C
C             CALL APXMKA (MSGUN, EPOCH, GPLAT,GPLON,GPALT,NLAT,NLON,NALT,
C            +            WK,LWK, IST)
C
C             CALL APXWRA (MSGUN, FILNAM,IUN, EPOCH,NEPOCH,
C            +            GPLAT,GPLON,GPALT,NLAT,NLON,NALT, WK,LWK, IST)
C
C          INPUTS:
C            MSGUN  = Fortran unit number to write diagnostic messages.
C            EPOCH  = Time formatted as UT year and fraction; e.g., 1990.0
C                     is 0 UT 1 Jan 1990.  For APXMKA EPOCH is single valued;
C                     for APXWRA EPOCH contains NEPOCH values.
C            GPLAT  = Grid point latitudes, an array of NLAT real values in
C                     ascending numerical order with each value in the range
C                     -90. to +90. degrees with positive north.
C            GPLON  = Grid point longitudes, an array of NLON real values in
C                     ascending numerical order with each value in the range
C                     -270. to 270. degrees with positive east.
C            GPALT  = Grid point altitudes, an array of NALT real values in
C                     ascending numerical order with none less than zero in
C                     units of km above ground.
C            NLAT   = Number of latitudes  in GPLAT.
C            NLON   = Number of longitudes in GPLON.
C            NALT   = Number of altitudes  in GPALT.
C            WK     = Work array which should not be altered between
C                     initialization (APXMKA, APXWRA or APXRDA) and use (APXGGC,
C                     APXALL, APXMALL, APXQ2G, APXA2G, APXM2G).
C            LWK    = Dimension of WK >=  NLAT*NLON*NALT*5 + NLAT+NLON+NALT
C
C          Additional INPUTS for APXWRA only:
C            FILNAM = file name where arrays are stored.
C            IUN    = Fortran unit number to be associated with FILNAM.
C            NEPOCH = Number of times in EPOCH.
C
C          RETURNS:
C            IST    = Return status: 0 (okay) or > 0 (failed)
C
C          Declarations for APXMKA arguments:
C
C            DIMENSION GPLAT(MLAT), GPLON(MLON), GPALT(MALT), WK(LWK)
C
C          where MLAT, MLON and MALT are, respectively, at least NLAT, NLON and
C          NALT.
C
C          Additional declarations for APXWRA arguments:
C
C            DIMENSION EPOCH(NEPOCH)
C            CHARACTER FILNAM*(LNAM)
C
C          where LNAM is the number of characters in FILNAM.
C------------------------------------------------------------------------------
C          APXRDA:
C
C          Read back tables previously created by calling APXWRA.
C
C            CALL APXRDA (MSGUN, FILNAM,IUN, DATE, WK,LWK, IST)
C
C          INPUTS:
C            MSGUN  = Fortran unit number to write diagnostic messages.
C            FILNAM = file name where arrays are stored.
C            IUN    = Fortran unit number to be associated with FILNAM.
C            DATE   = Time formatted as UT year and fraction; e.g., 1990.0
C                     is 0 UT 1 Jan 1990.  DATE is only used when FILNAM
C                     contains tables for more than one time.
C            WK,LWK = Same as APXMKA
C
C          RETURNS:
C            IST    = status: 0 (okay), -1 (minor date problem) or > 0 (failed)
C
C          Declarations for APXRDA formal arguments:
C
C            CHARACTER FILNAM*(LNAM)
C            DIMENSION WK(LWK)
C
C          where LNAM is the number of characters in FILNAM and WK, LWK are same
C          as APXMKA.
C------------------------------------------------------------------------------
C          APXGGC:
C
C          Get grid coordinates for tables currently in memory; i.e., after
C          calling APXMKA, APXWRA or APXRDA.
C
C            CALL APXGGC (MSGUN, WK,LWK, GPLAT,GPLON,GPALT,NLAT,NLON,NALT, IST)
C
C          INPUTS:
C            MSGUN  = Fortran unit number to write diagnostic messages.
C            WK,LWK = Same as APXMKA
C
C          RETURNS:
C            GPLAT  = Real array of grid point latitudes  in degrees.
C            GPLON  = Real array of grid point longitudes in degrees.
C            GPLAT  = Real array of grid point altitudes  in km.
C            NLAT   = Number of latitudes  in GPLAT
C            NLON   = Number of longitudes in GPLON.
C            NALT   = Number of altitudes  in GPALT.
C            IST    = Return status: 0 (okay).
C
C          Declarations for APXGGC formal arguments:
C
C            DIMENSION GPLAT(MLAT), GPLON(MLON), GPALT(MALT), WK(LWK)
C
C          where MLAT, MLON and MALT are, respectively, at least NLAT, NLON and
C          NALT.
C------------------------------------------------------------------------------
C          APXALL:
C
C          Determine Apex coordinates using tables currently in memory; i.e.,
C          after calling APXMKA, APXWRA, or APXRDA.
C
C            CALL APXALL (GLAT,GLON,ALT, WK, A,ALAT,ALON, IST)
C
C          INPUTS:
C            GLAT = Geographic latitude (degrees) must be within grid latitude
C                   limits.
C            GLON = Geographic longitude in degrees; must be within one
C                   revolution of the grid longitude limits.
C            ALT  = Altitude, km.
C            WK   = same as entry APXMKA
C          RETURNS:
C            A    = Apex radius, normalized by Req
C            ALAT = Apex latitude, degrees
C            ALON = Apex longitude, degrees
C            IST  = Return status:  okay (0); or failure (1).
C
C          Declarations for APXALL formal arguments:
C
C            DIMENSION WK(LWK)
C
C------------------------------------------------------------------------------
C          APXMALL:
C
C          Determine Modified Apex and Quasi-dipole coordinates using tables
C          currently in memory; i.e. after calling APXMKA, APXWRA, or APXRDA.
C          The above reference (Richmond, 1995) describes HR and all returned
C          variables.
C
C             CALL APXMALL (GLAT,GLON,ALT,HR, WK,                               !Inputs
C            +             B,BHAT,BMAG,SI,                                      !Mag Fld
C            +             ALON,                                                !Apx Lon
C            +             XLATM,VMP,WM,D,BE3,SIM,D1,D2,D3,E1,E2,E3,            !Mod Apx
C            +             XLATQD,F,F1,F2 , IST)                                !Qsi-Dpl
C
C          INPUTS:
C            GLAT   = Geographic latitude (degrees); must be within the grid
C                     latitude limits.
C            GLON   = Geographic longitude in degrees; must be within one
C                     revolution of the grid longitude limits.
C            ALT    = Altitude, km above ground.
C            HR     = Modified Apex coordinates reference altitude, km.
C            WK     = same as entry APXMKA
C          RETURNS:
C            B      = magnetic field components (east, north, up), in nT
C            BHAT   = components (east, north, up) of unit vector along the
C                     geomagnetic field direction
C            BMAG   = magnetic field magnitude, nT
C            SI     = sin(I) where I is the angle of inclination of the field
C                     line below the horizontal
C            ALON   = Apex longitude = Modified Apex longitude = Quasi-Dipole
C                     longitude, degrees
C            XLATM  = Modified Apex latitude, degrees
C            VMP    = magnetic potential, in T.m
C            WM     = Wm of reference above, in km**2 /nT; i.e., 10**15 m**2 /T)
C            D      = D of reference above
C            BE3    = B_e3 of reference above (= Bmag/D), in nT
C            SIM    = sin(I_m) of reference above
C            D1,D2,D3,E1,E2,E3 = Modified Apex coordinates base vectors, each with
C                     three components (east, north, up) as described in reference
C                     above
C            XLATQD = Quasi-Dipole latitude, degrees
C            F      = F for Quasi-Dipole coordinates described in reference above
C            F1,F2  = Quasi-Dipole coordinates base vectors with components
C                     (east, north) described in reference above
C            IST    = Return status: 0 (okay) or 1 (failed).
C
C          Declarations for APXMALL formal arguments:
C
C            DIMENSION WK(LWK), B(3),BHAT(3),D1(3),D2(3),D3(3), E1(3),E2(3),
C           +          E3(3), F1(2),F2(2)
C------------------------------------------------------------------------------
C          APXQ2G:
C
C          Convert from Quasi-Dipole to geographic coordinates using tables
C          currently in memory; i.e. after calling APXMKA, APXWRA, or APXRDA.
C
C            CALL APXQ2G (QDLAT,QDLON,ALT, WK, GDLAT,GDLON, IST)
C
C          INPUTS:
C            QDLAT = Quasi-Dipole latitude in degrees
C            QDLON = Quasi-Dipole longitude in degrees
C            ALT   = altitude in km
C            WK    = same as entry APXMKA
C          RETURNS:
C            GDLAT = geodetic latitude in degrees
C            GDLON = geodetic longitude in degrees
C            IST   = status: 0 (okay), -1 (results are not as close as PRECISE),
C                    or > 0 (failed)
C
C          Declarations for APXQ2G formal arguments:
C
C            DIMENSION WK(LWK)
C------------------------------------------------------------------------------
C          APXA2G:
C
C          Convert from Quasi-Dipole to geographic coordinates using tables
C          currently in memory; i.e. after calling APXMKA, APXWRA, or APXRDA.
C
C            CALL APXA2G (ALAT,ALON,ALT, WK, GDLAT,GDLON, IST)
C
C          INPUTS:
C            ALAT  = Apex latitude in degrees
C            ALON  = Apex longitude in degrees
C            ALT   = altitude in km
C            WK    = same as entry APXMKA
C          RETURNS:
C            GDLAT = geodetic latitude in degrees
C            GDLON = geodetic longitude in degrees
C            IST   = status: 0 (okay), -1 (results are not as close as PRECISE),
C                    or > 0 (failed)
C
C          Declarations for APXA2G formal arguments:
C
C            DIMENSION WK(LWK)
C------------------------------------------------------------------------------
C          APXM2G:
C
C          Convert from Modified Apex to geographic coordinates using tables
C          currently in memory; i.e. after calling APXMKA, APXWRA, or APXRDA.
C
C            CALL APXM2G (XLATM,ALON,ALT,HR, WK, GDLAT,GDLON, IST)
C
C          INPUTS:
C            XLATM = Modified Apex latitude in degrees
C            ALON  = Modified Apex longitude (same as Apex longitude) in degrees
C            ALT   = altitude in km
C            HR   = Reference altitude, km (used only for Modified Apex coords)
C            WK    = same as entry APXMKA
C          RETURNS:
C            GDLAT = geodetic latitude in degrees
C            GDLON = geodetic longitude in degrees
C            IST   = status: 0 (okay), -1 (results are not as close as PRECISE),
C                    or > 0 (failed)
C
C          Declarations for APXM2G formal arguments:
C            DIMENSION WK(LWK)
C
C------------------------------------------------------------------------------
C          EXTERNALS and related files
C
C          An asterisk (*) indicates those required when compiling this file
C
C          * apex.f     - Apex model subroutines
C            cossza.f   - subroutine which determines the cosine of the solar
C                         zenith angle
C            ggrid.f    - subroutine to produce grids for apxntrp tables
C          * magfld.f   - subroutines defining the International Geomagnetic
C                         Reference Field (IGRF) and geographic-geocentric
C                         coordinate conversion
C            magloctm.f - magnetic local time subroutine
C            Makefile   - 'make' file for Unix systems which defines source code
C                         dependencies; i.e. instructions to make executables
C            mkglob.f   - example program to create global look-up tables file
C            subsol.f   - subroutine to determine the sub-solar point
C            xapxntrp.f - example program to exercise apxntrp
C            xglob.f    - example program to exercise apxntrp using a previously
C                         created look-up table file
C
C------------------------------------------------------------------------------
C          INSTALLATION SPECIFICS:
C
C          (1) Set IRLF, defined below.
C          (2) Initialize memory to 0 when compiling; i.e. avoid options which
C              set to indefinite or 'NaN' or change KGMA initialization below.
C          (3) To make executables it is not necessary to use the Makefile or
C              prepare an object library; e.g., without these conveniences:
C
C                f77 myprog.f apex.f apxntrp.f magfld.f
C
C              where myprog.f is a program with some of the above CALLs
C
C
C------------------------------------------------------------------------------
C          HISTORY:
C
C          Aug 1994:
C          Initial version completed by A. D. Richmond, NCAR.
C
C          Sep 1995:
C          Changes were made with the objective to allow the user to completely
C          control definition of the interpolation grid.  While doing this the
C          order of the ENTRYs in the first subroutine and arguments lists were
C          changed:  APXMKA, APXWRA, APXRDA (formerly GETARR) and the other
C          ENTRYs now include a work array (WK) which holds arrays X,Y,Z,V,
C          GPLAT,GPLON and GPALT.  Subroutine SETLIM was removed and the grid
C          creation algorithm based on NVERT originally integral to GETARR and
C          SETLIM has been extracted, but still available in ggrid.f with examples
C          xapxntrp.f and mkglob.f.  Subroutine TSTDIM has a different role, so it
C          is now named CKGP (check grid points).  MAKEXYZV was also changed to
C          accomodate explicit grid point arrays instead of computed values from
C          an index, minimum and delta.  Only one format is written now, so that
C          it is possible to concatenate files for different epochs.  This
C          required changing delta time logic from fixed 5 yr epochs to computed
C          differences which may vary.  Also, the ties to DGRF/IGRF dates have
C          been removed.  R. Barnes.
C
C          Sep 1996:
C          Corrected bug in APXQ2G longitude iteration.  The latitude iteration
C          is now constrained to not step beyond the (partial global)
C          interpolation grid.  NITER was increased from 9 to 14 and code to
C          determine cos(lambda') (CLP) was revised by Art Richmond to reduce
C          truncation errors near the geographic poles.
C
C          Sep 1997:
C          Corrected comments, definition of COLAT.
C
C          Dec 1998:
C          Change GLON input to try +/-360 values before rejecting when out of
C          defined grid (GDLON) range; affects INTRP.
C
C          Feb-Mar 1999:
C          (1) Corrected a typo (bad variable name) in diagnostic print in
C          INTRP:  GLO -> GLON.  This error was probably introduced Dec 98, so
C          no-one had access to the bad copy and, therefore, no announcement is
C          needed. (2) Also modified APXMALL and APXQ2G:  When very close to a
C          geographic pole, gradients are recalculated after stepping away; in
C          this situation, the latitude input to INTRP was changed from GLAT to
C          GLATX.  This is affects gradients when within 0.1 degrees of a pole
C          (defined as the larger of GLATLIM, 89.9 deg, or the second largest
C          grid point). (3) Changed MAKEXYZV to make X,Y,Z, V constant for all
C          longitudes (at each alt) when poleward of POLA; POLA was added to
C          /APXCON/. (4) Replaced definition of DYLON in APXQ2G. (5) Reduced
C          NITER to 10 from 14.  Changes 3-5 fix a problem where APXQ2G
C          calculations were failing to satisify PRECISE within NITER iterations
C          at the pole. (6) Replace XYZ2APX with revised trigonometry to avoid a
C          truncation problem at the magnetic equator.  Most changes were
C          devised by Art Richmond and installed by Roy Barnes.
C
C          May 2000:
C          Relaxed acceptable input longitude check in CKGP to be +/- 270
C          degrees.  Also updated IGRF coefficients in magfld.f; now DGRF epochs
C          are 1900, 1905, ... 1995 and IGRF is for 2000-2005.
C
C          May 2004:
C          (1) Change definition of earth's equatorial radius (REQ) from the
C          IAU-1966 spheroid (6378.160 km) to the WGS-1984 spheroid (6378.137
C          km) in accordance with recent IGRF release; see file
C          $APXROOT/docs/igrf.2004.spheroid.desc. (2) Add APXA2G and APXM2G.
C
C          Please direct questions to Roy Barnes, NCAR
C          email: bozo@ucar.edu
C          phone: 303-497-1572
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
C                    transform from Quasi-Dipole to geodetic coordinates.
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
      CHARACTER FILNAM*(*)

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
C            RE, REQ    = 6371.2, 6378.137 (Mean and equatorial Earth radius)
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
      REQ  = 6378.137
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
	OPEN (IUN,FILE=FILNAM,ACCESS='direct',RECL=LDR,STATUS='unknown',
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

      ENTRY APXMALL (GLATM,GLONM,ALTM,HR, WK,                           !Inputs
     +              B,BHAT,BMAG,SI,                                     !Mag Fld
     +              ALON,                                               !Apx Lon
     +              XLATM,VMP,WM,D,BE3,SIM,D1,D2,D3,E1,E2,E3,           !Mod Apx
     +              XLATQD,F,F1,F2 , IST)                               !Qsi-Dpl
C          940822 A. D. Richmond, Sep 95 R. Barnes

C          Test to see if WK has been initialized
      CALNM = 'APXMALL'
      LCN   = 7
      IF (KGMA .LT. 1) GO TO 9300

C          Alias input variables to avoid same name in multiple entries
      GLAT = GLATM
      GLON = GLONM  ! also needed in case INTRP adjusts it by +/-360
      ALT  = ALTM

      CALL INTRP (GLAT,GLON,ALT ,WK(LBX),WK(LBY),WK(LBZ),WK(LBV),
     +           NLA,NLO,NAL,    WK(LLA),WK(LLO),WK(LAL),
     +           FX,FY,FZ,FV,
     +           DFXDTH,DFYDTH,DFZDTH,DFVDTH,DFXDLN,DFYDLN,DFZDLN,
     +           DFVDLN,DFXDH,DFYDH,DFZDH,DFVDH, CALNM(1:LCN),IST)

      IF (IST .NE. 0) THEN
	CALL SETMISS (XMISS, XLATM,ALON,VMP,B,BMAG,BE3,SIM,SI,F,D,WM,
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
     +            BMAG,SIM,SI,F,D,WM,BHAT,D1,D2,D3,E1,E2,E3,F1,F2)
      BE3 = BMAG/D

      IST = 0
      RETURN


C*******************************************************************************

      ENTRY APXALL (GLATA,GLONA,ALTA, WK, A,ALAT,ALON, IST)
C          940802 A. D. Richmond, Sep 95 R. Barnes

C          Test to see if WK has been initialized
      CALNM = 'APXALL'
      LCN   = 6
      IF (KGMA .LT. 1) GO TO 9300

C          Alias input variables to avoid same name in multiple entries
      GLAT = GLATA
      GLON = GLONA   ! also needed in case INTRPSC adjusts it by +/-360
      ALT  = ALTA

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

      ENTRY APXQ2G (QDLAT,QDLON,ALTQ, WK, GDLAT,GDLON, IST)
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

C          Alias input variables to avoid same name in multiple entries
      QDLA = QDLAT
      QDLO = QDLON
      QALT = ALTQ

C          Determine quasi-cartesian coordinates on a unit sphere of the
C          desired magnetic lat,lon in Quasi-Dipole coordinates.
  700 X0 = COS (QDLA*DTOR) * COS (QDLO*DTOR)
      Y0 = COS (QDLA*DTOR) * SIN (QDLO*DTOR)
      Z0 = SIN (QDLA*DTOR)

C          Initial guess:  use centered dipole, convert to geocentric coords
      CALL GM2GC (QDLA,QDLO,YLAT,YLON)

C          Iterate until (angular distance)**2 (units: radians) is within
C          PRECISE of location (QDLA,QDLO) on a unit sphere.
      NITER = 10
      DO 710 ITER=1,NITER
      CALL INTRP (YLAT,YLON,QALT, WK(LBX),WK(LBY),WK(LBZ),WK(LBV),
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
	CALL INTRP (GLATX,YLON,QALT, WK(LBX),WK(LBY),WK(LBZ),WK(LBV),
     +             NLA,NLO,NAL,     WK(LLA),WK(LLO),WK(LAL),
     +             FXDUM,FYDUM,FZDUM,FVDUM,
     +             DMXDTH,DMYDTH,DMZDTH,DMVDTH,DFXDLN,DFYDLN,DFZDLN,
     +             DFVDLN,DMXDH,DMYDH,DMZDH,DMVDH, CALNM(1:LCN),IST)
	CALL ADPL (GLATX,YLON,CTH,STH,FXDUM,FYDUM,FZDUM,FVDUM,
     +          DMXDTH,DMYDTH,DMZDTH,DMVDTH,DFXDLN,DFYDLN,DFZDLN,DFVDLN)
      ENDIF

C          At this point, FX,FY,FZ are approximate quasi-cartesian
C          coordinates on a unit sphere for the Quasi-Dipole coordinates
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

      IF (DIST2 .LE. PRECISE) GO TO 720
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

      YLON  = YLON + DYLON
      IF (YLON .GT.  WK(LLO+NLO-1)) YLON = YLON - 360.
      IF (YLON .LT.  WK(LLO)      ) YLON = YLON + 360.
  710 CONTINUE

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
      GO TO 730

  720 IST = 0

  730 GDLAT = YLAT
      GDLON = YLON
      IF (CALNM .EQ. 'APXQ2G') GO TO 800
      IF (CALNM .EQ. 'APXM2G') GO TO 900

      RETURN

C*******************************************************************************

      ENTRY APXA2G (ALAA2,ALOA2,ALTA2, WK, GDLAT,GDLON, IST)
C          Given Apex coordinates (ALAT, ALON, ALT) determine geographic
C          coordinates (GDLAT, GDLON).  May 2004, R.Barnes, NCAR

C          Test to see if WK has been initialized
      CALNM = 'APXA2G'
      LCN   = 6
      IF (KGMA .LT. 1) GO TO 9300

      HA   = ( (1./COS(ALAA2*DTOR))**2 - 1.)*REQ                ! apex altitude from equations 3.1 and 3.2 (Richmond, 1995)
      QDLA = SIGN (RTOD*ACOS( SQRT((RE+ALTA2)/(RE+HA))), ALTA2) ! Quasi-Dipole lat. from equation 6.1 (Richmond, 1995)
      QDLO = ALOA2
      QALT = ALTA2
      GO TO 700

  800 CONTINUE
      RETURN


C*******************************************************************************

      ENTRY APXM2G (ALAM2,ALOM2,ALTM2,HR, WK, GDLAT,GDLON, IST)
C          Given Modified Apex coordinates (XLATM, XLON, ALT,HR) determine geographic
C          coordinates (GDLAT, GDLON). May 2004, R.Barnes, NCAR

C          Test to see if WK has been initialized
      CALNM = 'APXM2G'
      LCN   = 6
      IF (KGMA .LT. 1) GO TO 9300

      QDLA = ACOS ( SQRT((RE+ALTM2)/(RE+HR)) * COS(ALAM2*DTOR) )*RTOD ! Quasi-Dipole lat. equation 6.2 (Richmond, 1995)
      QDLO = ALOM2
      QALT = ALTM2
      GO TO 700

  900 CONTINUE
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
     +    is',F6.2,' years')

 9190 WRITE (MSGU,'(A,'': Dimensions (lat,lon,alt) read from "'',A,''" f
     +or EPOCH '',F7.2,/,''        ('',I5,'','',I5,'','',I5,'') do not m
     +atch ('',I5,'','',I5,'','',I5,'') from the'',/,''        first EPO
     +CH read ('',F7.2,'')'')') CALNM(1:LCN), FILNAM(1:LFN), TI ,
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
C            QDLAT  = Quasi-Dipole latitude, degrees
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
C       re, req    = 6371.2, 6378.137
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
      SUBROUTINE GGRID (NVERT,GLAMN,GLAMX,GLOMN,GLOMX,ALTMN,ALTMX,
     +                 GPLAT,GPLON,GPALT,NLAT,NLON,NALT)

C          Define grid points suitable for preparing interpolation tables for
C          Apex coordinates, Modified Apex Coordinates and Quasi-Dipole
C          coordinates.
C
C          INPUTS:
C            NVERT = spatial resolution parameter must be at least 2 but 30 or
C                    or more is recommended.  Interpolation accuracy increases
C                    with larger NVERT until about 100 when points near the
C                    poles become too close for accurate east-west gradient
C                    determination affecting B(1), G, H, W, Bhat(1), D1(1),
C                    D2(1), D3, E1, E2, E3, F1(2), and F2(2)); this effects all
C                    latitudes when NVERT is 1000 or more.
C            GLAMN = minimum latitude (degrees) must be in the range -90 to 90.
C            GLAMX = maximum latitude (degrees) must be in the range GLAMN to 90.
C            GLOMN = minimum longitude (degrees) must be in the range -180 to 180.
C            GLOMX = maximum longitude (degrees) must be in the range GLOMN to
C                    the smaller of 540 or GLOMN+360.
C            ALTMN = minimum altitude (km) must be at least 0.
C            ALTMX = maximum altitude (km) must be at least ALTMN.
C          RETURNS:
C            GPLAT = Grid point latitudes (degrees north).
C            GPLON = Grid point longitudes (degrees east).
C            GPALT = Grid point altitudes (km above earth).
C            NLAT  = Number of assigned values in GPLAT.
C            NLON  = Number of assigned values in GPLON.
C            NALT  = Number of assigned values in GPALT.
C
C          This is the suggested method for defining grid points used to
C          prepare interpolation tables (by entry APXWRA of apxntrp.f).
C
C          When determining geomagnetic coordinates, the dipole component is
C          separated so it may be computed analytically (to improve accuracy)
C          and the non-dipole component is interpolated from the tables.  As
C          altitude increases, the non-dipole component becomes spatially
C          smoother and relatively less important, so vertical spacing can
C          increase with altitude; hence, the quantity
C
C            Re/(Re+height)
C
C          is divided into NVERT segments, where Re is the mean earth radius
C          (6371.2 km).  The latitudinal grid spacing is chosen to be
C          approximately the same at the Earth's surface as the vertical
C          distance between the two lowest grid points in altitude or
C          180/(3*NVERT).  The longitude spacing is chosen to be roughly
C          equivalent to low and middle latitude spacing or 360/(5*NVERT);
C          however longitude spacing decreases when approaching the poles.
C          Grid values returned are the minimal subset encompassing the
C          specified latitude, longitude and altitude extremes.
C
C          HISTORY:
C          Aug 1995: Adapted from the initial version of apxntrp.f written by
C                    A. D. Richmond.  This grid definition code was separated from
C                    apxntrp.f to allow arbitrary grid definition.  R. Barnes, NCAR.
C          May 2004: Clean up comments; change print diagnostics to unit 0 (stderr);
C                    alias GLOMX to GLONMX s.t. it can be changed w/o affecting the
C                    input argument.

C          Formal argument declarations
      DIMENSION GPLAT(*), GPLON(*), GPALT(*)

C          Local declarations:
      DOUBLE PRECISION DLON , DLAT , DIHT , DNV , DRE
      PARAMETER (RE = 6371.2 ,  DRE = 6371.2D0)
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

      LAMN =  MAX0 (INT((DBLE(GLAMN) +90.D0 )/DLAT)     , 0)
      LAMX =  MIN0 (INT((DBLE(GLAMX) +90.D0 )/DLAT+1.D0), 3*NVERT)

      LOMN  = MAX0 (INT((DBLE(GLOMN) +180.D0)/DLON)  , 0)
      GLOMXL= AMIN1 (GLOMX,GLOMN+360.)
      LOMX  = MIN0 (INT((DBLE(GLOMXL)+180.D0)/DLON+1.D0),10*NVERT)

      X = RE/(RE+ALTMX)/DIHT - 1.E-5
      IHTMN = AMAX1 (X,1.)
      IHTMN = MIN0 (IHTMN,NVERT-1)
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
      IF (GPLON(NLON-1) .GE. GLOMXL) NLON = NLON - 1
      GPALT(1) = AMAX1 (GPALT(1),0.)

      RETURN

 9100 WRITE (0,'(''GGRID:  GLAMX < GLAMN'')')
      CALL EXIT (1)
 9200 WRITE (0,'(''GGRID:  GLOMX < GLOMN'')')
      CALL EXIT (1)
 9300 WRITE (0,'(''GGRID:  ALTMX < ALTMN'')')
      CALL EXIT (1)
 9400 WRITE (0,'(''GGRID:  |GLAMN| > 90.'')')
      CALL EXIT (1)
 9500 WRITE (0,'(''GGRID:  |GLOMN| > 180.'')')
      CALL EXIT (1)
 9600 WRITE (0,'(''GGRID:  ALTMN < 0.'')')
      CALL EXIT (1)
      END
      SUBROUTINE APEX (DATE,DLAT,DLON,ALT,
     +                 A,ALAT,ALON,BMAG,XMAG,YMAG,ZMAG,V)
C          Calculate apex radius, latitude, longitude; and magnetic field and
C          scalar magnetic potential.
C
C          INPUTS:
C            DATE = Year and fraction (1990.0 = 1990 January 1, 0 UT)
C            DLAT = Geodetic latitude in degrees
C            DLON = Geodetic longitude in degrees
C            ALT = Altitude in km
C
C          RETURNS:
C            A    = (Apex height + REQ)/REQ, where REQ = equatorial Earth radius.
C                   A is analogous to the L value in invariant coordinates.
C            ALAT = Apex latitude in degrees (negative in S. magnetic hemisphere)
C            ALON = Apex longitude (geomagnetic longitude of apex) in degrees
C            BMAG = geomagnetic field magnitude (nT)
C            XMAG = geomagnetic field component (nT): north
C            YMAG = geomagnetic field component (nT): east
C            ZMAG = geomagnetic field component (nT): downward
C            V    = geomagnetic potential (T-m)
C
C          COMMON BLOCKS:
C            COMMON /DIPOLE/ COLAT,ELON,VP,CTP,STP
C
C          DIPOLE has IGRF variables obtained from routines in magfld.f:
C            COLAT = Geocentric colatitude of geomagnetic dipole north pole (deg)
C            ELON  = East longitude of geomagnetic dipole north pole (deg)
C            VP    = Magnitude (T-m) of dipole component of magnetic potential at
C                    geomagnetic pole and geocentric radius of 6371.2 km
C            CTP   = cosine of COLAT
C            STP   = sine   of COLAT
C
C------------------------------------------------------------------------------
C          HISTORY:
C          Aug 1994: First version completed on the 22nd by A.D. Richmond.
C          May 1999: Revise DS calculation in LINAPX to avoid divide by zero.
C          Apr 2004: - Change definition of earth's equatorial radius (REQ)
C                      from the IAU-1966 spheroid (6378.160 km) to the WGS-1984
C                      spheroid (6378.137 km); see description in
C                      '$APXROOT/docs/igrf.2004.spheroid.desc'.
C                    - Revise comments toward a consistent format so they are
C                      easy to read.
C                    - Replace computed GO TO in ITRACE with IF blocks.
C                    - Refine FNDAPX to insure |Bdown/Btot| < 1.E-6 at apex

      PARAMETER (RE = 6371.2, DTOR = .01745329251994330)
      COMMON /DIPOLE/ COLAT,ELON,VP,CTP,STP

      CALL COFRM (DATE)
      CALL DYPOL (CLATP,POLON,VPOL)
      COLAT = CLATP
      CTP   = COS(CLATP*DTOR)
      STP   = SQRT(1. - CTP*CTP)
      ELON  = POLON
      VP    = VPOL
      CALL LINAPX (DLAT,DLON,ALT, A,ALAT,ALON,XMAG,YMAG,ZMAG,BMAG)
      XMAG = XMAG*1.E5                 ! convert from gauss to nT
      YMAG = YMAG*1.E5
      ZMAG = ZMAG*1.E5
      BMAG = BMAG*1.E5
      CALL GD2CART (DLAT,DLON,ALT,X,Y,Z)
      CALL FELDG (3, X/RE,Y/RE,Z/RE, BX,BY,BZ,V)
      RETURN
      END

      SUBROUTINE LINAPX (GDLAT,GLON,ALT, A,ALAT,ALON,XMAG,YMAG,ZMAG,F)

C          Transform geographic coordinates to Apex coordinates.
C
C          INPUTS:
C            GDLAT = Latitude  (degrees, positive northward)
C            GLON  = Longitude (degrees, positive eastward)
C            ALT   = Height of starting point (km above mean sea level)
C
C          OUTPUTS:
C            A     = (Apex height + REQ)/REQ, where REQ = equatorial Earth radius.
C                    A is analogous to the L value in invariant coordinates.
C            ALAT  = Apex Lat. (deg)
C            ALON  = Apex Lon. (deg)
C            XMAG  = Geomagnetic field component (gauss): north
C            YMAG  = Geomagnetic field component (gauss): east
C            ZMAG  = Geomagnetic field component (gauss): down
C            F     = Geomagnetic field magnitude (gauss)
C
C          Trace the geomagnetic field line from the given location to find the
C          apex of the field line.  Before starting iterations to trace along
C          the field line: (1) Establish a step size (DS, arc length in km)
C          based on the geomagnetic dipole latitude; (2) determine the step
C          direction from the sign of the vertical component of the geomagnetic
C          field; and (3) convert to geocentric cartesian coordinates.  Each
C          iteration increments a step count (NSTP) and calls ITRACE to move
C          along the the field line until reaching the iteration count limit
C          (MAXS) or passing the apex (IAPX=2) and then calling FNDAPX to
C          determine the apex location from the last three step locations
C          (YAPX); however, if reaching the iteration limit, apex coordinates
C          are calculated by DIPAPX which assumes a simplified dipole field.
C
C          COMMON BLOCKS:
C            COMMON /APXIN/   YAPX(3,3)
C            COMMON /DIPOLE/  COLAT,ELON,VP,CTP,STP
C            COMMON /FLDCOMD/ BX, BY, BZ, BB
C            COMMON /ITRA/    NSTP, Y(3), YOLD(3), SGN, DS
C
C          APXIN has step locations determined in ITRACE:
C            YAPX  = Matrix of cartesian coordinates (loaded columnwise) of the
C                    three points about the apex.  Set in subroutine ITRACE.
C                                                                               
C          DIPOLE has IGRF variables obtained from routines in magfld.f:
C            COLAT = Geocentric colatitude of geomagnetic dipole north pole (deg)
C            ELON  = East longitude of geomagnetic dipole north pole (deg)
C            VP    = Magnitude (T-m) of dipole component of magnetic potential at
C                    geomagnetic pole and geocentric radius of 6371.2 km
C            CTP   = cosine of COLAT
C            STP   = sine   of COLAT
C                                                                               
C          FLDCOMD has geomagnetic field at current trace point:
C            BX    = X component (Gauss)
C            BY    = Y component (Gauss)
C            BZ    = Z component (Gauss)
C            BB    = Magnitude   (Gauss)
C
C          ITRA has field line tracing variables determined in LINAPX:
C            NSTP  = Step count.
C            Y     = Array containing current tracing point cartesian coordinates.
C            YOLD  = Array containing previous tracing point cartesian coordinates.
C            SGN   = Determines direction of trace.
C            DS    = Step size (arc length in km).
C                                                                               
C          REFERENCES:
C            Stassinopoulos E. G. , Mead Gilbert D., X-841-72-17 (1971) GSFC,
C            Greenbelt, Maryland
C                                                                               
C          EXTERNALS:
C            GD2CART = Convert geodetic to geocentric cartesian coordinates (in magfld.f)
C            CONVRT  = Convert geodetic to geocentric cylindrical or geocentric spherical
C                      and back (in magfld.f).
C            FELDG   = Obtain IGRF magnetic field components (in magfld.f).
C            ITRACE  = Follow a geomagnetic field line
C            DIPAPX  = Compute apex coordinates assuming a geomagnetic dipole field
C            FNDAPX  = Compute apex coordinates from the last three traced field line points
C
C------------------------------------------------------------------------------
C          HISTORY:
C          Oct 1973: Initial version completed on the 29th by Wally Clark, NOAA
C                    ERL Lab.
C          Feb 1988: Revised on the 1st by Harsh Anand Passi, NCAR.
C          Aug 1994: Revision by A. D. Richmond, NCAR.

      PARAMETER (MAXS = 200, RTOD = 57.2957795130823,   RE =6371.2,
     +                       DTOR = .01745329251994330, REQ=6378.137)
      COMMON /FLDCOMD/ BX, BY, BZ, BB
      COMMON /APXIN/   YAPX(3,3)
      COMMON /DIPOLE/  COLAT,ELON,VP,CTP,STP
      COMMON /ITRA/    NSTP, Y(3), YP(3), SGN, DS

C          Set step size based on the geomagnetic dipole latitude of the starting point
      CALL CONVRT (2,GDLAT,ALT,GCLAT,R)
      SINGML = CTP*SIN(GCLAT*DTOR) + STP*COS(GCLAT*DTOR)*
     +                                             COS((GLON-ELON)*DTOR)
C          May 1999: avoid possible divide by zero (when SINGML = 1.): the old version
C          limited DS to its value at 60 deg GM latitude with: DS = .06*R/(1.-SINGML*SINGML) - 370.
C                                                              IF (DS .GT. 1186.) DS = 1186.
      CGML2 = AMAX1 (0.25,1.-SINGML*SINGML)
      DS = .06*R/CGML2 - 370.

C          Initialize YAPX array
      DO 4 J=1,3
      DO 4 I=1,3
    4 YAPX(I,J) = 0.

C          Convert from geodetic to earth centered cartesian coordinates
      CALL GD2CART (GDLAT,GLON,ALT,Y(1),Y(2),Y(3))
      NSTP = 0

C          Get magnetic field components to determine the direction for
C          tracing the field line
      CALL FELDG (1,GDLAT,GLON,ALT,XMAG,YMAG,ZMAG,F)
      SGN = SIGN (1.,-ZMAG)

C          Use cartesian coordinates to get magnetic field components
C          (from which gradients steer the tracing)
   10 CALL FELDG (2, Y(1)/RE,Y(2)/RE,Y(3)/RE, BX,BY,BZ,BB)
      NSTP = NSTP + 1

      IF (NSTP .LT. MAXS) THEN
	CALL ITRACE (IAPX)                               ! trace along field line
	IF (IAPX .EQ. 1) GO TO 10
	CALL FNDAPX (ALT,ZMAG,A,ALAT,ALON)               ! (IAPX=2) => passed max radius; find its coordinates
      ELSE
	RHO = SQRT (Y(1)*Y(1) + Y(2)*Y(2))               ! too many steps; get apex from dipole approximation
	CALL CONVRT (3,XLAT,HT,RHO,Y(3))
	XLON = RTOD*ATAN2 (Y(2),Y(1))
	CALL FELDG (1,XLAT,XLON,HT,BNRTH,BEAST,BDOWN,BABS)
	CALL DIPAPX  (XLAT,XLON,HT,BNRTH,BEAST,BDOWN,A,ALON)
	ALAT = -SGN*RTOD*ACOS (SQRT(1./A))
      ENDIF

      RETURN
      END

      SUBROUTINE ITRACE (IAPX)

C          Follow a geomagnetic field line until passing its apex
C
C          INPUTS:
C            (all are in common blocks)
C          OUTPUTS:
C            IAPX = 2 (when apex passed) or 1 (not)
C                                                                               
C          This uses the 4-point Adams formula after initialization.
C          First 7 iterations advance point by 3 steps.
C                                                                               
C          COMMON BLOCKS:
C            COMMON /APXIN/   YAPX(3,3)
C            COMMON /FLDCOMD/ BX, BY, BZ, BB
C            COMMON /ITRA/    NSTP, Y(3), YOLD(3), SGN, DS
C
C          APXIN has step locations determined in ITRACE:
C            YAPX  = Matrix of cartesian coordinates (loaded columnwise) of the
C                    three points about the apex.  Set in subroutine ITRACE.
C                                                                               
C          FLDCOMD has geomagnetic field at current trace point:
C            BX    = X component (Gauss)
C            BY    = Y component (Gauss)
C            BZ    = Z component (Gauss)
C            BB    = Magnitude   (Gauss)
C                                                                               
C          ITRA has field line tracing variables determined in LINAPX:
C            NSTP  = Step count.
C            Y     = Array containing current tracing point cartesian coordinates.
C            YOLD  = Array containing previous tracing point cartesian coordinates.
C            SGN   = Determines direction of trace.
C            DS    = Step size (arc length in km).
C
C          REFERENCES:
C            Stassinopoulos E. G. , Mead Gilbert D., X-841-72-17 (1971) GSFC,
C            Greenbelt, Maryland
C------------------------------------------------------------------------------
C          HISTORY:
C          Oct 1973: Initial version completed on the 29th by W. Clark, NOAA ERL
C                    Laboratory.
C          Feb 1988: Revised by H. Passi, NCAR.
C          Apr 2004: Replace computed GO TO with IF blocks because some compilers
C                    are threatening to remove this old feature
C

      COMMON /APXIN/   YAPX(3,3)
      COMMON /FLDCOMD/ BX, BY, BZ, BB
      COMMON /ITRA/    NSTP, Y(3), YOLD(3), SGN, DS
      DIMENSION YP(3,4)
      SAVE

C          Statement function
      RDUS(D,E,F) = SQRT (D**2 + E**2 + F**2)

      IAPX = 1

C          Cartesian component magnetic field (partial) derivitives steer the trace
      YP(1,4) = SGN*BX/BB
      YP(2,4) = SGN*BY/BB
      YP(3,4) = SGN*BZ/BB

      IF (NSTP .LE. 7) THEN
	DO 10 I=1,3
	IF (NSTP .EQ. 1) THEN
	  D2        = DS/2.
	  D6        = DS/6.
	  D12       = DS/12.
	  D24       = DS/24.
	  YP(I,1)   = YP(I,4)
	  YOLD(I)   = Y(I)
	  YAPX(I,1) = Y(I)
	  Y(I)      = YOLD(I) + DS*YP(I,1)

	ELSE IF (NSTP .EQ. 2) THEN
	  YP(I,2) = YP(I,4)
	  Y(I)    = YOLD(I) + D2*(YP(I,2)+YP(I,1))

	ELSE IF (NSTP .EQ. 3) THEN
	  Y(I) = YOLD(I) + D6*(2.*YP(I,4)+YP(I,2)+3.*YP(I,1))

	ELSE IF (NSTP .EQ. 4) THEN
	  YP(I,2)   = YP(I,4)
	  YAPX(I,2) = Y(I)
	  YOLD(I)   = Y(I)
	  Y(I)      = YOLD(I) + D2*(3.*YP(I,2)-YP(I,1))

	ELSE IF (NSTP .EQ. 5) THEN
	  Y(I) = YOLD(I) + D12*(5.*YP(I,4)+8.*YP(I,2)-YP(I,1))

	ELSE IF (NSTP .EQ. 6) THEN
	  YP(I,3)   = YP(I,4)
	  YOLD(I)   = Y(I)
	  YAPX(I,3) = Y(I)
	  Y(I)      = YOLD(I) + D12*(23.*YP(I,3)-16.*YP(I,2)+5.*YP(I,1))

	ELSE IF (NSTP .EQ. 7) THEN
	  YAPX(I,1) = YAPX(I, 2)
	  YAPX(I,2) = YAPX(I, 3)
	  Y(I)      = YOLD(I) + D24*(9.*YP(I,4) + 19.*YP(I,3) -
     +                               5.*YP(I,2) +     YP(I,1))
	  YAPX(I,3) = Y(I)
	ENDIF
   10   CONTINUE
	IF (NSTP .EQ. 6 .OR. NSTP .EQ. 7) THEN        ! signal if apex passed
	  RC = RDUS (YAPX(1,3), YAPX(2,3), YAPX(3,3))
	  RP = RDUS (YAPX(1,2), YAPX(2,2), YAPX(3,2))
	  IF (RC .LT. RP) IAPX = 2
	ENDIF

      ELSE                 ! NSTP > 7

	DO 30 I=1,3
	YAPX(I,1) = YAPX(I,2)
	YAPX(I,2) = Y(I)
	YOLD(I)   = Y(I)
	Y(I)      = YOLD(I) + D24*(55.*YP(I,4) - 59.*YP(I,3) +
     +                             37.*YP(I,2) -  9.*YP(I,1))
	YAPX(I,3) = Y(I)

	DO 20 J=1,3
   20   YP(I,J) = YP(I,J+1)
   30   CONTINUE
	RC = RDUS (   Y(1),    Y(2),    Y(3))
	RP = RDUS (YOLD(1), YOLD(2), YOLD(3))
	IF (RC .LT. RP) IAPX = 2
      ENDIF

      RETURN
      END

       SUBROUTINE FNDAPX (ALT,ZMAG,A,ALAT,ALON)

C          Find apex coordinates once tracing (in subroutine ITRACE) has
C          signalled that the apex has been passed.
C          INPUTS:
C            ALT  = Altitude of starting point
C            ZMAG = Downward component of geomagnetic field at starting point
C          OUTPUT
C            A    = Apex radius, defined as (Apex height + Req)/Req, where
C                   Req = equatorial Earth radius.
C                   A is analogous to the L value in invariant coordinates.
C            ALAT = Apex Lat. (deg)
C            ALON = Apex Lon. (deg)
C
C          COMMON BLOCKS:
C            COMMON /APXIN/  YAPX(3,3)
C            COMMON /DIPOLE/ COLAT,ELON,VP,CTP,STP
C
C          APXIN has step locations determined in ITRACE:
C            YAPX  = Matrix of cartesian coordinates (loaded columnwise) of the
C                    three points about the apex.  Set in subroutine ITRACE.
C                                                                               
C          DIPOLE has IGRF variables obtained from routines in magfld.f:
C            COLAT = Geocentric colatitude of geomagnetic dipole north pole (deg)
C            ELON  = East longitude of geomagnetic dipole north pole (deg)
C            VP    = Magnitude (T-m) of dipole component of magnetic potential at
C                    geomagnetic pole and geocentric radius of 6371.2 km
C            CTP   = cosine of COLAT
C            STP   = sine   of COLAT
C                                                                               
C          EXTERNALS:
C            FINT = Second degree interpolation routine
C------------------------------------------------------------------------------
C          HISTORY:
C          Oct 1973: Initial version completed on the 23rd by Clark, W., NOAA
C                    Boulder.
C          Aug 1994: Revision on the 3rd by A.D. Richmond, NCAR
C          Apr 2004: Repair problem noted by Dan Weimer where the apex location
C                    produced by FINT may still have a non-zero vertical magnetic
C                    field component.

      PARAMETER (RTOD = 57.2957795130823,
     +           DTOR = .01745329251994330, REQ=6378.137)
      COMMON /APXIN/  YAPX(3,3)
      COMMON /DIPOLE/ COLAT,ELON,VP,CTP,STP
      DIMENSION BD(3), Y(3)

C          Get geodetic height and vertical (downward) component of the magnetic
C          field at last three points found by ITRACE
      DO 10 I=1,3
      RHO  = SQRT (YAPX(1,I)**2 + YAPX(2,I)**2)
      CALL CONVRT (3,GDLT,HT, RHO,YAPX(3,I))
      GDLN = RTOD*ATAN2 (YAPX(2,I),YAPX(1,I))
   10 CALL FELDG (1,GDLT,GDLN,HT, BN,BE,BD(I),BMAG)

C          Interpolate to where Bdown=0 to find cartesian coordinates at dip equator
      NITR = 0
   20 Y(1) = FINT (BD(1),BD(2),BD(3),YAPX(1,1),YAPX(1,2),YAPX(1,3), 0.)
      Y(2) = FINT (BD(1),BD(2),BD(3),YAPX(2,1),YAPX(2,2),YAPX(2,3), 0.)
      Y(3) = FINT (BD(1),BD(2),BD(3),YAPX(3,1),YAPX(3,2),YAPX(3,3), 0.)

C          Insure negligible Bdown or
C
C            |Bdown/Btot| < 2.E-6
C
C          For instance, Bdown must be less than 0.1 nT at low altitudes where
C          Btot ~ 50000 nT.  This ratio can be exceeded when interpolation is
C          not accurate; i.e., when the middle of the three points interpolated
C          is too far from the dip equator.  The three points were initially
C          defined with equal spacing by ITRACE, so replacing point 2 with the
C          most recently fit location will reduce the interpolation span.
      RHO  = SQRT (Y(1)**2 + Y(2)**2)
      GDLN = RTOD*ATAN2 (Y(2),Y(1))
      CALL CONVRT (3,GDLT,HTA, RHO,Y(3))
      CALL FELDG (1,GDLT,GDLN,HTA, BNA,BEA,BDA,BA)
      ABDOB = ABS(BDA/BA)

      IF (ABDOB .GT. 2.E-6) THEN
	IF (NITR .LT. 4) THEN        ! 4 was chosen because tests rarely required 2 iterations
	  NITR      = NITR + 1
	  YAPX(1,2) = Y(1)
	  YAPX(2,2) = Y(2)
	  YAPX(3,2) = Y(3)
	  BD(2)     = BDA
	  GO TO 20
	ELSE
	  WRITE (0,'(''APEX: Imprecise fit of apex: |Bdown/B| ='',1PE7.1
     +    )') ABDOB
	ENDIF
      ENDIF

C          Ensure altitude of the Apex is at least the initial altitude when
C          defining the Apex radius then use it to define the Apex latitude whose
C          hemisphere (sign) is inferred from the sign of the dip angle at the
C          starting point
      A = (REQ + AMAX1(ALT,HTA)) / REQ
      IF (A .LT. 1.) THEN
	WRITE (0,'(''APEX: A can not be less than 1; A, REQ, HTA: '',1P3
     +E15.7)') A,REQ,HTA
	CALL EXIT (1)
      ENDIF
      RASQ = ACOS (SQRT(1./A))*RTOD
      ALAT = SIGN (RASQ,ZMAG)

C          ALON is the dipole longitude of the apex and is defined using
C          spherical coordinates where
C            GP   = geographic pole.
C            GM   = geomagnetic pole (colatitude COLAT, east longitude ELON).
C            XLON = longitude of apex.
C            TE   = colatitude of apex.
C            ANG  = longitude angle from GM to apex.
C            TP   = colatitude of GM.
C            TF   = arc length between GM and apex.
C            PA   = ALON be geomagnetic longitude, i.e., Pi minus angle measured
C                   counterclockwise from arc GM-apex to arc GM-GP.
C          then, spherical-trigonometry formulas for the functions of the angles
C          are as shown below.  Notation uses C=cos, S=sin and STFCPA = sin(TF) * cos(PA),
C                                                              STFSPA = sin(TF) * sin(PA)
      XLON = ATAN2 (Y(2),Y(1))
      ANG  = XLON-ELON*DTOR
      CANG = COS (ANG)
      SANG = SIN (ANG)
      R    = SQRT (Y(1)**2 + Y(2)**2 + Y(3)**2)
      CTE  = Y(3)/R
      STE  = SQRT (1.-CTE*CTE)
      STFCPA = STE*CTP*CANG - CTE*STP
      STFSPA = SANG*STE
      ALON = ATAN2 (STFSPA,STFCPA)*RTOD
      RETURN
      END                                                                      

      SUBROUTINE DIPAPX (GDLAT,GDLON,ALT,BNORTH,BEAST,BDOWN, A,ALON)

C          Compute A, ALON from local magnetic field using dipole and spherical
C          approximation.
C
C          INPUTS:
C            GDLAT  = geodetic latitude, degrees
C            GDLON  = geodetic longitude, degrees
C            ALT    = altitude, km
C            BNORTH = geodetic northward magnetic field component (any units)
C            BEAST  = eastward magnetic field component
C            BDOWN  = geodetic downward magnetic field component
C          OUTPUTS:
C            A      = apex radius, 1 + h_A/R_eq
C            ALON   = apex longitude, degrees
C
C          Use spherical coordinates and define:
C            GP    = geographic pole.
C            GM    = geomagnetic pole (colatitude COLAT, east longitude ELON).
C            G     = point at GDLAT,GDLON.
C            E     = point on sphere below apex of dipolar field line passing
C                    through G.
C            TD    = dipole colatitude of point G, found by applying dipole
C                    formula for dip angle to actual dip angle.
C            B     = Pi plus local declination angle.  B is in the direction
C                    from G to E.
C            TG    = colatitude of G.
C            ANG   = longitude angle from GM to G.
C            TE    = colatitude of E.
C            TP    = colatitude of GM.
C            A     = longitude angle from G to E.
C            APANG = A + ANG
C            PA    = geomagnetic longitude, i.e., Pi minus angle measured
C                    counterclockwise from arc GM-E to arc GM-GP.
C            TF    = arc length between GM and E.
C          Then, using notation C=cos, S=sin, COT=cot, spherical-trigonometry
C          formulas for the functions of the angles are as shown below.  Note:
C            STFCPA = sin(TF) * cos(PA)
C            STFSPA = sin(TF) * sin(PA)
C
C          COMMON BLOCKS:
C            COMMON /DIPOLE/ COLAT,ELON,VP,CTP,STP
C
C          DIPOLE has IGRF variables obtained from routines in magfld.f:
C            COLAT = Geocentric colatitude of geomagnetic dipole north pole (deg)
C            ELON  = East longitude of geomagnetic dipole north pole (deg)
C            VP    = Magnitude (T-m) of dipole component of magnetic potential at
C                    geomagnetic pole and geocentric radius of 6371.2 km
C            CTP   = cosine of COLAT
C            STP   = sine   of COLAT
C------------------------------------------------------------------------------
C          HISTORY:
C          May 1994:  Completed on the 1st by A. D. Richmond

      PARAMETER (RTOD = 57.2957795130823,   RE =6371.2,
     +           DTOR = .01745329251994330, REQ=6378.137)
      COMMON /DIPOLE/ COLAT,ELON,VP,CTP,STP

      BHOR = SQRT(BNORTH*BNORTH + BEAST*BEAST)
      IF (BHOR .EQ. 0.) THEN
	ALON = 0.
	A    = 1.E34
	RETURN
      ENDIF
      COTTD  = BDOWN*.5/BHOR
      STD    = 1./SQRT(1.+COTTD*COTTD)
      CTD    = COTTD*STD
      SB     = -BEAST /BHOR
      CB     = -BNORTH/BHOR
      CTG    = SIN (GDLAT*DTOR)
      STG    = COS (GDLAT*DTOR)
      ANG    = (GDLON-ELON)*DTOR
      SANG   = SIN(ANG)
      CANG   = COS(ANG)
      CTE    = CTG*STD + STG*CTD*CB
      STE    = SQRT(1. - CTE*CTE)
      SA     = SB*CTD/STE
      CA     = (STD*STG - CTD*CTG*CB)/STE
      CAPANG = CA*CANG - SA*SANG
      SAPANG = CA*SANG + SA*CANG
      STFCPA = STE*CTP*CAPANG - CTE*STP
      STFSPA = SAPANG*STE
      ALON = ATAN2 (STFSPA,STFCPA)*RTOD
      R    = ALT + RE
      HA   = ALT + R*COTTD*COTTD
      A    = 1. + HA/REQ
      RETURN
      END

      FUNCTION FINT (X1,X2,X3,Y1,Y2,Y3, XFIT)
C          Second degree interpolation used by FNDAPX
C          INPUTS:
C            X1   = point 1 ordinate value
C            X2   = point 2 ordinate value
C            X3   = point 3 ordinate value
C            Y1   = point 1 abscissa value
C            Y2   = point 2 abscissa value
C            Y3   = point 3 abscissa value
C            XFIT = ordinate value to fit
C          RETURNS:
C            YFIT = abscissa value corresponding to XFIT
C
C          MODIFICATIONS:
C          Apr 2004: Change from subroutine to function, rename variables and
C                    add intermediates which are otherwise calculated twice
      X12 = X1-X2
      X13 = X1-X3
      X23 = X2-X3
      XF1 = XFIT-X1
      XF2 = XFIT-X2
      XF3 = XFIT-X3

      FINT = (Y1*X23*XF2*XF3 - Y2*X13*XF1*XF3 + Y3*X12*XF1*XF2) /
     +                                                     (X12*X13*X23)
      RETURN
      END
      SUBROUTINE COFRM (DATE)
C          Define the International Geomagnetic Reference Field (IGRF) as a
C          scalar potential field using a truncated series expansion with
C          Schmidt semi-normalized associated Legendre functions of degree n and
C          order m.  The polynomial coefficients are a function of time and are
C          interpolated between five year epochs or extrapolated at a constant
C          rate after the last epoch.
C
C          INPUTS:
C            DATE = yyyy.fraction (UT)
C          OUTPUTS (in comnon block MAGCOF):
C            NMAX = Maximum order of spherical harmonic coefficients used
C            GB   = Coefficients for magnetic field calculation
C            GV   = Coefficients for magnetic potential calculation
C            ICHG = Flag indicating when GB,GV have been changed in COFRM
C
C          It is fatal to supply a DATE before the first epoch.  A warning is
C          issued to Fortran unit 0 (stderr) if DATE is later than the
C          recommended limit, five years after the last epoch.
C
C          HISTORY (blame):
C          Apr 1983:  Written by Vincent B. Wickwar (Utah State Univ.) including
C          secular variation acceleration rate set to zero in case the IGRF
C          definition includes such second time derivitives.  The maximum degree
C          (n) defined was 10.
C
C          Jun 1986:  Updated coefficients adding Definitive Geomagnetic Reference
C          Field (DGRF) for 1980 and IGRF for 1985 (EOS Volume 7 Number 24).  The
C          designation DGRF means coefficients will not change in the future
C          whereas IGRF coefficients are interim pending incorporation of new
C          magnetometer data.  Common block MAG was replaced by MAGCOF, thus
C          removing variables not used in subroutine FELDG.  (Roy Barnes)
C
C          Apr 1992 (Barnes):  Added DGRF 1985 and IGRF 1990 as given in EOS
C          Volume 73 Number 16 April 21 1992.  Other changes were made so future
C          updates should:  (1) Increment NDGY; (2) Append to EPOCH the next IGRF
C          year; (3) Append the next DGRF coefficients to G1DIM and H1DIM; and (4)
C          replace the IGRF initial values (G0, GT) and rates of change indices
C          (H0, HT).
C
C          Apr 1994 (Art Richmond): Computation of GV added, for finding magnetic
C          potential.
C
C          Aug 1995 (Barnes):  Added DGRF for 1990 and IGRF for 1995, which were
C          obtained by anonymous ftp to geomag.gsfc.nasa.gov (cd pub, mget table*)
C          as per instructions from Bob Langel (langel@geomag.gsfc.nasa.gov) with
C          problems reported to baldwin@geomag.gsfc.nasa.gov.
C
C          Oct 1995 (Barnes):  Correct error in IGRF-95 G 7 6 and H 8 7 (see email
C          in folder).  Also found bug whereby coefficients were not being updated
C          in FELDG when IENTY did not change so ICHG was added to flag date
C          changes.  Also, a vestigial switch (IS) was removed from COFRM; it was
C          always zero and involved 3 branch if statements in the main polynomial
C          construction loop now numbered 200.
C
C          Feb 1999 (Barnes):  Explicitly initialize GV(1) in COFRM to avoid the
C          possibility of compiler or loader options initializing memory to
C          something else (e.g., indefinite).  Also simplify the algebra in COFRM
C          with no effect on results.
C
C          Mar 1999 (Barnes):  Removed three branch if's from FELDG and changed
C          statement labels to ascending order.
C
C          Jun 1999 (Barnes):  Corrected RTOD definition in GD2CART.
C
C          May 2000 (Barnes):  Replace IGRF 1995, add IGRF 2000, and extend the
C          earlier DGRF's back to 1900.  The coefficients came from an NGDC web
C          page.  Related documentation is in $APXROOT/docs/igrf.2000.*  where
C          $APXROOT, defined by 'source envapex', is traditionally ~bozo/apex).
C
C          Mar 2004 (Barnes):  Replace 1995 and 2000 coefficients; now both are
C          DGRF.  Coefficients for 2000 are degree 13 with precision increased to
C          tenths nT and accommodating this instigated changes:  (1) degree (NMAX)
C          is now a function of epoch (NMXE) to curtail irrelevant looping over
C          unused high order terms (n > 10 in epochs before 2000) when calculating
C          GB; (2) expand coefficients data statement layout for G1D and H1D,
C          formerly G1DIM and H1DIM; (3) omit secular variation acceleration terms
C          which were always zero; (4) increase array dimensions in common block
C          MAGCOF and associated arrays G and H in FELDG; (5) change earth's shape
C          in CONVRT from the IAU-1966 to the WGS-1984 spheroid; (6) eliminate
C          reference to 'definitive' in variables in COFRM which were not always
C          definitive; (7) change G to GB in COFRM s.t. arrays GB and GV in common
C          block MAGCOF are consistently named in all subroutines; (8) remove
C          unused constants in all five subroutines.  See EOS Volume 84 Number 46
C          November 18 2003, www.ngdc.noaa.gov/IAGA/vmod/igrf.html or local files
C          $APXROOT/docs/igrf.2004.*
C
C          Sept. 2005 (Maute): update with IGRF10 from 
C          http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html use script 
C          ~maute/apex.d/apex_update/igrf2f Note that the downloaded file the start 
C          column of the year in the first line has to be before the start of each 
C          number in the same column

      DOUBLE PRECISION F,F0
      COMMON /MAGCOF/ NMAX,GB(255),GV(225),ICHG
      DATA ICHG /-99999/
 
C          NEPO = Number of epochs
C          NGH  = Single dimensioned array size of 2D version (GYR or HYR)
C          NGHT = Single dimensioned array size of 2D version (GT  or HT)
      PARAMETER (NEPO = 22, NGH = 225*NEPO, NGHT = 225)
      DIMENSION GYR(15,15,NEPO), HYR(15,15,NEPO), EPOCH(NEPO),
     +          GT (15,15),      HT (15,15),       NMXE(NEPO),
     +          GY1D(NGH),       HY1D(NGH),
     +          GT1D(NGHT),      HT1D(NGHT)
      EQUIVALENCE (GYR(1,1,1),GY1D(1)), (HYR(1,1,1),HY1D(1)),
     +            (GT (1,1),  GT1D(1)), (HT (1,1),  HT1D(1))

      SAVE DATEL, EPOCH, NMXE, GYR, HYR, GT, HT, GY1D, HY1D, GT1D, HT1D
      DATA DATEL /-999./,
     +     EPOCH / 1900, 1905, 1910, 1915, 1920, 1925, 1930, 1935, 1940,
     +             1945, 1950, 1955, 1960, 1965, 1970, 1975, 1980, 1985,
     +             1990, 1995, 2000, 2005/,
     +     NMXE  /   10,   10,   10,   10,   10,   10,   10,   10,   10,
     +               10,   10,   10,   10,   10,   10,   10,   10,   10,
     +               10,   13,   13,   13/

C          g(n,m) for 1900
C          Fields across a line are (degree) n=1,13; lines are (order) m=0,13 as indicated
C          in column 6; e.g., for 1965 g(n=3,m=0) = 1297 or g(n=6,m=6) = -111
C
C           1       2       3      4      5      6      7      8      9
C                                        10     11     12     13          (n)
      DATA (GY1D(I),I=1,145) /0,
     O  -31543,   -677,  1022,   876,  -184,    63,    70,    11,     8,
     +                                   -3,     0,     0,     0,   2*0,
     1   -2298,   2905, -1469,   628,   328,    61,   -55,     8,    10,
     +                                   -4,     0,     0,     0,   3*0,
     2             924,  1256,   660,   264,   -11,     0,    -4,     1,
     +                                    2,     0,     0,     0,   4*0,
     3                    572,  -361,     5,  -217,    34,    -9,   -11,
     +                                   -5,     0,     0,     0,   5*0,
     4                           134,   -86,   -58,   -41,     1,    12,
     +                                   -2,     0,     0,     0,   6*0,
     5                                  -16,    59,   -21,     2,     1,
     +                                    6,     0,     0,     0,   7*0,
     6                                         -90,    18,    -9,    -2,
     +                                    4,     0,     0,     0,   8*0,
     7                                                  6,     5,     2,
     +                                    0,     0,     0,     0,   9*0,
     8                                                         8,    -1,
     +                                    2,     0,     0,     0,  10*0,
     9                                                               -1/
      DATA (GY1D(I),I=146,225) /
     +                                    2,     0,     0,     0,  11*0,
     O                                    0,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          g(n,m) for 1905
      DATA (GY1D(I),I=226,370) /0,
     O  -31464,   -728,  1037,   880,  -192,    62,    70,    11,     8,
     +                                   -3,     0,     0,     0,   2*0,
     1   -2298,   2928, -1494,   643,   328,    60,   -54,     8,    10,
     +                                   -4,     0,     0,     0,   3*0,
     2            1041,  1239,   653,   259,   -11,     0,    -4,     1,
     +                                    2,     0,     0,     0,   4*0,
     3                    635,  -380,    -1,  -221,    33,    -9,   -11,
     +                                   -5,     0,     0,     0,   5*0,
     4                           146,   -93,   -57,   -41,     1,    12,
     +                                   -2,     0,     0,     0,   6*0,
     5                                  -26,    57,   -20,     2,     1,
     +                                    6,     0,     0,     0,   7*0,
     6                                         -92,    18,    -8,    -2,
     +                                    4,     0,     0,     0,   8*0,
     7                                                  6,     5,     2,
     +                                    0,     0,     0,     0,   9*0,
     8                                                         8,     0,
     +                                    2,     0,     0,     0,  10*0,
     9                                                               -1/
      DATA (GY1D(I),I=371,450) /
     +                                    2,     0,     0,     0,  11*0,
     O                                    0,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          g(n,m) for 1910
      DATA (GY1D(I),I=451,595) /0,
     O  -31354,   -769,  1058,   884,  -201,    62,    71,    11,     8,
     +                                   -3,     0,     0,     0,   2*0,
     1   -2297,   2948, -1524,   660,   327,    58,   -54,     8,    10,
     +                                   -4,     0,     0,     0,   3*0,
     2            1176,  1223,   644,   253,   -11,     1,    -4,     1,
     +                                    2,     0,     0,     0,   4*0,
     3                    705,  -400,    -9,  -224,    32,    -9,   -11,
     +                                   -5,     0,     0,     0,   5*0,
     4                           160,  -102,   -54,   -40,     1,    12,
     +                                   -2,     0,     0,     0,   6*0,
     5                                  -38,    54,   -19,     2,     1,
     +                                    6,     0,     0,     0,   7*0,
     6                                         -95,    18,    -8,    -2,
     +                                    4,     0,     0,     0,   8*0,
     7                                                  6,     5,     2,
     +                                    0,     0,     0,     0,   9*0,
     8                                                         8,     0,
     +                                    2,     0,     0,     0,  10*0,
     9                                                               -1/
      DATA (GY1D(I),I=596,675) /
     +                                    2,     0,     0,     0,  11*0,
     O                                    0,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          g(n,m) for 1915
      DATA (GY1D(I),I=676,820) /0,
     O  -31212,   -802,  1084,   887,  -211,    61,    72,    11,     8,
     +                                   -3,     0,     0,     0,   2*0,
     1   -2306,   2956, -1559,   678,   327,    57,   -54,     8,    10,
     +                                   -4,     0,     0,     0,   3*0,
     2            1309,  1212,   631,   245,   -10,     2,    -4,     1,
     +                                    2,     0,     0,     0,   4*0,
     3                    778,  -416,   -16,  -228,    31,    -9,   -11,
     +                                   -5,     0,     0,     0,   5*0,
     4                           178,  -111,   -51,   -38,     2,    12,
     +                                   -2,     0,     0,     0,   6*0,
     5                                  -51,    49,   -18,     3,     1,
     +                                    6,     0,     0,     0,   7*0,
     6                                         -98,    19,    -8,    -2,
     +                                    4,     0,     0,     0,   8*0,
     7                                                  6,     6,     2,
     +                                    0,     0,     0,     0,   9*0,
     8                                                         8,     0,
     +                                    1,     0,     0,     0,  10*0,
     9                                                               -1/
      DATA (GY1D(I),I=821,900) /
     +                                    2,     0,     0,     0,  11*0,
     O                                    0,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          g(n,m) for 1920
      DATA (GY1D(I),I=901,1045) /0,
     O  -31060,   -839,  1111,   889,  -221,    61,    73,    11,     8,
     +                                   -3,     0,     0,     0,   2*0,
     1   -2317,   2959, -1600,   695,   326,    55,   -54,     7,    10,
     +                                   -4,     0,     0,     0,   3*0,
     2            1407,  1205,   616,   236,   -10,     2,    -3,     1,
     +                                    2,     0,     0,     0,   4*0,
     3                    839,  -424,   -23,  -233,    29,    -9,   -11,
     +                                   -5,     0,     0,     0,   5*0,
     4                           199,  -119,   -46,   -37,     2,    12,
     +                                   -2,     0,     0,     0,   6*0,
     5                                  -62,    44,   -16,     4,     1,
     +                                    6,     0,     0,     0,   7*0,
     6                                        -101,    19,    -7,    -2,
     +                                    4,     0,     0,     0,   8*0,
     7                                                  6,     6,     2,
     +                                    0,     0,     0,     0,   9*0,
     8                                                         8,     0,
     +                                    1,     0,     0,     0,  10*0,
     9                                                               -1/
      DATA (GY1D(I),I=1046,1125) /
     +                                    3,     0,     0,     0,  11*0,
     O                                    0,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          g(n,m) for 1925
      DATA (GY1D(I),I=1126,1270) /0,
     O  -30926,   -893,  1140,   891,  -230,    61,    73,    11,     8,
     +                                   -3,     0,     0,     0,   2*0,
     1   -2318,   2969, -1645,   711,   326,    54,   -54,     7,    10,
     +                                   -4,     0,     0,     0,   3*0,
     2            1471,  1202,   601,   226,    -9,     3,    -3,     1,
     +                                    2,     0,     0,     0,   4*0,
     3                    881,  -426,   -28,  -238,    27,    -9,   -11,
     +                                   -5,     0,     0,     0,   5*0,
     4                           217,  -125,   -40,   -35,     2,    12,
     +                                   -2,     0,     0,     0,   6*0,
     5                                  -69,    39,   -14,     4,     1,
     +                                    6,     0,     0,     0,   7*0,
     6                                        -103,    19,    -7,    -2,
     +                                    4,     0,     0,     0,   8*0,
     7                                                  6,     7,     2,
     +                                    0,     0,     0,     0,   9*0,
     8                                                         8,     0,
     +                                    1,     0,     0,     0,  10*0,
     9                                                               -1/
      DATA (GY1D(I),I=1271,1350) /
     +                                    3,     0,     0,     0,  11*0,
     O                                    0,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          g(n,m) for 1930
      DATA (GY1D(I),I=1351,1495) /0,
     O  -30805,   -951,  1172,   896,  -237,    60,    74,    11,     8,
     +                                   -3,     0,     0,     0,   2*0,
     1   -2316,   2980, -1692,   727,   327,    53,   -54,     7,    10,
     +                                   -4,     0,     0,     0,   3*0,
     2            1517,  1205,   584,   218,    -9,     4,    -3,     1,
     +                                    2,     0,     0,     0,   4*0,
     3                    907,  -422,   -32,  -242,    25,    -9,   -12,
     +                                   -5,     0,     0,     0,   5*0,
     4                           234,  -131,   -32,   -34,     2,    12,
     +                                   -2,     0,     0,     0,   6*0,
     5                                  -74,    32,   -12,     5,     1,
     +                                    6,     0,     0,     0,   7*0,
     6                                        -104,    18,    -6,    -2,
     +                                    4,     0,     0,     0,   8*0,
     7                                                  6,     8,     3,
     +                                    0,     0,     0,     0,   9*0,
     8                                                         8,     0,
     +                                    1,     0,     0,     0,  10*0,
     9                                                               -2/
      DATA (GY1D(I),I=1496,1575) /
     +                                    3,     0,     0,     0,  11*0,
     O                                    0,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          g(n,m) for 1935
      DATA (GY1D(I),I=1576,1720) /0,
     O  -30715,  -1018,  1206,   903,  -241,    59,    74,    11,     8,
     +                                   -3,     0,     0,     0,   2*0,
     1   -2306,   2984, -1740,   744,   329,    53,   -53,     7,    10,
     +                                   -4,     0,     0,     0,   3*0,
     2            1550,  1215,   565,   211,    -8,     4,    -3,     1,
     +                                    2,     0,     0,     0,   4*0,
     3                    918,  -415,   -33,  -246,    23,    -9,   -12,
     +                                   -5,     0,     0,     0,   5*0,
     4                           249,  -136,   -25,   -33,     1,    11,
     +                                   -2,     0,     0,     0,   6*0,
     5                                  -76,    25,   -11,     6,     1,
     +                                    6,     0,     0,     0,   7*0,
     6                                        -106,    18,    -6,    -2,
     +                                    4,     0,     0,     0,   8*0,
     7                                                  6,     8,     3,
     +                                    0,     0,     0,     0,   9*0,
     8                                                         7,     0,
     +                                    2,     0,     0,     0,  10*0,
     9                                                               -2/
      DATA (GY1D(I),I=1721,1800) /
     +                                    3,     0,     0,     0,  11*0,
     O                                    0,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          g(n,m) for 1940
      DATA (GY1D(I),I=1801,1945) /0,
     O  -30654,  -1106,  1240,   914,  -241,    57,    74,    11,     8,
     +                                   -3,     0,     0,     0,   2*0,
     1   -2292,   2981, -1790,   762,   334,    54,   -53,     7,    10,
     +                                   -4,     0,     0,     0,   3*0,
     2            1566,  1232,   550,   208,    -7,     4,    -3,     1,
     +                                    2,     0,     0,     0,   4*0,
     3                    916,  -405,   -33,  -249,    20,   -10,   -12,
     +                                   -5,     0,     0,     0,   5*0,
     4                           265,  -141,   -18,   -31,     1,    11,
     +                                   -2,     0,     0,     0,   6*0,
     5                                  -76,    18,    -9,     6,     1,
     +                                    6,     0,     0,     0,   7*0,
     6                                        -107,    17,    -5,    -2,
     +                                    4,     0,     0,     0,   8*0,
     7                                                  5,     9,     3,
     +                                    0,     0,     0,     0,   9*0,
     8                                                         7,     1,
     +                                    2,     0,     0,     0,  10*0,
     9                                                               -2/
      DATA (GY1D(I),I=1946,2025) /
     +                                    3,     0,     0,     0,  11*0,
     O                                    0,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          g(n,m) for 1945
      DATA (GY1D(I),I=2026,2170) /0,
     O  -30594,  -1244,  1282,   944,  -253,    59,    70,    13,     5,
     +                                   -3,     0,     0,     0,   2*0,
     1   -2285,   2990, -1834,   776,   346,    57,   -40,     7,   -21,
     +                                   11,     0,     0,     0,   3*0,
     2            1578,  1255,   544,   194,     6,     0,    -8,     1,
     +                                    1,     0,     0,     0,   4*0,
     3                    913,  -421,   -20,  -246,     0,    -5,   -11,
     +                                    2,     0,     0,     0,   5*0,
     4                           304,  -142,   -25,   -29,     9,     3,
     +                                   -5,     0,     0,     0,   6*0,
     5                                  -82,    21,   -10,     7,    16,
     +                                   -1,     0,     0,     0,   7*0,
     6                                        -104,    15,   -10,    -3,
     +                                    8,     0,     0,     0,   8*0,
     7                                                 29,     7,    -4,
     +                                   -1,     0,     0,     0,   9*0,
     8                                                         2,    -3,
     +                                   -3,     0,     0,     0,  10*0,
     9                                                               -4/
      DATA (GY1D(I),I=2171,2250) /
     +                                    5,     0,     0,     0,  11*0,
     O                                   -2,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          g(n,m) for 1950
      DATA (GY1D(I),I=2251,2395) /0,
     O  -30554,  -1341,  1297,   954,  -240,    54,    65,    22,     3,
     +                                   -8,     0,     0,     0,   2*0,
     1   -2250,   2998, -1889,   792,   349,    57,   -55,    15,    -7,
     +                                    4,     0,     0,     0,   3*0,
     2            1576,  1274,   528,   211,     4,     2,    -4,    -1,
     +                                   -1,     0,     0,     0,   4*0,
     3                    896,  -408,   -20,  -247,     1,    -1,   -25,
     +                                   13,     0,     0,     0,   5*0,
     4                           303,  -147,   -16,   -40,    11,    10,
     +                                   -4,     0,     0,     0,   6*0,
     5                                  -76,    12,    -7,    15,     5,
     +                                    4,     0,     0,     0,   7*0,
     6                                        -105,     5,   -13,    -5,
     +                                   12,     0,     0,     0,   8*0,
     7                                                 19,     5,    -2,
     +                                    3,     0,     0,     0,   9*0,
     8                                                        -1,     3,
     +                                    2,     0,     0,     0,  10*0,
     9                                                                8/
      DATA (GY1D(I),I=2396,2475) /
     +                                   10,     0,     0,     0,  11*0,
     O                                    3,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          g(n,m) for 1955
      DATA (GY1D(I),I=2476,2620) /0,
     O  -30500,  -1440,  1302,   958,  -229,    47,    65,    11,     4,
     +                                   -3,     0,     0,     0,   2*0,
     1   -2215,   3003, -1944,   796,   360,    57,   -56,     9,     9,
     +                                   -5,     0,     0,     0,   3*0,
     2            1581,  1288,   510,   230,     3,     2,    -6,    -4,
     +                                   -1,     0,     0,     0,   4*0,
     3                    882,  -397,   -23,  -247,    10,   -14,    -5,
     +                                    2,     0,     0,     0,   5*0,
     4                           290,  -152,    -8,   -32,     6,     2,
     +                                   -3,     0,     0,     0,   6*0,
     5                                  -69,     7,   -11,    10,     4,
     +                                    7,     0,     0,     0,   7*0,
     6                                        -107,     9,    -7,     1,
     +                                    4,     0,     0,     0,   8*0,
     7                                                 18,     6,     2,
     +                                   -2,     0,     0,     0,   9*0,
     8                                                         9,     2,
     +                                    6,     0,     0,     0,  10*0,
     9                                                                5/
      DATA (GY1D(I),I=2621,2700) /
     +                                   -2,     0,     0,     0,  11*0,
     O                                    0,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          g(n,m) for 1960
      DATA (GY1D(I),I=2701,2845) /0,
     O  -30421,  -1555,  1302,   957,  -222,    46,    67,    15,     4,
     +                                    1,     0,     0,     0,   2*0,
     1   -2169,   3002, -1992,   800,   362,    58,   -56,     6,     6,
     +                                   -3,     0,     0,     0,   3*0,
     2            1590,  1289,   504,   242,     1,     5,    -4,     0,
     +                                    4,     0,     0,     0,   4*0,
     3                    878,  -394,   -26,  -237,    15,   -11,    -9,
     +                                    0,     0,     0,     0,   5*0,
     4                           269,  -156,    -1,   -32,     2,     1,
     +                                   -1,     0,     0,     0,   6*0,
     5                                  -63,    -2,    -7,    10,     4,
     +                                    4,     0,     0,     0,   7*0,
     6                                        -113,    17,    -5,    -1,
     +                                    6,     0,     0,     0,   8*0,
     7                                                  8,    10,    -2,
     +                                    1,     0,     0,     0,   9*0,
     8                                                         8,     3,
     +                                   -1,     0,     0,     0,  10*0,
     9                                                               -1/
      DATA (GY1D(I),I=2846,2925) /
     +                                    2,     0,     0,     0,  11*0,
     O                                    0,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          g(n,m) for 1965
      DATA (GY1D(I),I=2926,3070) /0,
     O  -30334,  -1662,  1297,   957,  -219,    45,    75,    13,     8,
     +                                   -2,     0,     0,     0,   2*0,
     1   -2119,   2997, -2038,   804,   358,    61,   -57,     5,    10,
     +                                   -3,     0,     0,     0,   3*0,
     2            1594,  1292,   479,   254,     8,     4,    -4,     2,
     +                                    2,     0,     0,     0,   4*0,
     3                    856,  -390,   -31,  -228,    13,   -14,   -13,
     +                                   -5,     0,     0,     0,   5*0,
     4                           252,  -157,     4,   -26,     0,    10,
     +                                   -2,     0,     0,     0,   6*0,
     5                                  -62,     1,    -6,     8,    -1,
     +                                    4,     0,     0,     0,   7*0,
     6                                        -111,    13,    -1,    -1,
     +                                    4,     0,     0,     0,   8*0,
     7                                                  1,    11,     5,
     +                                    0,     0,     0,     0,   9*0,
     8                                                         4,     1,
     +                                    2,     0,     0,     0,  10*0,
     9                                                               -2/
      DATA (GY1D(I),I=3071,3150) /
     +                                    2,     0,     0,     0,  11*0,
     O                                    0,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          g(n,m) for 1970
      DATA (GY1D(I),I=3151,3295) /0,
     O  -30220,  -1781,  1287,   952,  -216,    43,    72,    14,     8,
     +                                   -3,     0,     0,     0,   2*0,
     1   -2068,   3000, -2091,   800,   359,    64,   -57,     6,    10,
     +                                   -3,     0,     0,     0,   3*0,
     2            1611,  1278,   461,   262,    15,     1,    -2,     2,
     +                                    2,     0,     0,     0,   4*0,
     3                    838,  -395,   -42,  -212,    14,   -13,   -12,
     +                                   -5,     0,     0,     0,   5*0,
     4                           234,  -160,     2,   -22,    -3,    10,
     +                                   -1,     0,     0,     0,   6*0,
     5                                  -56,     3,    -2,     5,    -1,
     +                                    6,     0,     0,     0,   7*0,
     6                                        -112,    13,     0,     0,
     +                                    4,     0,     0,     0,   8*0,
     7                                                 -2,    11,     3,
     +                                    1,     0,     0,     0,   9*0,
     8                                                         3,     1,
     +                                    0,     0,     0,     0,  10*0,
     9                                                               -1/
      DATA (GY1D(I),I=3296,3375) /
     +                                    3,     0,     0,     0,  11*0,
     O                                   -1,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          g(n,m) for 1975
      DATA (GY1D(I),I=3376,3520) /0,
     O  -30100,  -1902,  1276,   946,  -218,    45,    71,    14,     7,
     +                                   -3,     0,     0,     0,   2*0,
     1   -2013,   3010, -2144,   791,   356,    66,   -56,     6,    10,
     +                                   -3,     0,     0,     0,   3*0,
     2            1632,  1260,   438,   264,    28,     1,    -1,     2,
     +                                    2,     0,     0,     0,   4*0,
     3                    830,  -405,   -59,  -198,    16,   -12,   -12,
     +                                   -5,     0,     0,     0,   5*0,
     4                           216,  -159,     1,   -14,    -8,    10,
     +                                   -2,     0,     0,     0,   6*0,
     5                                  -49,     6,     0,     4,    -1,
     +                                    5,     0,     0,     0,   7*0,
     6                                        -111,    12,     0,    -1,
     +                                    4,     0,     0,     0,   8*0,
     7                                                 -5,    10,     4,
     +                                    1,     0,     0,     0,   9*0,
     8                                                         1,     1,
     +                                    0,     0,     0,     0,  10*0,
     9                                                               -2/
      DATA (GY1D(I),I=3521,3600) /
     +                                    3,     0,     0,     0,  11*0,
     O                                   -1,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          g(n,m) for 1980
      DATA (GY1D(I),I=3601,3745) /0,
     O  -29992,  -1997,  1281,   938,  -218,    48,    72,    18,     5,
     +                                   -4,     0,     0,     0,   2*0,
     1   -1956,   3027, -2180,   782,   357,    66,   -59,     6,    10,
     +                                   -4,     0,     0,     0,   3*0,
     2            1663,  1251,   398,   261,    42,     2,     0,     1,
     +                                    2,     0,     0,     0,   4*0,
     3                    833,  -419,   -74,  -192,    21,   -11,   -12,
     +                                   -5,     0,     0,     0,   5*0,
     4                           199,  -162,     4,   -12,    -7,     9,
     +                                   -2,     0,     0,     0,   6*0,
     5                                  -48,    14,     1,     4,    -3,
     +                                    5,     0,     0,     0,   7*0,
     6                                        -108,    11,     3,    -1,
     +                                    3,     0,     0,     0,   8*0,
     7                                                 -2,     6,     7,
     +                                    1,     0,     0,     0,   9*0,
     8                                                        -1,     2,
     +                                    2,     0,     0,     0,  10*0,
     9                                                               -5/
      DATA (GY1D(I),I=3746,3825) /
     +                                    3,     0,     0,     0,  11*0,
     O                                    0,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          g(n,m) for 1985
      DATA (GY1D(I),I=3826,3970) /0,
     O  -29873,  -2072,  1296,   936,  -214,    53,    74,    21,     5,
     +                                   -4,     0,     0,     0,   2*0,
     1   -1905,   3044, -2208,   780,   355,    65,   -62,     6,    10,
     +                                   -4,     0,     0,     0,   3*0,
     2            1687,  1247,   361,   253,    51,     3,     0,     1,
     +                                    3,     0,     0,     0,   4*0,
     3                    829,  -424,   -93,  -185,    24,   -11,   -12,
     +                                   -5,     0,     0,     0,   5*0,
     4                           170,  -164,     4,    -6,    -9,     9,
     +                                   -2,     0,     0,     0,   6*0,
     5                                  -46,    16,     4,     4,    -3,
     +                                    5,     0,     0,     0,   7*0,
     6                                        -102,    10,     4,    -1,
     +                                    3,     0,     0,     0,   8*0,
     7                                                  0,     4,     7,
     +                                    1,     0,     0,     0,   9*0,
     8                                                        -4,     1,
     +                                    2,     0,     0,     0,  10*0,
     9                                                               -5/
      DATA (GY1D(I),I=3971,4050) /
     +                                    3,     0,     0,     0,  11*0,
     O                                    0,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          g(n,m) for 1990
      DATA (GY1D(I),I=4051,4195) /0,
     O  -29775,  -2131,  1314,   939,  -214,    61,    77,    23,     4,
     +                                   -3,     0,     0,     0,   2*0,
     1   -1848,   3059, -2239,   780,   353,    65,   -64,     5,     9,
     +                                   -4,     0,     0,     0,   3*0,
     2            1686,  1248,   325,   245,    59,     2,    -1,     1,
     +                                    2,     0,     0,     0,   4*0,
     3                    802,  -423,  -109,  -178,    26,   -10,   -12,
     +                                   -5,     0,     0,     0,   5*0,
     4                           141,  -165,     3,    -1,   -12,     9,
     +                                   -2,     0,     0,     0,   6*0,
     5                                  -36,    18,     5,     3,    -4,
     +                                    4,     0,     0,     0,   7*0,
     6                                         -96,     9,     4,    -2,
     +                                    3,     0,     0,     0,   8*0,
     7                                                  0,     2,     7,
     +                                    1,     0,     0,     0,   9*0,
     8                                                        -6,     1,
     +                                    3,     0,     0,     0,  10*0,
     9                                                               -6/
      DATA (GY1D(I),I=4196,4275) /
     +                                    3,     0,     0,     0,  11*0,
     O                                    0,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          g(n,m) for 1995
      DATA (GY1D(I),I=4276,4420) /0,
     O  -29692,  -2200,  1335,   940,  -214,    68,    77,    25,     4,
     +                                   -3,     0,     0,     0,   2*0,
     1   -1784,   3070, -2267,   780,   352,    67,   -72,     6,     9,
     +                                   -6,     0,     0,     0,   3*0,
     2            1681,  1249,   290,   235,    68,     1,    -6,     3,
     +                                    2,     0,     0,     0,   4*0,
     3                    759,  -418,  -118,  -170,    28,    -9,   -10,
     +                                   -4,     0,     0,     0,   5*0,
     4                           122,  -166,    -1,     5,   -14,     8,
     +                                   -1,     0,     0,     0,   6*0,
     5                                  -17,    19,     4,     9,    -8,
     +                                    4,     0,     0,     0,   7*0,
     6                                         -93,     8,     6,    -1,
     +                                    2,     0,     0,     0,   8*0,
     7                                                 -2,    -5,    10,
     +                                    2,     0,     0,     0,   9*0,
     8                                                        -7,    -2,
     +                                    5,     0,     0,     0,  10*0,
     9                                                               -8/
      DATA (GY1D(I),I=4421,4500) /
     +                                    1,     0,     0,     0,  11*0,
     O                                    0,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          g(n,m) for 2000
      DATA (GY1D(I),I=4501,4645) /0,
     O-29619.4,-2267.7,1339.6, 932.3,-218.8,  72.3,  79.0,  24.4,   5.0,
     +                                 -2.6,   2.7,  -2.2,  -0.2,   2*0,
     1 -1728.2, 3068.4,-2288.0, 786.8, 351.4,  68.2, -74.0,   6.6,  9.4,
     +                                 -6.0,  -1.7,  -0.3,  -0.9,   3*0,
     2          1670.9,1252.1, 250.0, 222.3,  74.2,   0.0,  -9.2,   3.0,
     +                                  1.7,  -1.9,   0.2,   0.3,   4*0,
     3                  714.5,-403.0,-130.4,-160.9,  33.3,  -7.9,  -8.4,
     +                                 -3.1,   1.5,   0.9,   0.1,   5*0,
     4                         111.3,-168.6,  -5.9,   9.1, -16.6,   6.3,
     +                                 -0.5,  -0.1,  -0.2,  -0.4,   6*0,
     5                                -12.9,  16.9,   6.9,   9.1,  -8.9,
     +                                  3.7,   0.1,   0.9,   1.3,   7*0,
     6                                       -90.4,   7.3,   7.0,  -1.5,
     +                                  1.0,  -0.7,  -0.5,  -0.4,   8*0,
     7                                               -1.2,  -7.9,   9.3,
     +                                  2.0,   0.7,   0.3,   0.7,   9*0,
     8                                                      -7.0,  -4.3,
     +                                  4.2,   1.7,  -0.3,  -0.4,  10*0,
     9                                                            -8.2/
      DATA (GY1D(I),I=4646,4725) /
     +                                  0.3,   0.1,  -0.4,   0.3,  11*0,
     O                                 -1.1,   1.2,  -0.1,  -0.1,  12*0,
     1                                         4.0,  -0.2,   0.4,  13*0,
     2                                               -0.4,   0.0,  14*0,
     3                                                       0.1,  16*0/
C          g(n,m) for 2005
      DATA (GY1D(I),I=4726,4870) /0,
     O-29556.8,-2340.5,1335.7, 919.8,-227.6,  72.9,  79.8,  24.8,   5.6,
     +                                 -2.2,   2.9,  -2.2,  -0.2,   2*0,
     1 -1671.8, 3047.0,-2305.3, 798.2, 354.4,  69.6, -74.4,   7.7,  9.8,
     +                                 -6.3,  -1.6,  -0.3,  -0.9,   3*0,
     2          1656.9,1246.8, 211.5, 208.8,  76.6,  -1.4, -11.4,   3.6,
     +                                  1.6,  -1.7,   0.3,   0.3,   4*0,
     3                  674.4,-379.5,-136.6,-151.1,  38.6,  -6.8,  -7.0,
     +                                 -2.5,   1.5,   0.9,   0.3,   5*0,
     4                         100.2,-168.3, -15.0,  12.3, -18.0,   5.0,
     +                                 -0.1,  -0.2,  -0.4,  -0.4,   6*0,
     5                                -14.1,  14.7,   9.4,  10.0, -10.8,
     +                                  3.0,   0.2,   1.0,   1.2,   7*0,
     6                                       -86.4,   5.5,   9.4,  -1.3,
     +                                  0.3,  -0.7,  -0.4,  -0.4,   8*0,
     7                                                2.0, -11.4,   8.7,
     +                                  2.1,   0.5,   0.5,   0.7,   9*0,
     8                                                      -5.0,  -6.7,
     +                                  3.9,   1.8,  -0.3,  -0.3,  10*0,
     9                                                            -9.2/
      DATA (GY1D(I),I=4871,4950) /
     +                                 -0.1,   0.1,  -0.4,   0.4,  11*0,
     O                                 -2.2,   1.0,   0.0,  -0.1,  12*0,
     1                                         4.1,  -0.4,   0.4,  13*0,
     2                                                0.0,  -0.1,  14*0,
     3                                                      -0.3,  16*0/
C          h(n,m) for 1900
      DATA (HY1D(I),I=1,145) /16*0,
     1    5922,  -1061,  -330,   195,  -210,    -9,   -45,     8,   -20,
     +                                    2,     0,     0,     0,   3*0,
     2            1121,     3,   -69,    53,    83,   -13,   -14,    14,
     +                                    1,     0,     0,     0,   4*0,
     3                    523,  -210,   -33,     2,   -10,     7,     5,
     +                                    2,     0,     0,     0,   5*0,
     4                           -75,  -124,   -35,    -1,   -13,    -3,
     +                                    6,     0,     0,     0,   6*0,
     5                                    3,    36,    28,     5,    -2,
     +                                   -4,     0,     0,     0,   7*0,
     6                                         -69,   -12,    16,     8,
     +                                    0,     0,     0,     0,   8*0,
     7                                                -22,    -5,    10,
     +                                   -2,     0,     0,     0,   9*0,
     8                                                       -18,    -2,
     +                                    4,     0,     0,     0,  10*0,
     9                                                                2/
      DATA (HY1D(I),I=146,225) /
     +                                    0,     0,     0,     0,  11*0,
     O                                   -6,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          h(n,m) for 1905
      DATA (HY1D(I),I=226,370) /16*0,
     1    5909,  -1086,  -357,   203,  -193,    -7,   -46,     8,   -20,
     +                                    2,     0,     0,     0,   3*0,
     2            1065,    34,   -77,    56,    86,   -14,   -15,    14,
     +                                    1,     0,     0,     0,   4*0,
     3                    480,  -201,   -32,     4,   -11,     7,     5,
     +                                    2,     0,     0,     0,   5*0,
     4                           -65,  -125,   -32,     0,   -13,    -3,
     +                                    6,     0,     0,     0,   6*0,
     5                                   11,    32,    28,     5,    -2,
     +                                   -4,     0,     0,     0,   7*0,
     6                                         -67,   -12,    16,     8,
     +                                    0,     0,     0,     0,   8*0,
     7                                                -22,    -5,    10,
     +                                   -2,     0,     0,     0,   9*0,
     8                                                       -18,    -2,
     +                                    4,     0,     0,     0,  10*0,
     9                                                                2/
      DATA (HY1D(I),I=371,450) /
     +                                    0,     0,     0,     0,  11*0,
     O                                   -6,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          h(n,m) for 1910
      DATA (HY1D(I),I=451,595) /16*0,
     1    5898,  -1128,  -389,   211,  -172,    -5,   -47,     8,   -20,
     +                                    2,     0,     0,     0,   3*0,
     2            1000,    62,   -90,    57,    89,   -14,   -15,    14,
     +                                    1,     0,     0,     0,   4*0,
     3                    425,  -189,   -33,     5,   -12,     6,     5,
     +                                    2,     0,     0,     0,   5*0,
     4                           -55,  -126,   -29,     1,   -13,    -3,
     +                                    6,     0,     0,     0,   6*0,
     5                                   21,    28,    28,     5,    -2,
     +                                   -4,     0,     0,     0,   7*0,
     6                                         -65,   -13,    16,     8,
     +                                    0,     0,     0,     0,   8*0,
     7                                                -22,    -5,    10,
     +                                   -2,     0,     0,     0,   9*0,
     8                                                       -18,    -2,
     +                                    4,     0,     0,     0,  10*0,
     9                                                                2/
      DATA (HY1D(I),I=596,675) /
     +                                    0,     0,     0,     0,  11*0,
     O                                   -6,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          h(n,m) for 1915
      DATA (HY1D(I),I=676,820) /16*0,
     1    5875,  -1191,  -421,   218,  -148,    -2,   -48,     8,   -20,
     +                                    2,     0,     0,     0,   3*0,
     2             917,    84,  -109,    58,    93,   -14,   -15,    14,
     +                                    1,     0,     0,     0,   4*0,
     3                    360,  -173,   -34,     8,   -12,     6,     5,
     +                                    2,     0,     0,     0,   5*0,
     4                           -51,  -126,   -26,     2,   -13,    -3,
     +                                    6,     0,     0,     0,   6*0,
     5                                   32,    23,    28,     5,    -2,
     +                                   -4,     0,     0,     0,   7*0,
     6                                         -62,   -15,    16,     8,
     +                                    0,     0,     0,     0,   8*0,
     7                                                -22,    -5,    10,
     +                                   -2,     0,     0,     0,   9*0,
     8                                                       -18,    -2,
     +                                    4,     0,     0,     0,  10*0,
     9                                                                2/
      DATA (HY1D(I),I=821,900) /
     +                                    0,     0,     0,     0,  11*0,
     O                                   -6,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          h(n,m) for 1920
      DATA (HY1D(I),I=901,1045) /16*0,
     1    5845,  -1259,  -445,   220,  -122,     0,   -49,     8,   -20,
     +                                    2,     0,     0,     0,   3*0,
     2             823,   103,  -134,    58,    96,   -14,   -15,    14,
     +                                    1,     0,     0,     0,   4*0,
     3                    293,  -153,   -38,    11,   -13,     6,     5,
     +                                    2,     0,     0,     0,   5*0,
     4                           -57,  -125,   -22,     4,   -14,    -3,
     +                                    6,     0,     0,     0,   6*0,
     5                                   43,    18,    28,     5,    -2,
     +                                   -4,     0,     0,     0,   7*0,
     6                                         -57,   -16,    17,     9,
     +                                    0,     0,     0,     0,   8*0,
     7                                                -22,    -5,    10,
     +                                   -2,     0,     0,     0,   9*0,
     8                                                       -19,    -2,
     +                                    4,     0,     0,     0,  10*0,
     9                                                                2/
      DATA (HY1D(I),I=1046,1125) /
     +                                    0,     0,     0,     0,  11*0,
     O                                   -6,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          h(n,m) for 1925
      DATA (HY1D(I),I=1126,1270) /16*0,
     1    5817,  -1334,  -462,   216,   -96,     3,   -50,     8,   -20,
     +                                    2,     0,     0,     0,   3*0,
     2             728,   119,  -163,    58,    99,   -14,   -15,    14,
     +                                    1,     0,     0,     0,   4*0,
     3                    229,  -130,   -44,    14,   -14,     6,     5,
     +                                    2,     0,     0,     0,   5*0,
     4                           -70,  -122,   -18,     5,   -14,    -3,
     +                                    6,     0,     0,     0,   6*0,
     5                                   51,    13,    29,     5,    -2,
     +                                   -4,     0,     0,     0,   7*0,
     6                                         -52,   -17,    17,     9,
     +                                    0,     0,     0,     0,   8*0,
     7                                                -21,    -5,    10,
     +                                   -2,     0,     0,     0,   9*0,
     8                                                       -19,    -2,
     +                                    4,     0,     0,     0,  10*0,
     9                                                                2/
      DATA (HY1D(I),I=1271,1350) /
     +                                    0,     0,     0,     0,  11*0,
     O                                   -6,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          h(n,m) for 1930
      DATA (HY1D(I),I=1351,1495) /16*0,
     1    5808,  -1424,  -480,   205,   -72,     4,   -51,     8,   -20,
     +                                    2,     0,     0,     0,   3*0,
     2             644,   133,  -195,    60,   102,   -15,   -15,    14,
     +                                    1,     0,     0,     0,   4*0,
     3                    166,  -109,   -53,    19,   -14,     5,     5,
     +                                    2,     0,     0,     0,   5*0,
     4                           -90,  -118,   -16,     6,   -14,    -3,
     +                                    6,     0,     0,     0,   6*0,
     5                                   58,     8,    29,     5,    -2,
     +                                   -4,     0,     0,     0,   7*0,
     6                                         -46,   -18,    18,     9,
     +                                    0,     0,     0,     0,   8*0,
     7                                                -20,    -5,    10,
     +                                   -2,     0,     0,     0,   9*0,
     8                                                       -19,    -2,
     +                                    4,     0,     0,     0,  10*0,
     9                                                                2/
      DATA (HY1D(I),I=1496,1575) /
     +                                    0,     0,     0,     0,  11*0,
     O                                   -6,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          h(n,m) for 1935
      DATA (HY1D(I),I=1576,1720) /16*0,
     1    5812,  -1520,  -494,   188,   -51,     4,   -52,     8,   -20,
     +                                    2,     0,     0,     0,   3*0,
     2             586,   146,  -226,    64,   104,   -17,   -15,    15,
     +                                    1,     0,     0,     0,   4*0,
     3                    101,   -90,   -64,    25,   -14,     5,     5,
     +                                    2,     0,     0,     0,   5*0,
     4                          -114,  -115,   -15,     7,   -15,    -3,
     +                                    6,     0,     0,     0,   6*0,
     5                                   64,     4,    29,     5,    -3,
     +                                   -4,     0,     0,     0,   7*0,
     6                                         -40,   -19,    18,     9,
     +                                    0,     0,     0,     0,   8*0,
     7                                                -19,    -5,    11,
     +                                   -1,     0,     0,     0,   9*0,
     8                                                       -19,    -2,
     +                                    4,     0,     0,     0,  10*0,
     9                                                                2/
      DATA (HY1D(I),I=1721,1800) /
     +                                    0,     0,     0,     0,  11*0,
     O                                   -6,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          h(n,m) for 1940
      DATA (HY1D(I),I=1801,1945) /16*0,
     1    5821,  -1614,  -499,   169,   -33,     4,   -52,     8,   -21,
     +                                    2,     0,     0,     0,   3*0,
     2             528,   163,  -252,    71,   105,   -18,   -14,    15,
     +                                    1,     0,     0,     0,   4*0,
     3                     43,   -72,   -75,    33,   -14,     5,     5,
     +                                    2,     0,     0,     0,   5*0,
     4                          -141,  -113,   -15,     7,   -15,    -3,
     +                                    6,     0,     0,     0,   6*0,
     5                                   69,     0,    29,     5,    -3,
     +                                   -4,     0,     0,     0,   7*0,
     6                                         -33,   -20,    19,     9,
     +                                    0,     0,     0,     0,   8*0,
     7                                                -19,    -5,    11,
     +                                   -1,     0,     0,     0,   9*0,
     8                                                       -19,    -2,
     +                                    4,     0,     0,     0,  10*0,
     9                                                                2/
      DATA (HY1D(I),I=1946,2025) /
     +                                    0,     0,     0,     0,  11*0,
     O                                   -6,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          h(n,m) for 1945
      DATA (HY1D(I),I=2026,2170) /16*0,
     1    5810,  -1702,  -499,   144,   -12,     6,   -45,    12,   -27,
     +                                    5,     0,     0,     0,   3*0,
     2             477,   186,  -276,    95,   100,   -18,   -21,    17,
     +                                    1,     0,     0,     0,   4*0,
     3                    -11,   -55,   -67,    16,     2,   -12,    29,
     +                                  -20,     0,     0,     0,   5*0,
     4                          -178,  -119,    -9,     6,    -7,    -9,
     +                                   -1,     0,     0,     0,   6*0,
     5                                   82,   -16,    28,     2,     4,
     +                                   -6,     0,     0,     0,   7*0,
     6                                         -39,   -17,    18,     9,
     +                                    6,     0,     0,     0,   8*0,
     7                                                -22,     3,     6,
     +                                   -4,     0,     0,     0,   9*0,
     8                                                       -11,     1,
     +                                   -2,     0,     0,     0,  10*0,
     9                                                                8/
      DATA (HY1D(I),I=2171,2250) /
     +                                    0,     0,     0,     0,  11*0,
     O                                   -2,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          h(n,m) for 1950
      DATA (HY1D(I),I=2251,2395) /16*0,
     1    5815,  -1810,  -476,   136,     3,    -1,   -35,     5,   -24,
     +                                   13,     0,     0,     0,   3*0,
     2             381,   206,  -278,   103,    99,   -17,   -22,    19,
     +                                   -2,     0,     0,     0,   4*0,
     3                    -46,   -37,   -87,    33,     0,     0,    12,
     +                                  -10,     0,     0,     0,   5*0,
     4                          -210,  -122,   -12,    10,   -21,     2,
     +                                    2,     0,     0,     0,   6*0,
     5                                   80,   -12,    36,    -8,     2,
     +                                   -3,     0,     0,     0,   7*0,
     6                                         -30,   -18,    17,     8,
     +                                    6,     0,     0,     0,   8*0,
     7                                                -16,    -4,     8,
     +                                   -3,     0,     0,     0,   9*0,
     8                                                       -17,   -11,
     +                                    6,     0,     0,     0,  10*0,
     9                                                               -7/
      DATA (HY1D(I),I=2396,2475) /
     +                                   11,     0,     0,     0,  11*0,
     O                                    8,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          h(n,m) for 1955
      DATA (HY1D(I),I=2476,2620) /16*0,
     1    5820,  -1898,  -462,   133,    15,    -9,   -50,    10,   -11,
     +                                   -4,     0,     0,     0,   3*0,
     2             291,   216,  -274,   110,    96,   -24,   -15,    12,
     +                                    0,     0,     0,     0,   4*0,
     3                    -83,   -23,   -98,    48,    -4,     5,     7,
     +                                   -8,     0,     0,     0,   5*0,
     4                          -230,  -121,   -16,     8,   -23,     6,
     +                                   -2,     0,     0,     0,   6*0,
     5                                   78,   -12,    28,     3,    -2,
     +                                   -4,     0,     0,     0,   7*0,
     6                                         -24,   -20,    23,    10,
     +                                    1,     0,     0,     0,   8*0,
     7                                                -18,    -4,     7,
     +                                   -3,     0,     0,     0,   9*0,
     8                                                       -13,    -6,
     +                                    7,     0,     0,     0,  10*0,
     9                                                                5/
      DATA (HY1D(I),I=2621,2700) /
     +                                   -1,     0,     0,     0,  11*0,
     O                                   -3,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          h(n,m) for 1960
      DATA (HY1D(I),I=2701,2845) /16*0,
     1    5791,  -1967,  -414,   135,    16,   -10,   -55,    11,   -18,
     +                                    4,     0,     0,     0,   3*0,
     2             206,   224,  -278,   125,    99,   -28,   -14,    12,
     +                                    1,     0,     0,     0,   4*0,
     3                   -130,     3,  -117,    60,    -6,     7,     2,
     +                                    0,     0,     0,     0,   5*0,
     4                          -255,  -114,   -20,     7,   -18,     0,
     +                                    2,     0,     0,     0,   6*0,
     5                                   81,   -11,    23,     4,    -3,
     +                                   -5,     0,     0,     0,   7*0,
     6                                         -17,   -18,    23,     9,
     +                                    1,     0,     0,     0,   8*0,
     7                                                -17,     1,     8,
     +                                   -1,     0,     0,     0,   9*0,
     8                                                       -20,     0,
     +                                    6,     0,     0,     0,  10*0,
     9                                                                5/
      DATA (HY1D(I),I=2846,2925) /
     +                                    0,     0,     0,     0,  11*0,
     O                                   -7,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          h(n,m) for 1965
      DATA (HY1D(I),I=2926,3070) /16*0,
     1    5776,  -2016,  -404,   148,    19,   -11,   -61,     7,   -22,
     +                                    2,     0,     0,     0,   3*0,
     2             114,   240,  -269,   128,   100,   -27,   -12,    15,
     +                                    1,     0,     0,     0,   4*0,
     3                   -165,    13,  -126,    68,    -2,     9,     7,
     +                                    2,     0,     0,     0,   5*0,
     4                          -269,   -97,   -32,     6,   -16,    -4,
     +                                    6,     0,     0,     0,   6*0,
     5                                   81,    -8,    26,     4,    -5,
     +                                   -4,     0,     0,     0,   7*0,
     6                                          -7,   -23,    24,    10,
     +                                    0,     0,     0,     0,   8*0,
     7                                                -12,    -3,    10,
     +                                   -2,     0,     0,     0,   9*0,
     8                                                       -17,    -4,
     +                                    3,     0,     0,     0,  10*0,
     9                                                                1/
      DATA (HY1D(I),I=3071,3150) /
     +                                    0,     0,     0,     0,  11*0,
     O                                   -6,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          h(n,m) for 1970
      DATA (HY1D(I),I=3151,3295) /16*0,
     1    5737,  -2047,  -366,   167,    26,   -12,   -70,     7,   -21,
     +                                    1,     0,     0,     0,   3*0,
     2              25,   251,  -266,   139,   100,   -27,   -15,    16,
     +                                    1,     0,     0,     0,   4*0,
     3                   -196,    26,  -139,    72,    -4,     6,     6,
     +                                    3,     0,     0,     0,   5*0,
     4                          -279,   -91,   -37,     8,   -17,    -4,
     +                                    4,     0,     0,     0,   6*0,
     5                                   83,    -6,    23,     6,    -5,
     +                                   -4,     0,     0,     0,   7*0,
     6                                           1,   -23,    21,    10,
     +                                    0,     0,     0,     0,   8*0,
     7                                                -11,    -6,    11,
     +                                   -1,     0,     0,     0,   9*0,
     8                                                       -16,    -2,
     +                                    3,     0,     0,     0,  10*0,
     9                                                                1/
      DATA (HY1D(I),I=3296,3375) /
     +                                    1,     0,     0,     0,  11*0,
     O                                   -4,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          h(n,m) for 1975
      DATA (HY1D(I),I=3376,3520) /16*0,
     1    5675,  -2067,  -333,   191,    31,   -13,   -77,     6,   -21,
     +                                    1,     0,     0,     0,   3*0,
     2             -68,   262,  -265,   148,    99,   -26,   -16,    16,
     +                                    1,     0,     0,     0,   4*0,
     3                   -223,    39,  -152,    75,    -5,     4,     7,
     +                                    3,     0,     0,     0,   5*0,
     4                          -288,   -83,   -41,    10,   -19,    -4,
     +                                    4,     0,     0,     0,   6*0,
     5                                   88,    -4,    22,     6,    -5,
     +                                   -4,     0,     0,     0,   7*0,
     6                                          11,   -23,    18,    10,
     +                                   -1,     0,     0,     0,   8*0,
     7                                                -12,   -10,    11,
     +                                   -1,     0,     0,     0,   9*0,
     8                                                       -17,    -3,
     +                                    3,     0,     0,     0,  10*0,
     9                                                                1/
      DATA (HY1D(I),I=3521,3600) /
     +                                    1,     0,     0,     0,  11*0,
     O                                   -5,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          h(n,m) for 1980
      DATA (HY1D(I),I=3601,3745) /16*0,
     1    5604,  -2129,  -336,   212,    46,   -15,   -82,     7,   -21,
     +                                    1,     0,     0,     0,   3*0,
     2            -200,   271,  -257,   150,    93,   -27,   -18,    16,
     +                                    0,     0,     0,     0,   4*0,
     3                   -252,    53,  -151,    71,    -5,     4,     9,
     +                                    3,     0,     0,     0,   5*0,
     4                          -297,   -78,   -43,    16,   -22,    -5,
     +                                    6,     0,     0,     0,   6*0,
     5                                   92,    -2,    18,     9,    -6,
     +                                   -4,     0,     0,     0,   7*0,
     6                                          17,   -23,    16,     9,
     +                                    0,     0,     0,     0,   8*0,
     7                                                -10,   -13,    10,
     +                                   -1,     0,     0,     0,   9*0,
     8                                                       -15,    -6,
     +                                    4,     0,     0,     0,  10*0,
     9                                                                2/
      DATA (HY1D(I),I=3746,3825) /
     +                                    0,     0,     0,     0,  11*0,
     O                                   -6,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          h(n,m) for 1985
      DATA (HY1D(I),I=3826,3970) /16*0,
     1    5500,  -2197,  -310,   232,    47,   -16,   -83,     8,   -21,
     +                                    1,     0,     0,     0,   3*0,
     2            -306,   284,  -249,   150,    88,   -27,   -19,    15,
     +                                    0,     0,     0,     0,   4*0,
     3                   -297,    69,  -154,    69,    -2,     5,     9,
     +                                    3,     0,     0,     0,   5*0,
     4                          -297,   -75,   -48,    20,   -23,    -6,
     +                                    6,     0,     0,     0,   6*0,
     5                                   95,    -1,    17,    11,    -6,
     +                                   -4,     0,     0,     0,   7*0,
     6                                          21,   -23,    14,     9,
     +                                    0,     0,     0,     0,   8*0,
     7                                                 -7,   -15,     9,
     +                                   -1,     0,     0,     0,   9*0,
     8                                                       -11,    -7,
     +                                    4,     0,     0,     0,  10*0,
     9                                                                2/
      DATA (HY1D(I),I=3971,4050) /
     +                                    0,     0,     0,     0,  11*0,
     O                                   -6,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          h(n,m) for 1990
      DATA (HY1D(I),I=4051,4195) /16*0,
     1    5406,  -2279,  -284,   247,    46,   -16,   -80,    10,   -20,
     +                                    2,     0,     0,     0,   3*0,
     2            -373,   293,  -240,   154,    82,   -26,   -19,    15,
     +                                    1,     0,     0,     0,   4*0,
     3                   -352,    84,  -153,    69,     0,     6,    11,
     +                                    3,     0,     0,     0,   5*0,
     4                          -299,   -69,   -52,    21,   -22,    -7,
     +                                    6,     0,     0,     0,   6*0,
     5                                   97,     1,    17,    12,    -7,
     +                                   -4,     0,     0,     0,   7*0,
     6                                          24,   -23,    12,     9,
     +                                    0,     0,     0,     0,   8*0,
     7                                                 -4,   -16,     8,
     +                                   -2,     0,     0,     0,   9*0,
     8                                                       -10,    -7,
     +                                    3,     0,     0,     0,  10*0,
     9                                                                2/
      DATA (HY1D(I),I=4196,4275) /
     +                                   -1,     0,     0,     0,  11*0,
     O                                   -6,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          h(n,m) for 1995
      DATA (HY1D(I),I=4276,4420) /16*0,
     1    5306,  -2366,  -262,   262,    46,   -17,   -69,    11,   -20,
     +                                    1,     0,     0,     0,   3*0,
     2            -413,   302,  -236,   165,    72,   -25,   -21,    15,
     +                                    0,     0,     0,     0,   4*0,
     3                   -427,    97,  -143,    67,     4,     8,    12,
     +                                    4,     0,     0,     0,   5*0,
     4                          -306,   -55,   -58,    24,   -23,    -6,
     +                                    5,     0,     0,     0,   6*0,
     5                                  107,     1,    17,    15,    -8,
     +                                   -5,     0,     0,     0,   7*0,
     6                                          36,   -24,    11,     8,
     +                                   -1,     0,     0,     0,   8*0,
     7                                                 -6,   -16,     5,
     +                                   -2,     0,     0,     0,   9*0,
     8                                                        -4,    -8,
     +                                    1,     0,     0,     0,  10*0,
     9                                                                3/
      DATA (HY1D(I),I=4421,4500) /
     +                                   -2,     0,     0,     0,  11*0,
     O                                   -7,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          h(n,m) for 2000
      DATA (HY1D(I),I=4501,4645) /16*0,
     1  5186.1,-2481.6,-227.6, 272.6,  43.8, -17.4, -64.6,  11.9, -19.7,
     +                                  1.7,   0.1,  -0.4,  -0.9,   3*0,
     2          -458.0, 293.4,-231.9, 171.9,  63.7, -24.2, -21.5,  13.4,
     +                                  0.0,   1.3,   0.3,   0.2,   4*0,
     3                 -491.1, 119.8,-133.1,  65.1,   6.2,   8.5,  12.5,
     +                                  4.0,  -0.9,   2.5,   1.8,   5*0,
     4                        -303.8, -39.3, -61.2,  24.0, -21.5,  -6.2,
     +                                  4.9,  -2.6,  -2.6,  -0.4,   6*0,
     5                                106.3,   0.7,  14.8,  15.5,  -8.4,
     +                                 -5.9,   0.9,   0.7,  -1.0,   7*0,
     6                                        43.8, -25.4,   8.9,   8.4,
     +                                 -1.2,  -0.7,   0.3,  -0.1,   8*0,
     7                                               -5.8, -14.9,   3.8,
     +                                 -2.9,  -2.8,   0.0,   0.7,   9*0,
     8                                                      -2.1,  -8.2,
     +                                  0.2,  -0.9,   0.0,   0.3,  10*0,
     9                                                              4.8/
      DATA (HY1D(I),I=4646,4725) /
     +                                 -2.2,  -1.2,   0.3,   0.6,  11*0,
     O                                 -7.4,  -1.9,  -0.9,   0.3,  12*0,
     1                                        -0.9,  -0.4,  -0.2,  13*0,
     2                                                0.8,  -0.5,  14*0,
     3                                                      -0.9,  16*0/
C          h(n,m) for 2005
      DATA (HY1D(I),I=4726,4870) /16*0,
     1  5080.0,-2594.9,-200.4, 281.4,  42.7, -20.2, -61.4,  11.2, -20.1,
     +                                  2.4,   0.3,  -0.5,  -0.7,   3*0,
     2          -516.7, 269.3,-225.8, 179.8,  54.7, -22.5, -21.0,  12.9,
     +                                  0.2,   1.4,   0.3,   0.3,   4*0,
     3                 -524.5, 145.7,-123.0,  63.7,   6.9,   9.7,  12.7,
     +                                  4.4,  -0.7,   2.3,   1.7,   5*0,
     4                        -304.7, -19.5, -63.4,  25.4, -19.8,  -6.7,
     +                                  4.7,  -2.4,  -2.7,  -0.5,   6*0,
     5                                103.6,   0.0,  10.9,  16.1,  -8.1,
     +                                 -6.5,   0.9,   0.6,  -1.0,   7*0,
     6                                        50.3, -26.4,   7.7,   8.1,
     +                                 -1.0,  -0.6,   0.4,   0.0,   8*0,
     7                                               -4.8, -12.8,   2.9,
     +                                 -3.4,  -2.7,   0.0,   0.7,   9*0,
     8                                                      -0.1,  -7.9,
     +                                 -0.9,  -1.0,   0.0,   0.2,  10*0,
     9                                                              5.9/
      DATA (HY1D(I),I=4871,4950) /
     +                                 -2.3,  -1.5,   0.3,   0.6,  11*0,
     O                                 -8.0,  -2.0,  -0.8,   0.4,  12*0,
     1                                        -1.4,  -0.4,  -0.2,  13*0,
     2                                                1.0,  -0.5,  14*0,
     3                                                      -1.0,  16*0/
C          Secular variation rates are nominally okay through 2010
      DATA (GT1D(I),I=1,145) /0,
     O     8.8,  -15.0,  -0.3,  -2.5,  -2.6,  -0.8,  -0.4,  -0.2,     0,
     +                                    0,     0,     0,     0,   2*0,
     1    10.8,   -6.9,  -3.1,   2.8,   0.4,   0.2,   0.0,   0.2,     0,
     +                                    0,     0,     0,     0,   3*0,
     2            -1.0,  -0.9,  -7.1,  -3.0,  -0.2,  -0.2,  -0.2,     0,
     +                                    0,     0,     0,     0,   4*0,
     3                   -6.8,   5.9,  -1.2,   2.1,   1.1,   0.2,     0,
     +                                    0,     0,     0,     0,   5*0,
     4                          -3.2,   0.2,  -2.1,   0.6,  -0.2,     0,
     +                                    0,     0,     0,     0,   6*0,
     5                                 -0.6,  -0.4,   0.4,   0.2,     0,
     +                                    0,     0,     0,     0,   7*0,
     6                                         1.3,  -0.5,   0.5,     0,
     +                                    0,     0,     0,     0,   8*0,
     7                                                0.9,  -0.7,     0,
     +                                    0,     0,     0,     0,   9*0,
     8                                                       0.5,     0,
     +                                    0,     0,     0,     0,  10*0,
     9                                                                0/
      DATA (GT1D(I),I=146,225) /
     +                                    0,     0,     0,     0,  11*0,
     O                                    0,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
      DATA (HT1D(I),I=1,145) /16*0,
     1   -21.3,  -23.3,   5.4,   2.0,   0.1,  -0.4,   0.8,  -0.2,     0,
     +                                    0,     0,     0,     0,   3*0,
     2           -14.0,  -6.5,   1.8,   1.8,  -1.9,   0.4,   0.2,     0,
     +                                    0,     0,     0,     0,   4*0,
     3                   -2.0,   5.6,   2.0,  -0.4,   0.1,   0.2,     0,
     +                                    0,     0,     0,     0,   5*0,
     4                           0.0,   4.5,  -0.4,   0.2,   0.4,     0,
     +                                    0,     0,     0,     0,   6*0,
     5                                 -1.0,  -0.2,  -0.9,   0.2,     0,
     +                                    0,     0,     0,     0,   7*0,
     6                                         0.9,  -0.3,  -0.3,     0,
     +                                    0,     0,     0,     0,   8*0,
     7                                                0.3,   0.5,     0,
     +                                    0,     0,     0,     0,   9*0,
     8                                                       0.4,     0,
     +                                    0,     0,     0,     0,  10*0,
     9                                                                0/
      DATA (HT1D(I),I=146,225) /
     +                                    0,     0,     0,     0,  11*0,
     O                                    0,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
 
C          Do not need to load new coefficients if date has not changed
      ICHG = 0
      IF (DATE .EQ. DATEL) GO TO 300
      DATEL = DATE
      ICHG = 1
 
C          Trap out of range date:
      IF (DATE .LT. EPOCH(1)) GO TO 9100
      IF (DATE .GT. EPOCH(NEPO)+5.) WRITE(0,9200) DATE, EPOCH(NEPO) + 5.
 
      DO 100 I=1,NEPO
      IF (DATE .LT. EPOCH(I)) GO TO 110
      IY = I
  100 CONTINUE
  110 CONTINUE
 
      NMAX  = NMXE(IY)
      TIME  = DATE
      T     = TIME-EPOCH(IY)
      TO5   = T/5.
      IY1   = IY + 1
      GB(1) = 0.0
      GV(1) = 0.0
      I  = 2
      F0 = -1.0D-5
      DO 200 N=1,NMAX
      F0 = F0 * REAL(N)/2.
      F  = F0 / SQRT(2.0)
      NN = N+1
      MM = 1
      IF (IY .LT. NEPO) GB(I) = (GYR(NN,MM,IY) +                             ! interpolate (m=0 terms)
     +                          (GYR(NN,MM,IY1)-GYR(NN,MM,IY))*TO5) * F0
      IF (IY .EQ. NEPO) GB(I) = (GYR(NN,MM,IY) + GT(NN,MM)    *T  ) * F0     ! extrapolate (m=0 terms)
      GV(I) = GB(I) / REAL(NN)
      I = I+1
      DO 200 M=1,N
      F  = F / SQRT( REAL(N-M+1) / REAL(N+M) )
      NN = N+1
      MM = M+1
      I1 = I+1
      IF (IY .LT. NEPO) THEN                                                ! interpolate (m>0 terms)
	GB(I)  = (GYR(NN,MM,IY) +
     +           (GYR(NN,MM,IY1)-GYR(NN,MM,IY))*TO5) * F
	GB(I1) = (HYR(NN,MM,IY) +
     +           (HYR(NN,MM,IY1)-HYR(NN,MM,IY))*TO5) * F
      ELSE                                                                  ! extrapolate (m>0 terms)
	GB(I)  = (GYR(NN,MM,IY) +GT (NN,MM)    *T  ) * F
	GB(I1) = (HYR(NN,MM,IY) +HT (NN,MM)    *T  ) * F
      ENDIF
      RNN = REAL(NN)
      GV(I)  = GB(I)  / RNN
      GV(I1) = GB(I1) / RNN
  200 I = I+2
 
  300 CONTINUE

      RETURN
 
C          Error trap diagnostics:
 9100 WRITE (0,'(''COFRM:  DATE'',F9.3,'' preceeds earliest available ('
     +',F6.1,'')'')') DATE, EPOCH(1)
      CALL EXIT (1)
 9200 FORMAT('COFRM:  DATE',F9.3,' is after the last recommended for ext
     +rapolation (',F6.1,')')
      END
 
      SUBROUTINE DYPOL (COLAT,ELON,VP)
C          Computes parameters for dipole component of geomagnetic field.
C          COFRM must be called before calling DYPOL!
C          940504 A. D. Richmond
C
C          INPUT from COFRM through COMMON /MAGCOF/ NMAX,GB(255),GV(225),ICHG
C            NMAX = Maximum order of spherical harmonic coefficients used
C            GB   = Coefficients for magnetic field calculation
C            GV   = Coefficients for magnetic potential calculation
C            ICHG = Flag indicating when GB,GV have been changed
C
C          RETURNS:
C            COLAT = Geocentric colatitude of geomagnetic dipole north pole
C                    (deg)
C            ELON  = East longitude of geomagnetic dipole north pole (deg)
C            VP    = Magnitude, in T.m, of dipole component of magnetic
C                    potential at geomagnetic pole and geocentric radius
C                    of 6371.2 km
 
      PARAMETER (RTOD = 57.2957795130823, RE = 6371.2)
      COMMON /MAGCOF/ NMAX,GB(255),GV(225),ICHG
 
C          Compute geographic colatitude and longitude of the north pole of
C          earth centered dipole
      GPL   = SQRT (GB(2)**2 + GB(3)**2 + GB(4)**2)
      CTP   = GB(2) / GPL
      STP   = SQRT (1. - CTP*CTP)
      COLAT = ACOS (CTP) * RTOD
      ELON  = ATAN2 (GB(4),GB(3)) * RTOD
 
C          Compute magnitude of magnetic potential at pole, radius Re.
      VP = .2*GPL*RE
C          .2 = 2*(10**-4 T/gauss)*(1000 m/km) (2 comes through F0 in COFRM).
 
      RETURN
      END
 
      SUBROUTINE FELDG (IENTY,GLAT,GLON,ALT, BNRTH,BEAST,BDOWN,BABS)
C          Compute the DGRF/IGRF field components at the point GLAT,GLON,ALT.
C          COFRM must be called to establish coefficients for correct date
C          prior to calling FELDG.
C
C          IENTY is an input flag controlling the meaning and direction of the
C                remaining formal arguments:
C          IENTY = 1
C            INPUTS:
C              GLAT = Latitude of point (deg)
C              GLON = Longitude (east=+) of point (deg)
C              ALT  = Ht of point (km)
C            RETURNS:
C              BNRTH  north component of field vector (Gauss)
C              BEAST  east component of field vector  (Gauss)
C              BDOWN  downward component of field vector (Gauss)
C              BABS   magnitude of field vector (Gauss)
C
C          IENTY = 2
C            INPUTS:
C              GLAT = X coordinate (in units of earth radii 6371.2 km )
C              GLON = Y coordinate (in units of earth radii 6371.2 km )
C              ALT  = Z coordinate (in units of earth radii 6371.2 km )
C            RETURNS:
C              BNRTH = X component of field vector (Gauss)
C              BEAST = Y component of field vector (Gauss)
C              BDOWN = Z component of field vector (Gauss)
C              BABS  = Magnitude of field vector (Gauss)
C          IENTY = 3
C            INPUTS:
C              GLAT = X coordinate (in units of earth radii 6371.2 km )
C              GLON = Y coordinate (in units of earth radii 6371.2 km )
C              ALT  = Z coordinate (in units of earth radii 6371.2 km )
C            RETURNS:
C              BNRTH = Dummy variable
C              BEAST = Dummy variable
C              BDOWN = Dummy variable
C              BABS  = Magnetic potential (T.m)
C
C          INPUT from COFRM through COMMON /MAGCOF/ NMAX,GB(255),GV(225),ICHG
C            NMAX = Maximum order of spherical harmonic coefficients used
C            GB   = Coefficients for magnetic field calculation
C            GV   = Coefficients for magnetic potential calculation
C            ICHG = Flag indicating when GB,GV have been changed
C
C          HISTORY:
C          Apr 1983: written by Vincent B. Wickwar (Utah State Univ.).
C
C          May 1994 (A.D. Richmond): Added magnetic potential calculation
C
C          Oct 1995 (Barnes): Added ICHG
 
      PARAMETER (DTOR = 0.01745329251994330, RE = 6371.2)
      COMMON /MAGCOF/ NMAX,GB(255),GV(225),ICHG
      DIMENSION G(255), H(255), XI(3)
      SAVE IENTYP, G
      DATA IENTYP/-10000/
 
      IF (IENTY .EQ. 1) THEN
	IS   = 1
	RLAT = GLAT * DTOR
	CT   = SIN (RLAT)
	ST   = COS (RLAT)
	RLON = GLON * DTOR
	CP   = COS (RLON)
	SP   = SIN (RLON)
        CALL GD2CART (GLAT,GLON,ALT,XXX,YYY,ZZZ)
        XXX = XXX/RE
        YYY = YYY/RE
        ZZZ = ZZZ/RE
      ELSE
        IS   = 2
        XXX  = GLAT
        YYY  = GLON
        ZZZ  = ALT
      ENDIF
      RQ    = 1./(XXX**2+YYY**2+ZZZ**2)
      XI(1) = XXX*RQ
      XI(2) = YYY*RQ
      XI(3) = ZZZ*RQ
      IHMAX = NMAX*NMAX+1
      LAST  = IHMAX+NMAX+NMAX
      IMAX  = NMAX+NMAX-1
 
      IF (IENTY .NE. IENTYP .OR. ICHG .EQ. 1) THEN
        IENTYP = IENTY
	ICHG = 0
        IF (IENTY .NE. 3) THEN
	  DO 10 I=1,LAST
   10     G(I) = GB(I)
        ELSE
	  DO 20 I=1,LAST
   20     G(I) = GV(I)
        ENDIF
      ENDIF
 
      DO 30 I=IHMAX,LAST
   30 H(I) = G(I)

      MK = 3
      IF (IMAX .EQ. 1) MK=1

      DO 100 K=1,MK,2
      I  = IMAX
      IH = IHMAX

   60 IL = IH-I
      F = 2./FLOAT(I-K+2)
      X = XI(1)*F
      Y = XI(2)*F
      Z = XI(3)*(F+F)

      I = I-2
      IF (I .LT. 1) GO TO 90
      IF (I .EQ. 1) GO TO 80

      DO 70 M=3,I,2
      IHM = IH+M
      ILM = IL+M
      H(ILM+1) = G(ILM+1)+ Z*H(IHM+1) + X*(H(IHM+3)-H(IHM-1))
     +                                        -Y*(H(IHM+2)+H(IHM-2))
   70 H(ILM)   = G(ILM)  + Z*H(IHM)   + X*(H(IHM+2)-H(IHM-2))
     +                                        +Y*(H(IHM+3)+H(IHM-1))

   80 H(IL+2) = G(IL+2) + Z*H(IH+2) + X*H(IH+4) - Y*(H(IH+3)+H(IH))
      H(IL+1) = G(IL+1) + Z*H(IH+1) + Y*H(IH+4) + X*(H(IH+3)-H(IH))

   90 H(IL)   = G(IL)   + Z*H(IH)   + 2.*(X*H(IH+1)+Y*H(IH+2))
      IH = IL
      IF (I .GE. K) GO TO 60
  100 CONTINUE
 
      S = .5*H(1)+2.*(H(2)*XI(3)+H(3)*XI(1)+H(4)*XI(2))
      T = (RQ+RQ)*SQRT(RQ)
      BXXX = T*(H(3)-S*XXX)
      BYYY = T*(H(4)-S*YYY)
      BZZZ = T*(H(2)-S*ZZZ)
      BABS = SQRT(BXXX**2+BYYY**2+BZZZ**2)
      IF (IS .EQ. 1) THEN            ! (convert back to geodetic)
        BEAST = BYYY*CP-BXXX*SP
        BRHO  = BYYY*SP+BXXX*CP
        BNRTH =  BZZZ*ST-BRHO*CT
        BDOWN = -BZZZ*CT-BRHO*ST
      ELSEIF (IS .EQ. 2) THEN        ! (leave in earth centered cartesian)
        BNRTH = BXXX
        BEAST = BYYY
        BDOWN = BZZZ
      ENDIF
 
C          Magnetic potential computation makes use of the fact that the
C          calculation of V is identical to that for r*Br, if coefficients
C          in the latter calculation have been divided by (n+1) (coefficients
C          GV).  Factor .1 converts km to m and gauss to tesla.
      IF (IENTY.EQ.3) BABS = (BXXX*XXX + BYYY*YYY + BZZZ*ZZZ)*RE*.1
 
      RETURN
      END
 
      SUBROUTINE GD2CART (GDLAT,GLON,ALT,X,Y,Z)
C          Convert geodetic to cartesian coordinates by calling CONVRT
C          940503 A. D. Richmond
      PARAMETER (DTOR = 0.01745329251994330)
      CALL CONVRT (1,GDLAT,ALT,RHO,Z)
      ANG = GLON*DTOR
      X = RHO*COS(ANG)
      Y = RHO*SIN(ANG)
      RETURN
      END
 
      SUBROUTINE CONVRT (I,GDLAT,ALT,X1,X2)
C          Convert space point from geodetic to geocentric or vice versa.
C
C          I is an input flag controlling the meaning and direction of the
C            remaining formal arguments:
C
C          I = 1  (convert from geodetic to cylindrical geocentric)
C            INPUTS:
C              GDLAT = Geodetic latitude (deg)
C              ALT   = Altitude above reference ellipsoid (km)
C            RETURNS:
C              X1    = Distance from Earth's rotation axis (km)
C              X2    = Distance above (north of) Earth's equatorial plane (km)
C
C          I = 2  (convert from geodetic to spherical geocentric)
C            INPUTS:
C              GDLAT = Geodetic latitude (deg)
C              ALT   = Altitude above reference ellipsoid (km)
C            RETURNS:
C              X1    = Geocentric latitude (deg)
C              X2    = Geocentric distance (km)
C
C          I = 3  (convert from cylindrical geocentric to geodetic)
C            INPUTS:
C              X1    = Distance from Earth's rotation axis (km)
C              X2    = Distance from Earth's equatorial plane (km)
C            RETURNS:
C              GDLAT = Geodetic latitude (deg)
C              ALT   = Altitude above reference ellipsoid (km)
C
C          I = 4  (convert from spherical geocentric to geodetic)
C            INPUTS:
C              X1    = Geocentric latitude (deg)
C              X2    = Geocentric distance (km)
C            RETURNS:
C              GDLAT = Geodetic latitude (deg)
C              ALT   = Altitude above reference ellipsoid (km)
C
C
C          HISTORY:
C          940503 (A. D. Richmond):  Based on a routine originally written
C          by V. B. Wickwar.
C
C          Mar 2004: (Barnes) Revise spheroid definition to WGS-1984 to conform
C          with IGRF-9 release (EOS Volume 84 Number 46 November 18 2003).
C
C          REFERENCE: ASTRON. J. VOL. 66, p. 15-16, 1961
 
C          E2  = square of eccentricity of ellipse
C          REP = earth's polar      radius (km)
C          REQ = earth's equatorial radius (km)
      PARAMETER (RTOD = 57.2957795130823, DTOR = 0.01745329251994330,
     +           REP  = 6356.752, REQ = 6378.137, E2 = 1.-(REP/REQ)**2,
     +     E4 = E2*E2, E6 = E4*E2, E8 = E4*E4, OME2REQ = (1.-E2)*REQ,
     +     A21 =     (512.*E2 + 128.*E4 + 60.*E6 + 35.*E8)/1024. ,
     +     A22 =     (                        E6 +     E8)/  32. ,
     +     A23 = -3.*(                     4.*E6 +  3.*E8)/ 256. ,
     +     A41 =    -(           64.*E4 + 48.*E6 + 35.*E8)/1024. ,
     +     A42 =     (            4.*E4 +  2.*E6 +     E8)/  16. ,
     +     A43 =                                   15.*E8 / 256. ,
     +     A44 =                                      -E8 /  16. ,
     +     A61 =  3.*(                     4.*E6 +  5.*E8)/1024. ,
     +     A62 = -3.*(                        E6 +     E8)/  32. ,
     +     A63 = 35.*(                     4.*E6 +  3.*E8)/ 768. ,
     +     A81 =                                   -5.*E8 /2048. ,
     +     A82 =                                   64.*E8 /2048. ,
     +     A83 =                                 -252.*E8 /2048. ,
     +     A84 =                                  320.*E8 /2048. )
 
      IF (I .GE. 3) GO TO 300
 
C          Geodetic to geocentric
 
C          Compute RHO,Z
      SINLAT = SIN(GDLAT*DTOR)
      COSLAT = SQRT(1.-SINLAT*SINLAT)
      D      = SQRT(1.-E2*SINLAT*SINLAT)
      Z      = (ALT+OME2REQ/D)*SINLAT
      RHO    = (ALT+REQ/D)*COSLAT
      X1 = RHO
      X2 = Z
      IF (I .EQ. 1) RETURN
 
C          Compute GCLAT,RKM
      RKM   = SQRT(Z*Z + RHO*RHO)
      GCLAT = RTOD*ATAN2(Z,RHO)
      X1 = GCLAT
      X2 = RKM
      RETURN
 
C          Geocentric to geodetic
  300 IF (I .EQ. 3) THEN
         RHO = X1
         Z = X2
         RKM = SQRT(Z*Z+RHO*RHO)
         SCL = Z/RKM
         GCLAT = ASIN(SCL)*RTOD
      ELSEIF (I .EQ. 4) THEN
         GCLAT = X1
         RKM = X2
         SCL = SIN(GCLAT*DTOR)
      ELSE
         RETURN
      ENDIF
 
      RI = REQ/RKM
      A2 = RI*(A21+RI*(A22+RI* A23))
      A4 = RI*(A41+RI*(A42+RI*(A43+RI*A44)))
      A6 = RI*(A61+RI*(A62+RI* A63))
      A8 = RI*(A81+RI*(A82+RI*(A83+RI*A84)))
      CCL = SQRT(1.-SCL*SCL)
      S2CL = 2.*SCL*CCL
      C2CL = 2.*CCL*CCL-1.
      S4CL = 2.*S2CL*C2CL
      C4CL = 2.*C2CL*C2CL-1.
      S8CL = 2.*S4CL*C4CL
      S6CL = S2CL*C4CL+C2CL*S4CL
      DLTCL = S2CL*A2+S4CL*A4+S6CL*A6+S8CL*A8
      GDLAT = DLTCL*RTOD+GCLAT
      SGL = SIN(GDLAT*DTOR)
      ALT = RKM*COS(DLTCL)-REQ*SQRT(1.-E2*SGL*SGL)
      RETURN
      END

