      SUBROUTINE APEX (DATE,DLAT,DLON,ALT,
     +                A,ALAT,ALON,BMAG,XMAG,YMAG,ZMAG,V)
C Calculate apex radius, latitude, longitude; and magnetic field and potential
C 940822 A. D. Richmond
C
C          INPUTS:
C            DATE = Year and fraction (1990.0 = 1990 January 1, 0 UT)
C            DLAT = Geodetic latitude in degrees
C            DLON = Geodetic longitude in degrees
C            ALT = Altitude in km
C
C          RETURNS:
C            A    = (Apex height + REQ)/REQ, where REQ = equatorial Earth
C                    radius.  A is analogous to the L value in invariant
C                    coordinates.
C            ALAT = Apex latitude in degrees (negative in S. magnetic hemisphere)
C            ALON = Apex longitude (Geomagnetic longitude of apex)
C            BMAG = Magnitude of geomagnetic field in nT
C            XMAG,YMAG,ZMAG = North, east, and downward geomagnetic field components in nT
C            V    = Magnetic potential, in T.m
C
C     COMMON /DIPOLE/ COLAT,ELON,VP,CTP,STP
C          COLAT = Geocentric colatitude of geomagnetic dipole north pole (deg)
C          ELON  = East longitude of geomagnetic dipole north pole (deg)
C          VP    = Magnitude, in T.m, of dipole component of magnetic potential at
C                  geomagnetic pole and geocentric radius of 6371.2 km
C          CTP,STP = cosine, sine of COLAT
C
C          MODIFICATIONS:
C          May 1999:  Revise DS calculation in LINAPX to avoid divide by zero.

      PARAMETER (RTOD=5.72957795130823E1, DTOR=1.745329251994330E-2,
     +           RE=6371.2, REQ=6378.160)
      COMMON /DIPOLE/ COLAT,ELON,VP,CTP,STP

      CALL COFRM (DATE)
      CALL DYPOL (CLATP,POLON,VPOL)
      COLAT = CLATP
      CTP   = COS(CLATP*DTOR)
      STP   = SQRT(1. - CTP*CTP)
      ELON  = POLON
      VP    = VPOL
      CALL LINAPX (DLAT,DLON,ALT,A,ALAT,ALON,XMAG,YMAG,ZMAG,BMAG)
      XMAG = XMAG*1.E5
      YMAG = YMAG*1.E5
      ZMAG = ZMAG*1.E5
      BMAG = BMAG*1.E5
      CALL GD2CART (DLAT,DLON,ALT,X,Y,Z)
      CALL FELDG (3,X/RE,Y/RE,Z/RE,BX,BY,BZ,V)
      RETURN
      END

       SUBROUTINE LINAPX(GDLAT,GLON,ALT
     2 ,A,ALAT,ALON,XMAG,YMAG,ZMAG,F)
C***BEGIN PROLOGUE  LINAPX                                                      
C***DATE WRITTEN   731029   (YYMMDD)                                            
C***AUTHOR  CLARK, W., N.O.A.A. ERL LAB.                                        
C***REVISION DATE  880201   (YYMMDD) Harsh Anand Passi, NCAR
C***LATEST REVISION 940803 A. D. Richmond
C***PURPOSE  Transforms the geographic coordinates to apex coordinates.         
C***DESCRIPTION                                                                 
C     The method used is as follow:                                             
C       1. Calculates step size as a function of the geomagnetic                
C          dipole coordinates of the starting point.                            
C       2. Determine direction of trace                                         
C       3. Convert the geodetic coordinates of the starting point               
C          to the cartesian coordinates for tracing.                            
C       Loop:                                                                   
C       i)   Increment step count, if count > 200,                              
C            assume it is dipole field, call DIPAPX to                          
C            determine Apex coordinates else continue.                          
C       ii)  Get field components in X, Y and Z direction                       
C       iii) Call ITRACE to trace field line.                                   
C       iv)  Test if Apex passed call FNDAPX to determine Apex coordinates      
C            else loop:                                                         
C                                                                               
C   INPUT                                                                       
C     GDLAT  Latitude of starting point (Deg)                                   
C     GLON   Longitude (East=+) of starting point (Deg)                         
C     ALT    Ht. of starting point (Km)                                         
C                                                                               
C  OUTPUT
C          A  (Apex height + REQ)/REQ, where REQ = equatorial Earth radius.
C             A is analogous to the L value in invariant coordinates.
C       ALAT  Apex Lat. (deg)                                                   
C       ALON  Apex Lon. (deg)                                                   
C       XMAG  North component of magnetic field at starting point
C       YMAG  East component of magnetic field at starting point
C       ZMAG  Down component of magnetic field at starting point
C          F  Magnetic field magnitude at starting point
C
C     COMMON Blocks Used                                                        
C     /APXIN/ YAPX(3,3)                                                         
C       YAPX    Matrix of cartesian coordinates (loaded columnwise)             
C               of the 3 points about APEX. Set in subprogram ITRACE.           
C                                                                               
C   /DIPOLE/COLAT,ELON,VP,CTP,STP                           
C     COLAT   Geographic colatitude of the north pole of the                    
C             earth-centered dipole (Deg).                                      
C     ELON    Geographic longitude of the north pole of the                     
C             earth-centered dipole (Deg).                                      
C     VP      Magnetic potential magnitude at geomagnetic pole (T.m)
C     CTP     cos(COLAT*DTOR)
C     STP     sin(COLAT*DTOR)
C                                                                               
C     /ITRA/ NSTP, Y(3), YOLD(3), SGN, DS                                       
C     NSTP      Step count. Incremented in subprogram LINAPX.                   
C     Y         Array containing current tracing point cartesian coordinates.   
C     YOLD      Array containing previous tracing point cartesian coordinates.
C     SGN       Determines direction of trace. Set in subprogram LINAPX         
C     DS        Step size (Km) Computed in subprogram LINAPX.                   
C                                                                               
C     /FLDCOMD/ BX, BY, BZ, BB                                                  
C     BX        X comp. of field vector at the current tracing point (Gauss)    
C     BY        Y comp. of field vector at the current tracing point (Gauss)    
C     BZ        Z comp. of field vector at the current tracing point (Gauss)    
C     BB        Magnitude of field vector at the current tracing point (Gauss)  
C                                                                               
C***REFERENCES  Stassinopoulos E. G. , Mead Gilbert D., X-841-72-17             
C                 (1971) GSFC, Greenbelt, Maryland                              
C                                                                               
C***ROUTINES CALLED  GD2CART,CONVRT,FELDG,
C                    ITRACE,DIPAPX,FNDAPX                                       
C                                                                               
C***END PROLOGUE  LINAPX                                                        
      PARAMETER (MAXS = 200,
     +           RTOD=5.72957795130823E1, DTOR=1.745329251994330E-2,
     +           RE=6371.2, REQ=6378.160)
      COMMON /FLDCOMD/ BX, BY, BZ, BB
      COMMON /APXIN/   YAPX(3,3)
      COMMON /DIPOLE/  COLAT,ELON,VP,CTP,STP
      COMMON /ITRA/    NSTP, Y(3), YP(3),  SGN, DS

C          Determine step size as a function of geomagnetic dipole
C          coordinates of the starting point
      CALL CONVRT (2,GDLAT,ALT,GCLAT,R)
C     SINGML = .98*SIN(GCLAT*DTOR) + .199*COS(GCLAT*DTOR)*
      SINGML = CTP*SIN(GCLAT*DTOR) +  STP*COS(GCLAT*DTOR)*
     +                                             COS((GLON-ELON)*DTOR)

C          May 99: avoid possible divide by zero (when SINGML = 1.)
      CGML2 = AMAX1 (0.25,1.-SINGML*SINGML)
      DS = .06*R/CGML2 - 370.
Cold:      Limit DS to its value at 60 deg gm latitude.
Cold: DS = .06*R/(1.-SINGML*SINGML) - 370.
Cold: IF (DS .GT. 1186.) DS = 1186.

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
   10 IENTRY=2

      CALL FELDG (IENTRY,Y(1)/RE,Y(2)/RE,Y(3)/RE,BX,BY,BZ,BB)

      NSTP = NSTP + 1

C          Quit if too many steps
      IF (NSTP .GE. MAXS) THEN                                                 
        RHO = SQRT(Y(1)*Y(1) + Y(2)*Y(2))
        CALL CONVRT(3,XLAT,HT,RHO,Y(3))
        XLON = RTOD*ATAN2(Y(2),Y(1))
	CALL FELDG (1,XLAT,XLON,HT,BNRTH,BEAST,BDOWN,BABS)
	CALL DIPAPX (XLAT,XLON,HT,BNRTH,BEAST,BDOWN,A,ALON)
        ALAT = -SGN*RTOD*ACOS(SQRT(1./A))
        RETURN                                                                
      END IF                                                                   

C          Find next point using adams algorithm after 7 points
       CALL ITRACE (IAPX)
       IF (IAPX .EQ. 1) GO TO 10

C          Maximum radius just passed.  Find apex coords
      CALL FNDAPX (ALT,ZMAG,A,ALAT,ALON)
      RETURN
      END

      SUBROUTINE ITRACE (IAPX)
      SAVE
C***BEGIN PROLOGUE  ITRACE                                                      
C***DATE WRITTEN   731029   (YYMMDD)                                            
C***AUTHOR  CLARK, W. N.O.A.A. ERL LAB.                                         
C***REVISION DATE  880201, H. Passi 
C***PURPOSE  Field line integration routine.                                    
C***DESCRIPTION                                                                 
C     It uses 4-point ADAMS formula after initialization.                       
C     First 7 iterations advance point by 3 steps.                              
C                                                                               
C   INPUT                                                                       
C            Passed in through the common blocks ITRA, FLDCOMD.                 
C            See the description below.                                         
C   OUTPUT                                                                      
C    IAPX    Flag set to 2 when APEX passed, otherwise set to 1.                
C                                                                               
C            Passed out through the common block APXIN.                         
C            See the description below.                                         
C                                                                               
C     COMMON Blocks Used                                                        
C     /APXIN/ YAPX(3,3)                                                         
C     YAPX      Matrix of cartesian coordinates (loaded columnwise)             
C               of the 3 points about APEX. Set in subprogram ITRACE.           
C                                                                               
C     /FLDCOMD/ BX, BY, BZ, BB                                                  
C     BX        X comp. of field vector at the current tracing point (Gauss)    
C     BY        Y comp. of field vector at the current tracing point (Gauss)    
C     BZ        Z comp. of field vector at the current tracing point (Gauss)    
C     BB        Magnitude of field vector at the current tracing point (Gauss)  
C                                                                               
C     /ITRA/ NSTP, Y(3), YOLD(3), SGN, DS                                       
C     NSTP      Step count for line integration.                                
C               Incremented in subprogram LINAPX.                               
C     Y         Array containing current tracing point cartesian coordinates.   
C     YOLD      Array containing previous tracing point cartesian coordinates.  
C     SGN       Determines direction of trace. Set in subprogram LINAPX         
C     DS        Integration step size (arc length Km)                           
C               Computed in subprogram LINAPX.                                  
C                                                                               
C***REFERENCES  reference 1                                                     
C                                                                               
C***ROUTINES CALLED  None                                                       
C***COMMON BLOCKS    APXIN,FLDCOMD,ITRA                                         
C***END PROLOGUE  ITRACE                                                        
C **                                                                            
      COMMON /ITRA/    NSTP, Y(3), YOLD(3), SGN, DS
      COMMON /FLDCOMD/ BX, BY, BZ, BB
      COMMON /APXIN/   YAPX(3,3)
      DIMENSION  YP(3, 4)

C          Statement function
      RDUS(D,E,F) = SQRT( D**2 + E**2 + F**2 )
C                                                                               
      IAPX = 1
C          Field line is defined by the following differential equations
C          in cartesian coordinates:
      YP(1, 4) = SGN*BX/BB
      YP(2, 4) = SGN*BY/BB
      YP(3, 4) = SGN*BZ/BB
      IF (NSTP .GT. 7) GO TO 90

C          First seven steps use this block
      DO 80 I = 1, 3
      GO TO (10, 20, 30, 40, 50, 60, 70) NSTP
   10 D2  = DS/2.
      D6  = DS/6.
      D12 = DS/12.
      D24 = DS/24.
      YP(I,1)   = YP(I,4)
      YOLD(I)   = Y(I)
      YAPX(I,1) = Y(I)
      Y(I) = YOLD(I) + DS*YP(I, 1)
      GO TO 80

   20 YP(I, 2) = YP(I,4)
      Y(I) = YOLD(I) + D2*(YP(I,2)+YP(I,1))
      GO TO 80

   30 Y(I) = YOLD(I) + D6*(2.*YP(I,4)+YP(I,2)+3.*YP(I,1))
      GO TO 80

   40 YP(I,2)  = YP(I,4)
      YAPX(I,2)= Y(I)
      YOLD(I)  = Y(I)
      Y(I)     = YOLD(I) + D2*(3.*YP(I,2)-YP(I,1))
      GO TO 80

   50 Y(I) = YOLD(I) + D12*(5.*YP(I,4)+8.*YP(I,2)-YP(I,1))
      GO TO 80

   60 YP(I,3)  = YP(I,4)
      YOLD(I)  = Y(I)
      YAPX(I,3)= Y(I)
      Y(I)     = YOLD(I) + D12*(23.*YP(I,3)-16.*YP(I,2)+5.*YP(I,1))
      GO TO 80

   70 YAPX(I, 1) = YAPX(I, 2)
      YAPX(I, 2) = YAPX(I, 3)
      Y(I) = YOLD(I) + D24*(9.*YP(I,4)+19.*YP(I,3)-5.*YP(I,2)+YP(I,1))
      YAPX(I, 3) = Y(I)
   80 CONTINUE

C          Signal if apex passed
      IF ( NSTP .EQ. 6 .OR. NSTP .EQ. 7) THEN
	RC = RDUS( YAPX(1,3), YAPX(2,3), YAPX(3,3))
        RP = RDUS( YAPX(1,2), YAPX(2,2), YAPX(3,2))                             
	IF (RC .LT. RP) IAPX=2
      ENDIF
      RETURN

C          Stepping block for NSTP .gt. 7
   90 DO 110 I = 1, 3                                                           
      YAPX(I, 1) = YAPX(I, 2)
      YAPX(I, 2) = Y(I)
      YOLD(I) = Y(I)
      TERM = 55.*YP(I, 4) - 59.*YP(I, 3) + 37.*YP(I, 2) - 9.*YP(I, 1)
      Y(I) = YOLD(I) + D24*TERM
      YAPX(I, 3) = Y(I)

      DO 100 J = 1, 3
      YP(I, J) = YP(I, J+1)
  100 CONTINUE
  110 CONTINUE                                                                  

      RC = RDUS (   Y(1),    Y(2),    Y(3))
      RP = RDUS (YOLD(1), YOLD(2), YOLD(3))
      IF (RC .LT. RP) IAPX=2

      RETURN
      END

       SUBROUTINE FNDAPX (ALT,ZMAG,A,ALAT,ALON)
C***BEGIN PROLOGUE  FNDAPX                                                      
C***DATE WRITTEN   731023   (YYMMDD)                                            
C***AUTHOR  CLARK, W., NOAA BOULDER                                             
C***REVISION DATE  940803, A. D. Richmond, NCAR 
C***PURPOSE  Finds apex coords once tracing has signalled that the apex         
C            has been passed.  
C***DESCRIPTION                                                                 
C                                                                               
C     It uses second degree interpolation routine, FINT, to find                
C     apex latitude and apex longtitude.                                        
C   INPUT                                                                       
C     ALT    Altitude of starting point
C     ZMAG   Downward component of geomagnetic field at starting point
C     NMAX   Order of IGRF coefficients being used
C     G      Array of coefficients from COFRM
C   OUTPUT                                                                      
C          A  (Apex height + REQ)/REQ, where REQ = equatorial Earth radius.
C             A is analogous to the L value in invariant coordinates.
C       ALAT  Apex Lat. (deg)                                                   
C       ALON  Apex Lon. (deg)                                                   
C                                                                               
C***LONG DESCRIPTION                                                            
C                                                                               
C     COMMON Blocks Used                                                        
C     /APXIN/ YAPX(3,3)                                                         
C     YAPX      Matrix of cartesian coordinates (loaded columnwise)             
C               of the 3 points about APEX. Set in subprogram ITRACE.           
C                                                                               
C   /DIPOLE/COLAT,ELON,VP,CTP,STP                           
C     COLAT   Geocentric colatitude of the north pole of the                    
C             earth-centered dipole (Deg).                                      
C     ELON    Geographic longitude of the north pole of the                     
C             earth-centered dipole (Deg).                                      
C     CTP     cos(COLAT*DTOR)
C     STP     sin(COLAT*DTOR)
C                                                                               
C***ROUTINES CALLED  FINT                                                       
C***COMMON BLOCKS    APXIN,DIPOLE                                         
C***END PROLOGUE  FNDAPX                                                        
C                                                                               
      PARAMETER (RTOD=5.72957795130823E1,DTOR=1.745329251994330E-2)
      PARAMETER (RE=6371.2,REQ=6378.160)
      COMMON /APXIN/  YAPX(3,3)
      COMMON /DIPOLE/ COLAT,ELON,VP,CTP,STP
      DIMENSION Z(3), HT(3), Y(3)

C       ****
C       ****     GET GEODETIC FIELD COMPONENTS
C       ****
      IENTRY = 1
      DO 2 I = 1,3
	RHO = SQRT(YAPX(1,I)**2+YAPX(2,I)**2)
	CALL CONVRT (3,GDLT,HT(I),RHO,YAPX(3,I))
	GDLN = RTOD*ATAN2 (YAPX(2,I),YAPX(1,I))
	CALL FELDG (IENTRY,GDLT,GDLN,HT(I),X,YDUM,Z(I),F)
    2 CONTINUE
C     ****
C     ****     FIND CARTESIAN COORDINATES AT DIP EQUATOR BY INTERPOLATION
C     ****
      DO 3 I = 1,3
	CALL FINT(Z(1),Z(2),Z(3),YAPX(I,1),YAPX(I,2),YAPX(I,3),0.,Y(I))
    3 CONTINUE
C     ****
C     ****     FIND APEX HEIGHT BY INTERPOLATION
C     ****
      CALL FINT(Z(1),Z(2),Z(3),HT(1),HT(2),HT(3),0.,XINTER)
C          Ensure that XINTER is not less than original starting altitude (ALT)
      XINTER = AMAX1(ALT,XINTER)
      A = (REQ+XINTER)/(REQ)
C     ****
C     ****     FIND APEX COORDINATES , GIVING ALAT SIGN OF DIP AT
C     ****       STARTING POINT.  ALON IS THE VALUE OF THE GEOMAGNETIC
C     ****       LONGITUDE AT THE APEX.
C     ****
      IF (A.LT.1.) THEN
	WRITE(6,20)
  20    FORMAT (' BOMBED! THIS MAKES A LESS THAN ONE')
	STOP
      ENDIF
      RASQ = RTOD*ACOS(SQRT(1./A))
      ALAT = SIGN(RASQ,ZMAG)

C Algorithm for ALON:
C   Use spherical coordinates.
C   Let GP be geographic pole.
C   Let GM be geomagnetic pole (colatitude COLAT, east longitude ELON).
C   Let XLON be longitude of apex.
C   Let TE be colatitude of apex.
C   Let ANG be longitude angle from GM to apex.
C   Let TP be colatitude of GM.
C   Let TF be arc length between GM and apex.
C   Let PA = ALON be geomagnetic longitude, i.e., Pi minus angle measured 
C     counterclockwise from arc GM-apex to arc GM-GP.
C   Then, using notation C=cos, S=sin, spherical-trigonometry formulas 
C     for the functions of the angles are as shown below.  Note: STFCPA,
C     STFSPA are sin(TF) times cos(PA), sin(PA), respectively.

      XLON = ATAN2(Y(2),Y(1))
      ANG  = XLON-ELON*DTOR
      CANG = COS(ANG)
      SANG = SIN(ANG)
      R    = SQRT(Y(1)**2+Y(2)**2+Y(3)**2)
      CTE  = Y(3)/R
      STE  = SQRT(1.-CTE*CTE)
      STFCPA = STE*CTP*CANG - CTE*STP
      STFSPA = SANG*STE
      ALON = ATAN2(STFSPA,STFCPA)*RTOD
      RETURN
      END                                                                      

      SUBROUTINE DIPAPX(GDLAT,GDLON,ALT,BNORTH,BEAST,BDOWN,A,ALON)
C Compute a, alon from local magnetic field using dipole and spherical approx.
C 940501 A. D. Richmond
C Input:
C   GDLAT  = geodetic latitude, degrees
C   GDLON  = geodetic longitude, degrees
C   ALT    = altitude, km
C   BNORTH = geodetic northward magnetic field component (any units)
C   BEAST  = eastward magnetic field component
C   BDOWN  = geodetic downward magnetic field component
C Output:
C   A      = apex radius, 1 + h_A/R_eq
C   ALON   = apex longitude, degrees
C
C Algorithm:
C   Use spherical coordinates.
C   Let GP be geographic pole.
C   Let GM be geomagnetic pole (colatitude COLAT, east longitude ELON).
C   Let G be point at GDLAT,GDLON.
C   Let E be point on sphere below apex of dipolar field line passing through G.
C   Let TD be dipole colatitude of point G, found by applying dipole formula
C     for dip angle to actual dip angle.
C   Let B be Pi plus local declination angle.  B is in the direction 
C     from G to E.
C   Let TG be colatitude of G.
C   Let ANG be longitude angle from GM to G.
C   Let TE be colatitude of E.
C   Let TP be colatitude of GM.
C   Let A be longitude angle from G to E.
C   Let APANG = A + ANG
C   Let PA be geomagnetic longitude, i.e., Pi minus angle measured 
C     counterclockwise from arc GM-E to arc GM-GP.
C   Let TF be arc length between GM and E.
C   Then, using notation C=cos, S=sin, COT=cot, spherical-trigonometry formulas 
C     for the functions of the angles are as shown below.  Note: STFCPA,
C     STFSPA are sin(TF) times cos(PA), sin(PA), respectively.

      PARAMETER (RTOD=5.72957795130823E1,DTOR=1.745329251994330E-2)
      PARAMETER (RE=6371.2,REQ=6378.160)
      COMMON/DIPOLE/COLAT,ELON,VP,CTP,STP
      BHOR = SQRT(BNORTH*BNORTH + BEAST*BEAST)
      IF (BHOR.EQ.0.) THEN
	ALON = 0.
	A = 1.E34
	RETURN
      ENDIF
      COTTD = BDOWN*.5/BHOR
      STD = 1./SQRT(1.+COTTD*COTTD)
      CTD = COTTD*STD
      SB = -BEAST/BHOR
      CB = -BNORTH/BHOR
      CTG = SIN(GDLAT*DTOR) 
      STG = COS(GDLAT*DTOR)
      ANG = (GDLON-ELON)*DTOR 
      SANG = SIN(ANG)
      CANG = COS(ANG)
      CTE = CTG*STD + STG*CTD*CB
      STE = SQRT(1. - CTE*CTE)
      SA = SB*CTD/STE
      CA = (STD*STG - CTD*CTG*CB)/STE
      CAPANG = CA*CANG - SA*SANG
      SAPANG = CA*SANG + SA*CANG
      STFCPA = STE*CTP*CAPANG - CTE*STP
      STFSPA = SAPANG*STE
      ALON = ATAN2(STFSPA,STFCPA)*RTOD
      R = ALT + RE
      HA = ALT + R*COTTD*COTTD
      A = 1. + HA/REQ
      RETURN
      END

       SUBROUTINE FINT (A1, A2, A3, A4, A5, A6, A7, RESULT)
C***PURPOSE  Second degree interpolation routine                                
C***REFER TO  FNDAPX                                                            
       RESULT = ((A2-A3)*(A7-A2)*(A7-A3)*A4-(A1-A3)*(A7-A1)*(A7-A3)*A5+         
     +  (A1-A2)*(A7-A1)*(A7-A2)*A6)/((A1-A2)*(A1-A3)*(A2-A3))
       RETURN                                                                   
       END                                                                      
