      subroutine divve(glat,glon,alt,HR,WK,ed1,ed2,divvem)
C          INPUTS:
C            GLAT = Geographic (geodetic) latitude, degrees
C            GLON = Geographic (geodetic) longitude, degrees
C            ALT  = Altitude, km
C            HR   = Reference altitude, km 
C            WK   = Work array used internally, i.e., there is no need
C                     to access the contents.  WK should not be altered
C                     between initialization (by APXRDA) and use (by 
C                     APXMALL or APXQ2G).
C            ed1  = E_(d1), as defined by Richmond (1995), in V/m 
C            ed2  = E_(d2), as defined by Richmond (1995), in V/m 
C          RETURNS:
C            divvem = divergence of ExB/B^2 velocity, in s^{-1}
C Reference:  Richmond, A. D., Ionospheric Electrodynamics Using
C Magnetic Apex Coordinates, J. Geomag. Geoelectr., 47, 191-212, 1995.
      PARAMETER (RTOD=57.2957795130823, DTOR=0.01745329251994330,
     +           RE=6371.2, REQ=6378.160)
C WK is an array; divve doesn't need to know its dimension, so set to 1.
      DIMENSION WK(1)
C          Dimensions of non-scalar arguments returned by APXMALL:
      DIMENSION B(3),BHAT(3),
     +          D1(3),D2(3),D3(3), E1(3),E2(3),E3(3), F1(2),F2(2)
      dimension grdlbm2(3)
C Calculate geocentric cartesian coordinates x,y,z (in km)
      call GD2CART(GLAT,GLON,ALT,X,Y,Z)
C Calculate gradient of ALOG(B0**(-2)) [= -2*ALOG(B0)] in Cartesian 
C   coordinates, in units of m^{-1}
      call FELDG (2,(X+10.)/RE,Y/RE,Z/RE,BNRTH,BEAST,BDOWN,BPX)
      call FELDG (2,(X-10.)/RE,Y/RE,Z/RE,BNRTH,BEAST,BDOWN,BMX)
      call FELDG (2,X/RE,(Y+10.)/RE,Z/RE,BNRTH,BEAST,BDOWN,BPY)
      call FELDG (2,X/RE,(Y-10.)/RE,Z/RE,BNRTH,BEAST,BDOWN,BMY)
      call FELDG (2,X/RE,Y/RE,(Z+10.)/RE,BNRTH,BEAST,BDOWN,BPZ)
      call FELDG (2,X/RE,Y/RE,(Z-10.)/RE,BNRTH,BEAST,BDOWN,BMZ)
      dlbm2dx = (alog(BMX) - alog(BPX))/1.E4
      dlbm2dy = (alog(BMY) - alog(BPY))/1.E4
      dlbm2dz = (alog(BMZ) - alog(BPZ))/1.E4
C Rotate gradient to local (east,north,up) coordinates
C (1=east, 2=north, 3=up)
      coslat = cos(glat*DTOR)
      sinlat = sin(glat*DTOR)
      coslon = cos(glon*DTOR)
      sinlon = sin(glon*DTOR)
      grdlbm2(1) = -dlbm2dx*sinlon + dlbm2dy*coslon
      grdlbm2(2) = -(dlbm2dx*coslon + dlbm2dy*sinlon)*sinlat
     1  + dlbm2dz*coslat
      grdlbm2(3) = (dlbm2dx*coslon + dlbm2dy*sinlon)*coslat
     1  + dlbm2dz*sinlat
      CALL APXMALL (GLAT,GLON,ALT,HR, WK,
     +              B,BHAT,BABS,SI, QDLON,
     +              XLATM,VMP,W,D,BE3,SIM,D1,D2,D3,E1,E2,E3,
     +              QDLAT,F,F1,F2, IST)
C BE3 is in nT; need to convert to T
      ve1 =  ed2/(BE3*1.e-9)
      ve2 = -ed1/(BE3*1.e-9)
C divvem is dot product of ExB/B^2 velocity [ve1*E1(i) + ve2*E2(i)]
C   with gradient of ALOG(B0**(-2)) 
      divvem = 0.
      do i=1,3
	divvem = divvem + (ve1*E1(i) + ve2*E2(i))*grdlbm2(i)
      enddo
      return
      end
