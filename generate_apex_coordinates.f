      program generate_apex_coordinates
c
CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG
CG                                                CG
CG    GENERATES A 'GLOBAL' FLUX_TUBE IONOSPHERE   CG
CG          USING APEX COORDINATES                CG
CG                                                CG
CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG
c
      implicit none
      integer npts, nmp, nlp
      include 'npts.h'
      CHARACTER*6 date_string
      CHARACTER*33 output_file
      real l_value,re(nmp,nlp)
      real vmp_tube(npts,nmp,nlp)
      real vmp_apex(nmp,nlp)
      real*8 re_dble(nmp,nlp)
      real*8 q_dble(npts,nmp,nlp)
      real*8 blon_dble(nmp,nlp)
      real*8 gr_dble(npts,nmp,nlp)
      real*8 gcol_dble(npts,nmp,nlp)
      real*8 glon_dble(npts,nmp,nlp)
      real*8 bcol_dble(npts,nmp,nlp)
      real*8 vmp_tube_dble(npts,nmp,nlp)
      real*8 vmp_apex_dble(nmp,nlp)
      real*8 vmp_south_pole_dble
      real*8 B_magnitude_dble(npts,nmp,nlp)
      real*8 D1_dble(3,npts,nmp,nlp)
      real*8 D2_dble(3,npts,nmp,nlp)
      real*8 D3_dble(3,npts,nmp,nlp)
      real*8 E1_dble(3,npts,nmp,nlp)
      real*8 E2_dble(3,npts,nmp,nlp)
      real*8 E3_dble(3,npts,nmp,nlp)
      real*8 BE3_dble(npts,nmp,nlp)
      real*8 lat_base_90
      real*8 apex_height
      real*8 apex_radius
      real*8 dummy
      real*8 Apex_L_value(200) 
      real*8 Largest_Apex_L_value
      real*8 Smallest_Apex_L_value
      real*8 tiegcm_L_value(nlp)
      real*8 tiegcm_L_value2(nlp)
      real*8 date_dble
      real fix_h2(npts),this_q(npts)
      real glat_fix(npts),glon_fix(npts)
      real gr(npts,nmp,nlp),gcol(npts,nmp,nlp),glon(npts,nmp,nlp)
      real bcol(npts,nmp,nlp),blon(nmp,nlp)
      real r0,dtr
      real alt
      real ha
      real hr
      real alon
      real blonb
      real qdlat
      integer in(nmp,nlp),is(nmp,nlp),n_mid_point,mp,lp
      integer icount
      integer mgtype
      integer ndx
      integer i_top_point
      integer iht
      integer j_direction
      integer i
      integer i_first_index_tube_less_than_L4
      integer idummy
      integer i_smallest
      integer i_largest
      integer iNumber_of_tubes
      integer i_num_tubes_for_tiegcm
      integer i_doubled_tube
      integer i_new_number_of_tubes
      integer istop
      integer ndate
      integer ist
      real babs, si
      real xlonm , xlatm
      real vmp
      real w, d, be3, sim, f
      real pi
      real date
      real xlatqd
      real height
      real GLAT_north90
      real GLAT_south90
      real GLON_north90
      real GLON_south90
      real GLAT_sp
      real GLON_sp
      real GLAT_apex
      real GLON_apex
      real vmp_north
      real vmp_south
      real vmp_south_pole
      real q(npts,nmp,nlp)
      real rtod, dtor, r_earth
      integer MSGUN, IUN, MLAT, MLON, MALT, NGRF, LWK
cga
cga Declarations from Richmond.....
cga
      PARAMETER (RTOD=57.2957795130823, DTOR=0.01745329251994330,
     &           R_EARTH=6371.2)
      real B(3),BHAT(3),
     &     D1(3),D2(3),D3(3), 
     &     E1(3),E2(3),E3(3), 
     &     F1(2),F2(2)

C          FILNAM contain gridded arrays for the epoch given by "DATE."
      CHARACTER*80 FILNAM

C          Declarations needed for APXMKA, APXWRA, or APXRDA
C            MSGUN = Fortran unit number for diagnostics
C            IUN   = Fortran unit number for I/O
C            MLAT,MLON,MALT = Maximum number of grid latitudes, longitudes
C                    and altitudes.
C            NGRF  = Number of epochs in the current DGRF/IGRF; see COFRM in
C                    file magfld.f
      PARAMETER (MSGUN=6, IUN=12, MLAT=91,MLON=151,MALT=24, NGRF=8,      ! for 1965-2005
     +           LWK= MLAT*MLON*MALT*5 + MLAT+MLON+MALT)
      real WK(LWK)
cga
cga ..end of declarations from Richmond.
cga
cg .. The flux tubes are defined to be the same as the tubes
cg    used for the TIEGCM electrodynamics....
cg    The following is a file containing the required tiegcm apex heights....
      open(3,file='tiegcm_defined_apex_heights',
     &             form='formatted',status='old')
cg 
cg  We are searching for the largest and Smallest tubes within our chosen group
cg  So we first need to initialise the Largest and Smallest L-values....\
cg
      i_first_index_tube_less_than_L4 = 0
      Largest_Apex_L_value = 0.0
      Smallest_Apex_L_value = 1000.0
      do i=1,48
        read(3,*) idummy,lat_base_90,apex_height,apex_radius,dummy,dummy,dummy
        Apex_L_value(i) = apex_radius / R_EARTH
cg
cg  We are only interested in tubes which have 'L values' of less than 4.
cg  the .gt.0.0 part below is to filter out the negative values which are
cg  used at the poles (instead of infinity - look at the file - unit 3 above.
cg  You'll see what I mean).
cg
         if(Apex_L_value(i).le.4.0.and.Apex_L_value(i).gt.0.0) then
         if(i_first_index_tube_less_than_L4.eq.0) then
         i_first_index_tube_less_than_L4 = i
         endif
c        write(6,*) idummy,lat_base_90,Apex_L_value(i)
         if(Apex_L_value(i).lt.Smallest_Apex_L_value) then
           Smallest_Apex_L_value = Apex_L_value(i)
           i_smallest = i
         endif
         if(Apex_L_value(i).gt.Largest_Apex_L_value) then
           Largest_Apex_L_value = Apex_L_value(i)
           i_largest = i
         endif
         endif
      enddo
cg
      iNumber_of_tubes = i_smallest - i_largest + 1
cg
cg  Check that iNumber_of_tubes is equal to NLP.  If it isn't then
cg  the code needs to be recompiled with NLP set to iNumber_of_tubes....
cg
c     if(NLP.NE.iNumber_of_tubes) then 
c       write(6,*) '****************************************'
c       write(6,*) '* The number of tubes we are using is'
c       write(6,*) '* different from NLP as defined in npts.h'
c       write(6,*) '* Please edit npts.h setting NLP to ',
c    &              iNumber_of_tubes
c       write(6,*) '* then recompile the code and try again.'
c       write(6,*) '* STOPPED'
c       write(6,*) '****************************************'
c       stop
c     endif
c
      write(6,*) i_smallest,Smallest_Apex_L_value,i_largest,
     &           Largest_Apex_L_value,iNumber_of_tubes
      write(6,*) i_first_index_tube_less_than_L4
      i_num_tubes_for_tiegcm = ((i_first_index_tube_less_than_L4 - 1)*2)
     &  + (iNumber_of_tubes*2) + 1   ! the +1 is the central tube/point with apex height of 90km
      write(6,*) 'i_num_tubes_for_tiegcm =',i_num_tubes_for_tiegcm
cg
cg  Now we loop again from 1 to iNumber_of_tubes to get
cg  our final required L-values in the right order.....
cg
      do i=1,iNumber_of_tubes
        tiegcm_L_value(i) = Apex_L_value(i+i_largest-1)
        write(6,*) i,tiegcm_L_value(i)
      enddo
cg
cg  Double the number of tubes.......
cg
       do i = 1, iNumber_of_tubes
       i_doubled_tube = (i * 2) - 1
       tiegcm_L_value2(i_doubled_tube) = tiegcm_L_value(i)
       enddo
       do i = 1, iNumber_of_tubes - 1
       i_doubled_tube = (i * 2) 
       tiegcm_L_value2(i_doubled_tube) = (tiegcm_L_value(i)
     &                                   +tiegcm_L_value(i+1))/2.0
       enddo
cg
       i_new_number_of_tubes = ( 2 * iNumber_of_tubes ) - 1
       do i=1,i_new_number_of_tubes
         write(6,*) 'new ',i,tiegcm_L_value2(i)
       enddo
cg
cg  Now we have our Apex L values ordered correctly
cg  in the nlp grid...
cg
      istop = 0
      if(istop.eq.1) stop
c
c     open(8,file='plasma_tube_apex_coords.2000')
c     &             form='unformatted')
c
      data pi/3.141592654e0/
      dtr=pi/180.
      r0=6.370e06
c
      n_mid_point=(npts+1)/2
c
cga
cga ...from Richmond. Read in the Apex-grid dataset....
cga
C          Sample date (corresponding to 2000 July 2)
c     DATE  = 1980.
      read(5,*) date_string
      read(5,*) date

      write(6,*) 'Input Date = ',date

      DATE_dble  = dble(date)
      NDATE = 1
C          Set up IGRF coefficients for the date desired
c     CALL COFRM (DATE)

C          Name of file containing gridded arrays
c     FILNAM = 'Apex_grid-2000.5'
      FILNAM = 'Apex_grid_data'

      CALL APXRDA (MSGUN, FILNAM,IUN, DATE, WK,LWK, IST)
      IF (IST .NE. 0) then
        write(6,*) 'APXRDA error: IST= ',IST
        stop
      endif

cga ... end of from Richmond.
cga
c
c  loop over all the flux tubes required to set up initial field line
c  coordinates for each one........
c
c       mp = 1
      do 1000  mp=1,nmp
 3285 format('   initial setup for group no.',i3)
        write(6,3285) mp
cg
cg  firstly define the magnetic longitude of each group of tubes.....
cg
      blonb=(float(mp-1))*360./(float(nmp))
cg
cg  then the l values...
cg
cga
cga ..The L-values to use have been read in at the top - taken from the present CTIP
cga centred dipole code.  The outer tube is at L=3.5, going in to smaller numbers
cga (around 1).
cga
c       lp = 1
      do 1500 lp=1,nlp
cg
cg  Now we are using the tiegcm defined L values
cg  - not CTIP dipole ones....
cg
         l_value = tiegcm_L_value2(lp)
cg looping
cg
c         write(6,*) 'lp',lp,l_value
c
cga
cga alon and HA used by richmond...
cga ...note, alon in degrees, HA in km....
cga
      alon = blonb
      HA = (l_value - 1.0) * R_EARTH
c     write(6,*) 'Apex Height = ',HA
cga
        re(mp,lp)=l_value*r0
        blon(mp,lp)=blonb*dtr
C Set reference height to 90 km:
      HR = 90.
C Find geographic coordinates of apex of this field line (QDLAT=0.,
C   ALT=HA):
      QDLAT = 0.
      ALT = HA
      CALL APXQ2G (QDLAT,ALON,ALT, WK, GLAT_apex,GLON_apex, IST)
c     write (6,*) 'APEX GLAT= ',GLAT_apex,'  GLON= ',GLON_apex

cg
cg *********************** AT THE APEX ****************************
cg
C Get modified-apex and quasi-dipole coordinates and associated
C   parameters for apex:
      CALL APXMALL (GLAT_apex,GLON_apex,ALT,HR, WK,
     +            B,BHAT,BABS,SI,XLONM,
     +             XLATM,VMP,W,D,BE3,SIM,D1,D2,D3,E1,E2,E3,
     +             XLATQD,F,F1,F2 , IST)
c     write(6,*) 'Magnetic potential =',vmp
      vmp_apex(mp,lp) = vmp
cg
cg magnitude of B at the Apex.....
cg     
       B_magnitude_dble(n_mid_point,mp,lp) = dble(BABS*1.e-9)
       BE3_dble(n_mid_point,mp,lp) = dble(BE3)
cg
cg vectors needed to be output also....
cg
       do j_direction=1,3
       D1_dble(j_direction,n_mid_point,mp,lp) = dble(D1(j_direction))
       D2_dble(j_direction,n_mid_point,mp,lp) = dble(D2(j_direction))
       D3_dble(j_direction,n_mid_point,mp,lp) = dble(D3(j_direction))
       E1_dble(j_direction,n_mid_point,mp,lp) = dble(E1(j_direction))
       E2_dble(j_direction,n_mid_point,mp,lp) = dble(E2(j_direction))
       E3_dble(j_direction,n_mid_point,mp,lp) = dble(E3(j_direction))
       enddo
cg
cg
cga
cga Find the magnetic potential for the south pole sea level....
cga (QDLAT = -90., ALT = 0.)
      QDLAT = -90.
      ALT = 0.0
      CALL APXQ2G (QDLAT,ALON,ALT, WK, GLAT_sp,GLON_sp, IST)
      CALL APXMALL (GLAT_sp,GLON_sp,ALT,HR, WK,
     +            B,BHAT,BABS,SI,XLONM,
     +             XLATM,VMP,W,D,BE3,SIM,D1,D2,D3,E1,E2,E3,
     +             XLATQD,F,F1,F2 , IST)
c     write (6,*) 'SP GLAT= ',GLAT_sp,'  GLON= ',GLON_sp
c     write(6,*) 'Magnetic potential (sp) =',vmp
      vmp_south_pole = vmp
      vmp_south_pole_dble = dble(vmp_south_pole)
cga
cga calculate the potential for the northen and southern ends of the field line
cga (i.e. where the field line crosses 90km).....
cga
       ALT = 90.
       QDLAT = ACOS(SQRT((R_EARTH+ALT)/(R_EARTH+HA)))*RTOD
c      write (6,*) 'QDLAT= ',QDLAT
       CALL APXQ2G (QDLAT,ALON,ALT, WK, GLAT_north90,GLON_north90, IST)
      CALL APXMALL (GLAT_north90,GLON_north90,ALT,HR, WK,
     +            B,BHAT,BABS,SI,XLONM,
     +             XLATM,VMP,W,D,BE3,SIM,D1,D2,D3,E1,E2,E3,
     +             XLATQD,F,F1,F2 , IST)
c     write (6,*) 'GLAT_N ',GLAT_north90,'  GLON_N ',GLON_north90
c     write(6,*) 'Magnetic potential north =',vmp
      vmp_north = vmp
cga
cga .... and for the southern
cga
       ALT = 90.
       QDLAT = -ACOS(SQRT((R_EARTH+ALT)/(R_EARTH+HA)))*RTOD
c      write (6,*) 'QDLAT= ',QDLAT
       CALL APXQ2G (QDLAT,ALON,ALT, WK, GLAT_south90,GLON_south90, IST)
      CALL APXMALL (GLAT_south90,GLON_south90,ALT,HR, WK,
     +            B,BHAT,BABS,SI,XLONM,
     +             XLATM,VMP,W,D,BE3,SIM,D1,D2,D3,E1,E2,E3,
     +             XLATQD,F,F1,F2 , IST)
c      write (6,*) 'GLAT_S ',GLAT_south90,'  GLON_S ',GLON_south90
c     write(6,*) 'Magnetic potential south =',vmp
      vmp_south = vmp

c     write(6,*) 'Apex height ',mp,lp,HA
c     write(6,*) 'LAT ',glat_north90,glat_apex,glat_south90
c     write(6,*) 'LON',glon_north90,glon_apex,glon_south90
c     write(6,*) 'VMP',vmp_north,vmp_apex,vmp_south
cga
cga writing out
cga
cga
cga for plotting purposes make all longitudes +ve...
cga
c     if (glon_north90.lt.0.0) glon_north90 = glon_north90 + 360.
c     if (glon_apex.lt.0.0) glon_apex = glon_apex + 360.
c     if (glon_south90.lt.0.0) glon_south90 = glon_south90 + 360.
c     write(4,8979) glat_north90,glon_north90,
c    &          glat_apex,glon_apex,
c    &          glat_south90,glon_south90
 8979 format(6f7.1)
cga
cga  For each flux-tube we now have the geographic end points and the apex point.
cg  New technique is to use the Sqrt(ha - 90) formula from Art to calculate the
cg  positions of the points and the number of points per tube....
cg
       height = 90000.
       do iht = 1,1000
cg2007       height = height + float(iht-1) * sqrt((HA*1000.)-height)
c       height = height + (float(iht-1) * sqrt((HA*1000.)-height))*0.4
        if (lp.lt.11) then
        height = height + (float(iht-1) * sqrt((HA*1000.)-height))*0.2
        else
        height = height + (float(iht-1) * sqrt((HA*1000.)-height))*0.4
        endif
!       write(6,*) iht,height/1000.
        if(height.gt.(HA*1000.)) then
         i_top_point = iht - 1
         goto 1622
        endif
       enddo
 1622  continue
       ndx = i_top_point
       in(mp,lp)=n_mid_point-ndx
       is(mp,lp)=n_mid_point+ndx
       write(6,*) 'Number of points per hemi ',ndx,in(mp,lp),is(mp,lp)
       istop = 0
       if(istop.eq.1) stop

cg
cg ***************** NORTHERN HEMISPHERE POINTS *******************
cg 
       height = 90000.
       icount = 0
       do i = in(mp,lp),n_mid_point-1
       icount = icount + 1
cg2007       height = height + float(icount - 1) * sqrt((HA*1000.)-height)
c      height = height + (float(icount-1) * sqrt((HA*1000.)-height))*0.4
        if (lp.lt.11) then
         height = height + (float(icount-1) * sqrt((HA*1000.)-height))*0.2
        else
         height = height + (float(icount-1) * sqrt((HA*1000.)-height))*0.4
        endif
       fix_h2(i)=height/1000.
       ALT = fix_h2(i)
       QDLAT = ACOS(SQRT((R_EARTH+fix_h2(i))/(R_EARTH+HA)))*RTOD
       bcol(i,mp,lp) = (90. - qdlat)*DTOR
       CALL APXQ2G (QDLAT,ALON,ALT,WK,GLAT_fix(i),GLON_fix(i),IST)
       CALL APXMALL (GLAT_fix(i),GLON_fix(i),ALT,HR, WK,
     +            B,BHAT,BABS,SI,XLONM,
     +             XLATM,VMP,W,D,BE3,SIM,D1,D2,D3,E1,E2,E3,
     +             XLATQD,F,F1,F2 , IST)
cg
cg We need the magnitude of B as a global parameter (in Tesla)
cg     
       B_magnitude_dble(i,mp,lp) = dble(BABS*1.e-9)
       BE3_dble(i,mp,lp) = dble(BE3)
cg
cg vectors needed to be output also....
cg
       do j_direction=1,3
       D1_dble(j_direction,i,mp,lp) = dble(D1(j_direction))
       D2_dble(j_direction,i,mp,lp) = dble(D2(j_direction))
       D3_dble(j_direction,i,mp,lp) = dble(D3(j_direction))
       E1_dble(j_direction,i,mp,lp) = dble(E1(j_direction))
       E2_dble(j_direction,i,mp,lp) = dble(E2(j_direction))
       E3_dble(j_direction,i,mp,lp) = dble(E3(j_direction))
       enddo
cg
       this_q(i) = 0.0 - ((vmp - vmp_apex(mp,lp)) / vmp_south_pole)
       vmp_tube(i,mp,lp) = vmp
       enddo
cga
       this_q(n_mid_point) = 0.0 
       vmp_tube(n_mid_point,mp,lp) = vmp_apex(mp,lp)
       bcol(n_mid_point,mp,lp) = 90.0*DTOR
cga
cg
cg ***************** SOUTHERN HEMISPHERE POINTS *******************
cg 
       height = 90000.
       icount = 0
c       do i = n_mid_point+1, is(mp,lp)
       do i = is(mp,lp) , n_mid_point+1, - 1
       icount = icount + 1
cg2007       height = height + float(icount - 1) * sqrt((HA*1000.)-height)
c      height = height + (float(icount-1) * sqrt((HA*1000.)-height))*0.4
        if (lp.lt.11) then
         height = height + (float(icount-1) * sqrt((HA*1000.)-height))*0.2
        else
         height = height + (float(icount-1) * sqrt((HA*1000.)-height))*0.4
        endif
       fix_h2(i)=height/1000.
       ALT = fix_h2(i)
       QDLAT = -ACOS(SQRT((R_EARTH+fix_h2(i))/(R_EARTH+HA)))*RTOD
       bcol(i,mp,lp) = (90. - qdlat)*DTOR
       CALL APXQ2G (QDLAT,ALON,ALT,WK,GLAT_fix(i),GLON_fix(i),IST)
       CALL APXMALL (GLAT_fix(i),GLON_fix(i),ALT,HR, WK,
     +            B,BHAT,BABS,SI,XLONM,
     +             XLATM,VMP,W,D,BE3,SIM,D1,D2,D3,E1,E2,E3,
     +             XLATQD,F,F1,F2 , IST)
cg
cg magnitude of B as a global parameter for Southern hemisphere...
cg     
       B_magnitude_dble(i,mp,lp) = dble(BABS*1.e-9)
       BE3_dble(i,mp,lp) = dble(BE3)
cg
cg vectors needed to be output also....
cg
       do j_direction=1,3
       D1_dble(j_direction,i,mp,lp) = dble(D1(j_direction))
       D2_dble(j_direction,i,mp,lp) = dble(D2(j_direction))
       D3_dble(j_direction,i,mp,lp) = dble(D3(j_direction))
       E1_dble(j_direction,i,mp,lp) = dble(E1(j_direction))
       E2_dble(j_direction,i,mp,lp) = dble(E2(j_direction))
       E3_dble(j_direction,i,mp,lp) = dble(E3(j_direction))
       enddo
cg
       this_q(i) = 0.0 - ((vmp - vmp_apex(mp,lp)) / vmp_south_pole)
       vmp_tube(i,mp,lp) = vmp
       enddo
cga
cga  Create our double precision variables for output...
cga
       do i=in(mp,lp),is(mp,lp)
         gr(i,mp,lp) = (fix_h2(i)*1000.) + r0
         gcol(i,mp,lp) = (90. - glat_fix(i)) * dtr
         glon(i,mp,lp) = glon_fix(i)
         if(glon(i,mp,lp).lt.0.0) glon(i,mp,lp)=
     &      glon(i,mp,lp) + 360.
         glon(i,mp,lp) = glon(i,mp,lp) * dtr
         q(i,mp,lp) = this_q(i)
       enddo
cga
cga  Make sure we have our mid point (n_mid_point) properly defined...
cga
         gr(n_mid_point,mp,lp) = (HA*1000.) + r0
         gcol(n_mid_point,mp,lp) = (90. - glat_apex) * dtr
         glon(n_mid_point,mp,lp) = glon_apex
         if(glon(n_mid_point,mp,lp).lt.0.0) glon(n_mid_point,mp,lp)=
     &      glon(n_mid_point,mp,lp) + 360.
         glon(n_mid_point,mp,lp) = glon(n_mid_point,mp,lp) * dtr
         q(n_mid_point,mp,lp) = 0.0
cga
cga  Convert to double precision for output...
cga
       do i=in(mp,lp),is(mp,lp)
         gr_dble(i,mp,lp) = dble(gr(i,mp,lp))
         gcol_dble(i,mp,lp) = dble(gcol(i,mp,lp))
         glon_dble(i,mp,lp) = dble(glon(i,mp,lp))
         q_dble(i,mp,lp) = dble(q(i,mp,lp))
         bcol_dble(i,mp,lp) = dble(bcol(i,mp,lp))
         vmp_tube_dble(i,mp,lp) = dble(vmp_tube(i,mp,lp))
       enddo
         vmp_apex_dble(mp,lp) = dble(vmp_apex(mp,lp))
         re_dble(mp,lp) = dble(re(mp,lp))
         blon_dble(mp,lp) = dble(blon(mp,lp))
cga
cg
 1500 continue
 1000 continue
c
       istop = 0
       if(istop.eq.1) stop
       mgtype = 4

       output_file = 
     &     trim('GIP_apex_coords_etc.'//date_string//'.format')

       open(114,file=output_file,
     &    form='formatted',status='unknown')

         call calc_apex_params_2d_2(
     &           in,is,
     &           gr_dble,gcol_dble,glon_dble,
     &           re_dble,q_dble,blon_dble,
     &           bcol_dble,vmp_tube_dble,
     &           vmp_apex_dble,vmp_south_pole_dble,
     &           B_magnitude_dble,
     &           D1_dble,
     &           D2_dble,
     &           D3_dble,
     &           E1_dble,
     &           E2_dble,
     &           E3_dble,
     &           BE3_dble,date_dble)

c      write(8,1234) npts,nlp,nmp,mgtype
c      write(8,1234) in,is
c      write(8,5678) gr_dble,gcol_dble,glon_dble
c      write(8,5678) re_dble,q_dble,blon_dble
c      write(8,5678) bcol_dble,vmp_tube_dble
c      write(8,5678) vmp_apex_dble,vmp_south_pole_dble
c      write(8,5678) B_magnitude_dble
c      write(8,5678) D1_dble
c      write(8,5678) D2_dble
c      write(8,5678) D3_dble
c      write(8,5678) E1_dble
c      write(8,5678) E2_dble
c      write(8,5678) E3_dble
c      write(8,5678) BE3_dble

c1234  format(10i5)
c5678  format(10e12.4)

c      do lp = 1 , nlp
c        midpoint(lp) = (in(1,lp) + is(1,lp))/2
c      enddo

c      call newinterp_high_res_TEC(MLOw,MHIgh,
c    &                             IN,IS,
c    &                             midpoint,
c    &    ii1,ii2,ii3,ii4,facfac,dd,gr_dble,gcol_dble,glon_dble)


c      close(8)
c
      stop 'normal end'
      end
