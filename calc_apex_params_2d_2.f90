    subroutine calc_apex_params_2d_2( &
                  IN , IS, &
                  gr, gcol, glon, &
                  RE_apex , q_coordinate , BLOn, &
                  bcol,vmp_tube, &
                  vmp_apex,vmp_south_pole, &
                  B_magnitude_on_tubes_3D, &
                  Apex_D1, &
                  Apex_D2, &
                  Apex_D3, &
                  Apex_E1, &
                  Apex_E2, &
                  Apex_E3, &
                  Apex_BE3, &
                  DATE)

    IMPLICIT NONE

    INTEGER :: NPTS , NMP , NLP
    INCLUDE 'npts.h'
    INTEGER :: IN(NMP,NLP) , IS(NMP,NLP)
    INTEGER :: i , j
    INTEGER :: mp , lp
    REAL(kind=8) R0
    REAL(kind=8) PI
    REAL(kind=8) DTR
    REAL(kind=8) RTD
    REAL(kind=8) Date
    REAL(kind=8) RE_apex(NMP,NLP)
    REAL(kind=8) Apex_D1(3,NPTS,NMP,NLP)
    REAL(kind=8) Apex_D2(3,NPTS,NMP,NLP)
    REAL(kind=8) Apex_D3(3,NPTS,NMP,NLP)
    REAL(kind=8) Apex_E1(3,NPTS,NMP,NLP)
    REAL(kind=8) Apex_E2(3,NPTS,NMP,NLP)
    REAL(kind=8) Apex_E3(3,NPTS,NMP,NLP)
    REAL(kind=8) Apex_BE3(NPTS,NMP,NLP)
    REAL(kind=8) Apex_grdlbm2(3,NPTS,NMP,NLP)
    REAL(kind=8) ETA_APEX_3D(NPTS,NMP,NLP)
    REAL(kind=8) e3e3
    REAL(kind=8) Apex_BHAT(3,npts,nmp,nlp)
    REAL(kind=8) Apex_D(npts,nmp,nlp)
    REAL(kind=8) Apex_d1d1(npts,nmp,nlp)
    REAL(kind=8) Apex_d1d2(npts,nmp,nlp)
    REAL(kind=8) Apex_d2d2(npts,nmp,nlp)
    REAL(kind=8) Apex_BMAG(npts,nmp,nlp)
    REAL(kind=8) vmp_tube(npts,nmp,nlp)
    REAL(kind=8) vmp_apex(nmp,nlp)
    REAL(kind=8) vmp_south_pole
    REAL(kind=8) B_magnitude_on_tubes_3D(npts,nmp,nlp)
    REAL(kind=8) gr(NPTS,nmp,nlp)
    REAL(kind=8) gcol(NPTS,nmp,nlp)
    REAL(kind=8) glon(NPTS,nmp,nlp)
    REAL(kind=8) bcol(NPTS,NMP,NLP)
    REAL(kind=8) BLOn(NMP,NLP)
    REAL(kind=8) q_coordinate(NPTS,NMP,NLP)
    INTEGER :: i_count
    INTEGER :: IN_2d(NLP) , IS_2d(NLP)
    INTEGER :: IN_2d_3d(NMP,NLP) , IS_2d_3d(NMP,NLP)
    integer npts2
!   parameter (npts2=12279)
    parameter (npts2=13813)
    real(kind=8) gr_2d(npts2,nmp)
    real(kind=8) gcol_2d(npts2,nmp)
    real(kind=8) glon_2d(npts2,nmp)
    real(kind=8) q_coordinate_2d(npts2,nmp)
    real(kind=8) bcol_2d(npts2,nmp)
    real(kind=8) vmp_tube_2d(npts2,nmp)
    real(kind=8) B_magnitude_on_tubes_3D_2d(npts2,nmp)
    real(kind=8) Apex_D1_2d(3,npts2,nmp)
    real(kind=8) Apex_D2_2d(3,npts2,nmp)
    real(kind=8) Apex_D3_2d(3,npts2,nmp)
    real(kind=8) Apex_E1_2d(3,npts2,nmp)
    real(kind=8) Apex_E2_2d(3,npts2,nmp)
    real(kind=8) Apex_E3_2d(3,npts2,nmp)
    real(kind=8) Apex_BE3_2d(npts2,nmp)
    real(kind=8) Apex_grdlbm2_2d(3,npts2,nmp)
    real(kind=8) Eta_Apex_3d_2d(npts2,nmp)
    real(kind=8) apex_D_2d(npts2,nmp)
    real(kind=8) apex_BHAT_2d(3,npts2,nmp)
    real(kind=8) apex_d1d1_2d(npts2,nmp)
    real(kind=8) apex_d1d2_2d(npts2,nmp)
    real(kind=8) apex_d2d2_2d(npts2,nmp)
    real(kind=8) apex_BMAG_2d(npts2,nmp)

    integer midpoint(nlp)

    PARAMETER (PI=3.141592654,DTR=PI/180.0,RTD=180.0/PI,R0=6.370E06)


        call cofrm2(date)
        call get_apex2(gr, gcol, glon,npts,nmp,nlp, &
        Apex_grdlbm2,in,is)


!g
!g  For our new Apex code we calculate magnetic field parameters as 3D
!g  variables just once - at the start....
!g
    do mp=1,nmp
        do lp=1,nlp
            do i=in(mp,lp),is(mp,lp)
            !g
            !g Be3 as read in above is in nT - we need T therefore need to multiply by 1.e-09 here.
            !g This also feeds through to BMAG below (which is therfore automatically correctly in Tesla).....
            !g
                Apex_Be3(i,mp,lp) = Apex_Be3(i,mp,lp) * 1.e-09
            !g
            !g Added a minus sign here to ETA - otherwise DS comes out negative.....
            !g
                Eta_Apex_3d(i,mp,lp) = 0.0 - ( B_magnitude_on_tubes_3D(i,mp,lp) &
                / vmp_south_pole )
            enddo
        enddo
    enddo
!g
!g  We also need Bmag, D, bHat and Dhat as 3D variables.....
!g
! get D
    lp_loop1: do lp = 1,nlp
    mp_loop1: do mp = 1,nmp
    ipts_loop1: do i = in(mp,lp),is(mp,lp)
    e3e3=dot_product( Apex_E3(1:3,i,mp,lp), Apex_E3(1:3,i,mp,lp) )
    apex_D(i,mp,lp)=SQRT( e3e3 ) !(3.13) D * D = e3 * e3
enddo  ipts_loop1
enddo  mp_loop1
enddo  lp_loop1

! get BHAT
    lp_loop2: do lp = 1,nlp
    mp_loop2: do mp = 1,nmp
    do J=1,3
        apex_BHAT(J,in(mp,lp):is(mp,lp),mp,lp)=apex_D3(J,in(mp,lp):is(mp,lp),mp,lp)*apex_D(in(mp,lp):is(mp,lp),mp,lp)
    enddo
enddo mp_loop2
enddo lp_loop2

! get d1d2 etc
    lp_loop3: do lp = 1,nlp
    mp_loop3: do mp = 1,nmp
    ipts_loop3: do i = in(mp,lp),is(mp,lp)
    apex_d1d1(i,mp,lp)=dot_product(Apex_d1(1:3,i,mp,lp),Apex_d1(1:3,i,mp,lp))
    apex_d1d2(i,mp,lp)=dot_product(Apex_d1(1:3,i,mp,lp),Apex_d2(1:3,i,mp,lp))
    apex_d2d2(i,mp,lp)=dot_product(Apex_d2(1:3,i,mp,lp),Apex_d2(1:3,i,mp,lp))
enddo  ipts_loop3 !: do i = in(mp,lp),is(mp,lp)
enddo  mp_loop3   !: do mp = 1,nmp
enddo  lp_loop3   !: do lp = 1,nlp

! get BMAG
    lp_loop4: do lp = 1,nlp
    mp_loop4: do mp = 1,nmp
    apex_BMAG(in(mp,lp):is(mp,lp),mp,lp) = Apex_BE3(in(mp,lp):is(mp,lp),mp,lp)*apex_D(in(mp,lp):is(mp,lp),mp,lp)    !(4.13)
enddo  mp_loop4  !: do mp = 1,nmp
enddo  lp_loop4  !: do lp = 1,nlp

!   open(114,file='plasma_tube_apex_coords_2.2005.format_latest',form='formatted',status='unknown')

    i_count = 0
    do lp = 1 , nlp
    in_2d(lp) = i_count + 1
    is_2d(lp) = in_2d(lp) + is(1,lp) - in(1,lp)
    i_count = is_2d(lp)
    enddo
    
    do lp = 1 , nlp
      write(6,*) 'points ',lp,in(1,lp),is(1,lp),in_2d(lp),is_2d(lp)
    enddo

    do mp = 1 , nmp
    do lp = 1 , nlp
      
    gr_2d(in_2d(lp):is_2d(lp),mp) = gr(in(mp,lp):is(mp,lp),mp,lp)
    gcol_2d(in_2d(lp):is_2d(lp),mp) = gcol(in(mp,lp):is(mp,lp),mp,lp)
    glon_2d(in_2d(lp):is_2d(lp),mp) = glon(in(mp,lp):is(mp,lp),mp,lp)
    q_coordinate_2d(in_2d(lp):is_2d(lp),mp) = q_coordinate(in(mp,lp):is(mp,lp),mp,lp)
    bcol_2d(in_2d(lp):is_2d(lp),mp) = bcol(in(mp,lp):is(mp,lp),mp,lp)
    vmp_tube_2d(in_2d(lp):is_2d(lp),mp) = vmp_tube(in(mp,lp):is(mp,lp),mp,lp)
    B_magnitude_on_tubes_3D_2d(in_2d(lp):is_2d(lp),mp) = B_magnitude_on_tubes_3D(in(mp,lp):is(mp,lp),mp,lp)
    Apex_BE3_2d(in_2d(lp):is_2d(lp),mp) = Apex_BE3(in(mp,lp):is(mp,lp),mp,lp)
    Eta_Apex_3d_2d(in_2d(lp):is_2d(lp),mp) = Eta_Apex_3d(in(mp,lp):is(mp,lp),mp,lp)
    apex_D_2d(in_2d(lp):is_2d(lp),mp) = apex_D(in(mp,lp):is(mp,lp),mp,lp)
    apex_d1d1_2d(in_2d(lp):is_2d(lp),mp) = apex_d1d1(in(mp,lp):is(mp,lp),mp,lp)
    apex_d2d2_2d(in_2d(lp):is_2d(lp),mp) = apex_d2d2(in(mp,lp):is(mp,lp),mp,lp)
    apex_BMAG_2d(in_2d(lp):is_2d(lp),mp) = apex_BMAG(in(mp,lp):is(mp,lp),mp,lp)

    do j = 1 , 3
    Apex_D1_2d(j,in_2d(lp):is_2d(lp),mp) = Apex_D1(j,in(mp,lp):is(mp,lp),mp,lp)
    Apex_D2_2d(j,in_2d(lp):is_2d(lp),mp) = Apex_D2(j,in(mp,lp):is(mp,lp),mp,lp)
    Apex_D3_2d(j,in_2d(lp):is_2d(lp),mp) = Apex_D3(j,in(mp,lp):is(mp,lp),mp,lp)
    Apex_E1_2d(j,in_2d(lp):is_2d(lp),mp) = Apex_E1(j,in(mp,lp):is(mp,lp),mp,lp)
    Apex_E2_2d(j,in_2d(lp):is_2d(lp),mp) = Apex_E2(j,in(mp,lp):is(mp,lp),mp,lp)
    Apex_E3_2d(j,in_2d(lp):is_2d(lp),mp) = Apex_E3(j,in(mp,lp):is(mp,lp),mp,lp)
    Apex_grdlbm2_2d(j,in_2d(lp):is_2d(lp),mp) = Apex_grdlbm2(j,in(mp,lp):is(mp,lp),mp,lp)
    apex_BHAT_2d(j,in_2d(lp):is_2d(lp),mp) = apex_BHAT(j,in(mp,lp):is(mp,lp),mp,lp)
    enddo

    enddo
    enddo

    do mp = 1 , nmp
    do lp = 1 , nlp
       in_2d_3d(mp,lp) = in_2d(lp)
       is_2d_3d(mp,lp) = is_2d(lp)
    enddo
    enddo

        write(114,*) IN_2d_3d , IS_2d_3d
        write(114,*) gr_2d, gcol_2d, glon_2d
!       write(115,*) IN_2d_3d , IS_2d_3d
!       write(115,*) gr_2d, gcol_2d, glon_2d
        write(114,*) RE_apex , q_coordinate_2d , BLOn
        write(114,*) bcol_2d,vmp_tube_2d
        write(114,*) vmp_apex,vmp_south_pole
        write(114,*) B_magnitude_on_tubes_3D_2d
        write(114,*) Apex_D1_2d
        write(114,*) Apex_D2_2d
        write(114,*) Apex_D3_2d
        write(114,*) Apex_E1_2d
        write(114,*) Apex_E2_2d
        write(114,*) Apex_E3_2d
        write(114,*) Apex_BE3_2d
        write(114,*) Apex_grdlbm2_2d
        write(114,*) Eta_Apex_3d_2d
        write(114,*) e3e3
        write(114,*) apex_D_2d
        write(114,*) apex_BHAT_2d
        write(114,*) apex_d1d1_2d
        write(114,*) apex_d1d2_2d
        write(114,*) apex_d2d2_2d
        write(114,*) apex_BMAG_2d


       do lp = 1 , nlp
         midpoint(lp) = (in_2d_3d(1,lp) + is_2d_3d(1,lp))/2
       enddo

       call newinterp_high_res_TEC(IN_2d_3d,IS_2d_3d,midpoint,gr_2d,gcol_2d,glon_2d)

!   stop
    end subroutine calc_apex_params_2d_2


    subroutine get_apex2(gr, gcol, glon,npts,nmp,nlp, &
    apex_grdlbm2,in, is)

! Reference:  Richmond, A. D., Ionospheric Electrodynamics Using
! Magnetic Apex Coordinates, J. Geomag. Geoelectr., 47, 191-212, 1995.

! INPUTS:
! GLAT = Geographic (geodetic) latitude, degrees, must be within
! the grid domain (GPLAT(1) <= GLAT <= GPLAT(NLAT)).
! GLON = Geographic (geodetic) longitude, degrees, must be within
! one revolution of the grid domain:
! GPLON(1) <= GLON-360.,GLON, or GLON+360. <= GPLON(NLON))
! ALT  = Altitude, km
! WK   = same as entry APXMKA
! RETURNS:
! BMAG   = magnitude of magnetic field, in nT
! BE3    = B_e3 of reference above (= Bmag/D), in nT
! D1,D2,D3,E1,E2,E3 = components (east, north, up) of base vectors
! described in reference above
! IST    = Return status:  okay (0); or failure (1).

! Dimensions of non-scalar arguments to APXMALL:
! GPLAT(NLAT),GPLON(NLON),GPALT(NALT),WK(LWK),
! B(3),BHAT(3),D1(3),D2(3),D3(3), E1(3),E2(3),E3(3), F1(2),F2(2)

    implicit none

    integer, intent(in) :: npts,nmp,nlp,in(nmp,nlp),is(nmp,nlp)
    real(kind = 8), intent(in) :: &
    gr(npts,nmp,nlp) ,    & ! Earth radii[m]
    gcol(npts,nmp,nlp) ,  & ! geog. colatitude [rad]
    glon(npts,nmp,nlp)    ! geog. longitude [rad]
    real(kind = 8), intent(out) :: apex_grdlbm2(3,npts,nmp,nlp)
    REAL(kind=8) PI
    REAL(kind=8) DTR
    REAL(kind=8) RTD
    PARAMETER (PI=3.141592654,DTR=PI/180.0,RTD=180.0/PI)

! local:
    integer:: ipts,mp,lp
    real(kind = 8) :: glon_deg,glat_deg,z_km
! m
    real(kind = 8) :: x,y,z, &
    bnrth,beast,bdown,bpx,bmx,bpy,bmy,bpz,bmz, &
    re_plasma2,amount
    real(kind = 8) :: dlbm2dx,dlbm2dy,dlbm2dz, &
    coslat,sinlat,coslon,sinlon


    bnrth = 0.
    beast = 0.
    bdown = 0.
    amount = 10.
    re_plasma2 = 6360.

    do mp = 1,nmp     ! loop over # of flux tubes in longitude
        do lp = 1,nlp  ! loop over # of flux tubes in latitude
            do ipts = in(mp,lp),is(mp,lp)    ! loop over points along flux tube
                glat_deg = 90. - gcol(ipts,mp,lp)*rtd  ! geo.colat[rad] -> geo.lat[deg]
                glon_deg = glon(ipts,mp,lp)*rtd      ! geo.lon[rad] -> geo.lon[deg]
                z_km = gr(ipts,mp,lp)*1e-3  - re_plasma2      ! earth radii[m] -> altitude [km]
            
            !C Calculate geocentric cartesian coordinates x,y,z (in km)
                call GD2CART2(GLAT_deg,GLON_deg,z_km,X,Y,Z)
            !C Calculate gradient of ALOG(B0**(-2)) [= -2*ALOG(B0)] in Cartesian
            !C   coordinates, in units of m^{-1}
            !g
            !g
            !g
                call FELDG2 (2,(X+amount)/re_plasma2,Y/re_plasma2,Z/re_plasma2,BNRTH,BEAST,BDOWN,BPX)
                call FELDG2 (2,(X-amount)/re_plasma2,Y/re_plasma2,Z/re_plasma2,BNRTH,BEAST,BDOWN,BMX)
                call FELDG2 (2,X/re_plasma2,(Y+amount)/re_plasma2,Z/re_plasma2,BNRTH,BEAST,BDOWN,BPY)
                call FELDG2 (2,X/re_plasma2,(Y-amount)/re_plasma2,Z/re_plasma2,BNRTH,BEAST,BDOWN,BMY)
                call FELDG2 (2,X/re_plasma2,Y/re_plasma2,(Z+amount)/re_plasma2,BNRTH,BEAST,BDOWN,BPZ)
                call FELDG2 (2,X/re_plasma2,Y/re_plasma2,(Z-amount)/re_plasma2,BNRTH,BEAST,BDOWN,BMZ)
                ! Converted these logs to use LOG instead of ALOG. mjh
                dlbm2dx = (log(BMX) - log(BPX))/1.E4
                dlbm2dy = (log(BMY) - log(BPY))/1.E4
                dlbm2dz = (log(BMZ) - log(BPZ))/1.E4
            !C Rotate gradient to local (east,north,up) coordinates
            !C (1=east, 2=north, 3=up)
                coslat = cos(glat_deg*DTR)
                sinlat = sin(glat_deg*DTR)
                coslon = cos(glon_deg*DTR)
                sinlon = sin(glon_deg*DTR)
                apex_grdlbm2(1,ipts,mp,lp) = -dlbm2dx*sinlon + dlbm2dy*coslon
                apex_grdlbm2(2,ipts,mp,lp) = -(dlbm2dx*coslon + dlbm2dy*sinlon)*sinlat+dlbm2dz*coslat
                apex_grdlbm2(3,ipts,mp,lp) = (dlbm2dx*coslon + dlbm2dy*sinlon)*coslat+dlbm2dz*sinlat
            ! g
            enddo! loop over points along flux tube
        enddo  ! loop over # of flux tubes in latitude
    enddo	! loop over # of flux tubes in longitude

    return
    end subroutine get_apex2






    SUBROUTINE newinterp_high_res_TEC(IN,IS,midpoint,gr,gcol,glon)

    IMPLICIT NONE
    integer npts, nmp, nlp 
    include 'npts.h'
    REAL(KIND=8) :: DTR
    REAL(KIND=8) :: factor
    REAL(KIND=8) :: gcol
    REAL(KIND=8) :: geol1
    REAL(KIND=8) :: geol2
    REAL(KIND=8) :: &
    geolat , glon , gr , grin , height , &
    phngd , plat , plon
    REAL(KIND=8) :: R0 , thngd
    REAL(KIND=8) :: x1 , xn , y1 , yn , z1 , zn , PI
    REAL(KIND=8) :: &
    facc , facfac , dd
    INTEGER :: i , ic , ifailed , iheight , &
    ihem , ilon , IN , in1 , in2 , &
    IS
    INTEGER :: l , lp , m , mp ,  midpoint(nlp) , &
    inearst , i1 , i2 , i3 , i4 , ii1 , ii2 , ii3 , &
    ii4
    REAL(KIND=8) :: maxlat , minlat 
    REAL(KIND=8) :: original_distance, &
    sorted_distance
    INTEGER :: MHIgh(90) , MLOw(90) , max_ic, minny1(1), &
    minny(3)
    PARAMETER (PI=3.141592654,DTR=PI/180.0,R0=6.370E06)
    allocatable :: original_distance(:),sorted_distance(:)

    integer npts2
!   parameter (npts2=12279)
    parameter (npts2=13813)
    parameter (max_ic = 2*nmp*nlp)
    DIMENSION gr(NPTS2,NMP) , gcol(NPTS2,NMP) , &
    glon(NPTS2,NMP) , &
    grin(NPTS2) , IN(NMP,NLP) , IS(NMP,NLP)
    DIMENSION plat(max_ic) , plon(max_ic) , &
    facc(max_ic) ,  &
    i1(max_ic), i2(max_ic), i3(max_ic) , i4(max_ic) , &
    ii1(3,183,91,90), ii2(3,183,91,90), &
    ii3(3,183,91,90) , ii4(3,183,91,90) , &
    facfac(3,183,91,90) , &
    dd(3,183,91,90)
    REAL(KIND=8) :: glon1 , glon2
!g
    DO 300 iheight = 183 , 1 , -1
        height = FLOAT(iheight-1)*5000. + 90000. + R0
        write(6,*) iheight,height
    !g
    !g Loop over all the flux tubes......
    !g
        ic = 0
        DO 150 mp = 1 , NMP
            DO 140 lp = 1 , NLP
                DO 110 i = IN(mp,lp) , IS(mp,lp)
                    grin(i) = gr(i,mp)
                110 ENDDO
                DO 120 ihem = 1 , 2
                    IF ( ihem == 1 ) THEN
                        call FASTNEARHTPLA(grin,height,in(mp,lp),midpoint(lp),in1,in2, &
                        inearst,ifailed)
                    ELSE
                        call FASTNEARHTPLA(grin,height,is(mp,lp),midpoint(lp),in1,in2, &
                        inearst,ifailed)
                    ENDIF

                    IF ( ifailed == 0 ) THEN
                        ic = ic + 1
                    !g  interpolation factor for each point....
                        factor = (height-grin(in1))/(grin(in2)-grin(in1))
                        facc(ic)=factor
                    !g  position of this point in latitude and longitude....
                        plat(ic) = (((gcol(in2,mp)-gcol(in1,mp))* &
                        factor)+gcol(in1,mp))/DTR

                        glon2 = glon(in2,mp)
                        glon1 = glon(in1,mp)

                        if((glon2-glon1)/dtr > 10.) glon1 = glon1 + (360.*dtr)
                        if((glon2-glon1)/dtr < -10.) glon2 = glon2 + (360.*dtr)

                        plon(ic) = (((glon2 - glon1) * factor)+glon1)/DTR

                        IF ( plon(ic) >= 360. ) plon(ic) = plon(ic) - 360.
                        IF ( plon(ic) < 0. ) plon(ic) = plon(ic) + 360.

                    !g  relevent parameters at the point.....
                    !g
                    !g keep the flux-tube indexes for each 'ic' point in 4 arrays ......
                    !g
                        i1(ic)=mp
                        i2(ic)=lp
                        i3(ic)=in2
                        i4(ic)=in1
                    ENDIF
                120 ENDDO
            140 ENDDO
        150 ENDDO
    !g
    !g  Find out the max and min latitudes covered at each geo
    !g  longitude by points at 1000 km....
    !g
        IF ( iheight == 183 ) THEN
            DO 180 ilon = 1 , 90
                geol2 = FLOAT(ilon)*4.
                geol1 = FLOAT(ilon-2)*4.
                minlat = 90.
                maxlat = -90.
                DO 160 i = 1 , ic
                    IF ( ilon == 1 ) THEN
                        IF ( plon(i) < 4. .OR. plon(i) > 356. ) THEN
                            geolat = 90. - plat(i)
                            IF ( geolat > maxlat ) maxlat = geolat
                            IF ( geolat < minlat ) minlat = geolat
                        ENDIF
                    ELSEIF ( plon(i) <= geol2 .AND. plon(i) >= geol1 ) &
                        THEN
                        geolat = 90. - plat(i)
                        IF ( geolat > maxlat ) maxlat = geolat
                        IF ( geolat < minlat ) minlat = geolat
                    ENDIF
                160 ENDDO
                MHIgh(ilon) = INT(maxlat/2.) + 45
                MLOw(ilon) = INT(minlat/2.) + 47
                write(6,5760) ilon,maxlat,mhigh(ilon),minlat,mlow(ilon)
                5760 format('mhigh mlow ',i4,f5.1,i4,f6.1,i4)
            180 ENDDO
        ENDIF
    !g
    !g  Now interpolate all the values at a given height
    !g   onto the Geo grid
        DO 250 l = 1 , 90
            phngd = FLOAT(l-1)*4.
            DO 220 m = MLOw(l) , MHIgh(l)
                thngd = 180. - (2.*FLOAT(m-1))
                xn = height*SIN(thngd*DTR)*COS(phngd*DTR)
                yn = height*SIN(thngd*DTR)*SIN(phngd*DTR)
                zn = height*COS(thngd*DTR)
            !g
            !g  loop over all points
            !g
                if(allocated(original_distance)) deallocate(original_distance)
                if(allocated(sorted_distance)) deallocate(sorted_distance)
                allocate (original_distance(ic),sorted_distance(ic))
                DO 200 i = 1 , ic
                    x1 = height*SIN((plat(i))*DTR)*COS((plon(i))*DTR)
                    y1 = height*SIN((plat(i))*DTR)*SIN((plon(i))*DTR)
                    z1 = height*COS((plat(i))*DTR)
                !g Cartesian distance to neutral point....
                    original_distance(i) = SQRT((x1-xn)*(x1-xn)+(y1-yn)*(y1-yn)+(z1-zn) &
                    *(z1-zn))
                    sorted_distance(i) = original_distance(i)
                200 ENDDO
            !g
                minny1 = minloc(sorted_distance)
                minny(1) = minny1(1)
                sorted_distance(minny(1)) = 1.d50
                minny1 = minloc(sorted_distance)
                minny(2) = minny1(1)
                sorted_distance(minny(2)) = 1.d50
                minny1 = minloc(sorted_distance)
                minny(3) = minny1(1)
                sorted_distance(minny(3)) = 1.d50
            !g
                do i=1,3
                !g
                !g  The 3 closest points are defined by 4 integers and a factor (mp,lp,in2,in1,facfac)......
                !g
                    ii1(i,iheight,m,l)=i1(minny(i))
                    ii2(i,iheight,m,l)=i2(minny(i))
                    ii3(i,iheight,m,l)=i3(minny(i))
                    ii4(i,iheight,m,l)=i4(minny(i))
                    facfac(i,iheight,m,l)=facc(minny(i))
                    dd(i,iheight,m,l)=original_distance(minny(i))
                enddo

            !g lat,lon loop.....
            220 ENDDO
        250 ENDDO
    !g  height loop....
    300 ENDDO
!g
    write(114,*) facfac
    write(114,*) dd
    write(114,*) ii1
    write(114,*) ii2
    write(114,*) ii3
    write(114,*) ii4
    write(114,*) mlow
    write(114,*) mhigh
    close(114)
!g
    RETURN









    end SUBROUTINE newinterp_high_res_TEC

    SUBROUTINE FASTNEARHTPLA(ARR,VALue,M_IN,N_IN,M,N,nearst,IFAiled)

! this routine finds which two elements of ARR
! surround Value.  Arr is a monotonically increasing array
! The start and end points for the search are M_IN and N_IN.
! The exit points surrounding Value are N and M
    IMPLICIT NONE
    integer npts2
!   parameter (npts2=12279)
    parameter (npts2=13813)
    REAL(kind=8) :: ARR , VALue , dist1 , dist2
    INTEGER :: IFAiled , M , N ,  &
    M_IN , N_IN , n1 , nearst
    DIMENSION ARR(NPTS2)

    N=N_IN
    M=M_IN
    IFAiled = 0
    if((value <= arr(n) .AND. value >= arr(m)) .OR. &
    (value >= arr(n) .AND. value <= arr(m))) then
        234 n1 =(n - m)/2 + m
        IF ( ARR(n1) > VALue ) then
            n=n1
        else
            m=n1
        ENDIF
        if(abs(n-m) > 1) goto 234
        dist1=abs(arr(n)-value)
        dist2=abs(arr(m)-value)
        if(dist1 < dist2) then
            nearst=n
        else
            nearst=m
        endif
    else
        ifailed=1
    endif

    RETURN



    end SUBROUTINE FASTNEARHTPLA

      SUBROUTINE COFRM2 (DATE)
!          Assign DGRF/IGRF spherical harmonic coefficients, to degree and
!          order NMAX, for DATE, yyyy.fraction, into array G.  Coefficients
!          are interpolated from the DGRF dates through the current IGRF year.
!          Coefficients for a later DATE are extrapolated using the IGRF
!          initial value and the secular change coefficients.  A warning
!          message is issued if DATE is later than the last recommended
!          (5 yrs later than the IGRF).  An DATE input earlier than the
!          first DGRF (EPOCH(1)), results in a diagnostic and a STOP.
!
!          Output in COMMON /MAGCOF2/ NMAX,GB(144),GV(144),ICHG
!             NMAX = Maximum order of spherical harmonic coefficients used
!             GB   = Coefficients for magnetic field calculation
!             GV   = Coefficients for magnetic potential calculation
!             ICHG = Flag indicating when GB,GV have been changed in COFRM
!
!          HISTORY (blame):
!          COFRM and FELDG originated 15 Apr 83 by Vincent B. Wickwar
!          (formerly at SRI. Int., currently at Utah State).  Although set
!          up to accomodate second order time derivitives, the IGRF
!          (GTT, HTT) have been zero.  The spherical harmonic coefficients
!          degree and order is defined by NMAX (currently 10).
!
!          Jun 86:  Updated coefficients adding DGRF 1980 & IGRF 1985, which
!          were obtained from Eos Vol. 7, No. 24.  Common block MAG was
!          replaced by MAGCOF2, thus removing variables not used in subroutine
!          FELDG.  (Roy Barnes)
!
!          Apr 1992 (Barnes):  Added DGRF 1985 and IGRF 1990 as described
!          in EOS Vol 73 Number 16 Apr 21 1992.  Other changes were made so
!          future updates should:
!            (1) Increment NDGY;
!            (2) Append to EPOCH the next IGRF year;
!            (3) Append the next DGRF coefficients to G1DIM and H1DIM; and
!            (4) Replace the IGRF initial values (G0, GT) and rates of
!                change indices (H0, HT).
!
!          Apr 94 (Art Richmond): Computation of GV added, for finding
!          magnetic potential.
!
!          Aug 95 (Barnes):  Added DGRF for 1990 and IGRF for 1995, which were
!          obtained by anonymous ftp geomag.gsfc.nasa.gov (cd pub, mget table*)
!          as per instructions from Bob Langel (langel@geomag.gsfc.nasa.gov),
!          but, problems are to be reported to baldwin@geomag.gsfc.nasa.gov
 
!          Oct 95 (Barnes):  Correct error in IGRF-95 G 7 6 and H 8 7 (see
!          email in folder).  Also found bug whereby coefficients were not being
!          updated in FELDG when IENTY did not change. ICHG was added to flag
!          date changes.  Also, a vestigial switch (IS) was removed from COFRM:
!          It was always 0 and involved 3 branch if statements in the main
!          polynomial construction loop (now numbered 200).
 
!          Feb 99 (Barnes):  Explicitly initialize GV(1) in COFRM to avoid
!          possibility of compiler or loader options initializing memory
!          to something else (e.g., indefinite).  Also simplify the algebra
!          in COFRM; this does not effect results.

!          Mar 99 (Barnes):  Removed three branch if's from FELDG and changed
!          statement labels to ascending order

!          Jun 99 (Barnes):  Corrected RTOD definition in GD2CART.

!          May 00 (Barnes):  Replace IGRF 1995 with DGRF 1995, add IGRF
!          2000, and extend the earlier DGRF's backward to 1900.  A complete
!          set of coefficients came from a NGDC web page

      implicit none

      integer ICHG,nmax , I , IY , I1
      integer N , M , NN , MM
      integer NDGY , NYT , NGH
      REAL(kind=8) RNN
      REAL(kind=8) RTOD
      REAL(kind=8) DTOR
      REAL(kind=8) F , F0
      REAL(kind=8) G(144) , GV(144)
      REAL(kind=8) date , time , t

      PARAMETER (RTOD=57.2957795130823, DTOR=0.01745329251994330)
      COMMON /MAGCOF2/G,GV,ICHG,nmax
      DATA NMAX,ICHG /10,-99999/
 
      PARAMETER (NDGY=20 , NYT = NDGY+1 , NGH = 144*NDGY)
!          NDGY = Number of DGRF years of sets of coefficients
!          NYT  = Add one for the IGRF set (and point to it).
!          NGH  = Dimension of the equivalenced arrays
      real(kind=8) GYR(12,12,NYT) , HYR(12,12,NYT), EPOCH(NYT) , &
                G1DIM(NGH)     , H1DIM(NGH) , &
                G0(12,12) , GT(12,12) , GTT(12,12) , &
                H0(12,12) , HT(12,12) ,  HTT(12,12)

      EQUIVALENCE (GYR(1,1,1),G1DIM(1))  , (HYR(1,1,1),H1DIM(1)) , &
                  (GYR(1,1,NYT),G0(1,1)) , (HYR(1,1,NYT),H0(1,1))
 
      DATA EPOCH /1900, 1905, 1910, 1915, 1920, 1925, 1930, 1935, 1940, &
      1945, 1950, 1955, 1960, 1965, 1970, 1975, 1980, 1985, 1990, 1995, &
      2000/
!          D_/Dtime2 coefficients are 0
      DATA GTT/144*0./,HTT/144*0./

!          DGRF g(n,m) for 1900:
!          The "column" corresponds to "n" and
!          the "line" corresponds to "m" as indicated in column 6;
!          e.g., for 1965 g(0,3) = 1297. or g(6,6) = -111.
      DATA (G1DIM(I),I=1,144) /0, &
        -31543,  -677,  1022,  876, -184,   63,  70,  11,   8, -3,  2*0, &
         -2298,  2905, -1469,  628,  328,   61, -55,   8,  10, -4,  3*0, &
                  924,  1256,  660,  264,  -11,   0,  -4,   1,  2,  4*0, &
                         572, -361,    5, -217,  34,  -9, -11, -5,  5*0, &
                               134,  -86,  -58, -41,   1,  12, -2,  6*0, &
                                     -16,   59, -21,   2,   1,  6,  7*0, &
                                           -90,  18,  -9,  -2,  4,  8*0, &
                                                  6,   5,   2,  0,  9*0, &
                                                       8,  -1,  2, 10*0, &
                                                           -1,  2, 11*0, &
                                                                0, 13*0/
!          DGRF g(n,m) for 1905:
      DATA (G1DIM(I),I=145,288) /0, &
        -31464,  -728,  1037,  880, -192,   62,  70,  11,   8, -3,  2*0, &
         -2298,  2928, -1494,  643,  328,   60, -54,   8,  10, -4,  3*0, &
                 1041,  1239,  653,  259,  -11,   0,  -4,   1,  2,  4*0, &
                         635, -380,   -1, -221,  33,  -9, -11, -5,  5*0, &
                               146,  -93,  -57, -41,   1,  12, -2,  6*0, &
                                     -26,   57, -20,   2,   1,  6,  7*0, &
                                           -92,  18,  -8,  -2,  4,  8*0, &
                                                  6,   5,   2,  0,  9*0, &
                                                       8,   0,  2, 10*0, &
                                                           -1,  2, 11*0, &
                                                                0, 13*0/
!          DGRF g(n,m) for 1910:
      DATA (G1DIM(I),I=289,432) /0, &
        -31354,  -769,  1058,  884, -201,   62,  71,  11,   8, -3,  2*0, &
         -2297,  2948, -1524,  660,  327,   58, -54,   8,  10, -4,  3*0, &
                 1176,  1223,  644,  253,  -11,   1,  -4,   1,  2,  4*0, &
                         705, -400,   -9, -224,  32,  -9, -11, -5,  5*0, &
                               160, -102,  -54, -40,   1,  12, -2,  6*0, &
                                     -38,   54, -19,   2,   1,  6,  7*0, &
                                           -95,  18,  -8,  -2,  4,  8*0, &
                                                  6,   5,   2,  0,  9*0, &
                                                       8,   0,  2, 10*0, &
                                                           -1,  2, 11*0, &
                                                                0, 13*0/
!          DGRF g(n,m) for 1915:
      DATA (G1DIM(I),I=433,576) /0, &
        -31212,  -802,  1084,  887, -211,   61,  72,  11,   8, -3,  2*0, &
         -2306,  2956, -1559,  678,  327,   57, -54,   8,  10, -4,  3*0, &
                 1309,  1212,  631,  245,  -10,   2,  -4,   1,  2,  4*0, &
                         778, -416,  -16, -228,  31,  -9, -11, -5,  5*0, &
                               178, -111,  -51, -38,   2,  12, -2,  6*0, &
                                     -51,   49, -18,   3,   1,  6,  7*0, &
                                           -98,  19,  -8,  -2,  4,  8*0, &
                                                  6,   6,   2,  0,  9*0, &
                                                       8,   0,  1, 10*0, &
                                                           -1,  2, 11*0, &
                                                                0, 13*0/
!          DGRF g(n,m) for 1920:
      DATA (G1DIM(I),I=577,720) /0, &
        -31060,  -839,  1111,  889, -221,   61,  73,  11,   8, -3,  2*0, &
         -2317,  2959, -1600,  695,  326,   55, -54,   7,  10, -4,  3*0, &
                 1407,  1205,  616,  236,  -10,   2,  -3,   1,  2,  4*0, &
                         839, -424,  -23, -233,  29,  -9, -11, -5,  5*0, &
                               199, -119,  -46, -37,   2,  12, -2,  6*0, &
                                     -62,   44, -16,   4,   1,  6,  7*0, &
                                          -101,  19,  -7,  -2,  4,  8*0, &
                                                  6,   6,   2,  0,  9*0, &
                                                       8,   0,  1, 10*0, &
                                                           -1,  3, 11*0, &
                                                                0, 13*0/
!          DGRF g(n,m) for 1925:
      DATA (G1DIM(I),I=721,864) /0, &
        -30926,  -893,  1140,  891, -230,   61,  73,  11,   8, -3,  2*0, &
         -2318,  2969, -1645,  711,  326,   54, -54,   7,  10, -4,  3*0, &
                 1471,  1202,  601,  226,   -9,   3,  -3,   1,  2,  4*0, &
                         881, -426,  -28, -238,  27,  -9, -11, -5,  5*0, &
                               217, -125,  -40, -35,   2,  12, -2,  6*0, &
                                     -69,   39, -14,   4,   1,  6,  7*0, &
                                          -103,  19,  -7,  -2,  4,  8*0, &
                                                  6,   7,   2,  0,  9*0, &
                                                       8,   0,  1, 10*0, &
                                                           -1,  3, 11*0, &
                                                                0, 13*0/
!          DGRF g(n,m) for 1930:
      DATA (G1DIM(I),I=865,1008) /0, &
        -30805,  -951,  1172,  896, -237,   60,  74,  11,   8, -3,  2*0, &
         -2316,  2980, -1692,  727,  327,   53, -54,   7,  10, -4,  3*0, &
                 1517,  1205,  584,  218,   -9,   4,  -3,   1,  2,  4*0, &
                         907, -422,  -32, -242,  25,  -9, -12, -5,  5*0, &
                               234, -131,  -32, -34,   2,  12, -2,  6*0, &
                                     -74,   32, -12,   5,   1,  6,  7*0, &
                                          -104,  18,  -6,  -2,  4,  8*0, &
                                                  6,   8,   3,  0,  9*0, &
                                                       8,   0,  1, 10*0, &
                                                           -2,  3, 11*0, &
                                                                0, 13*0/ 
!          DGRF g(n,m) for 1935:
      DATA (G1DIM(I),I=1009,1152) /0, &
        -30715, -1018,  1206,  903, -241,   59,  74,  11,   8, -3,  2*0, &
         -2306,  2984, -1740,  744,  329,   53, -53,   7,  10, -4,  3*0, &
                 1550,  1215,  565,  211,   -8,   4,  -3,   1,  2,  4*0, &
                         918, -415,  -33, -246,  23,  -9, -12, -5,  5*0, &
                               249, -136,  -25, -33,   1,  11, -2,  6*0, &
                                     -76,   25, -11,   6,   1,  6,  7*0, &
                                          -106,  18,  -6,  -2,  4,  8*0, &
                                                  6,   8,   3,  0,  9*0, &
                                                       7,   0,  2, 10*0, &
                                                           -2,  3, 11*0, &
                                                                0, 13*0/
!          DGRF g(n,m) for 1940:
      DATA (G1DIM(I),I=1153,1296) /0, &
        -30654, -1106,  1240,  914, -241,   57,  74,  11,   8, -3,  2*0, &
         -2292,  2981, -1790,  762,  334,   54, -53,   7,  10, -4,  3*0, &
                 1566,  1232,  550,  208,   -7,   4,  -3,   1,  2,  4*0, &
                         916, -405,  -33, -249,  20, -10, -12, -5,  5*0, &
                               265, -141,  -18, -31,   1,  11, -2,  6*0, &
                                     -76,   18,  -9,   6,   1,  6,  7*0, &
                                          -107,  17,  -5,  -2,  4,  8*0, &
                                                  5,   9,   3,  0,  9*0, &
                                                       7,   1,  2, 10*0, &
                                                           -2,  3, 11*0, &
                                                                0, 13*0/
!          DGRF g(n,m) for 1945:
      DATA (G1DIM(I),I=1297,1440) /0, &
        -30594, -1244,  1282,  944, -253,   59,  70,  13,   5, -3,  2*0, &
         -2285,  2990, -1834,  776,  346,   57, -40,   7, -21, 11,  3*0, &
                 1578,  1255,  544,  194,    6,   0,  -8,   1,  1,  4*0, &
                         913, -421,  -20, -246,   0,  -5, -11,  2,  5*0, &
                               304, -142,  -25, -29,   9,   3, -5,  6*0, &
                                     -82,   21, -10,   7,  16, -1,  7*0, &
                                          -104,  15, -10,  -3,  8,  8*0, &
                                                 29,   7,  -4, -1,  9*0, &
                                                       2,  -3, -3, 10*0, &
                                                           -4,  5, 11*0, &
                                                               -2, 13*0/
!          DGRF g(n,m) for 1950:
      DATA (G1DIM(I),I=1441,1584) /0, &
        -30554, -1341,  1297,  954, -240,   54,  65,  22,   3, -8,  2*0, &
         -2250,  2998, -1889,  792,  349,   57, -55,  15,  -7,  4,  3*0, &
                 1576,  1274,  528,  211,    4,   2,  -4,  -1, -1,  4*0, &
                         896, -408,  -20, -247,   1,  -1, -25, 13,  5*0, &
                               303, -147,  -16, -40,  11,  10, -4,  6*0, &
                                     -76,   12,  -7,  15,   5,  4,  7*0, &
                                          -105,   5, -13,  -5, 12,  8*0, &
                                                 19,   5,  -2,  3,  9*0, &
                                                      -1,   3,  2, 10*0, &
                                                            8, 10, 11*0, &
                                                                3, 13*0/
!          DGRF g(n,m) for 1955:
      DATA (G1DIM(I),I=1585,1728) /0, &
        -30500, -1440,  1302,  958, -229,   47,  65,  11,   4, -3,  2*0, &
         -2215,  3003, -1944,  796,  360,   57, -56,   9,   9, -5,  3*0, &
                 1581,  1288,  510,  230,    3,   2,  -6,  -4, -1,  4*0, &
                         882, -397,  -23, -247,  10, -14,  -5,  2,  5*0, &
                               290, -152,   -8, -32,   6,   2, -3,  6*0, &
                                     -69,    7, -11,  10,   4,  7,  7*0, &
                                          -107,   9,  -7,   1,  4,  8*0, &
                                                 18,   6,   2, -2,  9*0, &
                                                       9,   2,  6, 10*0, &
                                                            5, -2, 11*0, &
                                                                0, 13*0/
!          DGRF g(n,m) for 1960:
      DATA (G1DIM(I),I=1729,1872) /0, &
        -30421, -1555,  1302,  957, -222,   46,  67,  15,   4,  1,  2*0, &
         -2169,  3002, -1992,  800,  362,   58, -56,   6,   6, -3,  3*0, &
                 1590,  1289,  504,  242,    1,   5,  -4,   0,  4,  4*0, &
                         878, -394,  -26, -237,  15, -11,  -9,  0,  5*0, &
                               269, -156,   -1, -32,   2,   1, -1,  6*0, &
                                     -63,   -2,  -7,  10,   4,  4,  7*0, &
                                          -113,  17,  -5,  -1,  6,  8*0, &
                                                  8,  10,  -2,  1,  9*0, &
                                                       8,   3, -1, 10*0, &
                                                           -1,  2, 11*0, &
                                                                0, 13*0/
!          DGRF g(n,m) for 1965:
      DATA (G1DIM(I),I=1873,2016) /0, &
        -30334, -1662,  1297,  957, -219,   45,  75,  13,   8, -2,  2*0, &
         -2119,  2997, -2038,  804,  358,   61, -57,   5,  10, -3,  3*0, &
                 1594,  1292,  479,  254,    8,   4,  -4,   2,  2,  4*0, &
                         856, -390,  -31, -228,  13, -14, -13, -5,  5*0, &
                               252, -157,    4, -26,   0,  10, -2,  6*0, &
                                     -62,    1,  -6,   8,  -1,  4,  7*0, &
                                          -111,  13,  -1,  -1,  4,  8*0, &
                                                  1,  11,   5,  0,  9*0, &
                                                       4,   1,  2, 10*0, &
                                                           -2,  2, 11*0, &
                                                                0, 13*0/
!          DGRF g(n,m) for 1970:
      DATA (G1DIM(I),I=2017,2160) /0, &
        -30220, -1781,  1287,  952, -216,   43,  72,  14,   8, -3,  2*0, &
         -2068,  3000, -2091,  800,  359,   64, -57,   6,  10, -3,  3*0, &
                 1611,  1278,  461,  262,   15,   1,  -2,   2,  2,  4*0, &
                         838, -395,  -42, -212,  14, -13, -12, -5,  5*0, &
                               234, -160,    2, -22,  -3,  10, -1,  6*0, &
                                     -56,    3,  -2,   5,  -1,  6,  7*0, &
                                          -112,  13,   0,   0,  4,  8*0, &
                                                 -2,  11,   3,  1,  9*0, &
                                                       3,   1,  0, 10*0, &
                                                           -1,  3, 11*0, &
                                                               -1, 13*0/
!          DGRF g(n,m) for 1975:
      DATA (G1DIM(I),I=2161,2304) /0, &
        -30100, -1902,  1276,  946, -218,   45,  71,  14,   7, -3,  2*0, &
         -2013,  3010, -2144,  791,  356,   66, -56,   6,  10, -3,  3*0, &
                 1632,  1260,  438,  264,   28,   1,  -1,   2,  2,  4*0, &
                         830, -405,  -59, -198,  16, -12, -12, -5,  5*0, &
                               216, -159,    1, -14,  -8,  10, -2,  6*0, &
                                     -49,    6,   0,   4,  -1,  5,  7*0, &
                                          -111,  12,   0,  -1,  4,  8*0, &
                                                 -5,  10,   4,  1,  9*0, &
                                                       1,   1,  0, 10*0, &
                                                           -2,  3, 11*0, &
                                                               -1, 13*0/
!          DGRF g(n,m) for 1980:
      DATA (G1DIM(I),I=2305,2448) /0, &
        -29992, -1997,  1281,  938, -218,   48,  72,  18,   5, -4,  2*0, &
         -1956,  3027, -2180,  782,  357,   66, -59,   6,  10, -4,  3*0, &
                 1663,  1251,  398,  261,   42,   2,   0,   1,  2,  4*0, &
                         833, -419,  -74, -192,  21, -11, -12, -5,  5*0, &
                               199, -162,    4, -12,  -7,   9, -2,  6*0, &
                                     -48,   14,   1,   4,  -3,  5,  7*0, &
                                          -108,  11,   3,  -1,  3,  8*0, &
                                                 -2,   6,   7,  1,  9*0, &
                                                      -1,   2,  2, 10*0, &
                                                           -5,  3, 11*0, &
                                                                0, 13*0/
!          DGRF g(n,m) for 1985:
      DATA (G1DIM(I),I=2449,2592) /0, &
        -29873, -2072,  1296,  936, -214,   53,  74,  21,   5, -4,  2*0, &
         -1905,  3044, -2208,  780,  355,   65, -62,   6,  10, -4,  3*0, &
                 1687,  1247,  361,  253,   51,   3,   0,   1,  3,  4*0, &
                         829, -424,  -93, -185,  24, -11, -12, -5,  5*0, &
                               170, -164,    4,  -6,  -9,   9, -2,  6*0, &
                                     -46,   16,   4,   4,  -3,  5,  7*0, &
                                          -102,  10,   4,  -1,  3,  8*0, &
                                                  0,   4,   7,  1,  9*0, &
                                                      -4,   1,  2, 10*0, &
                                                           -5,  3, 11*0, &
                                                                0, 13*0/
!          DGRF g(n,m) for 1990:
      DATA (G1DIM(I),I=2593,2736) /0, &
        -29775, -2131,  1314,  939, -214,   61,  77,  23,   4, -3,  2*0, &
         -1848,  3059, -2239,  780,  353,   65, -64,   5,   9, -4,  3*0, &
                 1686,  1248,  325,  245,   59,   2,  -1,   1,  2,  4*0, &
                         802, -423, -109, -178,  26, -10, -12, -5,  5*0, &
                               141, -165,    3,  -1, -12,   9, -2,  6*0, &
                                     -36,   18,   5,   3,  -4,  4,  7*0, &
                                           -96,   9,   4,  -2,  3,  8*0, &
                                                  0,   2,   7,  1,  9*0, &
                                                      -6,   1,  3, 10*0, &
                                                           -6,  3, 11*0, &
                                                                0, 13*0/
!          DGRF g(n,m) for 1995:
      DATA (G1DIM(I),I=2737,2880) /0, &
        -29682, -2197,  1329,  941, -210,   66,  78,  24,   4, -3,  2*0, &
         -1789,  3074, -2268,  782,  352,   64, -67,   4,   9, -4,  3*0, &
                 1685,  1249,  291,  237,   65,   1,  -1,   1,  2,  4*0, &
                         769, -421, -122, -172,  29,  -9, -12, -5,  5*0, &
                               116, -167,    2,   4, -14,   9, -2,  6*0, &
                                     -26,   17,   8,   4,  -4,  4,  7*0, &
                                           -94,  10,   5,  -2,  3,  8*0, &
                                                 -2,   0,   7,  1,  9*0, &
                                                      -7,   0,  3, 10*0, &
                                                           -6,  3, 11*0, &
                                                                0, 13*0/
!          DGRF h(n,m) for 1900:
      DATA (H1DIM(I),I=1,144) /13*0, &
          5922, -1061,  -330,  195, -210,   -9, -45,   8, -20,  2,  3*0, &
                 1121,     3,  -69,   53,   83, -13, -14,  14,  1,  4*0, &
                         523, -210,  -33,    2, -10,   7,   5,  2,  5*0, &
                               -75, -124,  -35,  -1, -13,  -3,  6,  6*0, &
                                       3,   36,  28,   5,  -2, -4,  7*0, &
                                           -69, -12,  16,   8,  0,  8*0, &
                                                -22,  -5,  10, -2,  9*0, &
                                                     -18,  -2,  4, 10*0, &
                                                            2,  0, 11*0, &
                                                               -6, 13*0/
!          DGRF h(n,m) for 1905:
      DATA (H1DIM(I),I=145,288) /13*0, &
          5909, -1086,  -357,  203, -193,   -7, -46,   8, -20,  2,  3*0, &
                 1065,    34,  -77,   56,   86, -14, -15,  14,  1,  4*0, &
                         480, -201,  -32,    4, -11,   7,   5,  2,  5*0, &
                               -65, -125,  -32,   0, -13,  -3,  6,  6*0, &
                                      11,   32,  28,   5,  -2, -4,  7*0, &
                                           -67, -12,  16,   8,  0,  8*0, &
                                                -22,  -5,  10, -2,  9*0, &
                                                     -18,  -2,  4, 10*0, &
                                                            2,  0, 11*0, &
                                                               -6, 13*0/
!          DGRF h(n,m) for 1910:
      DATA (H1DIM(I),I=289,432) /13*0, &
          5898, -1128,  -389,  211, -172,   -5, -47,   8, -20,  2,  3*0, &
                 1000,    62,  -90,   57,   89, -14, -15,  14,  1,  4*0, &
                         425, -189,  -33,    5, -12,   6,   5,  2,  5*0, &
                               -55, -126,  -29,   1, -13,  -3,  6,  6*0, &
                                      21,   28,  28,   5,  -2, -4,  7*0, &
                                           -65, -13,  16,   8,  0,  8*0, &
                                                -22,  -5,  10, -2,  9*0, &
                                                     -18,  -2,  4, 10*0, &
                                                            2,  0, 11*0, &
                                                               -6, 13*0/
!          DGRF h(n,m) for 1915:
      DATA (H1DIM(I),I=433,576) /13*0, &
          5875, -1191,  -421,  218, -148,   -2, -48,   8, -20,  2,  3*0, &
                  917,    84, -109,   58,   93, -14, -15,  14,  1,  4*0, &
                         360, -173,  -34,    8, -12,   6,   5,  2,  5*0, &
                               -51, -126,  -26,   2, -13,  -3,  6,  6*0, &
                                      32,   23,  28,   5,  -2, -4,  7*0, &
                                           -62, -15,  16,   8,  0,  8*0, &
                                                -22,  -5,  10, -2,  9*0, &
                                                     -18,  -2,  4, 10*0, &
                                                            2,  0, 11*0, &
                                                               -6, 13*0/
!          DGRF h(n,m) for 1920:
      DATA (H1DIM(I),I=577,720) /13*0, &
          5845, -1259,  -445,  220, -122,    0, -49,   8, -20,  2,  3*0, &
                  823,   103, -134,   58,   96, -14, -15,  14,  1,  4*0, &
                         293, -153,  -38,   11, -13,   6,   5,  2,  5*0, &
                               -57, -125,  -22,   4, -14,  -3,  6,  6*0, &
                                      43,   18,  28,   5,  -2, -4,  7*0, &
                                           -57, -16,  17,   9,  0,  8*0, &
                                                -22,  -5,  10, -2,  9*0, &
                                                     -19,  -2,  4, 10*0, &
                                                            2,  0, 11*0, &
                                                               -6, 13*0/
!          DGRF h(n,m) for 1925:
      DATA (H1DIM(I),I=721,864) /13*0, &
          5817, -1334,  -462,  216,  -96,    3, -50,   8, -20,  2,  3*0, &
                  728,   119, -163,   58,   99, -14, -15,  14,  1,  4*0, &
                         229, -130,  -44,   14, -14,   6,   5,  2,  5*0, &
                               -70, -122,  -18,   5, -14,  -3,  6,  6*0, &
                                      51,   13,  29,   5,  -2, -4,  7*0, &
                                           -52, -17,  17,   9,  0,  8*0, &
                                                -21,  -5,  10, -2,  9*0, &
                                                     -19,  -2,  4, 10*0, &
                                                            2,  0, 11*0, &
                                                               -6, 13*0/
!          DGRF h(n,m) for 1930:
      DATA (H1DIM(I),I=865,1008) /13*0, &
          5808, -1424,  -480,  205,  -72,    4, -51,   8, -20,  2,  3*0, &
                  644,   133, -195,   60,  102, -15, -15,  14,  1,  4*0, &
                         166, -109,  -53,   19, -14,   5,   5,  2,  5*0, &
                               -90, -118,  -16,   6, -14,  -3,  6,  6*0, &
                                      58,    8,  29,   5,  -2, -4,  7*0, &
                                           -46, -18,  18,   9,  0,  8*0, &
                                                -20,  -5,  10, -2,  9*0, &
                                                     -19,  -2,  4, 10*0, &
                                                            2,  0, 11*0, &
                                                               -6, 13*0/
!          DGRF h(n,m) for 1935:
      DATA (H1DIM(I),I=1009,1152) /13*0, &
          5812, -1520,  -494,  188,  -51,    4, -52,   8, -20,  2,  3*0, &
                  586,   146, -226,   64,  104, -17, -15,  15,  1,  4*0, &
                         101,  -90,  -64,   25, -14,   5,   5,  2,  5*0, &
                              -114, -115,  -15,   7, -15,  -3,  6,  6*0, &
                                      64,    4,  29,   5,  -3, -4,  7*0, &
                                           -40, -19,  18,   9,  0,  8*0, &
                                                -19,  -5,  11, -1,  9*0, &
                                                     -19,  -2,  4, 10*0, &
                                                            2,  0, 11*0, &
                                                               -6, 13*0/
!          DGRF h(n,m) for 1940:
      DATA (H1DIM(I),I=1153,1296) /13*0, &
          5821, -1614,  -499,  169,  -33,    4, -52,   8, -21,  2,  3*0, &
                  528,   163, -252,   71,  105, -18, -14,  15,  1,  4*0, &
                          43,  -72,  -75,   33, -14,   5,   5,  2,  5*0, &
                              -141, -113,  -15,   7, -15,  -3,  6,  6*0, &
                                      69,    0,  29,   5,  -3, -4,  7*0, &
                                           -33, -20,  19,   9,  0,  8*0, &
                                                -19,  -5,  11, -1,  9*0, &
                                                     -19,  -2,  4, 10*0, &
                                                            2,  0, 11*0, &
                                                               -6, 13*0/
!          DGRF h(n,m) for 1945:
      DATA (H1DIM(I),I=1297,1440) /13*0, &
          5810, -1702,  -499,  144,  -12,    6, -45,  12, -27,  5,  3*0, &
                  477,   186, -276,   95,  100, -18, -21,  17,  1,  4*0, &
                         -11,  -55,  -67,   16,   2, -12,  29,-20,  5*0, &
                              -178, -119,   -9,   6,  -7,  -9, -1,  6*0, &
                                      82,  -16,  28,   2,   4, -6,  7*0, &
                                           -39, -17,  18,   9,  6,  8*0, &
                                                -22,   3,   6, -4,  9*0, &
                                                     -11,   1, -2, 10*0, &
                                                            8,  0, 11*0, &
                                                               -2, 13*0/
!          DGRF h(n,m) for 1950:
      DATA (H1DIM(I),I=1441,1584) /13*0, &
          5815, -1810,  -476,  136,    3,   -1, -35,   5, -24, 13,  3*0, &
                  381,   206, -278,  103,   99, -17, -22,  19, -2,  4*0, &
                         -46,  -37,  -87,   33,   0,   0,  12,-10,  5*0, &
                              -210, -122,  -12,  10, -21,   2,  2,  6*0, &
                                      80,  -12,  36,  -8,   2, -3,  7*0, &
                                           -30, -18,  17,   8,  6,  8*0, &
                                                -16,  -4,   8, -3,  9*0, &
                                                     -17, -11,  6, 10*0, &
                                                           -7, 11, 11*0, &
                                                                8, 13*0/
!          DGRF h(n,m) for 1955:
      DATA (H1DIM(I),I=1585,1728) /13*0, &
          5820, -1898,  -462,  133,   15,   -9, -50,  10, -11, -4,  3*0, &
                  291,   216, -274,  110,   96, -24, -15,  12,  0,  4*0, &
                         -83,  -23,  -98,   48,  -4,   5,   7, -8,  5*0, &
                              -230, -121,  -16,   8, -23,   6, -2,  6*0, &
                                      78,  -12,  28,   3,  -2, -4,  7*0, &
                                           -24, -20,  23,  10,  1,  8*0, &
                                                -18,  -4,   7, -3,  9*0, &
                                                     -13,  -6,  7, 10*0, &
                                                            5, -1, 11*0, &
                                                               -3, 13*0/
!          DGRF h(n,m) for 1960:
      DATA (H1DIM(I),I=1729,1872) /13*0, &
          5791, -1967,  -414,  135,   16,  -10, -55,  11, -18,  4,  3*0, &
                  206,   224, -278,  125,   99, -28, -14,  12,  1,  4*0, &
                        -130,    3, -117,   60,  -6,   7,   2,  0,  5*0, &
                              -255, -114,  -20,   7, -18,   0,  2,  6*0, &
                                      81,  -11,  23,   4,  -3, -5,  7*0, &
                                           -17, -18,  23,   9,  1,  8*0, &
                                                -17,   1,   8, -1,  9*0, &
                                                     -20,   0,  6, 10*0, &
                                                            5,  0, 11*0, &
                                                               -7, 13*0/
!          DGRF h(n,m) for 1965:
      DATA (H1DIM(I),I=1873,2016) /13*0, &
          5776, -2016,  -404,  148,   19,  -11, -61,   7, -22,  2,  3*0, &
                  114,   240, -269,  128,  100, -27, -12,  15,  1,  4*0, &
                        -165,   13, -126,   68,  -2,   9,   7,  2,  5*0, &
                              -269,  -97,  -32,   6, -16,  -4,  6,  6*0, &
                                      81,   -8,  26,   4,  -5, -4,  7*0, &
                                            -7, -23,  24,  10,  0,  8*0, &
                                                -12,  -3,  10, -2,  9*0, &
                                                     -17,  -4,  3, 10*0, &
                                                            1,  0, 11*0, &
                                                               -6, 13*0/
!          DGRF h(n,m) for 1970:
      DATA (H1DIM(I),I=2017,2160) /13*0, &
          5737, -2047,  -366,  167,   26,  -12, -70,   7, -21,  1,  3*0, &
                   25,   251, -266,  139,  100, -27, -15,  16,  1,  4*0, &
                        -196,   26, -139,   72,  -4,   6,   6,  3,  5*0, &
                              -279,  -91,  -37,   8, -17,  -4,  4,  6*0, &
                                      83,   -6,  23,   6,  -5, -4,  7*0, &
                                             1, -23,  21,  10,  0,  8*0, &
                                                -11,  -6,  11, -1,  9*0, &
                                                     -16,  -2,  3, 10*0, &
                                                            1,  1, 11*0, &
                                                               -4, 13*0/
!          DGRF h(n,m) for 1975:
      DATA (H1DIM(I),I=2161,2304) /13*0, &
          5675, -2067,  -333,  191,   31,  -13, -77,   6, -21,  1,  3*0, &
                  -68,   262, -265,  148,   99, -26, -16,  16,  1,  4*0, &
                        -223,   39, -152,   75,  -5,   4,   7,  3,  5*0, &
                              -288,  -83,  -41,  10, -19,  -4,  4,  6*0, &
                                      88,   -4,  22,   6,  -5, -4,  7*0, &
                                            11, -23,  18,  10, -1,  8*0, &
                                                -12, -10,  11, -1,  9*0, &
                                                     -17,  -3,  3, 10*0, &
                                                            1,  1, 11*0, &
                                                               -5, 13*0/
!          DGRF h(n,m) for 1980:
      DATA (H1DIM(I),I=2305,2448) /13*0, &
          5604, -2129,  -336,  212,   46,  -15, -82,   7, -21,  1,  3*0, &
                 -200,   271, -257,  150,   93, -27, -18,  16,  0,  4*0, &
                        -252,   53, -151,   71,  -5,   4,   9,  3,  5*0, &
                              -297,  -78,  -43,  16, -22,  -5,  6,  6*0, &
                                      92,   -2,  18,   9,  -6, -4,  7*0, &
                                            17, -23,  16,   9,  0,  8*0, &
                                                -10, -13,  10, -1,  9*0, &
                                                     -15,  -6,  4, 10*0, &
                                                            2,  0, 11*0, &
                                                               -6, 13*0/
!          DGRF h(n,m) for 1985:
      DATA (H1DIM(I),I=2449,2592) /13*0, &
          5500, -2197,  -310,  232,   47,  -16, -83,   8, -21,  1,  3*0, &
                 -306,   284, -249,  150,   88, -27, -19,  15,  0,  4*0, &
                        -297,   69, -154,   69,  -2,   5,   9,  3,  5*0, &
                              -297,  -75,  -48,  20, -23,  -6,  6,  6*0, &
                                      95,   -1,  17,  11,  -6, -4,  7*0, &
                                            21, -23,  14,   9,  0,  8*0, &
                                                 -7, -15,   9, -1,  9*0, &
                                                     -11,  -7,  4, 10*0, &
                                                            2,  0, 11*0, &
                                                               -6, 13*0/
!          DGRF h(n,m) for 1990:
      DATA (H1DIM(I),I=2593,2736) /13*0, &
          5406, -2279,  -284,  247,   46,  -16, -80,  10, -20,  2,  3*0, &
                 -373,   293, -240,  154,   82, -26, -19,  15,  1,  4*0, &
                        -352,   84, -153,   69,   0,   6,  11,  3,  5*0, &
                              -299,  -69,  -52,  21, -22,  -7,  6,  6*0, &
                                      97,    1,  17,  12,  -7, -4,  7*0, &
                                            24, -23,  12,   9,  0,  8*0, &
                                                 -4, -16,   8, -2,  9*0, &
                                                     -10,  -7,  3, 10*0, &
                                                            2, -1, 11*0, &
                                                               -6, 13*0/
!          DGRF h(n,m) for 1995:
      DATA (H1DIM(I),I=2737,2880) /13*0, &
          5318, -2356,  -263,  262,   44,  -16, -77,  12, -19,  2,  3*0, &
                 -425,   302, -232,  157,   77, -25, -20,  15,  1,  4*0, &
                        -406,   98, -152,   67,   3,   7,  11,  3,  5*0, &
                              -301,  -64,  -57,  22, -21,  -7,  6,  6*0, &
                                      99,    4,  16,  12,  -7, -4,  7*0, &
                                            28, -23,  10,   9,  0,  8*0, &
                                                 -3, -17,   7, -2,  9*0, &
                                                     -10,  -8,  3, 10*0, &
                                                            1, -1, 11*0, &
                                                               -6, 13*0/
!          Initial coefficients g0 (IGRF for 2000):
      DATA G0 /0, &
        -29615, -2267,  1341,  935, -217,   72,  79,  25,   5, -2,  2*0, &
         -1728,  3072, -2290,  787,  351,   68, -74,   6,   9, -6,  3*0, &
                 1672,  1253,  251,  222,   74,   0,  -9,   3,  2,  4*0, &
                         715, -405, -131, -161,  33,  -8,  -8, -3,  5*0, &
                               110, -169,   -5,   9, -17,   6,  0,  6*0, &
                                     -12,   17,   7,   9,  -9,  4,  7*0, &
                                           -91,   8,   7,  -2,  1,  8*0, &
                                                 -2,  -8,   9,  2,  9*0, &
                                                      -7,  -4,  4, 10*0, &
                                                           -8,  0, 11*0, &
                                                               -1, 13*0/
!          D_/Dtime coefficients gt (IGRF for 2000-2005):
      DATA GT /0, &
          14.6, -12.4,   0.7, -1.3,  0.0,  1.0,-0.4,-0.3, 0.0,0.0,  2*0, &
          10.7,   1.1,  -5.4,  1.6, -0.7, -0.4,-0.4, 0.2, 0.0,0.0,  3*0, &
                 -1.1,   0.9, -7.3, -2.1,  0.9,-0.3,-0.3, 0.0,0.0,  4*0, &
                        -7.7,  2.9, -2.8,  2.0, 1.1, 0.4, 0.0,0.0,  5*0, &
                              -3.2, -0.8, -0.6, 1.1,-1.0, 0.0,0.0,  6*0, &
                                     2.5, -0.3,-0.2, 0.3, 0.0,0.0,  7*0, &
                                           1.2, 0.6,-0.5, 0.0,0.0,  8*0, &
                                               -0.9,-0.7, 0.0,0.0,  9*0, &
                                                    -0.4, 0.0,0.0, 10*0, &
                                                          0.0,0.0, 11*0, &
                                                              0.0, 13*0/
!          Initial coefficients h0 (IGRF for 2000):
      DATA H0 /13*0, &
          5186, -2478,  -227,  272,   44,  -17, -65,  12, -20,  1,  3*0, &
                 -458,   296, -232,  172,   64, -24, -22,  13,  0,  4*0, &
                        -492,  119, -134,   65,   6,   8,  12,  4,  5*0, &
                              -304,  -40,  -61,  24, -21,  -6,  5,  6*0, &
                                     107,    1,  15,  15,  -8, -6,  7*0, &
                                            44, -25,   9,   9, -1,  8*0, &
                                                 -6, -16,   4, -3,  9*0, &
                                                      -3,  -8,  0, 10*0, &
                                                            5, -2, 11*0, &
                                                               -8, 13*0/
!          D_/Dtime coefficients ht (IGRF for 2000-2005):
      DATA HT /13*0, &
         -22.5, -20.6,   6.0,  2.1, -0.1, -0.2, 1.1, 0.1, 0.0,0.0,  3*0, &
                 -9.6,  -0.1,  1.3,  0.6, -1.4, 0.0, 0.0, 0.0,0.0,  4*0, &
                       -14.2,  5.0,  1.7,  0.0, 0.3, 0.0, 0.0,0.0,  5*0, &
                               0.3,  1.9, -0.8,-0.1, 0.3, 0.0,0.0,  6*0, &
                                     0.1,  0.0,-0.6, 0.6, 0.0,0.0,  7*0, &
                                           0.9,-0.7,-0.4, 0.0,0.0,  8*0, &
                                                0.2, 0.3, 0.0,0.0,  9*0, &
                                                     0.7, 0.0,0.0, 10*0, &
                                                          0.0,0.0, 11*0, &
                                                              0.0, 13*0/

 
!          Trap out of range date:
      IF (DATE .LT. EPOCH(1)) GO TO 9100
      IF (DATE .GT. EPOCH(NYT)+5.) WRITE(6,9200) DATE
 
      DO 100 I=1,NYT
      IF (DATE .LT. EPOCH(I)) GO TO 110
      IY = I
  100 CONTINUE
  110 CONTINUE
 
      TIME = DATE
      T = TIME-EPOCH(IY)
      G(1)  = 0.0
      GV(1) = 0.0
      I = 2
      F0 = -1.0D-5
      DO 200 N=1,NMAX
      F0 = F0 * REAL(N)/2.
      F  = F0 / SQRT(2.0)
      NN = N+1
      MM = 1
      IF (IY .EQ. NYT) THEN
!          Extrapolate coefficients
        G(I) = ((GTT(NN,MM)*T + GT(NN,MM))*T + G0(NN,MM)) * F0
      ELSE
!          Interpolate coefficients
        G(I) = (GYR(NN,MM,IY) + &
               T/5.0 * (GYR(NN,MM,IY+1)-GYR(NN,MM,IY))) * F0
      ENDIF
      GV(I) = G(I) / REAL(NN)
      I = I+1
      DO 200 M=1,N
      F = F / SQRT( REAL(N-M+1) / REAL(N+M) )
      NN = N+1
      MM = M+1
      I1 = I+1
      IF (IY .EQ. NYT) THEN
!          Extrapolate coefficients
        G(I)  = ((GTT(NN,MM)*T + GT(NN,MM))*T + G0(NN,MM)) * F
        G(I1) = ((HTT(NN,MM)*T + HT(NN,MM))*T + H0(NN,MM)) * F
      ELSE
!          Interpolate coefficients
        G(I)  = (GYR(NN,MM,IY) + &
                T/5.0 * (GYR(NN,MM,IY+1)-GYR(NN,MM,IY))) * F
        G(I1) = (HYR(NN,MM,IY) + &
                T/5.0 * (HYR(NN,MM,IY+1)-HYR(NN,MM,IY))) * F
      ENDIF
      RNN = REAL(NN)
      GV(I)  = G(I)  / RNN
      GV(I1) = G(I1) / RNN
  200 I = I+2
 
      RETURN
 
!          Error trap diagnostics:
 9100 WRITE(6,'('' '',/, &
      & '' COFRM:  DATE'',F9.3,'' preceeds DGRF coefficients'', &
      &                         '' presently coded.'')') DATE
      STOP 'mor cod'
 9200 FORMAT(' ',/, &
      & ' COFRM:  DATE',F9.3,' is after the maximum', &
      &                         ' recommended for extrapolation.')
      END
 
 
      SUBROUTINE FELDG2 (IENTY,GLAT,GLON,ALT, BNRTH,BEAST,BDOWN,BABS)
!          Compute the DGRF/IGRF field components at the point GLAT,GLON,ALT.
!          COFRM must be called to establish coefficients for correct date
!          prior to calling FELDG.
!
!          IENTY is an input flag controlling the meaning and direction of the
!                remaining formal arguments:
!          IENTY = 1
!            INPUTS:
!              GLAT = Latitude of point (deg)
!              GLON = Longitude (east=+) of point (deg)
!              ALT  = Ht of point (km)
!            RETURNS:
!              BNRTH  north component of field vector (Gauss)
!              BEAST  east component of field vector  (Gauss)
!              BDOWN  downward component of field vector (Gauss)
!              BABS   magnitude of field vector (Gauss)
!
!          IENTY = 2
!            INPUTS:
!              GLAT = X coordinate (in units of earth radii 6371.2 km )
!              GLON = Y coordinate (in units of earth radii 6371.2 km )
!              ALT  = Z coordinate (in units of earth radii 6371.2 km )
!            RETURNS:
!              BNRTH = X component of field vector (Gauss)
!              BEAST = Y component of field vector (Gauss)
!              BDOWN = Z component of field vector (Gauss)
!              BABS  = Magnitude of field vector (Gauss)
!          IENTY = 3
!            INPUTS:
!              GLAT = X coordinate (in units of earth radii 6371.2 km )
!              GLON = Y coordinate (in units of earth radii 6371.2 km )
!              ALT  = Z coordinate (in units of earth radii 6371.2 km )
!            RETURNS:
!              BNRTH = Dummy variable
!              BEAST = Dummy variable
!              BDOWN = Dummy variable
!              BABS  = Magnetic potential (T.m)
!
!          INPUT from COFRM through COMMON /MAGCOF2/ NMAX,GB(144),GV(144),ICHG
!            NMAX = Maximum order of spherical harmonic coefficients used
!            GB   = Coefficients for magnetic field calculation
!            GV   = Coefficients for magnetic potential calculation
!            ICHG = Flag indicating when GB,GV have been changed
!
!          HISTORY:
!          COFRM and FELDG originated 15 Apr 83 by Vincent B. Wickwar
!          (formerly at SRI. Int., currently at Utah State).
!
!          May 94 (A.D. Richmond): Added magnetic potential calculation
!
!          Oct 95 (Barnes): Added ICHG
 
      implicit none
      integer ichg, IENTYP, is, nmax
      integer ihmax, last, imax , m
      integer MK, I ,IH , ienty, K, IL
      integer ihm , ilm 
      REAL(kind=8) RTOD
      REAL(kind=8) DTOR
      REAL(kind=8) RE
      REAL(kind=8) REQ
      REAL(kind=8) gb(144)
      REAL(kind=8) gv(144)
      REAL(kind=8) XI(3)
      REAL(kind=8) H(144)
      REAL(kind=8) G(144)
      REAL(kind=8) RLAT
      REAL(kind=8) GLAT
      REAL(kind=8) GLON
      REAL(kind=8) ALT
      REAL(kind=8) CT
      REAL(kind=8) ST
      REAL(kind=8) RLON
      REAL(kind=8) CP
      REAL(kind=8) SP
      REAL(kind=8) XXX
      REAL(kind=8) YYY
      REAL(kind=8) ZZZ
      REAL(kind=8) RQ
      REAL(kind=8) F
      REAL(kind=8) X , Y , Z
      REAL(kind=8) S , T
      REAL(kind=8) BXXX , BYYY , BZZZ
      REAL(kind=8) BABS
      REAL(kind=8) BEAST
      REAL(kind=8) BRHO
      REAL(kind=8) BNRTH
      REAL(kind=8) BDOWN

      PARAMETER (RTOD=57.2957795130823, DTOR=0.01745329251994330, &
                 RE=6371.2, REQ=6378.160)

      COMMON /MAGCOF2/GB,GV,ICHG,NMAX
      SAVE IENTYP, G
      DATA IENTYP/-10000/
 
      IF (IENTY .EQ. 1) THEN
        IS=1
        RLAT = GLAT*DTOR
        CT   = SIN(RLAT)
        ST   = COS(RLAT)
        RLON = GLON*DTOR
        CP   = COS(RLON)
        SP   = SIN(RLON)
        CALL GD2CART2 (GLAT,GLON,ALT,XXX,YYY,ZZZ)
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
      H(ILM+1) = G(ILM+1)+ Z*H(IHM+1) + X*(H(IHM+3)-H(IHM-1)) &
                                              -Y*(H(IHM+2)+H(IHM-2))
   70 H(ILM)   = G(ILM)  + Z*H(IHM)   + X*(H(IHM+2)-H(IHM-2)) &
                                              +Y*(H(IHM+3)+H(IHM-1))

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
 
!          Magnetic potential computation makes use of the fact that the
!          calculation of V is identical to that for r*Br, if coefficients
!          in the latter calculation have been divided by (n+1) (coefficients
!          GV).  Factor .1 converts km to m and gauss to tesla.
      IF (IENTY.EQ.3) BABS = (BXXX*XXX + BYYY*YYY + BZZZ*ZZZ)*RE*.1
 
      RETURN
      END
 
      SUBROUTINE GD2CART2 (GDLAT,GLON,ALT,X,Y,Z)
!          Convert geodetic to cartesian coordinates by calling CONVRT
!          940503 A. D. Richmond

      implicit none
      REAL(kind=8) dtor
      REAL(kind=8) GDLAT,GLON,ALT,RHO,Z
      REAL(kind=8) ANG,X,Y
      PARAMETER (DTOR=0.01745329251994330)

      CALL CONVRT2 (1,GDLAT,ALT,RHO,Z)

      ANG = GLON*DTOR
      X = RHO*COS(ANG)
      Y = RHO*SIN(ANG)

      RETURN
      END
 
      SUBROUTINE CONVRT2 (I,GDLAT,ALT,X1,X2)
!          Convert space point from geodetic to geocentric or vice versa.
!
!          I is an input flag controlling the meaning and direction of the
!            remaining formal arguments:
!
!          I = 1  (convert from geodetic to cylindrical)
!            INPUTS:
!              GDLAT = Geodetic latitude (deg)
!              ALT   = Altitude above reference ellipsoid (km)
!            RETURNS:
!              X1    = Distance from Earth's rotation axis (km)
!              X2    = Distance above (north of) Earth's equatorial plane (km)
!
!          I = 2  (convert from geodetic to geocentric spherical)
!            INPUTS:
!              GDLAT = Geodetic latitude (deg)
!              ALT   = Altitude above reference ellipsoid (km)
!            RETURNS:
!              X1    = Geocentric latitude (deg)
!              X2    = Geocentric distance (km)
!
!          I = 3  (convert from cylindrical to geodetic)
!            INPUTS:
!              X1    = Distance from Earth's rotation axis (km)
!              X2    = Distance from Earth's equatorial plane (km)
!            RETURNS:
!              GDLAT = Geodetic latitude (deg)
!              ALT   = Altitude above reference ellipsoid (km)
!
!          I = 4  (convert from geocentric spherical to geodetic)
!            INPUTS:
!              X1    = Geocentric latitude (deg)
!              X2    = Geocentric distance (km)
!            RETURNS:
!              GDLAT = Geodetic latitude (deg)
!              ALT   = Altitude above reference ellipsoid (km)
!
!
!          HISTORY:
!          940503 (A. D. Richmond):  Based on a routine originally written
!          by V. B. Wickwar.
!
!          REFERENCE:  ASTRON. J. VOL. 66, p. 15-16, 1961.

      implicit none
 
      integer I

      REAL(kind=8) rtod
      REAL(kind=8) dtor
      REAL(kind=8) RE
      REAL(kind=8) REQ
      REAL(kind=8) FLTNVRS
      REAL(kind=8) E2,E4,E6,E8
      REAL(kind=8) OME2REQ
      REAL(kind=8) A21,A22,A23,A41,A42,A43,A44
      REAL(kind=8) A61,A62,A63,A81,A82,A83,A84
      REAL(kind=8) SINLAT
      REAL(kind=8) COSLAT
      REAL(kind=8) D
      REAL(kind=8) Z
      REAL(kind=8) RHO
      REAL(kind=8) X1
      REAL(kind=8) X2
      REAL(kind=8) RKM
      REAL(kind=8) SCL
      REAL(kind=8) GCLAT
      REAL(kind=8) RI
      REAL(kind=8) A2
      REAL(kind=8) A4
      REAL(kind=8) A6
      REAL(kind=8) A8
      REAL(kind=8) CCL
      REAL(kind=8) S2CL
      REAL(kind=8) C2CL
      REAL(kind=8) S4CL
      REAL(kind=8) C4CL
      REAL(kind=8) S8CL
      REAL(kind=8) S6CL
      REAL(kind=8) DLTCL
      REAL(kind=8) GDLAT
      REAL(kind=8) SGL
      REAL(kind=8) ALT


      PARAMETER (RTOD=57.2957795130823, DTOR=0.01745329251994330 , &
        RE=6371.2 , REQ=6378.160 , FLTNVRS=298.25 , &
        E2=(2.-1./FLTNVRS)/FLTNVRS , E4=E2*E2 , E6=E4*E2 , E8=E4*E4 , &
        OME2REQ = (1.-E2)*REQ , &
           A21 =     (512.*E2 + 128.*E4 + 60.*E6 + 35.*E8)/1024. , &
           A22 =     (                        E6 +     E8)/  32. , &
           A23 = -3.*(                     4.*E6 +  3.*E8)/ 256. , &
           A41 =    -(           64.*E4 + 48.*E6 + 35.*E8)/1024. , &
           A42 =     (            4.*E4 +  2.*E6 +     E8)/  16. , &
           A43 =                                   15.*E8 / 256. , &
           A44 =                                      -E8 /  16. , &
           A61 =  3.*(                     4.*E6 +  5.*E8)/1024. , &
           A62 = -3.*(                        E6 +     E8)/  32. , &
           A63 = 35.*(                     4.*E6 +  3.*E8)/ 768. , &
           A81 =                                   -5.*E8 /2048. , &
           A82 =                                   64.*E8 /2048. , &
           A83 =                                 -252.*E8 /2048. , &
           A84 =                                  320.*E8 /2048. )
!          E2 = Square of eccentricity
 
      IF (I .GE. 3) GO TO 300
 
!          Geodetic to geocentric
 
!          Compute RHO,Z
      SINLAT = SIN(GDLAT*DTOR)
      COSLAT = SQRT(1.-SINLAT*SINLAT)
      D      = SQRT(1.-E2*SINLAT*SINLAT)
      Z      = (ALT+OME2REQ/D)*SINLAT
      RHO    = (ALT+REQ/D)*COSLAT
      X1 = RHO
      X2 = Z
      IF (I .EQ. 1) RETURN
 
!          Compute GCLAT,RKM
      RKM   = SQRT(Z*Z + RHO*RHO)
      GCLAT = RTOD*ATAN2(Z,RHO)
      X1 = GCLAT
      X2 = RKM
      RETURN
 
!          Geocentric to geodetic
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
