       program t_dat2qc_2
c
c=======================================================================
c
c  This program makes quality control on temperature data at the model
c  levels for SODA and converts T-insitu to T-potential
c
c=======================================================================
c
c
       parameter (idim=360, jdim=180, kdim=75, imn=12)    ! dimension of climatology
       parameter (kmb=55)                                 ! N of assimilation levels
       parameter (bmiss=-999.99)                          ! missing values
c
       character typ*3, yst*4, ynd*4, stamp*32
       real, dimension(jdim) :: glat
       real, dimension(idim) :: glon
       real, dimension(idim, jdim,kdim,imn) :: tclim, sclim
       real, dimension(kmb) :: lev, prof
c
c--- Levels for SODA-4 MOM5.1
       data lev/
     &  0.54128076539161, 1.68073467983143, 2.93995264891447, 
     &  4.33152148514951, 5.86935042405403, 7.56880992005022, 
     &  9.44688495964871, 11.5223443928033, 13.8159279323332, 
     &  16.3505526329353, 19.151540835806,  22.2468717521887, 
     &  25.66745905867,   29.4474570874552, 33.6245984138505,  
     &  38.2405658667745, 43.3414022110416, 48.9779609631741, 
     &  55.2064019936584, 62.0887357201759, 69.6934197853445, 
     &  78.0960121026315, 87.379884002135,  97.6369968516802, 
     &  108.968744882098, 121.48686590184,  135.314419992604, 
     &  150.586833952313, 167.453005943493, 186.076460191499, 
     &  206.636535284236, 229.329581189748, 254.370128943059, 
     &  281.991982454696, 312.449163513683, 346.01661821764,  
     &  382.990565399297, 423.688335745048, 468.447515403956, 
     &  517.624172536375, 571.589915322663, 630.727513161329, 
     &  695.424821496212, 766.066799209776, 843.025512765185, 
     &  926.648198508397, 1017.24370768399, 1115.06797498371, 
     &  1220.30949515272, 1333.07609257785, 1453.38443981473, 
     &  1581.15373742821, 1716.2046509593,  1858.26402924282, 
     &  2006.97519804695 / 
       
c 2161.9128824586  2322.60123413657 2488.5331512273  2659.18911348314 2834.05406813066 
c 3012.6313812063  3194.4533865926  3379.08852049285 3566.14536241277 3755.27409765924 
c 3946.16599026663 4138.55143999841 4332.19712659193 4526.90264998742 4722.49697666354 
c 4918.83491240816 5115.7937470637  5313.27015821824 5511.17741767395 5709.44291422576 

c--- Levels for SODA-MOM-5.0.2
c       data lev/5.03355, 15.10065, 25.21935, 35.35845, 45.57635,
c     &         55.86325, 66.26175, 76.80285, 87.57695, 98.62325,
c     &         110.0962, 122.1067, 134.9086, 148.7466, 164.0538,
c     &         181.3125, 201.2630, 224.7773, 253.0681, 287.5508,
c     &         330.0078, 382.3651, 446.7263, 524.9824, 618.7031,
c     &         728.6921, 854.9935, 996.7153, 1152.376, 1319.997,
c     &         1497.562, 1683.057, 1874.788, 2071.252/
c
c
c--- open temperature climatology file
       open(31, file = '../CLIMATOLOGY/WODB/tpot_clim_fill_75levs.grd',
     &      access = 'direct',
     &      form = 'unformatted',
     &      status = 'old',
     &      recl = idim*jdim*4)
c--- open salinity climatology file
       open(32, file = '../CLIMATOLOGY/WODB/s_clim_fill_75levs.grd',
     &      access = 'direct',
     &      form = 'unformatted',
     &      status = 'old',
     &      recl = idim*jdim*4)
c
c--- read model levels T-climatology
       do im = 1,imn
         do k = 1,kdim
	   irec = k+(im-1)*kdim
c	   write(6,*) 'rec= ', irec,'  len= ',idim*jdim
	   read(31,rec=irec) tclim(1:idim,1:jdim,k,im)
c	   write(6,*) k, tclim(183:185,90,k,im)
	 enddo
       enddo
       close (31)
c--- read model levels S-climatology
       do im = 1,imn
         do k = 1,kdim
	   irec = k+(im-1)*kdim
	   read(32,rec=irec) sclim(1:idim,1:jdim,k,im)
c	   write(6,*) k, sclim(183:185,90,k,im)
	 enddo
       enddo
       close (32)
c
c--- get 1x1 grid longitudes and latitudes
       do i = 1,idim
          glon(i) = 0.5+(i-1)*1.0
       enddo
       do j = 1,jdim
          glat(j) = -89.5+(j-1)*1.0
       enddo
c       
       do iyr = 2023, 2024                                ! years loop
         write(yst,'(i4)') iyr                            ! makes year as character
	 write(6,*) '--- For year:', iyr
	 idst = jday(1,1,iyr)                             ! Gregorian day to start
	 idnd = jday(12,31,iyr)                           ! Gregorian day to stop
c
c--- bring counts to zero
         ntr = 0                                          ! checked profiles
	 ntw = 0                                          ! writen good profiles
	 ndb = 0                                          ! bad profiles
	 nbd1 = 0                                         ! bad position
	 nbd2 = 0                                         ! 0.0E,0.0N point
	 nbd3 = 0                                         ! temperature > 33.0
	 nbd4 = 0                                         ! temperature < -4.0
	 nbd5 = 0                                         ! departure from climatology > 12.0
	 nbd6 = 0                                         ! all profile RMS from climatology > 7.5
	 nbd7 = 0                                         ! departure from climatology > 5.0 in deep ocean (>500m)
	 nbd8 = 0                                         ! vertical temperature gradient is > 0.5
c
c--- open temperature profiles interpolated to model levels file
         open(10, 
     &	 file = 'TMP/t_'//yst//'_55levs.dat',
c     &	 file = '../STEP_1/t_'//yst//'.dat',
     &        status = 'old')
c--- open output files for usable and bad data
c         open(20, file = '../STEP_QC/qc_tpot_'//yst//'.dat',
c     &        status = 'unknown')
c         open(21, file = '../STEP_QC/qc_tpot_'//yst//'.bad',
c     &        status = 'unknown')     
c
c--- open temperature profiles interpolated to model levels file
c         open(10, file = 'TMP/t_'//yst//'.dat',
c     &        status = 'old')
c--- open output files for usable and bad data
         open(20, file = 'TMP/qc_tpot_'//yst//'_55levs.dat',
     &        status = 'unknown')
         open(21, file = 'TMP/qc_tpot_'//yst//'_55levs.bad',
     &        status = 'unknown')
c
c--- main loop to read temperature profiles at model levels
         do ipr = 1,99999999             	 
           read(10,800,end=123) typ,istn,iyear,month,
     &                          iday,rlat,rlon,llv,
     &                          (prof(k), k=1,llv)
c           write(6,*) 'Typ - ',typ,',   stn# = ',istn
c
           if(llv.gt.kmb) llv=kmb                         ! control number of assimilation levels
c
           if(month.gt.12.or.month.lt.1) go to 50         ! control month number
           if(iday.gt.31.or.iday.lt.0) go to 50           ! control day number
           ind=jday(month,iday,iyear)                     ! convert day to Gregorian day number
           if(ind.lt.idst.or.ind.gt.idnd) go to 50        ! control the day is in current year
           ntr=ntr+1                                      ! count checked profiles
c
           if(rlon.lt.0.0) rlon=rlon+360.0                ! convert longitude to 0-360 values
c
           ii=indp(rlon,glon,idim)                         ! find climatology i-index for profile 
           jj=indp(rlat,glat,jdim)                         ! find climatology j-index for profile
c
           lmiss = 0
	   do k = 1,llv
	     if (prof(k).eq.bmiss) lmiss = lmiss+1
	   enddo
	   if (lmiss.eq.llv) goto 50                      ! check for any data in profile
c
  753      continue
           if(prof(llv).eq.bmiss) then                    ! find the deepest profile level
             llv=llv-1
             go to 753
           endif
c	   
           if (rlat.lt.-90.00.or.rlat.gt.90.00.
     &         or.rlon.lt.0.0.or.rlon.gt.360) then        ! check profile position
              nbd1=nbd1+1
              ibad=1
              go to 110
           endif
c
           if (rlon.eq.0.0.and.rlat.eq.0.0) then          ! remove 0.0E,0.0N point
             nbd2=nbd2+1
             ibad=2
             go to 110
           endif
c
c -- start profile checking
           rms=0.0
           kk=kmb
           if(llv.le.kmb) kk=llv
           kls=0
c           
           do 700 k=1,kk
             if(prof(k).eq.bmiss) go to 700               ! no obs. at this level
c
             if(prof(k).gt.33.0) then                     ! temperature is too high 
               nbd3=nbd3+1
               ibad=3
               go to 110
             endif
c
             if(prof(k).le.-4.0.or.
     &          tclim(ii,jj,k,month).le.-4.0) then         ! temperature is too low  
               nbd4=nbd4+1
               ibad=4
               go to 110
             endif
c
             atp=abs(prof(k)-tclim(ii,jj,k,month))
             rms=rms+atp**2
             kls=kls+1
c
             if(atp.gt.12.0) then                         ! obs. is too far from climatology
               nbd5=nbd5+1
               ibad=5
               go to 110
             endif
c
700        continue
c
           rms=rms/float(kls)
c
           if((sqrt(rms).gt.7.5).and.(kls.ge.5)) then     ! all profile is too far from climatology
             nbd6=nbd6+1
             ibad=6
             go to 110
           endif
c
           atp=abs(prof(kk)-tclim(ii,jj,kk,month))         ! obs. is too far from climatology under termocline
           if(kk.gt.40.and.atp.gt.5.0) then
             nbd7=nbd7+1
             ibad=7
             go to 110
           endif
c
           do 900 k=2,kk
             if((prof(k-1).eq.bmiss).or.(prof(k).eq.bmiss)) go to 900
             grad=(prof(k)-prof(k-1))/(lev(k)-lev(k-1))   ! temperature gradient is too high
             if(grad.gt.0.5) then                         
               nbd8=nbd8+1
               ibad=8
               go to 110
             endif
900        continue
c
c -- calculate T-potential
           do k = 1,llv
	     if(prof(k).ne.bmiss) then
 	       bar = lev(k)/9.93117
               prof(k) = tpot(sclim(ii,jj,k,month),prof(k),bar)
	     endif
	   enddo
c
c -- write out good data
           write(20,800) typ,istn,iyear,month,iday,rlat,rlon,
     &                   llv,(prof(kg), kg=1,llv)
           ntw=ntw+1
	   goto 50
c
c -- write out bad data
110        continue
           write(21,830) ibad,typ,istn,iyear,month,iday,rlat,rlon,
     &                   llv,(prof(kb), kb=1,llv)
           write(21,831) ii,jj,month,(tclim(ii,jj,kb,month), kb=1,llv)
c	   write(21,*) kk, k
           ndb=ndb+1
c
50         continue
c
         enddo
c
         close (10)
         close (20)
         close (21)
c
123      continue
c
         print *, '- total read & write & bad:', ntr, ntw, ndb
         print *, '- BAD data'
         print *, nbd1,nbd2,nbd3,nbd4,nbd5,nbd6,nbd7,nbd8
c
       enddo
c	
c       
c--- formats to read and write
800    format(a4,i8,i5,2i3,2f9.3,i3,55f8.3)
830    format('B',i1,a4,i8,i5,2i3,2f9.3,i3,55f8.3)
831    format('CL',2i4,i3,33x,55f8.3)
c       
       end
c
c
c=======================================================================
c
c
      function jday(mon,iday,iyr)
c
c=======================================================================
c     compute the julian day corresponding to the
c     day on the gregorian calender
c=======================================================================
c
      dimension dpm(12)
      data dpm /31.0, 28.0, 31.0, 30.0, 31.0, 30.0, 31.0, 31.0, 30.0,
     $          31.0, 30.0, 31.0/

        dpm(2) = 28.0
        if(mod(real(iyr),4.) .eq. 0.)dpm(2) = 29.0
c
c first calculate days without leap years
c
ccao??11/8/95  for 1950 XXX
ccao    iyrs = iyr-1950
ccao    days = 3282.0+real(iyrs)*365.0
c
        iyrs = iyr-1970
        days = 587.0+real(iyrs)*365.0
ccao?? need to think about iyrs+?? e.x.: iyr-1950, should iyrs+1
ccaojday num_leap = int(real(iyrs+1)/4.)
        num_leap = floor(real(iyrs+1)/4.)
        days = days + real(num_leap)
c
c now sum up for days this year
c
        sum = 0.
        if(mon .gt. 1)then
        do l = 1, mon-1
        sum = sum + dpm(l)
        end do
        days = days + sum
        end if
c
        jday = int(days) + iday
cassim
      return
      end
c
      function indp (value, array, ia)
c
c======================================================
c
c     indp = index of nearest data point within "array" corresponding to
c            "value".
c
c     inputs:
c
c     value  = arbitrary data...same units as elements in "array"
c     array  = array of data points  (must be monotonically increasing)
c     ia     = dimension of "array"
c
c     output:
c
c     indp =  index of nearest data point to "value"
c             if "value" is outside the domain of "array" then indp = 1
c             or "ia" depending on whether array(1) or array(ia) is
c             closest to "value"
c
c             note: if "array" is dimensioned array(0:ia) in the calling
c                   program, then the returned index should be reduced
c                   by one to account for the zero base.
c
c     example:
c
c     let model depths be defined by the following:
c     parameter (km=5)
c     dimension z(km)
c     data z /5.0, 10.0, 50.0, 100.0, 250.0/
c
c     k1 = indp (12.5, z, km)
c     k2 = indp (0.0, z, km)
c
c     k1 would be set to 2, & k2 would be set to 1 so that
c     z(k1) would be the nearest data point to 12.5 and z(k2) would
c     be the nearest data point to 0.0
c
c======================================================
c
      dimension array(ia)
c
      do i=2,ia
        if (array(i) .lt. array(i-1)) then
          print *,            
     &   ' => Error: array must be monotonically increasing in "indp"'
     &,  '           when searching for nearest element to value=',value
          print *, '           array(i) < array(i-1) for i=',i
          print *, '           array(i) for i=1..ia follows:'
          do ii=1,ia
            print *, 'i=',ii, ' array(i)=',array(ii)
          enddo
          stop '=>indp'
        endif
      enddo
      if (value .lt. array(1) .or. value .gt. array(ia)) then
        if (value .lt. array(1))  indp = 1
        if (value .gt. array(ia)) indp = ia
        return
      else
        do i=2,ia
          if (value .le. array(i)) then
            indp = i
            if (array(i)-value .gt. value-array(i-1)) indp = i-1
            go to 101
          endif
        enddo
101     continue
      endif
      return
      end
c
c
c======================================================
c
      function tpot(S,T,P)	    
c
c======================================================
c
c   This function calculates potential temparature
c   from:
c   S - salinity (PSU);
c   T - insitu temparature (degree C)
c   P - pressure (bar)
c
c======================================================
c
      a0 = 3.6504e-4
      a1 = 8.3198e-5
      a2 = 5.4065e-7
      a3 = 4.0274e-9

      b0 = 1.7439e-5
      b1 = 2.9778e-7

      c0 = 8.9309e-7
      c1 = 3.1628e-8
      c2 = 2.1987e-10

      d0 = 4.1057e-9
      d1 = 1.6056e-10
      d2 = 5.0484e-12
      
      tpot = T - P *(a0+a1* T -a2* T * T +a3* T * T * T ) 
      tpot = tpot - P *( S -35)*(b0-b1* T ) 
      tpot = tpot - P * P *(c0-c1* T +c2* T * T ) 
      tpot = tpot +d0*( S -35)* P * P - P * P * P *(-d1+d2* T ) 

      return
      end
c
c
c======================================================
c
      



 
