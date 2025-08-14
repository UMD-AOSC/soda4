      PROGRAM wodFOR
      
c
c  Gennady Chepurin 12/23/2004
c
c  NOTE - check iyst & typ first !!!!!!!!!!!
c
c     This main program (wodFOR)
c     calls the subroutine OCLread (WODREAD200X if the data file are in WOD01(05)
c     format or the subroutine OCLread1998 if the data are in WOD98 format).
c     These subroutines do the actual reading of the ASCII format, and load
c     the data into arrays which are passed back to the main program.
c     
c     It is intended that the subroutine OCLread provides an example of how
c     to extract the data and variables from the ASCII format,
c     whereas the main part of the wodFOR program provides an example of how
c     these data can be made accessible/workable as a series of arrays.
c
c     Missing values used in this dataset = bmiss = -999.99
c
c***********************************************************
c     
c   Parameters (constants):
c     
c     maxlevel  - maximum number of depth levels, also maximum
c                   number of all types of variables
c     maxcalc   - maximum number of measured and calculated
c                   depth dependent variables
c     kdim      - number of standard depth levels
c     bmiss     - binary missing value marker
c     maxtcode  - maximum number of different taxa variable codes
c     maxtax    - maximum number of taxa sets
c     maxpinf   - maximum number of distinct measured variable
c                 information codes
c
c******************************************************************
      
ccaostd
      parameter (iyst=2023, iynd=2024, kmb=55)
c      parameter (iyst=2018, iynd=2018, kmb=34)
c      parameter (iyst=2008, iynd=2009, kmb=23)
c      parameter (iyst=1951, iynd=2001, kmb=40)
c
      parameter (maxlevel=30000, maxcalc=100)
      parameter (kdim=40, bmiss=-999.99)
      parameter (maxtcode=25, maxtax=2000)
      parameter (maxpinf=25)
      
c******************************************************************
c
c   Character Arrays:
c
c     cc        - NODC country code
c     chars     - OCL character data: 1. originators cruise code,
c                                     2. originators station cod
c     filename  - file name
c
c*****************************************************************
      
      character*2  cc
      character*15 chars(2)
ccao
      character*80 filename
      character toutfile*21, soutfile*21, yr*4, typ*3
cgena  spr* files for tests
c      character tprfile*14, sprfile*14
c
      character dum*80
 
c******************************************************************
c
c   Arrays:
c
c     isig()    - number of significant figures in (1) latitude, (2) longitude
c                  and (3) time
c     iprec()   - precision of (1) latitude, (2) longitude, (3) time
c     ip2()     - variable codes for variables in profile
c     ierror()  - whole profile error codes for each variable
c     
c     ipi()     - primary investigators information
c                   1. primary investigators
c                   2. for which variable
c
c     jsig2()   - number of significant figures in each secondary header variable
c     jprec2()  - precision of each secondary header variable
c     sechead() - secondary header variables
c
c     jsigb()   - number of significant figures in each biological variable
c     jprecb()  - precision of each biological variable
c     bio()     - biological data
c
c     depth()   - depth of each measurement
c
c     jtot2()   - number of bytes in each secondary header variable
c     jtotb()   - number of bytes in each biological variable
c
c     msig()    - number of significant figures in each measured variable at
c                  each level of measurement
c     mprec()   - precision of each measured variable at each
c                  level of measurement
c
c     mtot()    - number of digits in each measured variable at
c                  each level of measurement
c
c     temp()    - variable data at each level
c     iderror() - error flags for each variable at each depth level
c     iorigflag()- originators flags for each variable and depth
c
c     isec()    - variable codes for secondary header data
c     ibio()    - variable codes for biological data
c     parminf()  - parameter specific information
c     jprecp()   - precision for parameter specific information
c     jsigp()    - number of significant figures for parameter specific
c                information
c     jtotp()    - number of digits in for parameter specific information
c
c     itaxnum() - different taxonomic and biomass variable
c                  codes found in data
c     vtax()    - value of taxonomic variables and biomass variables
c
c     jsigtax() - number of significant figures in taxon values and
c                  biomass variables
c     jprectax()- precision of taxon values and biomass variables
c
c     jtottax() - number of bytes in taxon values 
c     itaxerr() - error codes for taxon data
c     itaxorigerr() - originators error codes for taxon data
c
c     nbothtot()- total number of taxa variables
c     stdz(40) - standard depth levels
c
c*******************************************************************

      integer isig(3), iprec(3), ip2(0:maxlevel), ierror(maxlevel),
     &        ipi(maxlevel,2)
      dimension jsig2(maxlevel),jprec2(maxlevel),sechead(maxlevel)
      dimension jsigb(maxlevel),jprecb(maxlevel),bio(maxlevel)
      dimension depth(maxlevel)
      dimension jtot2(maxlevel),jtotb(maxlevel)
      dimension msig(maxlevel,maxcalc), mprec(maxlevel,maxcalc)
      dimension mtot(maxlevel,maxcalc)
      dimension temp(maxlevel,maxcalc),iderror(maxlevel,0:maxcalc)
      dimension isec(maxlevel),ibio(maxlevel)
      dimension parminf(maxpinf,0:maxcalc),jsigp(maxpinf,0:maxcalc)
      dimension jprecp(maxpinf,0:maxcalc),jtotp(maxpinf,0:maxcalc)
      dimension iorigflag(maxlevel,0:maxcalc)
      dimension itaxnum(maxtcode,maxtax),vtax(0:maxtcode,maxtax)
      dimension jsigtax(maxtcode,maxtax),jprectax(maxtcode,maxtax)
      dimension jtottax(maxtcode,maxtax),itaxerr(maxtcode,maxtax)
      dimension itaxorigerr(maxtcode,maxtax)
ccao
c     dimension stdz(kdim)
      dimension stdz(kdim), ff(maxlevel), dd(kmb), tt(kmb)

      common /thedata/ depth,temp
      common /flags/ ierror,iderror
      common /oflags/ iorigflag
      common /significant/ msig
      common /precision/ mprec
      common /totfigs/ mtot
      common /second/ jsig2,jprec2,jtot2,isec,sechead
      common /parminfo/ jsigp,jprecp,jtotp,parminf
      common /biology/ jsigb,jprecb,jtotb,ibio,bio
      common /taxon/ jsigtax,jprectax,jtottax,itaxerr,
     &     vtax,itaxnum,nbothtot,itaxorigerr
      
      data stdz/ 0., 10., 20., 30., 50., 75., 100., 125., 150.,
     &     200., 250., 300., 400., 500., 600., 700., 800., 900.,
     &     1000., 1100., 1200., 1300., 1400., 1500., 1750., 2000.,
     &     2500., 3000., 3500., 4000., 4500., 5000., 5500., 6000.,
     &     6500., 7000., 7500., 8000., 8500., 9000./
      
c**************************************************************
c
c     nf is the input file indentification number
c
c**************************************************************

      data nf/11/
ccao
      data ntf/20/
      data nsf/30/
cgena file numbers for  tests
c      data ntpr/21/
c      data nspr/31/
c      
      data typ/'XBT'/
c
c**************************************************************
c Levels for SODA-4 MOM5.1
c**************************************************************
       data dd/
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
c
c
c**************************************************************
c Levels for SODA-MOM5.0.2
c**************************************************************
c       data dd/5.03355, 15.10065, 25.21935, 35.35845, 45.57635,
c     &        55.86325, 66.26175, 76.80285, 87.57695, 98.62325,
c     &        110.0962, 122.1067, 134.9086, 148.7466, 164.0538,
c     &        181.3125, 201.2630, 224.7773, 253.0681, 287.5508,
c     &        330.0078, 382.3651, 446.7263, 524.9824, 618.7031,
c     &        728.6921, 854.9935, 996.7153, 1152.376, 1319.997,
c     &        1497.562, 1683.057, 1874.788, 2071.252/
c
c**************************************************************
c levels for SODA-POP
c**************************************************************
c      data dd/5., 15., 25., 35., 46., 57., 70., 82., 96.,
c     *        112., 129., 148., 171., 197., 229., 268., 317.,
c     *        381., 465., 579., 729., 918., 1139./
     
c    *        381., 465., 579., 729., 918., 1139., 1378., 1625.,
c    *        1875., 2125., 2375., 2624., 2874., 3124., 3374., 
c    *        3624., 3874., 4124., 4374., 4624., 4874., 5124., 
c    *        5374./
c**************************************************************
c NCEP's levels (D. Behringer 11.04.2007)
c**************************************************************
c      data dd/   5.,   15.,   25.,   35.,   45.,   55.,   65., 
c     &          75.,   86.,   95.,  105.,  115.,  125.,  135.,
c     &         145.,  155.,  165.,  175.,  185.,  195.,  205.,
c     &         215.,  225.,  
c     &     238.4779,  262.2945,  303.0287,  366.7978,  459.0910,
c     &     584.6193,  747.1870,  949.5881, 1193.5300, 1479.5880,
c     &    1807.1870, 2174.6190, 2579.0910, 3016.7980, 3483.0290,
c     &    3972.2940, 4478.4780/
      
c**************************************************************
c
c     Get user input file name from which profiles will be
c     taken.  Open this file.
c
c**************************************************************

c User can modify the next section to read from a text file listing
c different input data files as a do-loop, for example, as opposed
c a single data input file.

ccao
c
      do 100 iy=iyst, iynd
c      
c  itc1 - readed temperature profiles number
c  itc2 - writen temperature profiles number
c  isc1 - readed salinity profiles number
c  isc2 - writen salinity profiles number
        itc1 = 0
	itc2 = 0
	isc1 = 0
	isc2 = 0 
c      
        write(yr, '(i4)') iy
        call system('ls 
c     &	/aosc/okhotsk/chepurin/DATA4SODA/OBS_DATA/*'//yr//'* > tmp')
     
     &	/aosc/iceland2/chepurin/okhotsk/DATA4SODA/TMP/*'//yr//'*.gz 
     &  > tmp')    
c     &	/aosc/okhotsk/chepurin/DATA4SODA/TMP/*'//yr//'*.gz > tmp')    
c     &	/aosc/laptev/chepurin/DATA4SODA/OBS_DATA/*'//yr//'* > tmp')
c     &  /aosc/beaufort/arctic1/chepurin/DATA4SODA/OBS_DATA/*'//yr//'*
c     &  > tmp')
c     &  /data/arctic1/chepurin/DATA4SODA/OBS_DATA/*'//yr//'* > tmp')
c     &  /data/pacific3/chepurin/STD_DATA/ORIGINAL/*'//yr//'* > tmp')
     
        open (50,file='tmp', status='old')
     
        toutfile='TMP/t_'//yr//'_55levs.dat'
        soutfile='TMP/s_'//yr//'_55levs.dat'
        open(ntf, file=toutfile, status='unknown')
        open(nsf, file=soutfile, status='unknown')
cgena spr files for tests
c        tprfile='TMP/t_'//yr//'.tst'
c	sprfile='TMP/s_'//yr//'.tst'
c	open(ntpr,file=tprfile, status='unknown')
c	open(nspr,file=sprfile, status='unknown')
	
        do 101 ityp=1,11
	  ieof = 0
	  
	  read(50,'(a80)', end=102) dum
	  write(6,*) dum
	  read(dum,'(46x,a3)') typ	  
	  write(6,*) 'data type - '//typ
	  
          read(dum,'(a53)') filename
	  write(6,*) 'File => '//filename

          call system('gzip -d '//dum)
          open(nf,file=filename, status='old')
      
c**************************************************************
c
c   SUBROUTINE "OCLread":  READS IN A SINGLE PROFILE FROM THE ASCII 
c                          FILE AND STORES THE DATA INTO ARRAYS (PASSED
c                          OR SHARED BETWEEN OCLread AND OCLFOR).
c   -------------------------------------------------------------------
c
c   Passed Variables:
c     
c     nf      - file identification number for input file
c     jj      - OCL profile number
c     cc      - NODC country code
c     icruise - NODC cruise number
c     iyear   - year of profile
c     month   - month of profile
c     iday    - day of profile
c     time    - time of profile
c     rlat    - latitude of profile
c     rlon    - longitude of profile
c     levels  - number of depth levels of data
c     istdlev - observed (0) or standard (1) levels
c     nvar   - number of variables recorded in profile
c     ip2(i)  - variable codes of variables in profile
c     nsecond - number of secondary header variables
c     nbio    - number of biological variables
c     isig()  - number of significant figures in (1) latitude, (2) longitude
c                and (3) time
c     iprec() - precision of (1) latitude, (2) longitude, (3) time
c     bmiss   - missing value marker
c     ieof    - set to one if end of file has been encountered
c
c   Common/Shared Variables and Arrays (see COMMON area of program):
c
c     depth(x)   - depth in meters (x = depth level)
c     temp(x,y)  - variable data (x = depth level, y = variable ID = ip2(i))
c                ... see also nvar, ip2, istdlev, levels above ...
c     sechead(i) - secondary header data (i = secondary header ID = isec(j))
c     isec(j)   - secondary header ID (j = #sequence in profile (1st, 2nd, 3rd))
c                ... see also nsecond above ...
c     bio(i)     - biology header data (i = biol-header ID = ibio(j))
c     ibio(j)    - biology header ID (j = #sequence in profile (1st, 2nd, 3rd))
c                ... see also nbio above ...
c     nbothtot   - number of taxa set / biomass variables
c     vtax(i,j)  - taxonomic/biomass array, where j = (1..nbothtot)
c                   For each entry (j=1..nbothtot), there are vtax(0,j)
c                   sub-entries.  [Note:  The number of sub-entries is 
c                   variable for each main entry.]  vtax also holds the
c                   value of the sub-entries.
c     itaxnum(i,j)- taxonomic code or sub-code 
c     chars       - OCL character data: 1. originators cruise code,
c                                       2. originators station cod
c     npi        - number of PI codes
c     ipi        - Primary Investigator information
c                  1. primary investigator
c                  2. variable investigated
c     

c***************************************************************
      
      iVERSflag = 0             !- default is "WOD-2001"
      
c Modify this loop to increase or decrease the number of stations read
c from any particular WMO from 10 to any number (e.g., 500000)

ccao
c      do 50 ij=1,10          !- MAIN STATION LOOP 
       do 50 ij=1,1000000
         
c      chars(1)= '               '
c      chars(2)= '               '
       
       if (iVERSflag .ne. 1) then
          
c     .   Read in as "WOD-2001" format

          call WODREAD200X(nf,jj,cc,icruise,iyear,month,iday,
     &         time,rlat,rlon,levels,istdlev,nvar,ip2,nsecond,nbio,
     &         isig,iprec,bmiss,ieof,chars,ipi,npi,iVERSflag)
	  
c     .   ONLY happens if format rejected (rewind and try as WOD-1999)
          if (iVERSflag .eq. 1) then
             print *,' '
             print *,
     &            ' This data file is not in the WOD-2001 format.',
     &            '    Trying it in the WOD-1998 format.'
ccao
             stop 'wod98'
             print *,' '
             print *,' '
             rewind(nf) !- rewind file
          endif

       endif

       if (iVERSflag .eq. 1) then 
          
c     .   Read in as "WOD-1998" format          
c     .  OCLread2001 rejected format.  Must be WOD-1998 format.
          
          call OCLread1998(nf,jj,cc,icruise,iyear,month,iday,
     &         time,rlat,rlon,levels,istdlev,nvar,ip2,nsecond,nbio,
     &         isig,iprec,bmiss,ieof,chars,ipi,npi)
       endif
       
       if ( ieof.gt.0 ) goto 4  !- Exit
       
c***************************************************************
c
ccaostd
ccao want only standard level data
c
       if (istdlev .eq. 0) then
c         print *, 'observed level data?'
c         stop 'obd data'
       endif
c
c**************************************************************
c     
c     WRITE HEADER INFORMATION TO THE SCREEN
c
c     cc         - country code (a2)
c     icruise    - OCL internal cruise identifier (i8)
c     rlat       - latitude (f7.3)
c     rlon       - longitude (f7.8)
c     iyear      - year (i4)
c     month      - month (i2)
c     iday       - day (i2)
c     time       - time (GMT) (f5.2)
c     jj         - OCL internal profile identifier (i8)
c     levels     - number of depth levels measured (i4)
c     
800   format(1x,a2,i8,1x,f7.3,1x,f8.3,1x,i4,1x,i2,1x,i2,
     &      1x,f7.2,1x,i8,1x,i6)

c    &'CC  cruise Latitde Longitde YYYY MM DD    Time  Station #levels'

c      write(6,800) cc,icruise,rlat,rlon,iyear,month,iday,
c     &      time,jj,levels
cgena      write(6,807) ij,rlat,rlon,iyear,month,iday,levels
807   format(1x,i7,1x,f7.3,1x,f8.3,1x,i4,1x,i2,1x,i2, 1x,i6)     
c      
c      write(6,*) 'Number of variables in this profile: ',nvar
c      
c**************************************************************
c
c     WRITE VARIABLE-CODE (column headings) TO THE SCREEN
c     
c     nvar          - number of variables (1...nvar)
c     ip2(1..nvar)  - variable code for each variable present
c     
c     Example:  
c       For a profile with just Temperature[1], Oxygen[3], Pressure[25]:
c     
c       The variable sequence is ip2(1)=Temperature, ip2(2)=Oxygen, 
c          ip2(3)=Pressure
c     
c       nvar = 3
c
c       ip2(1) = 1, ip2(2) = 3, ip3(3) = 25
c     
c     Note:  If "nvar = 0", biology only station.
c     format(5x,1a,5x,10(3x,i2,11x))
c 801   format(5x,a5,3x,10(i2,8x,a2,3x))
c801   format(5x,a5,4x,10(i2,8x,a2,3x))
c      
c      if (nvar .gt. 0) then
c         
c         write(6,801) "z  fo",((ip2(n),"fo"),n=1,nvar)
c         write(6,801) "z  fo",(ip2(n),n=1,nvar)
c801    format(5x,a5,4x,10(i2,3x))

c**************************************************************
c
c     WRITE DEPTH-DEPENDENT VARIABLE DATA TO THE SCREEN
c     Print depth (depth(n)), error flags for depth (iderror(n,0)),
c     each variable (temp(n,1..nvar)), and error flags for each
c     variables (iderror(n,1..nvar))
c
c802      format(1x,f6.1,1x,i1,i1,14(f9.3,' (',i1,') ',i1,i1))
c         do 80 n=1,levels
c            write(6,802) depth(n),iderror(n,0),iorigflag(n,0),
c    &            (temp(n,ip2(j)),msig(n,ip2(j)),
c    &            iderror(n,ip2(j)),iorigflag(n,ip2(j)),j=1,nvar)
c80       continue
          
c***************************************************************
c
c     PRINT ENTIRE-PROFILE ERROR FLAGS
c8021     format('VarFlag:    ',11x,11(i1,14x))
c         write(6,8021)(ierror(ip2(j)),j=1,nvar)

c***************************************************************
c
ccao - write out data
c
        do ll=1,nvar
c
          nout=0
c
c -- temp or salt --
c
          if (ip2(ll).eq.1) then
	    nout = ntf
	    itc1 = itc1+1
	  endif
          if (ip2(ll).eq.2) then
	    nout = nsf
	    isc1 = isc1+1
	  endif
cgena  file numbers for tests
c          npr=0
c          if(ip2(ll).eq.1) npr=ntpr
c          if(ip2(ll).eq.2) npr=nspr
	  
c
          if(nout.ne.0) then
c
c -- quality control --
c
cgena          if (ierror(ip2(ll)).ne.0) goto 50
          if (ierror(ip2(ll)).ne.0) goto 351
c
c -- vertical interpolation --
c
          do kk=1, levels
            ff(kk) = temp(kk,ip2(ll))
          enddo
c	  
          do kk=1,kmb
            tt(kk)=bmiss
          enddo
c
          kmm=kmb
c--- interpolation from standart levels
          if (istdlev.eq.1) then
            if (levels.le.2) goto 50
            if(stdz(levels).lt.dd(kmb)) then
              do kk=kmb, 2, -1
                if(stdz(levels).lt.dd(kk).and.
     &             stdz(levels).ge.dd(kk-1)) then
                  kmm=kk-1
                  go to 345
                endif
              enddo
            endif
 345        continue
            call sintrp(ff, levels, tt, kmm, stdz, dd, bmiss)
	  endif
c
c--- interpolation from observed levels
          if (istdlev.eq.0) then
            iii=0
	    do kk=1, kmb
	      if(depth(1).le.dd(kk).and.
     &           depth(levels).ge.dd(kk)) then
                iii=1
		goto 347
	      endif
	    enddo
347         continue
            if(iii.eq.0) goto 50    	  
            if(depth(levels).lt.dd(kmb)) then
              do kk=kmb, 2, -1
                if(depth(levels).lt.dd(kk).and.
     &		   depth(levels).ge.dd(kk-1)) then
                  kmm=kk-1
                  go to 346
                endif
              enddo
            endif
 346        continue
            call sintrp(ff, levels, tt, kmm, depth, dd, bmiss)
	  endif
c
c--- write data on the model levels
          do n=1, kmm
	    if(tt(n).ne.bmiss) goto 350
	  enddo
	  goto 351 
350       continue	  
	  if(ip2(ll).eq.1) itc2 = itc2+1
	  if(ip2(ll).eq.2) isc2 = isc2+1
          write(nout,111) typ, jj, iyear, month, iday, rlat, rlon,
     *             kmm, (tt(n), n=1,kmm)
351       continue          
          endif
	  
        enddo
c
cgena write data on original levels for tests
c          write(npr,112) typ, jj, iyear, month, iday, rlat, rlon,
c     &             levels, (depth(n), n=1,levels)
c          write(npr,113) (ff(n),n=1,levels)
c112    format(a4,i8,i5,2i3,2f9.3,i3,1000f8.3)
c113    format(44x,1000f8.3)        
c
c 111   format(a4,i8,i5,2i3,2f9.3,i3,23f8.3)
c 111   format(a4,i8,i5,2i3,2f9.3,i3,40f8.3)
 111   format(a4,i8,i5,2i3,2f9.3,i3,55f8.3)
c
50    continue !- End of MAIN LOOP
      
 4    continue !- EXIT 
c
      close (nf)
      call system('gzip -9 '//filename)
101   continue
102   continue   
      close (ntf)
      close (nsf)
cgena close  files for test
c      close (ntps)
c      close (nspr)
c
      close (50)
      call system('rm -f tmp')
c
      write(6,*) 'Year - ', iy
      write(6,*)  'T-profiles:  read = ',itc1,',  write = ',itc2
      write(6,*)  'S-profiles:  read = ',isc1,',  write = ',isc2
c
100   continue !- End of All Years

9999  continue

      stop
      end
c=================================================================================
c----------------------------------------------------------------

      SUBROUTINE WODREAD200X(nf,jj,cc,icruise,iyear,month,iday,
     &     time,rlat,rlon,levels,isoor,nvar,ip2,nsecond,nbio,
     &     isig,iprec,bmiss,ieof,chars,ipi,npi,iVERSflag)
      
c     This subroutine reads in the WOD ASCII format and loads it
c     into arrays which are common/shared with the calling program.

c*****************************************************************
c
c   Passed Variables:
c
c     nf       - file identification number for input file
c     jj       - WOD cast number
c     cc       - NODC country code
c     icruise  - NODC cruise number
c     iyear    - year of cast
c     month    - month of cast
c     iday     - day of cast
c     time     - time of cast
c     rlat     - latitude of cast
c     rlon     - longitude of cast
c     levels   - number of depth levels of data
c     isoor    - observed (0) or standard (1) levels
c     nvar     - number of variables recorded in cast
c     ip2      - variable codes of variables in cast
c     nsecond  - number of secondary header variables
c     nbio     - number of biological variables
c     isig     - number of significant figures in (1) latitude, (2) longitude,
c                 and (3) time
c     iprec    - precision of (1) latitude, (2) longitude, (3) time
c     itotfig  - number of digits in (1) latitude, (2) longitude, (3) time
c     bmiss    - missing value marker
c     ieof     - set to one if end of file has been encountered
c     chars    - character data: 1=originators cruise code,
c                                2=originators station code
c     npi      - number of PI codes
c     ipi      - Primary Investigator information
c                  1. primary investigator
c                  2. variable investigated
c
c     iVERSflag  -  set to "1" if data are in WOD-1998 format. 
c                (subroutine exits so 1998 subroutine can be run)
c
c   Common/Shared Variables and Arrays (see COMMON area of program):
c
c     depth(x)   - depth in meters (x = depth level)
c     temp(x,y)  - variable data (x = depth level, y = variable ID = ip2(i))
c                ... see also nvar, ip2, istdlev, levels above ...
c     sechead(i) - secondary header data (i = secondary header ID = isec(j))
c     isec(j)    - secondary header ID (j = #sequence (1st, 2nd, 3rd))
c                ... see also nsecond above ...
c     bio(i)     - biology header data (i = biol-header ID = ibio(j))
c     ibio(j)    - biology header ID (j = #sequence (1st, 2nd, 3rd))
c                ... see also nbio above ...
c     nbothtot   - number of taxa set / biomass variables
c     vtax(i,j)  - taxonomic/biomass array, where j = (1..nbothtot)
c                   For each entry (j=1..nbothtot), there are vtax(0,j)
c                   sub-entries.  [Note:  The number of sub-entries is
c                   variable for each main entry.]  vtax also holds the
c                   value of the sub-entries.
c    itaxnum(i,j)- taxonomic code or sub-code
c    parminf(i,j)- variable specific information
c    origflag(i,j)- originators data flags
c
c***************************************************************


c******************************************************************
c
c   Parameters (constants):
c
c     maxlevel - maximum number of depth levels, also maximum
c                 number of all types of variables
c     maxcalc  - maximum number of measured and calculated
c                 depth dependent variables
c     maxtcode - maximum number of different taxa variable codes
c     maxtax   - maximum number of taxa sets
c     maxpinf - number of distinct variable specific information
c               variables
c
c******************************************************************

      parameter (maxlevel=30000, maxcalc=100)
      parameter (maxtcode=25, maxtax=2000, maxpinf=25)

c******************************************************************
c
c   Character Variables:
c
c     cc       - NODC country code
c     xchar    - dummy character array for reading in each 80
c                 character record
c     aout     - format specifier (used for FORTRAN I/O)
c     ichar    - cast character array
c     
c******************************************************************

      character*2  cc
      character*4  aout
      character*15 chars(2)
      character*80 xchar
      character*1500000 ichar
      
      data aout /'(iX)'/
      
c******************************************************************
c
c    Arrays:
c
c     isig     - number of significant figures in (1) latitude, (2) longitude,
c                 and (3) time
c     iprec    - precision of (1) latitude, (2) longitude, (3) time
c     itotfig  - number of digits in (1) latitude, (2) longitude, (3) time
c     ip2      - variable codes for variables in cast
c     ierror   - whole profile error codes for each variable
c     jsig2    - number of significant figures in each secondary header variable
c     jprec2   - precision of each secondary header variable
c     jtot2    - number of digits in each secondary header variable
c     sechead  - secondary header variables
c     jsigb    - number of significant figures in each biological variable
c     jprecb   - precision of each biological variable
c     jtotb    - number of digits in each biological variable
c     bio      - biological data
c     idsig    - number of significant figures in each depth measurement
c     idprec   - precision of each depth measurement
c     idtot    - number of digits in each depth measurement
c     depth    - depth of each measurement
c     msig     - number of significant figures in each measured variable at
c                 each level of measurement
c     mprec    - precision of each measured variable at each
c                 level of measurement
c     mtot     - number of digits in each measured variable at
c                 each level of measurement
c     temp     - variable data at each level
c     iderror  - error flags for each variable at each depth level
c     iorigflag- originators flags for each variable and depth
c     isec     - variable codes for secondary header data
c     ibio     - variable codes for biological data
c     parminf  - variable specific information
c     jprecp   - precision for variable specific information
c     jsigp    - number of significant figures for variable specific
c                information
c     jtotp    - number of digits in for variable specific information
c     itaxnum  - different taxonomic and biomass variable
c                 codes found in data
c     vtax     - value of taxonomic variables and biomass variables
c     jsigtax  - number of significant figures in taxon values and
c                 biomass variables
c     jprectax - precision of taxon values and biomass variables
c     jtottax  - number of digits in taxon values and biomass
c                 variables
c     itaxerr  - taxon variable error code
c     itaxorigerr - taxon originators variable error code
c     nbothtot - total number of taxa and biomass variables
c     ipi      - Primary investigator informationc
c                 1. primary investigator
c                 2. variable investigated
c
c*******************************************************************

      dimension isig(3), iprec(3), ip2(0:maxlevel), ierror(maxlevel)
      dimension itotfig(3),ipi(maxlevel,2)
      dimension jsig2(maxlevel), jprec2(maxlevel), sechead(maxlevel)
      dimension jsigb(maxlevel), jprecb(maxlevel), bio(maxlevel)
      dimension idsig(maxlevel),idprec(maxlevel), depth(maxlevel)
      dimension jtot2(maxlevel),jtotb(maxlevel),idtot(maxlevel)
      dimension msig(maxlevel,maxcalc), mprec(maxlevel,maxcalc)
      dimension mtot(maxlevel,maxcalc)
      dimension temp(maxlevel,maxcalc),iderror(maxlevel,0:maxcalc)
      dimension isec(maxlevel),ibio(maxlevel)
      dimension parminf(maxpinf,0:maxcalc),jsigp(maxpinf,0:maxcalc)
      dimension jprecp(maxpinf,0:maxcalc),jtotp(maxpinf,0:maxcalc)
      dimension iorigflag(maxlevel,0:maxcalc)
      dimension itaxnum(maxtcode,maxtax),vtax(0:maxtcode,maxtax)
      dimension jsigtax(maxtcode,maxtax),jprectax(maxtcode,maxtax)
      dimension jtottax(maxtcode,maxtax),itaxerr(maxtcode,maxtax)
      dimension itaxorigerr(maxtcode,maxtax)

c*******************************************************************
c     
c   Common Arrays and Variables:
c
c*******************************************************************
      
      common /thedata/ depth,temp
      common /flags/ ierror,iderror
      common /oflags/ iorigflag
      common /significant/ msig
      common /precision/ mprec
      common /totfigs/ mtot
      common /second/ jsig2,jprec2,jtot2,isec,sechead
      common /parminfo/ jsigp,jprecp,jtotp,parminf
      common /biology/ jsigb,jprecb,jtotb,ibio,bio
      common /taxon/ jsigtax,jprectax,jtottax,itaxerr,
     &     vtax,itaxnum,nbothtot,itaxorigerr
      
c******************************************************************
c     
c     Read in the first line of a cast into dummy character
c     variable xchar
c     
c
c     WOD-2005   First byte of each "cast record" is char "A".
c
c     WOD-1998   First byte of each "cast recond" is a number.
c
c******************************************************************

      read(nf,'(a80)',end=500) xchar
      
      if ( xchar(1:1) .ne. 'B' .and. xchar(1:1) .ne. 'A' .and.
     *     xchar(1:1) .ne. 'C' ) then

         iVERSflag = 1 !- not WOD-2005 format, must be WOD-1998
         return

      else
         if ( xchar(1:1) .eq. 'C' ) then
          iVERSflag=2   !- WOD-2013 format
         else
          iVERSflag = 0 !- WOD-2005 format
         endif
      endif
      
c******************************************************************
c
c     The first seven characters of a cast contain the
c     number of characters which make up the entire cast.  Read
c     this number into nchar
c     
c******************************************************************

      read(xchar(2:2),'(i1)') inc
      write(aout(3:3),'(i1)') inc
      read(xchar(3:inc+2),aout) nchar

c******************************************************************
c
c     Place the first line of the cast into the cast holder
c     character array (ichar)
c
c******************************************************************

      ichar(1:80) = xchar

c******************************************************************
c
c     Calculate the number of full (all 80 characters contain information)
c     lines in this cast.  Subtract one since the first line was
c     already read in.
c
c******************************************************************

      nlines = nchar/80

c*****************************************************************
c
c     Read each line into the dummy variable
c
c*****************************************************************

      do 49 n0 = 2,nlines

       read(nf,'(a80)') xchar

c*****************************************************************
c
c     Place the line into the whole cast array
c
c*****************************************************************

       n = 80*(n0-1)+1
       ichar(n:n+79)=xchar

49    continue

c*****************************************************************
c
c     If there is a last line with partial information, read in
c     this last line and place it into the whole cast array
c
c*****************************************************************

      if ( nlines*80 .lt. nchar .and. nlines .gt. 0) then

       read(nf,'(a80)') xchar

       n = 80*nlines+1
       ichar(n:nchar) = xchar

      endif
       
c*****************************************************************
c
c   Extract header information from the cast array
c
c     jj       - WOD cast number  
c     cc       - NODC country code  
c     icruise  - NODC cruise number
c     iyear    - year of cast
c     month    - month of cast
c     iday     - day of cast
c
c*****************************************************************

      istartc=inc+3
      read(ichar(istartc:istartc),'(i1)') inc
      write(aout(3:3),'(i1)') inc
      read(ichar(istartc+1:istartc+inc),aout) jj
      istartc=istartc+inc+1

      cc = ichar(istartc:istartc+1)
      istartc=istartc+2

      read(ichar(istartc:istartc),'(i1)') inc
      write(aout(3:3),'(i1)') inc
      read(ichar(istartc+1:istartc+inc),aout) icruise
      istartc=istartc+inc+1

      read(ichar(istartc:istartc+3),'(i4)') iyear
      istartc=istartc+4
      read(ichar(istartc:istartc+1),'(i2)') month
      istartc=istartc+2
      read(ichar(istartc:istartc+1),'(i2)') iday
      istartc=istartc+2

c*****************************************************************
c
c   SUBROUTINE "charout":  READS IN AN WOD ASCII FLOATING-POINT
c                          VALUE SEQUENCE (i.e. # sig-figs,
c                          # total figs, precision, value itself).
c                          * THIS WILL BE CALLED TO EXTRACT MOST 
c   Examples:              FLOATING POINT VALUES IN THE WOD ASCII.
c
c     VALUE  Precision    WOD ASCII
c     -----  ---------    ---------
c     5.35       2        332535
c     5.         0        1105
c     15.357     3        55315357
c    (missing)            -
c
c   ---------------------------------------------------------------
c
c  Read in time of cast (time) using CHAROUT subroutine:
c
c     istartc  - position in character array to begin to read
c                 in data
c     isig     - number of digits in data value
c     iprec    - precision of data value
c     ichar    - character array from which to read data
c     time     - data value
c     bmiss    - missing value marker
c
c*****************************************************************

      call charout(istartc,isig(3),iprec(3),itotfig(3),ichar,
     *             time,bmiss)

c*****************************************************************
c
c     Read in latitude (rlat) and longitude (rlon) using CHAROUT:
c     
c        Negative latitude is south.
c        Negative longitude is west.
c     
c*****************************************************************

      call charout(istartc,isig(1),iprec(1),itotfig(1),ichar,
     *             rlat,bmiss)
      call charout(istartc,isig(2),iprec(2),itotfig(2),ichar,
     *             rlon,bmiss)

c*****************************************************************
c     
c     Read in the number of depth levels (levels) using CHAROUT:
c
c*****************************************************************

      read(ichar(istartc:istartc),'(i1)') inc
      write(aout(3:3),'(i1)') inc
      read(ichar(istartc+1:istartc+inc),aout) levels
      istartc=istartc+inc+1

c*****************************************************************
c
c     Read in whether data is on observed levels (isoor=0) or
c     standard levels (isoor=1)
c
c*****************************************************************

      read(ichar(istartc:istartc),'(i1)') isoor
      istartc=istartc+1

c*****************************************************************
c
c     Read in number of variables in cast
c
c*****************************************************************

      read(ichar(istartc:istartc+1),'(i2)') nvar
      istartc=istartc+2

c*****************************************************************
c
c     Read in the variable codes (ip2()), the whole profile
c       error flags (ierror(ip2())), and variable specific
c       information (iorigflag(,ip2()))
c
c*****************************************************************

      do 30 n = 1,nvar

       read(ichar(istartc:istartc),'(i1)') inc
       write(aout(3:3),'(i1)') inc
       read(ichar(istartc+1:istartc+inc),aout) ip2(n)
       istartc=istartc+inc+1

       read(ichar(istartc:istartc),'(i1)') ierror(ip2(n))
       istartc=istartc+1

       read(ichar(istartc:istartc),'(i1)') inc
       write(aout(3:3),'(i1)') inc
       read(ichar(istartc+1:istartc+inc),aout) npinf
       istartc=istartc+inc+1

       do 305 n2=1,npinf

        read(ichar(istartc:istartc),'(i1)') inc
        write(aout(3:3),'(i1)') inc
        read(ichar(istartc+1:istartc+inc),aout) nn
        istartc=istartc+inc+1

        call charout(istartc,jsigp(nn,ip2(n)),jprecp(nn,ip2(n)),
     &  jtotp(nn,ip2(n)),ichar, parminf(nn,ip2(n)),bmiss) 

305    continue

30    continue

c****************************************************************
c
c     Read in number of bytes in character data
c
c****************************************************************

      read(ichar(istartc:istartc),'(i1)') inc
      istartc=istartc+1

      npi=0
      chars(1)(1:4)='NONE'
      chars(2)(1:4)='NONE'

      if ( inc .gt. 0 ) then
       write(aout(3:3),'(i1)') inc
       read(ichar(istartc+1:istartc+inc),aout) inchad
       istartc=istartc+inc

c****************************************************************
c
c    Read in number of character and primary investigator arrays
c
c****************************************************************

      read(ichar(istartc:istartc),'(i1)') ica
      istartc=istartc+1

c****************************************************************
c
c    Read in character and primary investigator data
c      1 - originators cruise code
c      2 - originators station code
c      3 - primary investigators information
c
c****************************************************************

      do 45 nn=1,ica

       read(ichar(istartc:istartc),'(i1)') icn
       istartc=istartc+1

       if ( icn .lt. 3 ) then
        read(ichar(istartc:istartc+1),'(i2)') ns
        istartc=istartc+2
        chars(icn)= '               '
        chars(icn)= ichar(istartc:istartc+ns-1)
        istartc= istartc+ns
       else
        read(ichar(istartc:istartc+1),'(i2)') npi
        istartc=istartc+2
        do 505 n=1,npi
         read(ichar(istartc:istartc),'(i1)') inc
         write(aout(3:3),'(i1)') inc
         read(ichar(istartc+1:istartc+inc),aout) ipi(n,2)
         istartc=istartc+inc+1

         read(ichar(istartc:istartc),'(i1)') inc
         write(aout(3:3),'(i1)') inc
         read(ichar(istartc+1:istartc+inc),aout) ipi(n,1)
         istartc=istartc+inc+1
505     continue
       endif

45    continue

      endif

c****************************************************************
c
c     Read in number of bytes in secondary header variables
c
c****************************************************************

      read(ichar(istartc:istartc),'(i1)') inc
      istartc=istartc+1
      if ( inc .gt. 0 ) then
       write(aout(3:3),'(i1)') inc
       read(ichar(istartc+1:istartc+inc),aout) insec
       istartc=istartc+inc

c****************************************************************
c
c     Read in number of secondary header variables (nsecond)
c
c****************************************************************

       read(ichar(istartc:istartc),'(i1)') inc
       write(aout(3:3),'(i1)') inc
       read(ichar(istartc+1:istartc+inc),aout) nsecond
       istartc=istartc+inc+1

c****************************************************************
c
c     Read in secondary header variables (sechead())
c
c****************************************************************

       do 35 n = 1,nsecond

        read(ichar(istartc:istartc),'(i1)') inc
        write(aout(3:3),'(i1)') inc
        read(ichar(istartc+1:istartc+inc),aout) nn
        istartc=istartc+inc+1

        call charout(istartc,jsig2(nn),jprec2(nn),jtot2(nn),ichar,
     &  sechead(nn),bmiss) 

        isec(n) = nn

35     continue

       endif

c****************************************************************
c
c     Read in number of bytes in biology variables 
c
c****************************************************************

      read(ichar(istartc:istartc),'(i1)') inc
      istartc=istartc+1

      nbio=0
      inbio=0
      if ( inc .gt. 0 ) then
       write(aout(3:3),'(i1)') inc
       read(ichar(istartc+1:istartc+inc),aout) inbio
       istartc=istartc+inc

c****************************************************************
c
c     Read in number of biological variables (nbio)
c
c****************************************************************

      read(ichar(istartc:istartc),'(i1)') inc
      write(aout(3:3),'(i1)') inc
      read(ichar(istartc+1:istartc+inc),aout) nbio
      istartc=istartc+inc+1

c****************************************************************
c
c     Read in biological variables (bio())
c
c****************************************************************

      do 40 n = 1,nbio

       read(ichar(istartc:istartc),'(i1)') inc
       write(aout(3:3),'(i1)') inc
       read(ichar(istartc+1:istartc+inc),aout) nn
       istartc=istartc+inc+1

       call charout(istartc,jsigb(nn),jprecb(nn),jtotb(nn),ichar,
     & bio(nn),bmiss)

       ibio(n) = nn

40    continue

c****************************************************************
c
c     Read in biomass and taxonomic variables
c
c****************************************************************

      read(ichar(istartc:istartc),'(i1)') inc
      write(aout(3:3),'(i1)') inc
      read(ichar(istartc+1:istartc+inc),aout) nbothtot
      istartc=istartc+inc+1

      do 41 n = 1,nbothtot

       itaxtot=0
       read(ichar(istartc:istartc),'(i1)') inc
       write(aout(3:3),'(i1)') inc
       read(ichar(istartc+1:istartc+inc),aout) nn
       istartc=istartc+inc+1

       vtax(0,n)=nn

       do 42 n2 =1,nn

        itaxtot=itaxtot+1

        read(ichar(istartc:istartc),'(i1)') inc
        write(aout(3:3),'(i1)') inc
        read(ichar(istartc+1:istartc+inc),aout) itaxnum(itaxtot,n)
        istartc=istartc+inc+1
        call charout(istartc,jsigtax(itaxtot,n),jprectax(itaxtot,n),
     &   jtottax(itaxtot,n),ichar,vtax(itaxtot,n),bmiss)

        read(ichar(istartc:istartc),'(i1)') itaxerr(itaxtot,n)
        istartc=istartc+1
        read(ichar(istartc:istartc),'(i1)') itaxorigerr(itaxtot,n)
        istartc=istartc+1

42     continue

41    continue
      endif

c****************************************************************
c
c     Read in measured and calculated depth dependent variables
c       along with their individual reading flags
c
c****************************************************************

      do 50 n = 1,levels

        if ( isoor.eq.0 .or. iVERSflag .eq. 2 ) then

        call charout(istartc,idsig(n),idprec(n),idtot(n),ichar,
     & depth(n),bmiss)

        read(ichar(istartc:istartc),'(i1)') iderror(n,0)
        istartc=istartc+1
        read(ichar(istartc:istartc),'(i1)') iorigflag(n,0)
        istartc=istartc+1

       endif

       do 55 i = 1,nvar
     
        call charout(istartc,msig(n,ip2(i)),mprec(n,ip2(i)),
     & mtot(n,ip2(i)),ichar,temp(n,ip2(i)),bmiss)

       if ( temp(n,ip2(i)) .gt. bmiss ) then

       read(ichar(istartc:istartc),'(i1)') iderror(n,ip2(i))
       istartc=istartc+1
       read(ichar(istartc:istartc),'(i1)') iorigflag(n,ip2(i))
       istartc=istartc+1

       else
    
        iderror(n,ip2(i))=0
        iorigflag(n,ip2(i))=0
        msig(n,ip2(i))=0
        mprec(n,ip2(i))=0
        mtot(n,ip2(i))=0

       endif

55     continue

50     continue

       return

500    ieof = 1

       return
       end
c=================================================================================

C-----------------------------------------------

      SUBROUTINE OCLREAD1998(nf,jj,cc,icruise,iyear,month,iday,
     &     time,rlat,rlon,levels,isoor,nvar,ip2,nsecond,nbio,
     &     isig,iprec,bmiss,ieof,chars,ipi,npi)
      
c     This subroutine reads in the OCL ASCII format and loads it
c     into arrays which are common/shared with the calling program.

c*****************************************************************
c
c   Passed Variables:
c
c     nf       - file identification number for input file
c     jj       - OCL profile number
c     cc       - NODC country code
c     icruise  - NODC cruise number
c     iyear    - year of profile
c     month    - month of profile
c     iday     - day of profile
c     time     - time of profile
c     rlat     - latitude of profile
c     rlon     - longitude of profile
c     levels   - number of depth levels of data
c     isoor    - observed (0) or standard (1) levels
c     nvar     - number of variables recorded in profile
c     ip2      - variable codes of variables in profile
c     nsecond  - number of secondary header variables
c     nbio     - number of biological variables
c     isig     - number of significant figures in (1) latitude, (2) longitude,
c                 and (3) time
c     iprec    - precision of (1) latitude, (2) longitude, (3) time
c     itotfig  - number of digits in (1) latitude, (2) longitude, (3) time
c     bmiss    - missing value marker
c     ieof     - set to one if end of file has been encountered
c     chars    - character data: 1=originators cruise code,
c                                2=originators station code
c     npi      - number of PI codes
c     ipi      - Primary Investigator information
c                  1. primary investigator
c                  2. variable investigated
c
c   Common/Shared Variables and Arrays (see COMMON area of program):
c
c     depth(x)   - depth in meters (x = depth level)
c     temp(x,y)  - variable data (x = depth level, y = variable ID = ip2(i))
c                ... see also nvar, ip2, istdlev, levels above ...
c     sechead(i) - secondary header data (i = secondary header ID = isec(j))
c     isec(j)    - secondary header ID (j = #sequence in profile (1st, 2nd, 3rd))
c                ... see also nsecond above ...
c     bio(i)     - biology header data (i = biol-header ID = ibio(j))
c     ibio(j)    - biology header ID (j = #sequence in profile (1st, 2nd, 3rd))
c                ... see also nbio above ...
c     nbothtot   - number of taxa set / biomass variables
c     vtax(i,j)  - taxonomic/biomass array, where j = (1..nbothtot)
c                   For each entry (j=1..nbothtot), there are vtax(0,j)
c                   sub-entries.  [Note:  The number of sub-entries is
c                   variable for each main entry.]  vtax also holds the
c                   value of the sub-entries.
c    itaxnum(i,j)- taxonomic code or sub-code
c
c***************************************************************


c******************************************************************
c
c   Parameters (constants):
c
c     maxlevel - maximum number of depth levels, also maximum
c                 number of all types of variables
c     maxcalc  - maximum number of measured and calculated
c                 depth dependent variables
c     maxtcode - maximum number of different taxa variable codes
c     maxtax   - maximum number of taxa sets
c
c******************************************************************

      parameter (maxlevel=10000, maxcalc=100)
c     parameter (maxlevel=6000, maxcalc=200)
      parameter (maxtcode=25, maxtax=2000)

c******************************************************************
c
c   Character Variables:
c
c     cc       - NODC country code
c     xchar    - dummy character array for reading in each 80
c                 character record
c     aout     - format specifier (used for FORTRAN I/O)
c     ichar    - profile character array
c     
c******************************************************************

      character*2  cc
      character*4  aout
      character*15 chars(2)
      character*80 xchar
      character*300000 ichar
      
      data aout /'(iX)'/
      
c******************************************************************
c
c    Arrays:
c
c     isig     - number of significant figures in (1) latitude, (2) longitude,
c                 and (3) time
c     iprec    - precision of (1) latitude, (2) longitude, (3) time
c     itotfig  - number of digits in (1) latitude, (2) longitude, (3) time
c     ip2      - variable codes for variables in profile
c     ierror   - whole profile error codes for each variable
c     jsig2    - number of significant figures in each secondary header variable
c     jprec2   - precision of each secondary header variable
c     jtot2    - number of digits in each secondary header variable
c     sechead  - secondary header variables
c     jsigb    - number of significant figures in each biological variable
c     jprecb   - precision of each biological variable
c     jtotb    - number of digits in each biological variable
c     bio      - biological data
c     idsig    - number of significant figures in each depth measurement
c     idprec   - precision of each depth measurement
c     idtot    - number of digits in each depth measurement
c     depth    - depth of each measurement
c     msig     - number of significant figures in each measured variable at
c                 each level of measurement
c     mprec    - precision of each measured variable at each
c                 level of measurement
c     mtot     - number of digits in each measured variable at
c                 each level of measurement
c     temp     - variable data at each level
c     iderror  - error flags for each variable at each depth level
c     isec     - variable codes for secondary header data
c     ibio     - variable codes for biological data
c     itaxnum  - different taxonomic and biomass variable
c                 codes found in data
c     vtax     - value of taxonomic variables and biomass variables
c     jsigtax  - number of significant figures in taxon values and
c                 biomass variables
c     jprectax - precision of taxon values and biomass variables
c     jtottax  - number of digits in taxon values and biomass
c                 variables
c     itaxerr  - taxon variable error code
c     nbothtot - total number of taxa and biomass variables
c     ipi      - Primary investigator informationc
c                 1. primary investigator
c                 2. variable investigated
c
c*******************************************************************

      dimension isig(3), iprec(3), ip2(0:maxlevel), ierror(maxlevel)
      dimension itotfig(3),ipi(maxlevel,2)
      dimension jsig2(maxlevel), jprec2(maxlevel), sechead(maxlevel)
      dimension jsigb(maxlevel), jprecb(maxlevel), bio(maxlevel)
      dimension idsig(maxlevel),idprec(maxlevel), depth(maxlevel)
      dimension jtot2(maxlevel),jtotb(maxlevel),idtot(maxlevel)
      dimension msig(maxlevel,maxcalc), mprec(maxlevel,maxcalc)
      dimension mtot(maxlevel,maxcalc)
      dimension temp(maxlevel,maxcalc),iderror(maxlevel,0:maxcalc)
      dimension isec(maxlevel),ibio(maxlevel)
      dimension itaxnum(maxtcode,maxtax),vtax(0:maxtcode,maxtax)
      dimension jsigtax(maxtcode,maxtax),jprectax(maxtcode,maxtax)
      dimension jtottax(maxtcode,maxtax),itaxerr(maxtcode,maxtax)

c*******************************************************************
c     
c   Common Arrays and Variables:
c
c*******************************************************************
      
      common /thedata/ depth,temp
      common /flags/ ierror,iderror
      common /significant/ msig
      common /precision/ mprec
      common /totfigs/ mtot
      common /second/ jsig2,jprec2,jtot2,isec,sechead
      common /biology/ jsigb,jprecb,jtotb,ibio,bio
      common /taxon/ jsigtax,jprectax,jtottax,itaxerr,
     &     vtax,itaxnum,nbothtot,itaxorigerr
      
c******************************************************************
c     
c     Read in the first line of a profile into dummy character
c     variable xchar
c     
c******************************************************************

      read(nf,'(a80)',end=500) xchar

c******************************************************************
c
c     The first seven characters of a profile contain the
c     number of characters which make up the entire profile.  Read
c     this number into nchar
c     
c******************************************************************

      read(xchar(1:1),'(i1)') inc
      write(aout(3:3),'(i1)') inc
      read(xchar(2:inc+1),aout) nchar

c******************************************************************
c
c     Place the first line of the profile into the profile holder
c     character array (ichar)
c
c******************************************************************

      ichar(1:80) = xchar

c******************************************************************
c
c     Calculate the number of full (all 80 characters contain information)
c     lines in this profile.  Subtract one since the first line was
c     already read in.
c
c******************************************************************

      nlines = nchar/80

c*****************************************************************
c
c     Read each line into the dummy variable
c
c*****************************************************************

      do 49 n0 = 2,nlines

       read(nf,'(a80)') xchar

c*****************************************************************
c
c     Place the line into the whole profile array
c
c*****************************************************************

       n = 80*(n0-1)+1
       ichar(n:n+79)=xchar

49    continue

c*****************************************************************
c
c     If there is a last line with partial information, read in
c     this last line and place it into the whole profile array
c
c*****************************************************************

      if ( nlines*80 .lt. nchar .and. nlines .gt. 0) then

       read(nf,'(a80)') xchar

       n = 80*nlines+1
       ichar(n:nchar) = xchar

      endif
       
c*****************************************************************
c
c   Extract header information from the profile array
c
c     jj       - OCL profile number  
c     cc       - NODC country code  
c     icruise  - NODC cruise number
c     iyear    - year of profile 
c     month    - month of profile
c     iday     - day of profile 
c
c*****************************************************************

      istartc=inc+2
      read(ichar(istartc:istartc),'(i1)') inc
      write(aout(3:3),'(i1)') inc
      read(ichar(istartc+1:istartc+inc),aout) jj
      istartc=istartc+inc+1

      cc = ichar(istartc:istartc+1)
      istartc=istartc+2

      read(ichar(istartc:istartc),'(i1)') inc
      write(aout(3:3),'(i1)') inc
      read(ichar(istartc+1:istartc+inc),aout) icruise
      istartc=istartc+inc+1

      read(ichar(istartc:istartc+3),'(i4)') iyear
      istartc=istartc+4
      read(ichar(istartc:istartc+1),'(i2)') month
      istartc=istartc+2
      read(ichar(istartc:istartc+1),'(i2)') iday
      istartc=istartc+2

c*****************************************************************
c
c   SUBROUTINE "charout":  READS IN AN OCL ASCII FLOATING-POINT
c                          VALUE SEQUENCE (i.e. # sig-figs,
c                          # total figs, precision, value itself).
c                          * THIS WILL BE CALLED TO EXTRACT MOST 
c   Examples:              FLOATING POINT VALUES IN THE OCL ASCII.
c
c     VALUE  Precision    OCL ASCII
c     -----  ---------    ---------
c     5.35       2        332535
c     5.         0        1105
c     15.357     3        55315357
c    (missing)            -
c
c   ---------------------------------------------------------------
c
c  Read in time of profile (time) using CHAROUT subroutine:
c
c     istartc  - position in character array to begin to read
c                 in data
c     isig     - number of digits in data value
c     iprec    - precision of data value
c     ichar    - character array from which to read data
c     time     - data value
c   -999.99    - missing value marker (bmiss)
c
c*****************************************************************

      call charout(istartc,isig(3),iprec(3),itotfig(3),ichar,
     *             time,bmiss)

c*****************************************************************
c
c     Read in latitude (rlat) and longitude (rlon) using CHAROUT:
c     
c        Negative latitude is south.
c        Negative longitude is west.
c     
c*****************************************************************

      call charout(istartc,isig(1),iprec(1),itotfig(3),ichar,
     *             rlat,bmiss)
      call charout(istartc,isig(2),iprec(2),itotfig(3),ichar,
     *             rlon,bmiss)

c*****************************************************************
c     
c     Read in the number of depth levels (levels) using CHAROUT:
c
c*****************************************************************

      read(ichar(istartc:istartc),'(i1)') inc
      write(aout(3:3),'(i1)') inc
      read(ichar(istartc+1:istartc+inc),aout) levels
      istartc=istartc+inc+1

c*****************************************************************
c
c     Read in whether data is on observed levels (isoor=0) or
c     standard levels (isoor=1)
c
c*****************************************************************

      read(ichar(istartc:istartc),'(i1)') isoor
      istartc=istartc+1

c*****************************************************************
c
c     Read in number of variables in profile
c
c*****************************************************************

      read(ichar(istartc:istartc+1),'(i2)') nvar
      istartc=istartc+2

c*****************************************************************
c
c     Read in the variable codes (ip2()) and the whole profile 
c       error flags (ierror(ip2()))
c
c*****************************************************************

      do 30 n = 1,nvar

       read(ichar(istartc:istartc),'(i1)') inc
       write(aout(3:3),'(i1)') inc
       read(ichar(istartc+1:istartc+inc),aout) ip2(n)
       istartc=istartc+inc+1

       read(ichar(istartc:istartc),'(i1)') ierror(ip2(n))
       istartc=istartc+1

30    continue

c****************************************************************
c
c     Read in number of bytes in character data
c
c****************************************************************

      read(ichar(istartc:istartc),'(i1)') inc
      istartc=istartc+1
      if ( inc .gt. 0 ) then
       write(aout(3:3),'(i1)') inc
       read(ichar(istartc+1:istartc+inc),aout) inchad
       istartc=istartc+inc

c****************************************************************
c
c    Read in number of character and primary investigator arrays
c
c****************************************************************

      npi=0
      chars(1)(1:4)='NONE'
      chars(2)(1:4)='NONE'
      read(ichar(istartc:istartc),'(i1)') ica
      istartc=istartc+1

c****************************************************************
c
c    Read in character and primary investigator data
c      1 - originators cruise code
c      2 - originators station code
c      3 - primary investigators information
c
c****************************************************************

      do 45 nn=1,ica

       read(ichar(istartc:istartc),'(i1)') icn
       istartc=istartc+1

       if ( icn .lt. 3 ) then
        read(ichar(istartc:istartc+1),'(i2)') ns
        istartc=istartc+2
        chars(icn)= '               '
        chars(icn)= ichar(istartc:istartc+ns-1)
        istartc= istartc+ns
       else
        read(ichar(istartc:istartc+1),'(i2)') npi
        istartc=istartc+2
        do 505 n=1,npi
         read(ichar(istartc:istartc),'(i1)') inc
         write(aout(3:3),'(i1)') inc
         read(ichar(istartc+1:istartc+inc),aout) ipi(n,2)
         istartc=istartc+inc+1

         read(ichar(istartc:istartc),'(i1)') inc
         write(aout(3:3),'(i1)') inc
         read(ichar(istartc+1:istartc+inc),aout) ipi(n,1)
         istartc=istartc+inc+1
505     continue
       endif

45    continue

      endif

c****************************************************************
c
c     Read in number of bytes in secondary header variables
c
c****************************************************************

      read(ichar(istartc:istartc),'(i1)') inc
      istartc=istartc+1
      if ( inc .gt. 0 ) then
       write(aout(3:3),'(i1)') inc
       read(ichar(istartc+1:istartc+inc),aout) insec
       istartc=istartc+inc

c****************************************************************
c
c     Read in number of secondary header variables (nsecond)
c
c****************************************************************

       read(ichar(istartc:istartc),'(i1)') inc
       write(aout(3:3),'(i1)') inc
       read(ichar(istartc+1:istartc+inc),aout) nsecond
       istartc=istartc+inc+1

c****************************************************************
c
c     Read in secondary header variables (sechead())
c
c****************************************************************

       do 35 n = 1,nsecond

        read(ichar(istartc:istartc),'(i1)') inc
        write(aout(3:3),'(i1)') inc
        read(ichar(istartc+1:istartc+inc),aout) nn
        istartc=istartc+inc+1

        call charout(istartc,jsig2(nn),jprec2(nn),jtot2(nn),ichar,
     &  sechead(nn),bmiss) 

        isec(n) = nn

35     continue

       endif

c****************************************************************
c
c     Read in number of bytes in biology variables 
c
c****************************************************************

      read(ichar(istartc:istartc),'(i1)') inc
      istartc=istartc+1

      if ( inc .gt. 0 ) then
       write(aout(3:3),'(i1)') inc
       read(ichar(istartc+1:istartc+inc),aout) inbio
       istartc=istartc+inc

c****************************************************************
c
c     Read in number of biological variables (nbio)
c
c****************************************************************

      read(ichar(istartc:istartc),'(i1)') inc
      write(aout(3:3),'(i1)') inc
      read(ichar(istartc+1:istartc+inc),aout) nbio
      istartc=istartc+inc+1

c****************************************************************
c
c     Read in biological variables (bio())
c
c****************************************************************

      do 40 n = 1,nbio

       read(ichar(istartc:istartc),'(i1)') inc
       write(aout(3:3),'(i1)') inc
       read(ichar(istartc+1:istartc+inc),aout) nn
       istartc=istartc+inc+1

       call charout(istartc,jsigb(nn),jprecb(nn),jtotb(nn),ichar,
     & bio(nn),bmiss)

       ibio(n) = nn

40    continue

c****************************************************************
c
c     Read in biomass and taxonomic variables
c
c****************************************************************

      read(ichar(istartc:istartc),'(i1)') inc
      write(aout(3:3),'(i1)') inc
      read(ichar(istartc+1:istartc+inc),aout) nbothtot
      istartc=istartc+inc+1

      do 41 n = 1,nbothtot

       itaxtot=0
       read(ichar(istartc:istartc),'(i1)') inc
       write(aout(3:3),'(i1)') inc
       read(ichar(istartc+1:istartc+inc),aout) nn
       istartc=istartc+inc+1

       vtax(0,n)=nn

       do 42 n2 =1,nn

        itaxtot=itaxtot+1

        read(ichar(istartc:istartc),'(i1)') inc
        write(aout(3:3),'(i1)') inc
        read(ichar(istartc+1:istartc+inc),aout) itaxnum(itaxtot,n)
        istartc=istartc+inc+1
        call charout(istartc,jsigtax(itaxtot,n),jprectax(itaxtot,n),
     &   jtottax(itaxtot,n),ichar,vtax(itaxtot,n),bmiss)

        read(ichar(istartc:istartc),'(i1)') itaxerr(itaxtot,n)
        istartc=istartc+1

42     continue

41    continue
      endif

c****************************************************************
c
c     Read in measured and calculated depth dependent variables
c       along with their individual reading flags
c
c****************************************************************

      do 50 n = 1,levels

       if ( isoor.eq.0 ) then

        call charout(istartc,idsig(n),idprec(n),idtot(n),ichar,
     & depth(n),bmiss)

        read(ichar(istartc:istartc),'(i1)') iderror(n,0)
        istartc=istartc+1

       endif

       do 55 i = 1,nvar
     
        call charout(istartc,msig(n,ip2(i)),mprec(n,ip2(i)),
     & mtot(n,ip2(i)),ichar,temp(n,ip2(i)),bmiss)

       if ( temp(n,ip2(i)) .gt. bmiss ) then

       read(ichar(istartc:istartc),'(i1)') iderror(n,ip2(i))
       istartc=istartc+1

       else
    
        iderror(n,ip2(i))=0

       endif

55     continue

50     continue

       return

500    ieof = 1

       return
       end

C------------------------------------------------------------------
      SUBROUTINE CHAROUT(istartc,jsig,jprec,jtot,ichar,value,bmiss)
      
c     This subroutine reads a single real value from the
c     OCL ASCII format.  This value consists of four
c     components:  # significant figures, # total figures,
c     precision, and the value. 
     
c   Examples:

c     VALUE  Precision    OCL ASCII
c     -----  ---------    ---------
c     5.35       2        332535
c     5.         0        1105
c     15.357     3        55315357
c    (missing)            -           
     
c******************************************************
c     
c   Passed Variables:
c
c     istartc    - starting point to read in data
c     jsig       - number of significant figures in data value
c     jprec      - precision of data value
c     jtot       - number of figures in data value
c     ichar      - character array from which to read data
c     value      - data value
c     bmiss      - missing value marker
c
c*****************************************************

c*****************************************************
c
c   Character Array:
c
c     cwriter    - format statement (FORTRAN I/O)
c
c****************************************************

      character*6 cwriter
      character*(*) ichar
      
      data cwriter /'(fX.X)'/
      
c****************************************************
c     
c     Check if this is a missing value (number of 
c       figures = '-')
c
c****************************************************

      if ( ichar(istartc:istartc) .eq. '-' ) then

       istartc = istartc+1
       value = bmiss
       return

      endif
       
c****************************************************
c
c     Read in number of significant figure, total
c       figures and precision of value
c
c****************************************************

      read(ichar(istartc:istartc),'(i1)') jsig
      read(ichar(istartc+1:istartc+1),'(i1)') jtot
      read(ichar(istartc+2:istartc+2),'(i1)') jprec
      istartc=istartc+3

c****************************************************
c
c     Write these values into a FORTRAN format statement
c
c       e.g. "553" --> '(f5.3)'
c            "332" --> '(f3.2)'
c
c****************************************************

      write(cwriter(3:3),'(i1)') jtot
      write(cwriter(5:5),'(i1)') jprec

c****************************************************
c
c     Read in the data value using thhe FORTRAN 
c       format statement created above (cwriter).
c
c****************************************************

      read(ichar(istartc:istartc+jtot-1),cwriter) value

c****************************************************
c
c     Update the character array position (pointer)
c       and send it back to the calling program.
c
c****************************************************

      istartc=istartc+jtot

      return
      end
c
        SUBROUTINE SINTRP(CN,IMC,FN,IMF,CX,FX,xmiss)                        
        DIMENSION CN(IMC),CX(1),FN(IMF),FX(1)                        
C                                                                   
C --- SINTRP INTERPOLATES DATA FROM GRID CN TO GRID FN WHERE       
C ---    THE GRIDSPACE CN COMPLETELY OVERLIES THE GRIDSPACE FN.   
C                                                                
C ---    IMC =  DIMENSION OF CN                                 
C ---    IMF =  DIMENSION OF FN ( IN CALLING ROUTINE)          
C ---    CX = DISTRIBUTION OF POINTS IN CN                    
C ---    FX = DISTRIBUTION OF POINTS IN FN                   
C ---    CX AND FX MUST BE MONOTONICALLY INCREASING         
c
c    Xianhe Cao modefied (adding xmiss) on 2/13/97
c ---    xmiss = missing value
c
C                                                          
        IC=IMC                                                    
        IF=IMF                                                   
c        IF(CX(1).LE.FX(1).AND.CX(IMC).GE.FX(IMF))GOTO 10        
c        PRINT 9999                                             
c        GO TO 20                                              
10      DO 12 I=2,IC                                         
        IF(CX(I).GE.CX(I-1))GO TO 12                        
        PRINT 9998                                         
        GO TO 20                                          
12      CONTINUE                                         
        DO 16 I=2,IF                                    
        IF(FX(I).GE.FX(I-1))GO TO 16                   
        PRINT 9998                                    
        GO TO 20                                     
16      CONTINUE                                    
        GO TO 30                                                    
20      CONTINUE                                                   
        PRINT 9911, CX(1),CX(IC),FX(1),FX(IF)                     
        PRINT 9994                                               
        PRINT 9993,(CX(I),I=1,IC)                               
        PRINT 9982                                             
        PRINT 9993,(FX(I),I=1,IF)                             
*        STOP
        goto 130                                                 
30      CONTINUE                                            
C                                                                
C --- COME HERE IF GRID IS MONOTONIC                            
C                                                              
        DO 100 I=1,IF                                         
        DO 120 L=1,IC                                        
        IF(CX(L).LT.FX(I))GO TO 120                         
        IE=L                                               
        IW=L-1                                            
        IF(IW.EQ.0)IE=2                                  
        IF(IW.EQ.0)IW=1                                 
        GO TO 121                                            
120     CONTINUE                                            
121     CONTINUE                                           
c
        if(cn(iw).eq.xmiss.or.cn(ie).eq.xmiss) then
          fn(i)=xmiss
          go to 100
        endif
c
        A=FX(I)-CX(IW)                                    
        B=CX(IE)-FX(I)                                   
        TX=CX(IE)-CX(IW)                                
        FN(I  )=(   A*CN(IE)+B*CN(IW))/TX              
100     continue   
130     continue                                     
        RETURN                                        
9911    FORMAT(' CX(1),CX(IC),FX(1),FX(IF)',         
     *  /,1X,4E12.4)                                
9999    FORMAT(43H0GRID INTERPOLATED TO EXCEEDS GRID INT FROM)    
9998    FORMAT(30H0CX OR FX            MONOTONIC)                
9994    FORMAT(19H0CX    FX    FOLLOW)                          
9993    FORMAT(1X,10E12.4)                                     
9982    FORMAT(' FX')                                                 
        END                                                          
