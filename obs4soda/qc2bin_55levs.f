      PROGRAM  DATA_BIN                                                    
c
c===========================================================
c X. Cao   12/2004 
c 
c !!! NOTE: there are MISSING values in the datasets !!!
C                                                                       
C  THIS IS TO BIN TEMPERATURE PROFILE DATA TO A                         
C    5-DAY AVERAGE AND 1 DEGREE SQUARE FILE                             
C                                                                       
c ------  !!!!  decide lgp & DATE first  !!!!  ------
c
c===========================================================
c
      PARAMETER (IMAX=360,JMAX=180,KMAX=55,lgp=1,bmiss=-999.99)
      DIMENSION IBIN(IMAX,JMAX,KMAX),TBIN(IMAX,JMAX,KMAX)       
      DIMENSION T(KMAX),ISTART(lgp),IFINISH(lgp)
      character yst*4, ynd*4        
      character(len=51), parameter ::   tmpdir=
     &         '/aosc/iceland2/chepurin/okhotsk/DATA4SODA/SOFT/TMP/' 
C     &         '/aosc/okhotsk/chepurin/DATA4SODA/SOFT/TMP/' 
C                                                                       
C == DECIDE THE STARTING & ENDING DATES FIRST 
c  !!!!!!!!!!!!!  date must with 2 or 7 !!!!!!!!!!!!!!!
C
c      DATA ISTART/10087/
c      DATA IFINISH/10447/
C
      iyrst = 2023
      iyrnd = 2024
c
      ISTART = jday(1,1,iyrst)+3
*+2
*+1
      IDAT=ISTART(1)   
      write(*,*) istart, idat                                                 
      do 234 lp = iyrst, iyrnd
c
        iyst=lp
	iynd=lp
c
        write(yst,'(i4)') iyst
c
        print *, ' -- For year:', iyst
c      
*        IDAT = IDAT-5
        IFINISH = jday(12,31,lp)-1
c      
        do l=1,lgp
          print *, idat-2, ifinish(l)+2, ifinish(l)-idat+4
        enddo
C
c        OPEN(UNIT=11, file=tmpdir//'qc_tpot_'//yst//'_55levs.dat',
c     &       FORM='FORMATTED', STATUS='OLD')

c        OPEN(UNIT=8, file=tmpdir//'tpot_'//yst//'_55levs.bin', 
c     &       FORM='FORMATTED', status='unknown')
C
        OPEN(UNIT=11, file=tmpdir//'qc_s_'//yst//'_55levs.dat',
     &       FORM='FORMATTED', STATUS='OLD')
C
        OPEN(UNIT=8, file=tmpdir//'s_'//yst//'_55levs.bin', 
     *       FORM='FORMATTED', status='unknown')
C
        DO 9999 in=1,lgp
c
        ib=0
c
        IDAT=ISTART(in)                                                      
        IEND=IFINISH(in)
C                                                                        
 999    continue
c
        REWIND (11)
C                                                                       
        DO I=1,IMAX                                                   
        DO J=1,JMAX                                                   
        DO K=1,KMAX                                                   
         IBIN(I,J,K)=0                                                    
         TBIN(I,J,K)=0.0                                                  
        enddo
        enddo
        enddo
C                                                                       
C == READ IN DATA, FILL IN TO APPROPRIATE LOCATION 
C == IN EACH ARRAY        
C
        DO 100 ICNT=1,5000000                                            
C                                                                       
          DO K=1,KMAX                                                  
            T(K)=0.                                                       
          enddo
C                                                                       
          READ(11,111,END=200) iy, im, id, XLA,XLO,N,(T(I),I=1,N)
          id=jday(im,id,iy)
C                                                                     
C -- CHECK THE DATE                                                 
C                                                                       
          IF(ID.LT.IDAT-2.OR.ID.GT.IDAT+2) GO TO 100                      
C                                                                       
C -- CHECK THE DATA LOCATION OF BIN BOX                             
C          !!!  ( LAT: -89.5->89.5; LON:1.0->360.0 ) !!!                    
C                                                                    
          if(xlo.lt.0.5) xlo=xlo+360.0
          ILO=NINT(XLO)                                               
          if(xla.eq.90.0) xla=89.9
          JLA=NINT(XLA+90.5)
          if(ilo.lt.1.or.ilo.gt.imax.or.jla.lt.1.or.jla.gt.jmax) then
            print *, ' ILO, JLA', ilo, jla
            stop 'ilo&jla'
          endif
C                                                                       
C -- FILL IN THE DATA                                                   
C                                                                       
          DO 40 K=1,N                                                     
           if(t(k).eq.bmiss) go to 40
           IBIN(ILO,JLA,K)=IBIN(ILO,JLA,K)+1                             
           TBIN(ILO,JLA,K)=TBIN(ILO,JLA,K)+T(K)                          
40        CONTINUE                                                        
C                                                                       
100     CONTINUE                                                          
200     CONTINUE
C                                                                       
C == AVERAGE AND WRITE OUT BIN DATA                                     
C                                                                       
        DO 500 I=1,IMAX                                                   
        DO 500 J=1,JMAX                                                   
C
        N=0
        DO K=KMAX,1,-1                                                  
          IF(IBIN(I,J,K).NE.0) then
            n=k
            go to 800
          endif
        enddo                                                           
800     CONTINUE                                                        
C                                                                       
        IF(N.EQ.0) GO TO 500                                            
C                                                                    
        XLO=FLOAT(I)                                                
        XLA=FLOAT(J)-90.5                                               
C                                                                       
        DO 60 K=1,N                                                     
          if(IBIN(I,J,K).eq.0) then
            TBIN(I,J,K)=bmiss
            go to 60
          endif
          TBIN(I,J,K)=TBIN(I,J,K)/float(IBIN(I,J,K))
60      CONTINUE                                                        
C                                                                       
        WRITE(8,112) IDAT,XLA,XLO,N,(TBIN(I,J,K),K=1,N)
        ib=ib+1
 500    CONTINUE                                                          
C
c       WRITE(6,*) 'DAY =>', IDAT                                   
        IDAT=IDAT+5                                                       
        IF(IDAT.GT.IEND) GO TO 555                                        
        GO TO 999                                                         
C                                                                       
 555  CONTINUE                                                          
C                                                                       
      print *, 'total binned data ', ib 
c
      CLOSE(11)
      CLOSE(8)                                                          
C
9999  continue
C
234   continue
c
111   format(12x,i5,2i3,2f9.3,i3,55f8.3)
112   format(i6,2f9.3,i3,55f8.3)
c
      STOP                                                              
      END                                                               
c
      function jday(mon,iday,iyr)
c
c======================================================
c
c     compute the julian day corresponding to the
c     day on the gregorian calender
c
c======================================================
c
        dimension dpm(12)
        data dpm /31.0, 28.0, 31.0, 30.0, 31.0, 30.0, 31.0, 31.0, 30.0,
     $          31.0, 30.0, 31.0/

        dpm(2) = 28.0
        if(mod(real(iyr),4.) .eq. 0.)dpm(2) = 29.0
c
c first calculate days without leap years
c
        iyrs = iyr-1970
        days = 587.0+real(iyrs)*365.0
ccao?? need to think about iyrs+?? e.x.: iyr-1950, should iyrs+1
        num_leap = int(real(iyrs+1)/4.)
        if(iyrs.lt.0)then
          days = days - real(num_leap)
        else
          days = days + real(num_leap)
        end if
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
      return
      end
