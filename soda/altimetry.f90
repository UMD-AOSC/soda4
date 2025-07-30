! caoalpha 1/31/06 monthly alpha

! subroutine read_altm(nid) 
subroutine   read_altm(nid, m) 
    use params
    use io_units
    use altimetry

    ! total record number of each altimetry data
    integer, parameter ::  &
        krgeo = 62 ,  &
        kre1c = 18 ,  &
        kre1g = 13 ,  &
        krer2 = 77 ,  &
        krjs1 = 82 ,  &
        krtpx = 413 

    ! arrays of start & ending date of each record, each data
    real, dimension (krgeo) :: geo_datst, geo_datnd
    real, dimension (kre1c) :: e1c_datst, e1c_datnd
    real, dimension (kre1g) :: e1g_datst, e1g_datnd
    real, dimension (krer2) :: er2_datst, er2_datnd
    real, dimension (krjs1) :: js1_datst, js1_datnd
    real, dimension (krtpx) :: tpx_datst, tpx_datnd

    ! caoalpha
  ! read (ioalpha) alpha 
    read (ioalpha, rec=m) alpha 
  ! read (iodtdzm) dtdzm
  ! read (iodsdzm) dsdzm
    print *, '   '
  ! print *, ' alpha, dtdzm & dsdzm', alpha(180, 64, 2),
    print *, ' alpha', alpha(180, 64, 2)
  ! & dtdzm(180, 64, 2), dsdzm(180, 64, 2)

    if (lgcl_geosat) then
        write (6, '(/, a)') ' processing GEOSAT data'
        read (io_geo_date) geo_datst, geo_datnd
        call read_ssh(nid, krgeo, geo_datst, geo_datnd, 'geo', io_geo_msfh, 1)
        close (io_geo_date)
    endif

    if (lgcl_ers1_c) then
        write (6, '(/, a)') ' processing ERS1_C data'
        read (io_es1_date) e1c_datst, e1c_datnd
        call read_ssh(nid, kre1c, e1c_datst, e1c_datnd, 'es1', io_es1_msfh, 1)
        close (io_e1c_date)
    endif

    if (lgcl_ers1_g) then
        write (6, '(/, a)') ' processing ERS1_G data'
        read (io_es1_date) e1g_datst, e1g_datnd
        call read_ssh(nid, kre1g, e1g_datst, e1g_datnd, 'es1', io_es1_msfh, 1)
        close (io_e1g_date)
    endif

    if (lgcl_ers2) then
        write (6, '(/, a)') ' processing ERS2 data'
        read (io_es2_date) er2_datst, er2_datnd
        call read_ssh(nid, ikrer2, er2_datst, er2_datnd, 'es2', io_es2_msfh, 1)
        close (io_es2_date)
    endif

    if (lgcl_jason1) then
        write (6, '(/, a)') ' processing JASON-1 data'
        read (io_js1_date) js1_datst, js1_datnd
        call read_ssh(nid, krjs1, js1_datst, js1_datnd, 'js1', io_js1_msfh, 2)
        close (io_js1_date)
    endif

    if (lgcl_topex) then
        write (6, '(/, a)') ' processing TOPEX data'
        read (io_tpx_date) tpx_datst, tpx_datnd
        call read_ssh(nid, krtpx, tpx_datst, tpx_datnd, 'tpx', io_tpx_msfh, 2)
        close (io_tpx_date)
    endif

    ! caoalpha
    close(ioalpha)
  ! close(iodtdzm)
  ! close(ioalpha_s)
  ! close(iodsdzm)

    return
end subroutine read_altm




! subroutine read_ssh(nid, kr, datst, datnd, typ, iomsfh)    
subroutine   read_ssh(nid, kr, datst, datnd, typ, iomsfh, jd)    

!-----------------------------------------------------------------------
!
! read in altimetry data (geosat, ers1_c, ers1_g, ers2, topex)
! from unit io_obs_unit; then write resulted temperature and
! salinity profiles out to io_obs_work_unit with type:
! (for temperature) 'geo', 'es1', 'es2', 'js1' or 'tpx',
! (for salinity)    'sgo', 'se1', 'se2', 'sj1' or 'stp'
! for use in temperature OI
!
!                data cover period          cycle     total records
!               -------------------        -------   ---------------
!    geosat       11/8/86 - 9/30/89         17 days        62
!    ers1_c       4/14/92 - 12/20/93        35 days        18
!    ers1_g       3/24/95 - 6/2/96          35 days        13
!    ers2         4/29/95 - 8/26/02         35 days        77
!    jason-1      1/15/02 - 4/7/04          10 days        82
!    topex        12/31/92-3/18/04          10 days       413
!
! each altimetry data has kr*4 records, contains 4 (only 3 here) 
! fields:
!
!    ssh - sea surface height (anomaly)
!    time - Julian date of the observation
!    disp _ rsm of ssh (if >0.3, ignore it)
!    n - number of observations in each 1*1 square (not in this 
!        datasets)
!
! undefined value (land) is -9999.0
!
! time in the dataset is different than the Julian date used in 
! this code, therefore that will be substracted from time
!
! generally the start time of each record is the same as the
! ending time of the previous record; there are exceptions
!
! ssh field has substracted its mean, therefore we need to substract
! first guess dynamic height mean (mean calculted from non-altimetry
! run done in advance) here from calculated first guess dynamic 
! height at this OI time. the mean is with different time duration 
! for different data type:
!                     geosat : 10/1986 - 10/1987
!                     ers1_c :  1/1992 - 12/1993
!                     ers1_g :  5/1995 -  4/1996
!                     ers2   :  5/1995 -  4/1999
!                     jason-1:  1/2002 - 12/2003
!                     topex  :  1/1993 - 12/2002
!
! caoalpha
! mean dtdz (dtdzm) is needed for calculating temperature profiles 
! from altimetry data, the time duration for the mean is the whole 
! time period we use altimetry data in OI, undefined value is 
! -1.e34.
! mean dsdz (dsdzm) is the same as mean dtdz (dtdzm)
!
! exclude data out of time window, out of value range and on land
!
!-----------------------------------------------------------------------

    use params
    use io_units
    use obs_types
    use model_grid
    use altimetry
    
    integer ::  &
        kr,  & ! total record number of the data
        iomsfh ! mean ssh unit

    character (len=3) ::  &
        typ,  &! data typ: 'geo', 'es1', 'es2' 'js1' or 'tpx'
        sty    ! 'sgo', 'se1', 'se2' 'sj1' or 'stp'

    ! start & ending dates of the data records
    real, dimension (kr) :: datst, datnd   

    integer, parameter :: ig = 360, jg = 180

    ! caonew
    real, parameter ::  &
        dbase1 = 39999.0,  &! Julian date base for 'geo', 'es1' and 'es2'         
        dbase2 = 3795.0 ,  &! Julian date base for 'js1', 'tpx'
        undef  = -9990.0    ! undefined value in data (actually it's -9999.0)

    real :: dbase  ! Julian date base
      
    integer ::  &
        jst, jnd,       &! date time window for choosing data
        irst, irnd,     &! start & ending data record number chosed to read data
        jd,             &! idicator to chose dbase
        ii, jj,         &! nearest model grid number of data
        nb5, nb9, nbs9, &! "bad" data counts
        ktemp,          &! resulted temperature profile levels
        ksalt,          &! resulted salinity profile levels
        nw_altm,        &! number of temp profiles written out
        nw_altm_s,      &! number of salt profiles written out
        s_obs_index,    &! obs_index for altm->salt
        l, k, ir, i, j, nrec

    real ::  &
        rjst, rjnd,  &! real values of jst, jnd 
        glat, glon    ! data latitudes & longitudes

    real, dimension (imt, jmt) :: mssh  ! mean ssh of each data

    real, dimension (ig, jg) ::  &
        ssh,   &! observed sea surface height (anomaly)
        time,  &! Julian date of the observation
        disp    ! rsm of ssh 

    real, dimension (ksh) :: temp  ! resulted temperature profile
    real, dimension (ksh) :: salt  ! resulted salinity profile


    if (typ .eq. 'geo') sty = 'sgo'
    if (typ .eq. 'es1') sty = 'se1'
    if (typ .eq. 'es2') sty = 'se2'
    if (typ .eq. 'js1') sty = 'sj1'
    if (typ .eq. 'tpx') sty = 'stp'
 
    call find_obs_index(sty)
    s_obs_index = obs_index

    call find_obs_index(typ)

    jst = nid - obs_intt(obs_index)
    jnd = nid + obs_intt(obs_index)

    write (6, '(a,a,2i8)') '   time window:', typ, jst, jnd

    ! caonew
    if (jd .eq. 1) dbase = dbase1
    if (jd .eq. 2) dbase = dbase2

    do l = 1, kr
      datst(l) = datst(l) - dbase
      datnd(l) = datnd(l) - dbase
    end do

    read (iomsfh) mssh

    nw_altm   = 0
    nw_altm_s = 0


    ! get start and ending record number for reading the data by selected time window
    irst = 0
    irnd = 0
    rjst = float(jst)
    rjnd = float(jnd)

    if (rjst .ge. datst(kr) .and. rjst .le. datnd(kr)) then
        irst = kr
    else
        ! cao lt or le ??
        do l = 1, kr - 1
            if (rjst .ge. datst(l) .and. rjst .lt. datst(l+1)) then
                irst = l
                exit 
            endif
        enddo
    endif


    if (rjnd .ge. datst(1) .and. rjnd .le. datnd(1)) then
        irnd = 1
    else
        do l = 2, kr
            if (rjnd .gt. datnd(l-1) .and. rjnd .le. datnd(l)) then
                irnd = l
                exit
            endif
        enddo
    endif


    if (irst .eq. 0 .and. irnd .eq. 0) then
        write(6,*) '  no data chosen'
        return
    endif

    if (irst .eq. 0 .and. irnd .ne. 0) irst = 1
    if (irst .ne. 0 .and. irnd .eq. 0) irnd = kr
 
    print *, '  records read from record ', irst, ' to', irnd
    print *, '         from date ', datst(irst), '  to', datnd(irnd)


    ! read altimetry data, do quality control and get temperature profiles
    nb5  = 0
    nb9  = 0
    nbs9 = 0

    LOOP_TIME: do ir = irst, irnd
        !cao !!
        if(typ .eq. 'tpx' .and. ir .eq. 108) cycle

        nrec = (ir - 1)*3 + 1
        read (io_obs_unit(obs_index), rec=nrec) ssh
        nrec = nrec + 1
        read (io_obs_unit(obs_index), rec=nrec) time
        nrec = nrec + 1
        read (io_obs_unit(obs_index), rec=nrec) disp

        LOOP_SPACE_I: do i = 1, ig
        LOOP_SPACE_J: do j = 1, jg
            time(i, j) = time(i, j) - dbase

            if (time(i, j) .lt. undef .or. ssh(i, j) .lt. undef .or. disp(i, j) .lt. undef) cycle
            if (int(time(i, j)) .lt. jst .or. int(time(i, j)) .gt. jnd) cycle
            if (abs(ssh(i, j)) .ge. 1.0) cycle
            if (disp(i, j) .gt. 0.3) cycle
          
            glat = -89.5 + (j - 1)*1.0
            glon = 1.5 + (i - 1)*1.0
            ! cao??  if (abs(glat) .ge. 50.0) cycle LOOP_SPACE_J
            ii = indp(glon, xt, imt)
            jj = indp(glat, yt, jmt)

            if (kmt(ii, jj) .eq. 0) cycle

            ! altm -> temperature
            temp = 0.0
            ! ccaosatm  !!!!!!!  if (lgcl_sst) temp(1) use OI(1)
            temp(1) = ctemp_first_guess(ii, jj, 1)

            ! cao 12/22/2003 Gena changed:
          ! if (temp(1) .gt. 33.0 .or. temp(1) .lt. -4.0) cycle
            if (temp(1) .gt. 30.0 .or. temp(1) .lt. -3.0) go to 555

            do k = 2, ksh
                ktemp = k
                if (kmt(ii, jj) .lt. k) then
                    ktemp = k - 1
                    exit
                endif

                ! dtdzm: lower level - upper level !!! (mostly negative)
                ! ccaoalpha

                ! if (alpha(ii, jj, k) .gt. 0.0 .or. alpha(ii, jj, k) .lt. -400.  &
                !     .or. dtdzm(ii, jj, k) .gt. 0.0 .or. dtdzm(ii, jj, k) .lt. -99999.9) then

                if (alpha(ii, jj, k) .gt. 9000.0) then
                    temp(k) = ctemp_first_guess(ii, jj, k)

                    ! cao ??? 35 -9 ?? why don't check guess_cor in mixlayer
                    ! cao 12/22/2003 Gena changed:
                  ! if (temp(k) .gt. 35.0 .or. temp(k) .lt. -9.0) then
                    if (temp(k) .gt. 30.0 .or. temp(k) .lt. -3.0) then
                        go to 555 
                    else
                        cycle
                    endif
                endif           

                ! caonew unit ?????????
                ! ccaoalpha
                ! temp(k) = alpha(ii, jj, k) * dtdzm(ii, jj, k) *
                temp(k) = alpha(ii, jj, k) * (ssh(i, j) - (sfh_first_guess(ii, jj) - mssh(ii, jj)))    

                if (abs(temp(k)) .gt. 9.0) then
                    nb9 = nb9 + 1
                    go to 555            
                endif
 
                temp(k) = temp_first_guess(ii, jj, k) + temp(k)
                ! cao ??
              ! if (k .ge. 3 .and. (temp(k) - temp(k - 1)) .gt. 5.0) then
                if (zt(k) .ge. 40.0 .and. (temp(k) - temp(k-1)) .gt. 5.0) then
                    nb5 = nb5 + 1
                    go to 555           
                endif

                if (temp(k) .gt. 30.0 .or. temp(k) .lt. -3.0) go to 555
            end do  ! end of "do k = 2, ksh"

            ! write out temperature profiles
            write (io_obs_work_unit(obs_index), 222) int(time(i, j)), glat, glon, ktemp, (temp(k), k=1, ktemp)
            nw_altm = nw_altm + 1

555         continue


            ! altm -> salinity
            salt    = 0.0
            salt(1) = csalt_first_guess(ii, jj, 1)

            if (salt(1) .gt. 38.0 .or. salt(1) .lt. 27.0) cycle

            do k = 2, ksh
                ksalt = k

                if (kmt(ii, jj) .lt. k) then
                    ksalt = k - 1
                    exit     
                endif

                ! dsdzm: lower level - upper level 

                ! caoalpha
              ! if (alpha(ii, jj, k) .gt. 0.0 .or. alpha(ii, jj, k) .lt. -400.  &
              !     .or. dsdzm(ii, jj, k) .lt. -99999.9) then

                if (alpha(ii, jj, k) .gt. 9000.0) then
                    salt(k) = csalt_first_guess(ii, jj, k)

                    if (salt(k) .gt. 37.0 .or. salt(k) .lt. 27.0 ) then
                        cycle LOOP_SPACE_J  ! an outter loop
                    else
                        cycle ! current loop
                    endif
                endif

              ! caoalpha
              ! salt(k) = alpha(ii, jj, k) * dsdzm(ii, jj, k) *
                salt(k) = alpha(ii, jj, k) * dsdt(ii, jj, k) * (ssh(i, j) - (sfh_first_guess(ii, jj) - mssh(ii, jj)))

                if (abs(salt(k)) .gt. 0.3) then
                    nbs9 = nbs9 + 1
                    cycle LOOP_SPACE_J
                endif

                salt(k) = salt_first_guess(ii, jj, k) + salt(k)

                if (salt(k) .gt. 37.0 .or. salt(k) .lt. 27.0) cycle LOOP_SPACE_J
            end do  ! end of "do k = 2, ksh

            ! write out temperature profiles
            write(io_obs_work_unit(s_obs_index), 222) int(time(i, j)), glat, glon, ksalt, (salt(k), k=1, ksalt)
            nw_altm_s = nw_altm_s + 1

        end do LOOP_SPACE_J
        end do LOOP_SPACE_I
    end do LOOP_TIME  ! end of "do ir = irst, irnd"


    write(6, '(3x,a,i6)'  ) 'number of temp large inversion (5.0)', nb5
    write(6, '(3x,a,i6)'  ) 'number of large delta temp (9.0)'    , nb9
    write(6, '(3x,a,i6)'  ) 'number of large delta salt (0.3)'    , nbs9
    write(6, '(/a,a3,i9)' ) ' total temp written of ', typ, nw_altm
    write(6, '(a,a3,i9,/)') ' total salt written of ', sty, nw_altm_s

222 format(i6, 2f9.3, i3, 50f8.3)

    close(iomsfh)

    return
end subroutine read_ssh
