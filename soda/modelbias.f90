! caosbias   11/05
! subroutine modelbias(typ, ntu, nut, mdt, biasc)
subroutine   modelbias(typ, nut,      mdt, biasc)
 
!=======================================================================
!   read binned data and correct model bias
!=======================================================================
!
!   ntu - the order of calling this subroutine (for sin & cos)
!   nut - unit to read binned data;
!   mdt - maximum number of observations
!   bias - corrected model bias
!       np - EOF's number
!	x,y - horizontal coordinates
!	k - depth
!
!   tau(np,k) - EOF's time coefficients
!   lam(np,k) - EOF's eigen values
!   G5(x,y,k,np) - EOFs over 5x5 degree grid
!   G1(x,y,k,np) - EOFs over model grid
!
!=======================================================================
 
    use bias
 
    real :: sum, sum1
    character(len=3) :: ntype, typ

    double precision, allocatable, dimension(:,:) :: bcov
    double precision, allocatable, dimension(:  ) :: zzz
    integer         , allocatable, dimension(:  ) :: xobs, yobs, ipvt

    double precision, dimension(2) :: det
    double precision               :: rcond, ttt

    real   , dimension(imt, jmt, kmc    ) :: biasc
    real   , dimension(nx , ny , kmc    ) :: bias5
    real   , dimension(kmb              ) :: stobs
    real   , dimension(nx , ny , kmc    ) :: pobs
    real   , dimension(nx , ny , kmc, np) :: G5
    real   , dimension(imt, jmt, kmc, np) :: G1
    real   , dimension(np , kmc         ) :: tau, lam
    integer, dimension(nx , ny, kmc     ) :: nobs
    integer, dimension(kmc              ) :: n5obs
 
    print *, ' model bias called -', typ 

    pai = 3.141592654

    ! caosbias tt
    ! -- for biasc 

    nbs = 1
    nbk = kmc
    if (typ .eq. 'sst') nbk = 1
    if (typ .eq. 'tmp') nbs = 2
 
    do i_p = 1, np; do i_z = 1, kmc; 
        do i_y = 1, ny; do i_x = 1, nx;
            ! caosbias  !! doesn't need G5, G1
            G5(i_x, i_y, i_z, i_p) = G5_t(i_x, i_y, i_z, i_p)

            ! if (typ .eq. 'sst' .or. typ .eq. 'tmp')
            !   * G5(i_x, i_y, i_z, i_p) = G5_t(i_x, i_y, i_z, i_p)
            ! if (typ .eq. 'slt') 
            !   * G5(i_x, i_y, i_z, i_p) = G5_s(i_x, i_y, i_z, i_p)
        enddo; enddo;
 
        do i_y = 1, jmt; do i_x = 1, imt;
            ! caosbias
            G1(i_x, i_y, i_z, i_p) = G1_t(i_x, i_y, i_z, i_p)

            ! if (typ .eq. 'sst' .or. typ .eq. 'tmp')
            !   * G1(i_x, i_y, i_z, i_p) = G1_t(i_x, i_y, i_z, i_p)
            ! if (typ .eq. 'slt') 
            !   * G1(i_x, i_y, i_z, i_p) = G1_s(i_x, i_y, i_z, i_p)
        enddo; enddo;
    enddo; enddo;
 
    ! make array tau
    do k = 1, kmc
        tau(1, k) = 1.0
        lam(1, k) = 1.0

        !caosbias
      ! tau(3, k) = sin(2.0*pai*10.0*ntu/360.0)
        tau(2, k) = cos(2.0*pai*(10.0*ntu - 4.0)/365.0)
        tau(3, k) = sin(2.0*pai*(10.0*ntu - 4.0)/365.0)

        do l = 2, 3
          lam(l, k) = 1.0
        enddo
    enddo
     
  ! print *, ' TAU -0'
  ! do k = 1, kmc
  !     print *, (tau(ip, k), ip = 1, 3)
  ! enddo

    !----------------------------------------------------------
    ! bias forecast
    !----------------------------------------------------------
    do i = 1, nx; do j = 1, ny; do k = 1, kmc;
        sum = 0.0
        do ip = 1, np
            sum = sum + tau(ip, k)*G5(i, j, k, ip)
        enddo
        bias5(i, j, k) = sum
    enddo; enddo; enddo
 
    !==========================================================
    ! start bias correction for the Indian Ocean
    !==========================================================
     
    ! all initial vlues set to zero
    do i = 1, nx; do j = 1, ny; do k = 1, kmc;
        pobs(i, j, k) = 0
        nobs(i, j, k) = 0
    enddo; enddo; enddo;

    nob  = 0
    nob1 = 0
    rewind nut
 
    !---------------------------------------------------------
    ! read observation - first guess data from the unit "nut"
    !---------------------------------------------------------
    do n = 1, 9999999
        ! exit loop when reaching end of file
        read(nut, 222, end=110) jd, xlat, xlng, ii, jj, ntype, ncnt, (stobs(i), i=1, ncnt)

        !caoslt 10/03
        if (typ .eq. 'slt') then
            do i = 1, ncnt
                stobs(i) = stobs(i)/1000.0
            enddo
        endif

        if (ncnt .gt. kmc) ncnt = kmc

        nob = nob + 1
        if (nob .gt. mdt) then
            print *, 'No. of observations is more than maximum'
            print *, 'nob', nob, mdt 
            stop 'modelbias'
        end if

        !------------------------------------------------------------------
        ! separate area 20E-125E; 40S-30N and bin data by 5x5 degree bins 
        !------------------------------------------------------------------
        if ((xlng .gt. 20.0) .and. (xlng .lt. 125.0)) then
            if ((xlat .gt. -40.0) .and. (xlat .lt. 30.0)) then
                ix = floor(xlng/5.0) + 1
                iy = floor((xlat + 60)/5.0) + 1
                ! test
                if (ix .lt. 0  .or. iy .lt. 0 ) print *, ' lt : ix,iy', ix, iy
                if (ix .gt. nx .or. iy .gt. ny) print *, ' gt :ix,iy' , ix, iy
 
                do k = 1, ncnt
                    pobs(ix, iy, k) = pobs(ix, iy, k) + stobs(k)
                    nobs(ix, iy, k) = nobs(ix, iy, k) + 1
                enddo
            endif
        endif  
    end do

110 continue


    ! get averages in 5x5 degree bins
    do i = 1, nx; do j = 1, ny; do k = 1, kmc;
        if (nobs(i, j, k) .ne. 0) then
            ! caopr
          ! pobs(i, j, k) = pobs(i, j, k) / float(nobs(i, j, k))
            pobs(i, j, k) = (pobs(i, j, k) / float(nobs(i, j, k))) - bias5(i, j, k)
        endif
    enddo; enddo; enddo;

    ! get number of 5x5 bins with data for each level
    do k = 1, kmc
        n5obs(k) = 0

        do i = 1, nx; do j = 1, ny;
            if (nobs(i, j, k) .ne. 0) then
                n5obs(k) = n5obs(k) + 1
            endif
        enddo; enddo;
    enddo

    !----------------------------------------------------
    ! start loop for bias correction for diffrent levels
    !----------------------------------------------------
    !caosbias tt
  ! do k = 1, kmc
    do k = nbs, nbk
        n5 = n5obs(k)
        if (n5 .gt. 0) then
            allocate(xobs(n5), yobs(n5), ipvt(n5), bcov(n5, n5), zzz(n5), stat=istat)
            if (istat .ne. 0) stop 'allocate arrays'
       
            ! get positions of 5x5 bins with data
            no = 1
            do i = 1, nx; do j = 1, ny;
                if (nobs(i, j, k) .ne. 0) then
                    xobs(no) = i
                    yobs(no) = j
                    no       = no + 1
                endif
            enddo; enddo;

            ! get bias covariance matrix
            do i = 1, n5; do j = 1, n5;
                bsum = 0
                do ip = 1, np
                    bsum = bsum + lam(ip, k) * G5(xobs(i), yobs(i), k, ip) * G5(xobs(j), yobs(j), k, ip)
                enddo
                bcov(i, j) = bsum
            enddo; enddo;

            ! add "observations" error covariance matix to bias covariance matrix
            do i = 1, n5
                ! caosbias 
                if (typ .eq. 'slt') then
                    bcov(i, i) = bcov(i, i) + 0.756/4.5
                else
                    bcov(i, i) = bcov(i, i) + 0.756
                endif

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! "0.756" should be calculated as a sum of residual EOF's eigenvalues
                ! plus observations error
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            enddo 

            !-------------------------------------------------------------------
            ! invert matrix bcov
            ! (follow by subroutine toa)
            ! sub. dchco & dchdi are dgeco & dgedi from oak ridge national lib.
            !-------------------------------------------------------------------
            call dchco(bcov, n5, n5, ipvt, rcond, zzz)
 
            ttt = 1.0 + rcond
            if (ttt .eq. 1.0) stop 'dchco'
 
            call dchdi(bcov, n5, n5, ipvt, det, zzz, 01)

            ! update bias coefficients
            do ip = 1, np
                sum = 0

                do j = 1, n5
                    sum1 = 0
                    do i = 1, n5
                        sum1 = sum1 + bcov(i, j) * pobs(xobs(i), yobs(i), k)
                    enddo
                    sum = sum + lam(ip, k) * G5(xobs(j), yobs(j), k, ip) * sum1
                enddo

                tau(ip, k) = tau(ip, k) + sum
            enddo

            ! get corrected bias
            do i = 1, 125; do j = 1, jmt;
                sum = 0
                do ip = 1, np
                    ! ttt
                    if (tau(ip, k) .gt. 2.0 ) tau(ip, k) = 2.0
                    if (tau(ip, k) .lt. -2.0) tau(ip, k) = -2.0

                    sum = sum + tau(ip, k) * G1(i, j, k, ip)
                enddo
                biasc(i, j, k) = sum
            enddo; enddo;
 
            deallocate(xobs, yobs, ipvt, bcov, zzz, stat=istat)
            if (istat .ne. 0) stop 'deallocate arrays'
        endif  ! end of "if (n5 .gt. 0) then"
    end do

  ! print *, ' TAU -1'
  ! do k=1,kmc
  !     print *, (tau(ip, k), ip=1,3)
  ! enddo
 
    !==========================================================
    ! bias has been corrected in the Indian Ocean
    !==========================================================
     
     
    !==========================================================
    ! start bias correction for the Pacific Ocean
    !==========================================================
 
    ! all initial vlues set to zero
    do i = 1, nx; do j = 1, ny; do k = 1, kmc;
        pobs(i, j, k) = 0
        nobs(i, j, k) = 0
    enddo; enddo; enddo;

    nob  = 0
    nob1 = 0
    rewind nut
 
    !---------------------------------------------------------
    ! read observation - first guess data from the unit "nut"
    !---------------------------------------------------------
    do n = 1, 9999999
        ! exit loop when reaching end of file
        read(nut, 222, end=111) jd, xlat, xlng, ii, jj, ntype, ncnt, (stobs(i), i=1, ncnt)

        ! caoslt 10/03
        if (typ .eq. 'slt') then
            do i = 1, ncnt
                stobs(i) = stobs(i)/1000.0
            enddo
        endif

        if (ncnt .gt. kmc) ncnt = kmc

        nob = nob + 1
        if (nob .gt. mdt) then
            print *, 'No. of observations is more than maximum'
            print *, 'nob', nob, mdt         
            stop 'modelbias'
        end if

        !------------------------------------------------------------------
        ! separate area 125E-285E; 40S-55N and bin data by 5x5 degree bins 
        !------------------------------------------------------------------
        if ((xlng .gt. 125.0) .and. (xlng .lt. 285.0)) then
            if ((xlat .gt. -40.0) .and. (xlat .lt. 55.0)) then
                ix = floor(xlng/5.0) + 1
                iy = floor((xlat+60)/5.0) + 1
                ! test
                if (ix .lt. 0  .or. iy .lt. 0 ) print *, ' lt : ix,iy', ix, iy
                if (ix .gt. nx .or. iy .gt. ny) print *, ' gt :ix,iy' , ix, iy
 
                do k = 1, ncnt
                    pobs(ix, iy, k) = pobs(ix, iy, k) + stobs(k)
                    nobs(ix, iy, k) = nobs(ix, iy, k) + 1
                enddo
            endif
        endif  
    end do

111 continue

    !---------------------------------
    ! get averages in 5x5 degree bins
    !---------------------------------
    do i = 1, nx; do j = 1, ny; do k = 1, kmc;
        if (nobs(i, j, k) .ne. 0) then
            ! caopr
          ! pobs(i, j, k) = pobs(i, j, k) / float(nobs(i, j, k))
            pobs(i, j, k) = (pobs(i, j, k) / float(nobs(i, j, k))) - bias5(i, j, k)
        endif
    enddo; enddo; enddo;

    !-------------------------------------------------
    ! get number of 5x5 bins with data for each level
    !-------------------------------------------------
    do k = 1, kmc
        n5obs(k) = 0

        do i = 1, nx; do j = 1, ny;
            if (nobs(i, j, k) .ne. 0) then
                n5obs(k) = n5obs(k) + 1
            endif
        enddo; enddo;
    enddo

    !----------------------------------------------------
    ! start loop for bias correction for diffrent levels
    !----------------------------------------------------
    ! caosbias tt
  ! do k = 1, kmc
    do k = nbs, nbk
        n5 = n5obs(k)
        if (n5 .gt. 0) then
            allocate(xobs(n5), yobs(n5), ipvt(n5), bcov(n5, n5), zzz(n5), stat=istat)
            if (istat .ne. 0) stop 'allocate arrays'
       
            !-------------------------------------
            ! get positions of 5x5 bins with data
            !-------------------------------------
            no = 1
            do i = 1, nx; do j = 1, ny;
                if (nobs(i, j, k) .ne. 0) then
                    xobs(no) = i
                    yobs(no) = j
                    no       = no + 1
                endif
            enddo; enddo;

            !----------------------------
            ! get bias covariance matrix
            !----------------------------
            do i = 1, n5; do j = 1, n5;
                bsum = 0
                do ip = 1, np
                    bsum = bsum + lam(ip, k) * G5(xobs(i), yobs(i), k, ip) * G5(xobs(j), yobs(j), k, ip)
                enddo
                bcov(i, j) = bsum
            enddo; enddo;

            !---------------------------------------------------------------------
            ! add "observations" error covariance matix to bias covariance matrix
            !---------------------------------------------------------------------
            do i = 1, n5
                ! caosbias
                if (typ .eq. 'slt') then
                    bcov(i, i) = bcov(i, i) + 0.756/4.5
                else
                    bcov(i, i) = bcov(i, i) + 0.756
                endif

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! "0.756" should be calculated as a sum of residual EOF's eigenvalues
                ! plus observations error
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            enddo 

            !-------------------------------------------------------------------
            ! invert matrix bcov
            ! (follow by subroutine toa)
            ! sub. dchco & dchdi are dgeco & dgedi from oak ridge national lib.
            !-------------------------------------------------------------------
            call dchco(bcov, n5, n5, ipvt, rcond, zzz)
 
            ttt = 1.0 + rcond
            if (ttt .eq. 1.0) stop 'dchco'
 
            call dchdi(bcov, n5, n5, ipvt, det, zzz, 01)

            !--------------------------
            ! update bias coefficients
            !--------------------------
            do ip = 1, np
                sum = 0

                do j = 1, n5
                    sum1 = 0
                    do i = 1, n5
                        sum1 = sum1 + bcov(i, j) * pobs(xobs(i), yobs(i), k)
                    enddo
                    sum = sum + lam(ip, k) * G5(xobs(j), yobs(j), k, ip) * sum1
                enddo

                tau(ip, k) = tau(ip, k) + sum
            enddo

            !--------------------
            ! get corrected bias
            !--------------------
            do i = 126, 285; do j = 1, jmt;
                sum = 0

                do ip = 1, np
                    ! ttt
                    if (tau(ip, k) .gt. 2.0 ) tau(ip, k) = 2.0
                    if (tau(ip, k) .lt. -2.0) tau(ip, k) = -2.0
              
                    sum = sum + tau(ip, k) * G1(i, j, k, ip)
                enddo

                biasc(i, j, k) = sum
            enddo; enddo;
 
            deallocate(xobs, yobs, ipvt, bcov, zzz, stat=istat)

            if (istat .ne. 0) stop 'deallocate arrays'
        endif  ! end of "if (n5 .gt. 0) then"
    end do


  ! print *, ' TAU -2'
  ! do k = 1, kmc
  !     print *, (tau(ip, k), ip=1,3)
  ! enddo
 
    !==========================================================
    ! bias has been corrected in the Pacific Ocean
    !==========================================================
     
      
    !==========================================================
    ! start bias correction for the Atlantic Ocean
    !==========================================================
     
    !-------------------------------
    ! all initial vlues set to zero
    !-------------------------------
    do i = 1, nx; do j = 1, ny; do k = 1, kmc;
        pobs(i, j, k) = 0
        nobs(i, j, k) = 0
    enddo; enddo; enddo;

    nob  = 0
    nob1 = 0
    rewind nut
 
    !---------------------------------------------------------
    ! read observation - first guess data from the unit "nut"
    !---------------------------------------------------------
    do n = 1, 9999999
        read(nut, 222, end=112) jd, xlat, xlng, ii, jj, ntype, ncnt, (stobs(i), i=1, ncnt)

        ! caoslt 10/03
        if (typ .eq. 'slt') then
            do i = 1, ncnt
                stobs(i) = stobs(i)/1000.0
            enddo
        endif

        if (ncnt .gt. kmc) ncnt = kmc

        nob = nob + 1
        if (nob .gt. mdt) then
            print *, 'No. of observations is more than maximum'
            print *, 'nob', nob, mdt         
            stop 'modelbias'
        end if

        !------------------------------------------------------------------
        ! separate area 285E-360E; 40S-55N and bin data by 5x5 degree bins 
        !------------------------------------------------------------------
        if ((xlng .gt. 285.0) .and. (xlng .lt. 360.0)) then
            if ((xlat .gt. -40.0) .and. (xlat .lt. 55.0)) then
                ix = floor(xlng/5.0) + 1
                iy = floor((xlat+60)/5.0) + 1
                ! test
                if (ix .lt. 0  .or. iy .lt. 0 ) print *, ' lt : ix,iy', ix, iy
                if (ix .gt. nx .or. iy .gt. ny) print *, ' gt :ix,iy' , ix, iy
           
                do k = 1, ncnt
                    pobs(ix, iy, k) = pobs(ix, iy, k) + stobs(k)
                    nobs(ix, iy, k) = nobs(ix, iy, k) + 1
                enddo
            endif
        endif  
    end do

112 continue

    !---------------------------------
    ! get averages in 5x5 degree bins
    !---------------------------------
    do i = 1, nx; do j = 1, ny; do k = 1, kmc;
        if (nobs(i, j, k) .ne. 0) then
            ! caopr
          ! pobs(i, j, k) = pobs(i, j, k) / float(nobs(i, j, k))
            pobs(i, j, k) = (pobs(i, j, k) / float(nobs(i, j, k))) - bias5(i, j, k)
        endif
    enddo; enddo; enddo;

    !-------------------------------------------------
    ! get number of 5x5 bins with data for each level
    !-------------------------------------------------
    do k = 1, kmc
        n5obs(k) = 0

        do i = 1, nx; do j = 1, ny;
            if (nobs(i, j, k) .ne. 0) then
                n5obs(k) = n5obs(k) + 1
            endif
        enddo; enddo;
    enddo

    !----------------------------------------------------
    ! start loop for bias correction for diffrent levels
    !----------------------------------------------------
    ! caosbias tt
  ! do k = 1, kmc
    do k = nbs, nbk
        n5 = n5obs(k)
        if (n5 .gt. 0) then
            allocate(xobs(n5), yobs(n5), ipvt(n5), bcov(n5, n5), zzz(n5), stat=istat)
            if (istat .ne. 0) stop 'allocate arrays'
       
            !-------------------------------------
            ! get positions of 5x5 bins with data
            !-------------------------------------
            no = 1
            do i = 1, nx; do j = 1, ny;
                if (nobs(i, j, k) .ne. 0) then
                    xobs(no) = i
                    yobs(no) = j
                    no       = no + 1
                endif
            enddo; enddo;

            !----------------------------
            ! get bias covariance matrix
            !----------------------------
            do i = 1, n5; do j = 1, n5;
                bsum = 0
                do ip = 1, np
                    bsum = bsum + lam(ip, k) * G5(xobs(i), yobs(i), k, ip) * G5(xobs(j), yobs(j), k, ip)
                enddo
                bcov(i, j) = bsum
            enddo; enddo;

            !---------------------------------------------------------------------
            ! add "observations" error covariance matix to bias covariance matrix
            !---------------------------------------------------------------------
            do i = 1, n5
                ! caosbias
                if (typ .eq. 'slt') then
                    bcov(i, i) = bcov(i, i) + 0.756/4.5
                else
                    bcov(i, i) = bcov(i, i) + 0.756
                endif

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! "0.756" should be calculated as a sum of residual EOF's eigenvalues
                ! plus observations error
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            enddo 

            !-------------------------------------------------------------------
            ! invert matrix bcov
            ! (follow by subroutine toa)
            ! sub. dchco & dchdi are dgeco & dgedi from oak ridge national lib.
            !-------------------------------------------------------------------
            call dchco(bcov, n5, n5, ipvt, rcond, zzz)
 
            ttt = 1.0 + rcond
            if (ttt .eq. 1.0) stop 'dchco'
 
            call dchdi(bcov, n5, n5, ipvt, det, zzz, 01)

            !--------------------------
            ! update bias coefficients
            !--------------------------
            do ip = 1, np
                sum = 0

                do j = 1, n5
                    sum1 = 0
                    do i = 1, n5
                        sum1 = sum1 + bcov(i, j) * pobs(xobs(i), yobs(i), k)
                    enddo
                    sum = sum + lam(ip, k) * G5(xobs(j), yobs(j), k, ip) * sum1
                enddo

                tau(ip, k) = tau(ip, k) + sum
            enddo

            !--------------------
            ! get corrected bias
            !--------------------
            ! cao !!!!!!!!!!!!!!
          ! do i = 286, 361
            do i = 286, 360; do j = 1, jmt;
                sum = 0

                do ip = 1, np
                    ! ttt
                    if (tau(ip, k) .gt. 2.0 ) tau(ip, k) = 2.0
                    if (tau(ip, k) .lt. -2.0) tau(ip, k) = -2.0
                    sum = sum + tau(ip, k) * G1(i, j, k, ip)
                enddo

                biasc(i, j, k) = sum
            enddo; enddo;
 
            deallocate(xobs, yobs, ipvt, bcov, zzz, stat=istat)
            if (istat .ne. 0) stop 'deallocate arrays'
        endif  ! end of "if (n5 .gt. 0) then"
    end do     ! end of "do 202 k = nbs, nbk"


  ! print *, ' TAU -3'
  ! do k=1,kmc
  !   print *, (tau(ip,k), ip=1,3)
  ! enddo

    ! caosbias    11/05 
    do j = 160, jmt
        biasc(:, j, :) = 0.0
    enddo
 
    !==========================================================
    ! bias has been corrected in the Atlantic Ocean
    !==========================================================
 
222 format(i6, 2f8.2, 2i4, a6, i4, 5f8.3/,5(10f8.3,/))
 
    return
end subroutine modelbias
