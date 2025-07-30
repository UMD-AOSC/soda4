subroutine tempoa(nid, ivar, array)
!=======================================================================
!
!     this is the main subroutine to do objective analysis of temp.
!     call tobj for each level to do objective analysis
!     write out analysis & error of temperature anomaly
!     add analysis of anomaly to corrected guess field to get intact
!     analysis result
!
!
!     inputs:
!
!     n2 = unit number for reading in obs. data after checks
!     n3 = unit number for writing out analysis of anomaly  
!     n4 = unit number for writing out error of t anomaly  
!     ts = array containing corrected guess temp field
!
!     outputs:
!    
!     ts = array containing final analysis of temperature
!
!=======================================================================
 
    use model_grid
    use obs_types
    use var_types
    use loadts

    use mpi

 
    real, dimension(imt, jmt, km) :: array
    real, dimension(imt, jmt    ) :: a1, er1
 
    integer :: nproc, OMP_GET_NUM_PROCS
 
    nproc = OMP_GET_NUM_PROCS()
  ! if (myrank == 0) then
  !     write(6, *) 'myrank == 0, Available ', nproc, ' processors'
  ! end if

  ! if (myrank == 100) then
  !     write(6, *) 'myrank == 100, Available ', nproc, ' processors'
  ! end if

    !----------------------------------------------------------------------
    !     do objective analysis for each level up to kmb levels
    !----------------------------------------------------------------------
    if (var_type(ivar) .eq. 'sst' ) then
        k_start = 1
        k_end   = 1
    else
        if (var_type(ivar) .eq. 'xbt' .and. lgcl_sst) then
            k_start = 2
        else
            k_start = 1
        endif
        k_end = kmb
    endif
 
    if (var_type(ivar) .eq. 'ctd') then
        k_start = 1
        k_end   = kmb
    endif
 
    !----------------------------------------------------------------------
    !     get observation profiles for analysis
    !----------------------------------------------------------------------
    call getts(nid, ivar, k_end)
 
    do lev = k_start, k_end
        call tobj(nid, ivar, lev, a1, er1)
 
        !----------------------------------------------------------------------
        !     smoothing the result
        !----------------------------------------------------------------------
 
      ! print *, 'going to smoothing'
      ! call smth5pt(a1,lev)
      ! print *, 'back from smoothing'
 
        !----------------------------------------------------------------------
        !     bring in topography to make t on land zero
        !     add analysis of anomaly to corrected guess field
        !----------------------------------------------------------------------

      ! if (myrank == 0) then 
            do i = 1, imt; do j = 1, jmt;
                if (lev .gt. kmt(i, j)) then
                    array(i, j, lev) = 0.0
                else
                    array(i, j, lev) = array(i, j, lev) + a1(i, j)
                endif
            end do; end do;
      ! end if
    end do


    deallocate(jdatg, xlatg, xlngg, jlatg, ilngg, typeg, tempg, stat=istat)
    if (istat .ne. 0) stop 'dealloc loadts'

    return
end subroutine tempoa

 
 
 
subroutine getts(nid, ivar, k_end)
!=======================================================================
!     read temperature binned data and store in common batches 
!=======================================================================
 
    use params
    use io_units
    use obs_types
    use var_types
    use batches
    use loadts 

    use mpi
 
    real, dimension(kmb) :: tempt
 
    maxd = num_var_used(ivar)
 
    allocate(jdatg(maxd), xlatg(maxd), xlngg(maxd), typeg(maxd),  &
             jlatg(maxd), ilngg(maxd), tempg(maxd,k_end),         &
             stat=istat)
    if (istat .ne. 0) stop 'allocate loadts'
 
    iobs  = 0
    ibch  = 0
    nbmax = 0
    tempg = tmiss
 
    rewind io_var_work_unit(ivar)

    if (myrank == 0) then
        print *, '  '
        write(6, *) 'myrank == 0, getts    ', var_type(ivar), num_var_used(ivar)
    end if
 
    do nob = 1, num_var_used(ivar)
        ! exit loop when reaching end of file
        read(io_var_work_unit(ivar), 200, end=111)     &
            jdatg(nob), xlatg(nob), xlngg(nob),        &
            ilngg(nob), jlatg(nob), typeg(nob), ncnt,  & 
            (tempt(i), i=1, ncnt)
 
        do k = 1, ncnt
            tempg(nob, k) = tempt(k)
            iobs(k)       = iobs(k) + 1
        enddo
 
        !----------------------------------------------------------------------
        !     count number of data in each batch 
        !----------------------------------------------------------------------
 
!$omp parallel do
        do nb = 1, nbatch
            if (xlatg(nob) .gt. y_north(nb, 2)) cycle
            if (xlatg(nob) .lt. y_south(nb, 2)) cycle

            if (x_west(nb, 2) .gt. x_east(nb, 2)) then
                if (xlngg(nob) .lt. 180.) then
                    if (xlngg(nob) .gt. x_east(nb, 2)) cycle
                else
                    if (xlngg(nob) .lt. x_west(nb, 2)) cycle
                endif
            else
                if (xlngg(nob) .lt. x_west(nb, 2) .or. xlngg(nob) .gt. x_east(nb, 2)) cycle
            endif
 
            do k = 1, ncnt
                ibch(nb, k) = ibch(nb, k) + 1
            enddo
        end do
!omp end parallel do
    end do

111 continue

    do k = 1, k_end
        nbmax(k) = maxval(ibch(:, k))
        if (myrank == 0) then
          write(6, '(" myrank == 0, total obs for lev ",i2," = ",i6," max is ",i4)') k, iobs(k), nbmax(k)
        end if
    enddo

    if (myrank == 0) then
        print *, '  '
    end if
 
200 format(i6, 2f8.2, 2i4, a6, i4, 5f8.3/,3(10f8.3,/))
 
    return
end subroutine getts

 
 
 
subroutine tobj(nid, ivar, lev, psi, err)
!=======================================================================
!  Put discription from another tobj variant
!=======================================================================
 
    use model_grid
    use batches
    use var_types
    use loadbh

    use mpi
 
    real :: a2, er2
    real, dimension(imt, jmt) :: psi, err
    integer :: ixp, iyp
    integer :: m, mttl, in2, jn2
    integer :: tid, omp_get_thread_num

    integer :: i, j, k
 
    save ixp, iyp, m, a2, er2, in2, jn2
 
!$omp threadprivate(ixp, iyp, m, a2, er2, in2, jn2)
    lv = 0
    if (lgcl_vert .and. var_type(ivar) .ne. 'sst') lv = 1
    lvt = 1 + 2*lv
    if (lev .le. 2) then
        if (myrank == 0) then
            print *, 'use obs data for levels of', lvt
            print *, '  '
        end if
    endif
 
    psi  = 0.
    err  = 0.
    mttl = 0
 
!---------------------------------------------------------------------
! == for each grid point
!---------------------------------------------------------------------

!$omp parallel shared(psi, err)                    &
!$omp private(ib1, lk,                             &
!$omp kgtm, xlim, lvd, bttl, innb, indc_get, tid)  &
!$omp reduction(+:mttl)
 
!$omp do schedule(static)
  ! do ib1 = 1, nbatch
    do ib1 = bs, be
 
        !-----------------------------------------------------------
        !    allocate space to store obs data in one patch
        !-----------------------------------------------------------
        tid  = omp_get_thread_num()
        maxb = 0
 
        do lk = lev - lv, lev + lv
            if (lk .lt. 1 .or. lk .gt. kmb) cycle
            maxb = maxb + ibch(ib1, lk)
        end do

        if (maxb .lt. 3) then
          m = 0
          go to 48
        endif
 
        if (ib1.eq.5207)then
	  write(*,*) "--- Allocate memory for the batch ib1 =",ib1 &
	            ,", maxb =",maxb
	endif
        allocate(bti(maxb), btp(maxb), bxg(maxb), byg(maxb), bxm(maxb),  &
                 nii(maxb), njj(maxb), bym(maxb), btt(maxb), nlv(maxb),  &
                 stat=istat)
 
        !-----------------------------------------------------------
        !     get data for each batch by criterion xlim
        !     no. of data is limited to mbx
        !------------------------------------------------------------
        kgtm = 0
        m    = 0
        xlim = 0.0075

65      xlim = xlim*10.

	btp = "-"
	bti = 0.0
	bxg = 0.0 
	byg = 0.0
	bxm = 0.0
	bym = 0.0
	btt = 0.0
	nii = 0
	njj = 0
	nlv = 0

!        if((ib1.eq.5207).or.(ib1.eq.8568).or.(ib1.eq.9764)) then
!	   write(*,*) "myrunk =",myrank,"ib1 =",ib1,"lev =",lev
!	   write(*,*) "nii ->",nii
!	   write(*,*) "njj ->",njj
!	 end if
	   
        call getstm(ib1, ivar, lev, nid, xlim, m, lvd,  &
                    maxb, bti, btp, bxg, byg, nii, njj, bxm, bym, btt, nlv)
 
!        if((ib1.eq.5207).or.(ib1.eq.8568).or.(ib1.eq.9764)) then
!	   write(*,*) "---After getstm, ib1 =",ib1,"lev =",lev &
!	             ,"maxb =",maxb,", m =",m,", kgtm =",kgtm
!	   write(*,*) "nii ->",nii
!	   write(*,*) "njj ->",njj
!	   write(*,*) "bti ->",bti
!	   write(*,*) "bxg ->",bxg
!	   write(*,*) "byg ->",byg
!	   write(*,*) "bxm ->",bxm
!	   write(*,*) "bym ->",bym
!	 end if
 
        kgtm = kgtm + 1
 
        !----------------------------------------------------------------------
        ! m was 15 35 30 35   ??????
        !----------------------------------------------------------------------
        if (m .gt. mbx) then          
            if (kgtm .lt. 5) then
                xlim = xlim/5.
                go to 65
            else
                bttl = 0.0
                do i = 1, mbx
                    bttl = bttl + btt(i)
                enddo
                if (bttl .ne. 0.0) then
                    m = mbx
                    go to 48
                else
                    m = 0
                    go to 941
                endif
            endif
        endif
 
941     continue
 
        if (m .lt. 3 .or. lvd .eq. 0) then
            m = 0
            if (allocated(bti)) then
                deallocate(bti, btp, bxg, byg, nii, njj, bxm, bym, btt, nlv, stat=istat)
            end if
          ! if(istat.ne.0) stop 'dealloc loadbh'
        endif
 
        !----------------------------------------------------------------------
        !     zero out the transfer arrays (between tobj and toa)
        !----------------------------------------------------------------------
 
48      continue
 
        a2  = 0.
        er2 = 0.
 
        if (.not. (m .eq. 0)) then
            ! compute psi
            ! get each grid's location (m) and time (s), where we evaluate psi
            tp = nid*86400.
 
            ixp = i_model_west(ib1, 1)
            if (ixp .gt. 360) ixp = ixp - 360
            iyp = j_model_south(ib1, 1)
 
            ! get objective analysis and error
            call toa(lev, m, tp, ixp, iyp, a2, er2, ib1,  &
                     maxb, bti, btp, bxg, byg, nii, njj, bxm, bym, btt, nlv)
 
            er2 = 1.0

            deallocate(bti, btp, bxg, byg, nii, njj, bxm, bym, btt, nlv, stat=istat)
        end if  
 
        jn2 = j_model_south(ib1, 1)
        in2 = i_model_west (ib1, 1)
        if (in2 .gt. 360) in2 = in2 - 360

        psi(in2, jn2) = a2
        err(in2, jn2) = er2
 
        mttl = mttl + m
    end do  ! end of "do ib1 = bs, be"
!$omp end do
!$omp end parallel

    call mpi_barrier(MPI_COMM_WORLD, ierr)
    do j = 1, jmt; do i = 1, imt;
        call mpi_allreduce(psi(i, j), psi(i, j), 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)  !Ligang Chen: need to figure out.
        call mpi_allreduce(err(i, j), err(i, j), 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
    end do; end do;  

    call mpi_allreduce(mttl, mttl, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    if (myrank == 0) then
       write(6, *) 'Total sum of obs. for all patches at level', lev, ' is', mttl
      ! write(6, *) "psi(*, 100) = ", psi(:, 100)
    end if
 
    !------------------------------------------------------------
    ! == end analysis for grid point
    !--------------------------------------------------------------
 
    return    
end subroutine tobj

 
 
 
subroutine getstm(nb, ivar, lev, nid, xlm, icnt1, lvd,  &
                  mx, bti1, btp1, bxg1, byg1, nii1, njj1, bxm1, bym1, btt1, nlv1)
!=======================================================================
! Put discription from another getstm's variant
!=======================================================================
 
    use model_grid
    use constants
    use batches
    use var_types
    use loadts
    use loadbh

    use mpi
 
    character(len=3) :: mtyp
 
    integer :: nb, mx
    integer, save :: jdag, iobm, jobm, icnt, kbch 
    real   , save :: xlt, xlg, xd, yd, dt, dx, dy, phim, phid, phia
 
    dimension lll(lvt), temp(lvt)
    real, dimension(lvt) :: cc5
 
    real            , dimension(mx) :: bti1, bxg1, byg1, bxm1, bym1, btt1
    integer         , dimension(mx) :: nlv1, nii1, njj1
    character(len=3), dimension(mx) :: btp1
 
    data indc_get/0/
    data innb/0/
 
    save indc_get, innb, kcnt

 
!$omp threadprivate(indc_get, innb, kcnt,              &
!$omp xlt, xlg, xd, yd, dt, dx, dy, phim, phid, phia,  &
!$omp jdag, iobm, jobm, icnt, kbch)

    if (nb .ne. innb) then    
        indc_get = 0
        innb     = nb
    endif
 
    !------------------------------------------------------------------
    ! get batch center point in m 
    !------------------------------------------------------------------
 
    xpm = x_center(nb, 1) 
    ypm = y_center(nb, 1)
 
    icnt = 0
    kbch = 0
    if (indc_get .eq. 0) kcnt = 0
!    if (indc_get .ne. 0 .and. kcnt .lt. mbx) icnt = kcnt
 
    do ib = 1, maxd
        mtyp = typeg(ib)
        xlt  = xlatg(ib)
        xlg  = xlngg(ib)
 
!        if (indc_get .ne. 0 .and. kcnt .lt. mbx .and. mtyp .eq. var_type(ivar)) cycle
 
        !----------------------------------------------------------------------
        ! keep data near our box
        ! eliminate all data outside 300 km beyond patch
        !----------------------------------------------------------------------
 
        if (xlt .gt. y_north(nb, 2)) cycle
        if (xlt .lt. y_south(nb, 2)) cycle
 
        if (x_west(nb, 2) .gt. x_east(nb, 2)) then
            if (xlg .lt. 180.) then
                if (xlg .gt. x_east(nb, 2)) cycle
            else
                if (xlg .lt. x_west(nb, 2)) cycle
            endif
        else
            if (xlg .lt. x_west(nb, 2) .or. xlg .gt. x_east(nb, 2)) cycle
        end if
 
        llev = 0
        do kk = 1, lvt
            temp(kk) = 0.0
            lll (kk) = 0
        end do
 
        ! not use upper & lower data of altimetry (and sst)
        ! when doing vertical correlation
        lvv = 0
        if (lgcl_vert .and. (mtyp .eq. 'xbt' .or. mtyp .eq. 'ctd' .or. mtyp .eq. 'pst')) lvv = 1
 
        do ll = lev - lvv, lev + lvv
            if (ll .lt. 1 .or. ll .gt. kmb) cycle
            if (tempg(ib, ll) .eq. tmiss) cycle

            llev       = llev + 1
            temp(llev) = tempg(ib, ll)
            lll(llev)  = ll
        end do
 
        if (llev .eq. 0) cycle
 
        !----------------------------------------------------------------------
        ! calculate correlation function & eliminate data
        !----------------------------------------------------------------------
 
        jdag = jdatg(ib)
        iobm = ilngg(ib)
        jobm = jlatg(ib)

        dt = (nid - jdag)*86400.
	
        xd = xlg
        yd = xlt	

        dx = (xd - xpm)*111000.*cos((xlt*deg_to_radians))
        dy = (yd - ypm)*111000.	
 
        do nn = 1, llev
            cc5(nn) = cov33(dx, dy, dt, xlt,                               &
                            y_center(nb, 1), iobm, i_model_center(nb, 1),  &
                            jobm, j_model_center(nb, 1),                   &
                            lev, lev, lll(nn), mtyp)
        enddo

        do nn = 1, llev
            if (cc5(nn) .le. xlm) cycle

            icnt = icnt + 1
            if (mtyp .eq. var_type(ivar)) kbch = kbch + 1

            bti1(icnt) = float(jdag)*86400.
            btp1(icnt) = mtyp
            bxg1(icnt) = xlg
            byg1(icnt) = xlt
            nii1(icnt) = iobm
            njj1(icnt) = jobm
            bxm1(icnt) = xd
            bym1(icnt) = yd
            nlv1(icnt) = lll(nn)
            btt1(icnt) = temp(nn)
        end do
    end do    ! end of "do ib = 1, maxd"

 
    if (indc_get .eq.0 .or. kcnt .ge. mbx) kcnt = kbch

    if (.not. (icnt .eq. 0)) then
        phim = 0.
        phid = 0.
        do k = 1, icnt
            phim = phim + btt1(k)
            phid = phid + btt1(k)**2
        enddo

        phia = phid/float(icnt) - (phim/float(icnt))**2
        if (phia .lt. 0.0) phia = 0.0
        phid = sqrt(phia)
 
        if (phid .eq. 0.) icnt = 0
    end if 
 
    if (indc_get .eq. 0) indc_get = 1
 
    lvd = 0
    do l = 1, icnt
        if (nlv1(l) .eq. lev) lvd = lvd + 1  
    enddo
 
      icnt1 = icnt
       
    return
end subroutine getstm

 
 
 
subroutine toa(lev, num, tp, ixp, iyp, a2, er2, itch,  &
               mx, bti1, btp1, bxg1, byg1, nii1, njj1, bxm1, bym1, btt1, nlv1)
!=======================================================================
! Put here discription from another toa's variant
!=======================================================================
 
    use model_grid
    use constants
    use loadbh

    use mpi
 
    double precision aa(num, num), zzz(num), det(2), rcond, ttt
    dimension c(num), ainv(num, num)
    integer, dimension(num) :: ipvt
 
    integer :: num, ixp, iyp, itch
    integer :: iobj, jobj, iobi, jobi, iip, jip
    real    :: x1, y1, bttm, bttd, b, amp, temp3, yj, yi
    real    :: theta, tp, err, a2, er2, covd, e, thetm, ainvm, cainm
    real, dimension(num, num) :: dtt, dxx, dyy
    real, dimension(num     ) :: ddt, ddx, ddy
 
    real   , dimension(mx) :: bti1, bxg1, byg1, bxm1, bym1, btt1
    integer, dimension(mx) :: nlv1, nii1, njj1
    character(len=3), dimension(mx) :: btp1
 
    save yj, yi, iobj, iobi, jobj, jobi, iip, jip

 
!$omp threadprivate(yj,yi,iobj,iobi,jobj,jobi,iip,jip)
    x1 = xt(ixp)
    y1 = yt(iyp)
 
    !----------------------------------------------------------------------
    ! estimate observation statistics
    !----------------------------------------------------------------------
 
    bttm = 0.
    bttd = 0.
    do j = 1, num
        bttm = bttm + btt1(j)
        bttd = bttd + btt1(j)**2
    enddo
 
    bttd = sqrt(bttd/float(num) - (bttm/float(num))**2)
    bttm = bttm/float(num)
 
    !----------------------------------------------------------------------
    ! compute aa - the observation position corrrelation matrix
    !----------------------------------------------------------------------
    b     = 1./160000.
    amp   = b**2/3.
    temp3 = (3./b**2)

    do j = 1, num
        yj   = byg1(j)
        iobj = nii1(j)
        jobj = njj1(j)
 
        do i = 1, num
            dxx(i, j) = (bxm1(j) - bxm1(i))*111000.*cos((yt(iyp)*deg_to_radians))
            dyy(i, j) = (bym1(j) - bym1(i))*111000.
            dtt(i, j) = bti1(j) - bti1(i)
 
            yi = byg1(i)
 
            iobi = nii1(i)
            jobi = njj1(i)
 
            covd = cov33(dxx(i, j), dyy(i, j), dtt(i, j), yj, yi, iobj, iobi,  &
                         jobj, jobi, lev, nlv1(j), nlv1(i), btp1(j))
 
            !----------------------------------------------------------------------
            !     the off-diagonal elements contain the covariance
            !----------------------------------------------------------------------
 
            aa(i, j) = covd * bttd**2
 
            !----------------------------------------------------------------------
            ! the diagonal elements contain the signal variance
            ! decide e (error covariance = % of ob covariance)
            !----------------------------------------------------------------------
 
            if (i .eq. j) then
              ! if (btp1(i) .eq. 'sst') e = 0.2
              ! if (btp1(i) .eq. 'nst') e = 0.2
              ! if (btp1(i) .eq. 'xbt') e = 0.2
              ! if (btp1(i) .eq. 'ctd') e = 0.4
              ! if (btp1(i) .eq. 'pst') e = 1.3    
                if (btp1(i) .eq. 'sst') e = 1.4
                if (btp1(i) .eq. 'nst') e = 1.4
                if (btp1(i) .eq. 'xbt') e = 1.4
                if (btp1(i) .eq. 'ctd') e = 1.4
              ! if (btp1(i) .eq. 'pst') e = 1.8
                if (btp1(i) .eq. 'pst') e = 4.9
 
                if (btp1(i) .eq. 'geo') e = 0.75
                if (btp1(i) .eq. 'es1' .or. btp1(i) .eq. 'es2') e = 0.5
                ! caonew
                if (btp1(i) .eq. 'js1') e = 0.5
                if (btp1(i) .eq. 'tpx') e = 0.5

                ! caoerr
              ! if (btp1(i) .eq. 'sgo') e = 0.85
              ! if (btp1(i) .eq. 'se1' .or. btp1(i) .eq. 'se2') e = 0.75
              ! if (btp1(i) .eq. 'sj1') e = 0.75
              ! if (btp1(i) .eq. 'stp') e = 0.75
                if (btp1(i) .eq. 'sgo') e = 2.0 
                if (btp1(i) .eq. 'se1' .or. btp1(i) .eq. 'se2') e = 2.0 
                if (btp1(i) .eq. 'sj1') e = 2.0 
                if (btp1(i) .eq. 'stp') e = 2.0 
 
                aa(i, j) = aa(i, j) + bttd**2*e
            endif
        end do  ! end of "do i = 1, num"
    end do      ! end of "do j = 1, num"
 
    !----------------------------------------------------------------------
    ! invert aa
    ! sub. dchco & dchdi are dgeco & dgedi from oak ridge national lib.
    !----------------------------------------------------------------------
 
    call dchco(aa, num, num, ipvt, rcond, zzz)
 
    ttt = 1.0 + rcond
    if (ttt .eq. 1.0) stop 'dchco'
 
    call dchdi(aa, num, num, ipvt, det, zzz, 01)
 
    do i = 1, num
    do j = 1, num
        ainv(i, j) = aa(i, j)
    enddo
    enddo
 
    !----------------------------------------------------------------------
    ! begin main computation
    !----------------------------------------------------------------------
    iip = ixp
    jip = iyp
    yi  = yt(iyp)
 
    !----------------------------------------------------------------------
    ! compute c
    !----------------------------------------------------------------------
    do j = 1, num
        ddx(j) = (x1 - bxm1(j))*111000.*cos((yt(iyp)*deg_to_radians))
        ddy(j) = (y1 - bym1(j))*111000.
        ddt(j) = tp - bti1(j)
 
        iob = nii1(j)
        job = njj1(j)
        yj  = byg1(j)

        covd = cov33(ddx(j), ddy(j), ddt(j), yi, yj, iip, iob,  &
                     iyp, job, lev, lev, nlv1(j), btp1(j))
 
        c(j) = covd * bttd**2
    end do
 
    !----------------------------------------------------------------------
    ! compute the summations
    !----------------------------------------------------------------------
 
    thetm = 0.
    ainvm = 0.
    cainm = 0.
 
    do j = 1, num; do i = 1, num;
        cainm = cainm + c(j)*ainv(j, i)
        ainvm = ainvm + ainv(i, j)
        thetm = thetm + ainv(i, j)*btt1(j)
    enddo; enddo;
 
    thetm = thetm/ainvm
    theta = thetm
 
    !----------------------------------------------------------------------
    ! compute theta and err. err is normalized by 2. times the var.
    !----------------------------------------------------------------------
 
    do j = 1, num; do i = 1, num;
        theta = theta + c(j) * ainv(j, i) * (btt1(i) - thetm)
    enddo; enddo;
 
    covd = amp*temp3
    covd = covd/(60.*86400.)
    err  = covd * bttd**2 + (1. - cainm)**2/ainvm
 
    do j = 1, num; do i = 1, num;
        err = err - c(i) * ainv(i, j) *c(j)
    enddo; enddo;
 
    if (err .lt. 0.0) err = .000001
    err = err/(bttd**2)
 
    a2  = theta
    er2 = err
 
    return
end subroutine toa
 
 

 
! caocov33
! function cov33(dx, dy, dt, lat1, lat2, i1, i2, j1, j2, lev, lev1, lev2)
function   cov33(dx, dy, dt, lat1, lat2, i1, i2, j1, j2, lev, lev1, lev2, typ)
!=======================================================================
!   lon1, lon2 = the longitudes (degree) of two points
!   lat1, lat2 = the latitudes (degree) of two points
!                use the middle latitude in cov33
!=======================================================================
 
    use model_grid
 
    character(len=3) :: typ

    integer :: i1, i2, j1, j2, lev, lev1, lev2
    integer, save :: numsfh1, numsfh2

    real :: lat1, lat2, dx, dy, dt
    real, save :: z, x, y, t, sfh1, sfh2, lat_t

 
!$omp threadprivate(z, x, y, t, sfh1, sfh2, lat_t, numsfh1, numsfh2)
    lat_t = abs((lat1 + lat2)/2.)
    z     = zt(lev)
    if (z.ge.700.0) z=700.0   ! below 700.0 meters covariance doesn't changes with depth
    if (lat_t.ge.60.0) lat_t=60.0  ! for latitudes higher than 60, covariance is as it is at 60 degree
 
    ! depth "z" in meter
    x = 1000.*(450. + (375. - 450.)/50.*lat_t)*(1200. - z*(1. - 0.02*lat_t))/1200.
    y = 1000.*(250. + (375. - 250.)/50.*lat_t)*(1200. + z*(1. - 0.02*lat_t))/1200.

    t = 60.*24.*3600.
 
  ! cov33 = exp(-abs(dx)/x - abs(dy)/y - abs(dt)/t)
    cov33 = exp(-sqrt((dx/x)**2 + (dy/y)**2) - abs(dt)/t)
    ! caovert
    if (lgcl_vert) cov33 = cov33*(cvrt(lev1, lev2))
    ! caocov33
    if (z .gt.500.0) return
 
    ! -- get first guess sfh on obs position
    call vobs_4grd(i1, j1, sfh_first_guess, sfh1, numsfh1)
    call vobs_4grd(i2, j2, sfh_first_guess, sfh2, numsfh2)
 
    if (numsfh1 .ne. 0 .and. numsfh2 .ne. 0) cov33 = cov33 * exp(-((sfh1 - sfh2)/rmssfh)**2)
 
    return
end function cov33

 
 
 
subroutine vobs_4grd(ii, jj, array, value, num)
!=======================================================================
!  X. Cao write on 11/2/99
!
!  for an obs position, get its value through interpolating from 
!  surrounding 4 model grid values
!=======================================================================
 
    use model_grid

    use mpi
 
    integer :: ii, jj, num
    real    :: value
    dimension array(imt, jmt)
 
    integer :: i, j, ii2, idi, jj2, idj
    ! caocovindp 2/01
  ! ii = indp(rlon, xt, imt)
  ! jj = indp(rlat, yt, jmt)
 
    num   = 0
    value = 0.0
 
    ii2 = ii - 1
    idi = -1
    jj2 = jj - 1

    idj = -1
    if (ii .eq. 1) then
        !??????
        ii2 = 360
        idi = 359
    endif
    if (jj .eq. 1) then
        jj2 = 2
        idj = 1
    endif
 
    do i = ii, ii2, idi; do j = jj, jj2, idj;
        if (kmt(i, j) .ne. 0) then
            num   = num + 1
            value = value + array(i, j)
        endif

        if (i < 1 .or. i > imt .or. j < 1 .or. j > jmt) then
            write(*, *) "ERROR: in vobs_4grd, after averaging!"
        end if
    enddo; enddo;

    if (num .ne. 0) value = value/float(num)
 
    return
end subroutine vobs_4grd
