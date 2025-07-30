! caonst 1/20/05 for adding 'nst' data
! caochg 1/5/05  change layout, print, no important things
 
subroutine write_var_anom(nid)
    use model_grid
    use io_units
    use obs_types
    use var_types
    use one_d

    use mpi
 
    character(len=3) :: sty
 
    ! kp     : writing out profile levels
    ! kp_slt : 'pst' (xbt-->salt) writing out profile levels
 
    integer kp, kp_slt
    ! G.Ch (06.03.08)
    integer, dimension(500000) :: iis, jjs, ids
    ! end G.Ch (06.03.08)

 
    ! xbt-->salt flag for "call psd_slt": 0-do; 1-not  
 
    integer salt_flag  
 
    num_obs_read = 0
    num_var_used = 0
 
    idst = nid
    idnd = nid
 
    rewind io_var_work_unit(1)
    rewind io_var_work_unit(2)
 
    ! caonst
  ! if (lgcl_sst)then
  !     k_obs_start = 2
  ! else
  !     k_obs_start = 1
  ! endif
    k_obs_start = 1
    if (lgcl_sst) k_obs_start = 2
    if (lgcl_nst) k_obs_start = 3
 
  ! do kr = k_obs_start, obs_read
    ! G.Ch (06.03.08)
    do kr = obs_read, k_obs_start, -1
    ! end G.Ch (06.03.08)
        rewind io_obs_work_unit(kr)
 
        do icnt = 1, 500000
            td = 0.
            read(io_obs_work_unit(kr), 831, end=111) id, xla, xlo, n, (td(i), i=1, n)
            num_obs_read(kr) = num_obs_read(kr) + 1
 
            ! -- get ii & jj from latitude & longitude
            ii = indp(xlo, xt, imt)
            jj = indp(xla, yt, jmt)
 
            salt_flag = 0
 
            ! G.Ch (06.03.08)
            if (obs_type(kr) .eq. 'xbt') then
                do is = 1, num_var_used(2)
                    if (ii .eq. iis(is) .and. jj .eq. jjs(is) .and. id .eq. ids(is)) then
                        salt_flag = 1
                        exit
                    endif
                enddo
            endif
            ! end G.Ch (06.03.08)
 
            ! -- get anomaly  
            do k = 1, n
                if (k .gt. kmt(ii, jj)) then
                    kp     = k - 1
                    kp_slt = 0
                    go to 180
                endif
 
                ! -- for salt data ...test for limits
 
                ! caosatm
                ! caochg Is sty needed??
                sty = obs_type(kr)

              ! caosatm if (obs_type(kr) .eq. 'ctd') then
                if (obs_type(kr) .eq. 'ctd' .or. sty(1:1) .eq. 's') then
                    tp(k) = td(k) - ss(ii, jj, k)
 
                    ! caonst
                  ! if(abs(tp(k)).gt.1.0) then
                    if (td(k) .lt. 0.0 .or. td(k) .gt. 38. .or. abs(tp(k)) .gt. 1.0) then
                        kp = k - 1
                        go to 180
                    endif
                else
                    ! -- otherwise test temp data for limits
                    tp(k) = td(k) - ts(ii, jj, k)
                    ! caonst
                  ! if (td(k) .le. -3.0 .or. ts(ii, jj, k) .le. -3.0) then
                    if (td(k) .le. -3.0 .or. td(k) .gt. 33.0 .or. ts(ii, jj, k) .le. -3.0 .or. abs(tp(k)) .gt. 6.0 ) then
                        kp = k - 1 
                        if (salt_flag .eq. 0) kp_slt = kp
                        go to 180
                    endif
                endif
 
                ! -- generate salt difference from xbt data by t-s relationship
                if (obs_type(kr) .eq. 'xbt' .and. salt_flag .eq. 0) then
                    inversion = 0 
                    tmp       = td(k)
                    tmp1      = ts(ii, jj, k)
                    slt       = 0.
                    slt1      = 0.
      
                  !  call psd_slt(ii, jj, tmp1, slt1, inversion)
                  !  inversion = 0 
                  !  call psd_slt(ii, jj, tmp, slt, inversion)

                  ! if (inversion .eq. 0) sp(k) = slt - ss(ii, jj, k)
                  !  if(inversion.eq.0) sp(k)=slt-slt1
		  
		    sp(k) = dsdt(ii,jj,k)*(tmp-tmp1)
		  
                    ! caonst
                  ! if (abs(sp(k)) .gt. 1.0 .or. inversion .eq. 1) then
                  !  if (abs(sp(k)) .gt. 1.0 .or. inversion .eq. 1 .or. slt .lt. 0.0 .or. slt .gt. 38.0 ) then
                    if (abs(sp(k)) .gt. 1.0) then
                        kp_slt    = k - 1
                        salt_flag = 1
                    endif
                endif
            end do    ! end of "do k = 1, n"
 
            kp = n
            if (salt_flag .eq. 0) kp_slt = kp
 
180         continue
 
            if (.not. (kp .eq. 0)) then  ! go to 190
                ! -- write out salt data
                if (obs_type(kr) .eq. 'ctd' .or. sty(1:1) .eq. 's') then
                    ivar = 2

                    if (myrank == 0) then
                        write(io_var_work_unit(2), 920)  &
                        !&  id, xla, xlo, kp, obs_type(kr), ii, jj, (tp(i), i=1, kp)
                            id, xla, xlo, ii, jj, obs_type(kr), kp, (tp(i), i=1, kp)
                    end if

                    num_var_used(ivar) = num_var_used(ivar) + 1
                    ! G.Ch (06.03.08)
                    iis(num_var_used(ivar)) = ii
                    jjs(num_var_used(ivar)) = jj
                    ids(num_var_used(ivar)) = id
                    ! end G.Ch (06.03.08)
                else
                    if (obs_type(kr) .eq. 'xbt') then
                        ! Gena suggested to add ".and. salt_flag .eq. 0" on 20141217
                        if (kp_slt .ne. 0 .and. salt_flag .eq. 0) then  
                            ivar = 2

                            if (myrank == 0) then
                                write(io_var_work_unit(2), 920)  &
                                ! caonst
                                !&  id, xla, xlo, kp_slt, 'pst', ii, jj, (sp(k), k=1, kp_slt)
                                    id, xla, xlo, ii, jj, 'pst', kp_slt, (sp(k), k=1, kp_slt)
                            end if

                            num_var_used(ivar) = num_var_used(ivar) + 1
                        endif
                    endif
 
                    ! -- write out temp data
                    ivar = 1

                    if (myrank == 0) then
                        write(io_var_work_unit(1), 920)  &
                             ! caonst
                         !&  id, xla, xlo, kp, obs_type(kr), ii, jj,
                             id, xla, xlo, ii, jj, obs_type(kr), kp,  &
                             (tp(k), k=1, kp)
                    end if

                    num_var_used(ivar) = num_var_used(ivar) + 1
                endif
 
                if (id .lt. idst) idst = id
                if (id .gt. idnd) idnd = id
            end if
        end do  ! end of "do icnt = 1, 500000"

111     continue

        if (myrank == 0) then 
            print *, '  '
            print *, ' Writing data anorm', kr, obs_type(kr)
 
            write(stdout, '(3x,a,i6,2x,a,i6)') 'time window is', nid - obs_intt(kr), ' to', nid + obs_intt(kr)
            write(stdout, '(3x,a,i8)') 'this data read in for analysis is', num_obs_read(kr)
            ! caoall
            write(stdout, '(3x,a,2i8)') 'total data saved for analysis - temp & salt :', num_var_used(1), num_var_used(2)
        end if
    end do  ! end of "do kr = obs_read, k_obs_start, -1"

    if (myrank == 0) then 
        endfile io_var_work_unit(1)
        endfile io_var_work_unit(2)
    end if

    if (myrank == 0) then 
        write(stdout, 432) 'saved obs time coverage:  ', idst,' - ', idnd, ',  for ', nid
    end if
 
432 format(/2x, a26, i5, a3, i5, a8, i5)
831   format(i6, 2f9.3, i3, 34f8.3)  ! lgchen: modified.
! 831 format(i6, 2f9.3, i3, 22f8.3)  ! lgchen: originally 22, should be 34
    ! caonst
! 920 format(2x,' day',i5,' lat, lng ',2f10.2,'nmes ',i5,a3,' ij',2i4/, 5(10f8.3,/))
  920 format(i6, 2f8.2, 2i4, a6, i4, 5f8.3/,5(10f8.3,/))  ! ligang: original format
! 920 format(i6, 2f8.2, 2i4, a6, i4, 34f8.3)  ! ligang: can not change to this, read part do not correspond!
! 929 format(2x, i5, 2f10.2, i5, a3, 2i4, 5(10f8.3))
 
    return
end subroutine write_var_anom

 

 
subroutine write_sst_anom(nid, typ)
    use model_grid
    use io_units
    use obs_types
    use var_types
    use one_d

    use mpi

    character(len=3) :: typ
    save num_sst
 
    ivar = 3
    call find_obs_index(typ)
 
    idst = nid
    idnd = nid
    num_obs_read(obs_index) = 0
    if (typ .eq. 'sst') then
        num_sst            = 0
        num_var_used(ivar) = 0
    endif
    if (typ .eq. 'nst') then
        num_var_used(ivar) = num_sst
        num_sst            = 0
    endif
 
    rewind io_obs_unit(obs_index)
    if (typ .eq. 'sst') rewind io_var_work_unit(ivar)
 
    do isst = 1, 999999999
        read(io_obs_unit(obs_index), 832, end=211) id, xla, xlo, td(1)
 
        if (id .lt. nid - obs_intt(obs_index)) cycle
        if (id .gt. nid + obs_intt(obs_index)) exit
        num_obs_read(obs_index) = num_obs_read(obs_index) + 1
 
        if (xlo .lt. 0.0 .or. xlo .gt. 360.0) cycle
 
        ii = indp(xlo, xt, imt)
        jj = indp(xla, yt, jmt)
        if (ii .eq. 0 .or. jj .eq. 0 .or. ii .gt. imt .or. jj .gt. jmt) stop 'ij-sst'
 
        if (kmt(ii, jj) .eq. 0) cycle
 
        tp(1) = td(1) - ts(ii, jj, 1)

        if (td(1) .lt. -3.0 .or. td(1) .gt. 33.0 .or. abs(tp(1)) .gt. 9.0 ) cycle

        if (myrank == 0) then 
            write(io_var_work_unit(ivar), 921) id, xla, xlo, ii, jj, typ, tp(1)
        end if

        num_sst = num_sst + 1
 
        if (id .lt. idst) idst = id
        if (id .gt. idnd) idnd = id
    end do
211 continue

 
    num_var_used(ivar) = num_var_used(ivar) + num_sst

    if (myrank == 0) then 
        write(stdout, '(/,a,a3,a,i4,a,i4,a,i6,a,i6)')   &
            ' read ', obs_type(obs_index),              &
            ' data from unit', io_obs_unit(obs_index),  &
            ', intt=', obs_intt(obs_index),             &
            ', days from', idst, ' to', idnd
 
        write(stdout, '(4x,a,i8)') 'total sst data read for analysis is', num_obs_read(obs_index)
        write(stdout, '(4x,a,i8)') 'total sst data used for analysis is', num_sst
    end if

    if (myrank == 0) then 
        if (lgcl_nst) then
            if (typ .eq. 'nst') endfile io_var_work_unit(ivar)
        else
            endfile io_var_work_unit(ivar)
        endif
    end if
 
    return
 
832 format(1x, i5, 1x, f5.1, 1x, f5.1, 2x, f5.2)
921 format(i6, 2f8.2, 2i4, a6, '   1', f8.3)
 
end subroutine write_sst_anom
