PROGRAM SODA
!=======================================================================
!
!     this routine performs ocean data assimilation of temperature.
!     it gets t field of time level tau from slab memory as first guess,
!     uses sst to correct mixing layer t, then performs oi analysis          
!     and puts analysed t back to slab memory for both tau & tau-1 time
!     levels, which makes t,s,u,v the same for both time levels.
!
!     the main procedures are:
!
!  --  replaces the mixed layer temperature of guess with sea surface temp.
!  --  gets anomaly of temp(tob-tguess) with respect to observation station
!      by using corrected guess & 60 day binned observation xbt data       
!  --  corrects temperature anomaly field and eleminates any value
!      greater than 4 sigma, twice by subrontine dtatpr
!  --  calculates the mean & variance of final temperature anomaly
!  --  performs objective analysis of temp anomaly field using dr. gandin's
!      version in 5-degree * 5-degree section; gets analysis & error fields
!      of temp anomaly; adds anomaly to corrected guess = final analysis
!      of temperature.
!
!=======================================================================
 
    use params
    use io_units
    use model_grid
    use obs_types
    use var_types
    use one_d
    use bias
    use mpi 
 
    namelist /soda_inputs/  &
        year, month, day,                               &
        do_deepocn_relaxation, deepocn_relaxation_time, &
        lgcl_vert, lgcl_bias,                           &
        lgcl_ctd, lgcl_sst, lgcl_nst, lgcl_altm,        &
        lgcl_geosat,lgcl_ers1_c,lgcl_ers1_g,lgcl_ers2,  &
        lgcl_jason1, lgcl_topex


    ! Initial settings
    call initialize_mpi

  ! write(6, *) "myrank = ", myrank
 
    ! read namelist
    open (99, file='soda_in', status='old', form='formatted')
    read (99, nml=soda_inputs)
    close(99)

  ! if (myrank == 0 .or. myrank == 511) then
  !     write(*, *) "myrank = ", myrank
  !     write(stdout, nml=soda_inputs)
  ! end if


    ! -- julian date of doing analysis  
    nid  = jday(month, day, year)
    idst = nid
    idnd = nid
 
    if (myrank == 0) then
        write(*, *) "myrank = ", myrank
        write(stdout, '(/,a,i6,i3,i3,3x,a,i8/)')  &
            ' analysis time - year, month, day & julian date:',  &
            year, month, day, '&', nid
    end if
 
    ! set the variables to be analyzed
    var_type(1) = 'xbt'
    var_type(2) = 'ctd'
    var_type(3) = 'sst'
 
    ! open all files needed
    obs_read = 0
    call open_files

    if (myrank == 0) then
        do kr = 1, obs_read
            write(stdout, '(2x,i3,2x,a3,2x,i3,2x,i3)') kr, obs_type(kr), io_obs_unit(kr), obs_intt(kr)
        end do
    end if
 
    !  either read ...or generate the grid
    call generate_grid_xy
 
    ! define batches for data assimilation
    call getbatch
 
    ! read in vertical correlation table
    if (lgcl_vert) call read_cvrt      
    ! ccaoalpha
    ! read in dsdt for bias & altimetry assim
 
    ! if (lgcl_bias .or. lgcl_altm) call read_dsdt
     call read_dsdt
     
    !----------------------------------------------------------------------
    !  get first guess field (model forecast)          
    !----------------------------------------------------------------------

  ! if (myrank == 0) then 
      ! write(*, *) "Before reading model direct, myrank = ", myrank
        call read_model_direct
      ! write(*, *) "Done reading model_direct, myrank = ", myrank
  ! end if
  ! call mpi_bcast(sfh_first_guess, imt*jmt, MPI_REAL, sfh_first_guess, imt*jmt  &
  !     , MPI_REAL, 0, MPI_COMM_WORLD, ierr)
  ! call mpi_bcast(temp_first_guess, imt*jmt*km, MPI_REAL, temp_first_guess, imt*jmt*km  &
  !     , MPI_REAL, 0, MPI_COMM_WORLD, ierr)
  ! call mpi_bcast(salt_first_guess, imt*jmt*km, MPI_REAL, salt_first_guess, imt*jmt*km  &
  !     , MPI_REAL, 0, MPI_COMM_WORLD, ierr)


    ! set corrected first guess to first guess
 
    ctemp_first_guess = temp_first_guess
    csalt_first_guess = salt_first_guess
 
    ! Read levitus temp and salt
  ! if (myrank == 0) write(*, *) "Before read_levitus ..."
    call read_levitus
 
    ! generate a kmt file - number of levels available
  ! if (myrank == 0) write(*, *) "Before make_kmt ..."
    call make_kmt
 
    if (lgcl_bias) then
        ntu = (month - 1)*3 + int(day/10) + 1
        print *, ' Bias NTU =', ntu
        call read_bias
    endif
 
    !----------------------------------------------------------------------
    !     do temp level 1 OI using sst obs, if lgcl_sst
    !     ts(:,:,1) changed, and put to sst_obs(:,:) for sub. sstcor
    !----------------------------------------------------------------------
 
  ! if (myrank == 0) write(*, *) "Before lgcl_sst ..."
    if (lgcl_sst) then
        ivar = 3
        ts   = temp_first_guess
        ss   = salt_first_guess
     
      ! if (myrank == 0) then  ! only on CPU will write out the anomaly work file, WRONG! it also involves some calculation for all CPUs!
            call write_sst_anom(nid, 'sst')
            if (lgcl_nst) call write_sst_anom(nid, 'nst')
      ! end if
        call mpi_barrier(MPI_COMM_WORLD, ierr)
 
        if (lgcl_bias) then
            call modelbias('sst', io_var_work_unit(ivar), num_var_used(ivar), temp_bias)
            ts(:, :, 1) = temp_first_guess(:, :, 1) + temp_bias(:, :, 1)

            call write_sst_anom(nid, 'sst')
            if (lgcl_nst) call write_sst_anom(nid, 'nst')
        endif
 
      ! get mixed layer depth from first guess
	call mixlayer_depth(nid)

        if (myrank == 0) then
            write(*, *) "myrank == 0, before tempoa(sst)!"
        end if
 
        call tempoa(nid, ivar, ts)
 
        sst_obs(:, :) = ts(:, :, 1)
        sss_obs(:, :) = ss(:, :, 1)
        !- mix correction      endif
 
        !----------------------------------------------------------------------
        !     replace the mixed layer temp of first guess field
        !     with sea surface temp at each latitude
        !----------------------------------------------------------------------
 
      !  call sstcor(nid)
	ctemp_first_guess(:, :, 1) = sst_obs(:, :)
      
    endif  ! end of "if (lgcl_sst) then"


  ! if (myrank == 0) then
  !     write(*, *) "before call find_obs_index('xbt'), myrank = ", myrank
  ! end if

    if (myrank == 0) then  ! be careful of this "if" statement 
        call find_obs_index('xbt')
        call read_xbt(nid, obs_index)
 
        if (lgcl_ctd) then
          ! if (myrank == 0) then
          !     write(*, *) "before call find_obs_index('ctd'), myrank = ", myrank
          ! end if  

            call find_obs_index('ctd')
            call read_xbt(nid, obs_index)
        endif
    end if

    call mpi_barrier(MPI_COMM_WORLD, ierr)
 
    if (lgcl_altm) call read_altm(nid, month)
 
!----------------------------------------------------------------------
!  objective analysis of temp 
!----------------------------------------------------------------------
    ts = ctemp_first_guess 
    ss = csalt_first_guess 

  ! if (myrank == 0) then  ! some calculation and assignment to global vars involved so all CPUs need to participate
      ! write(*, *) "before call write_var_anom(nid), myrank = ", myrank
        call write_var_anom(nid)
  ! end if
    call mpi_barrier(MPI_COMM_WORLD, ierr)

    if (myrank == 0) then
        write(*, *) "myrank == 0, after write_var_anom(temp, salt) and before tempoa(temp)!"
    end if   

    if (lgcl_bias) then
        ivar = 1
        call modelbias('tmp', io_var_work_unit(ivar), num_var_used(ivar), temp_bias)

        ! caosbias  k starts from 2 !!!! because of "write_var_anom"
        do k = 2, kmc
            ts(:, :, k) = ctemp_first_guess(:, :, k) + temp_bias(:, :, k)
        enddo
 
        ivar = 2
        ! ttt
        ! call modelbias('slt', io_var_work_unit(ivar), num_var_used(ivar), salt_bias)
        salt_bias = temp_bias * dsdt(:,:,1:kmc)

        do k = 1, kmc
            ss(:, :, k) = csalt_first_guess(:, :, k) + salt_bias(:, :, k)
        enddo

        call write_var_anom(nid)
    endif  ! end of "if (lgcl_bias) then"

  ! if (myrank == 0) then
  !     write(*, *) "before call tempoa(nid, ivar, ts) and after mpi_barrier, myrank = ", myrank
  ! end if
 
    ivar = 1
    call tempoa(nid, ivar, ts)

  ! write(*, *) "after call tempoa(nid, ivar, ts), myrank = ", myrank

    if (myrank == 0) then
        write(6, *) "myrank == 0, after tempoa(temp)"
        write(6, *) ts(180, 65, 1), ts(180, 65, 2)
    end if
 
    
    call mpi_barrier(MPI_COMM_WORLD, ierr)

  ! if (myrank == 0) then
      ! write(*, *) "before call tempoa(nid, ivar, ss) and after mpi_barrier, myrank = ", myrank  ! succeeded here
  ! end if

    if (myrank == 0) then
        write(*, *) "myrank == 0, before tempoa(salt)!"
    end if   

    ivar = 2 
    call tempoa(nid, ivar, ss)

  ! write(*, *) "after call tempoa(nid, ivar, ss), myrank = ", myrank
 
    if (lgcl_bias) ts(:, :, 1) = sst_obs(:, :)
 
! correct temp and salt in the mixed layer
    ctemp_first_guess = ts
    csalt_first_guess = ss
 
    call mixlayer_cor(nid)
 
    ts = ctemp_first_guess
    ss = csalt_first_guess

    if (myrank == 0) then 
        print *, "myrank == 0, after tempoa(salt)"
        print *, 'OI results of sst, temp & salt (10 day)'
        print *, ts(180, 65, 1), ts(180, 65, 2), ss(180, 65, 1)
        print *, '   '
    end if
 
    !----------------------------------------------------------------------
    !  write out incremental correctors (divided by timesteps of 2nd loop run)
    !  for msmth = 2 loop
    !----------------------------------------------------------------------

    call mpi_barrier(MPI_COMM_WORLD, ierr)
    if (myrank == 0) then 
        ! delta_t is the amount of time during which the correction is applied
        ! (currently 10(days)*24.(hours)*3600. (seconds)
         
        ! lgchen: previously used
      ! delta_t = 10.*24.*3600. 
        delta_t = 1.0
        ts      = (ts - temp_first_guess)/(delta_t)

        ! Gena 20150808: I set all big correction values to 0.0 in soda.f90 routine (temp > 10.0, salt >2.0). 
        ! It happens near a cost line, I think because first guess regriding. 
        where (abs(ts(:, :, 1:34)) .gt. 10.0) ts(:, :, 1:34) = 0.0 

        ! Ligang Chen: check if do deep ocean relaxation to climatology
        if (do_deepocn_relaxation) then
            ! Gena 20150808 email: We do not need divide deepocn_relaxation_time by 10 in soda.f90, I changed this routine.
          ! ts(:, :, 35:50) = (levitus_temp(:, :, 35:50) - temp_first_guess(:, :, 35:50)) / (deepocn_relaxation_time/10.0)
            ts(:, :, 35:50) = (levitus_temp(:, :, 35:50) - temp_first_guess(:, :, 35:50)) / (deepocn_relaxation_time     )
            where (abs(ts(:, :, 35:50)) .gt. 1.0) ts(:, :, 35:50) = 0.0

          ! write(*, *) "levitus_temp(:, 40, 40): ", levitus_temp(:, 40, 40)
          ! write(*, *) "temp_first_guess(:, 40, 40): ", temp_first_guess(:, 40, 40)
        end if

        do k = 1, km
            write(io_output, rec=k) ts(:, :, k)
        end do

 
        ! factor 1000 to convert psu to g/g
      ! ss = (ss - salt_first_guess)/(1000.*delta_t)  ! lgchen: ss too small
        ss = (ss - salt_first_guess)/delta_t

        ! Gena 20150808: I set all big correction values to 0.0 in soda.f90 routine (temp > 10.0, salt >2.0). 
        ! It happens near a cost line, I think because first guess regriding.   
        where (abs(ss(:, :, 1:34)) .gt. 2.0) ss(:, :, 1:34) = 0.0

        ! Ligang Chen: check if do deep ocean relaxation to climatology
        if (do_deepocn_relaxation) then
          ! ss(:, :, 35:50) = (levitus_salt(:, :, 35:50) - salt_first_guess(:, :, 35:50)) / (deepocn_relaxation_time/10.0)
            ss(:, :, 35:50) = (levitus_salt(:, :, 35:50) - salt_first_guess(:, :, 35:50)) / (deepocn_relaxation_time     )
            where (abs(ss(:, :, 35:50)) .gt. 0.2) ss(:, :, 35:50) = 0.0

          ! write(*, *) "levitus_salt(:, 40, 40): ", levitus_salt(:, 40, 40)
          ! write(*, *) "salt_first_guess(:, 40, 40): ", salt_first_guess(:, 40, 40)
        end if

        ! cao!!!!!!!!!
        do k = 1, km
      !      do j = 155, 180; do i = 1, 360;
      !          ss(i,j,k) = 0.
      !      end do; end do;
            write(io_output, rec=km+k) ss(:,:,k)
        end do
 
        print *, 'OI increments of sst, temp & salt (10 day)'
        print *, ts(180, 65, 1), ts(180, 65, 2), ss(180, 65, 1)
        print *, '   '
    end if

    ! ttt
    if (lgcl_bias) then
        do k = 1, kmc
            do i = 1, imt; do j = 1, jmt;
                if (dsdt(i, j, k) .eq. 0.0) temp_bias(i, j, k) = 0.0
            enddo; enddo;

            write(io_temp_bias, rec=k) temp_bias(:, :, k)
            write(io_salt_bias, rec=k) salt_bias(:, :, k)
        enddo
    endif

    call finalize_mpi()

    stop
end program soda
