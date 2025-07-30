! caoalpha 1/31/2006 monthly alpha    
! caosbias 11/2005   add salt bias
! caonst   1/19/05   for adding "nst" data 
! caovert  1/7/2005  for adding vertical corr.

SUBROUTINE OPEN_FILES
    ! open files for assimilation
    use params
    use io_units
    use obs_types
    use var_types
    use bias

    use mpi 

    character(len=80) :: filename
  ! write(stdout, '(/,12x,a,/)') ' open files for assimilation'
 
    ! initialize
    iunit = 10
  ! lrec_sgl = imt * jmt * 4  ! lgchen: for gfortran
   lrec_sgl = imt * jmt      ! lgchen: for intel ifort
 
    ! open a file (files.dat) that contains a listing of opened files
    filename = 'files.dat'
    io_list  = iunit
    open(unit=iunit, file=filename, status='unknown')
    write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
    iunit = iunit + 1
 
 
    ! open files that contain model temp and salinity
    filename = 'first_guess.dat' 
    io_model = iunit
    open(unit=iunit, file=filename, status='old', form='unformatted', access='direct', recl=lrec_sgl)
    write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
    iunit = iunit + 1

    ! caosfh
    ! filename = 'first_guess_sfh.dat'
    ! iosfh = iunit
    ! open(unit=iunit, file=filename, status='old', form='unformatted')
    ! write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
    ! iunit = iunit + 1
    
    ! open files to put the incremental updates fields
    filename  = 'correctors.dat'
    io_output = iunit
    open(unit=iunit, file=filename, status='unknown', form='unformatted', access='direct', recl=lrec_sgl)
    write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
    iunit = iunit + 1
 
    ! Open work files for temperature and salt data
    filename            = 'temp_var.wrk'
    io_var_work_unit(1) = iunit
    open(unit=iunit, file=filename, status='unknown')
    write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
    iunit = iunit + 1
    
    filename            = 'salt_var.wrk'
    io_var_work_unit(2) = iunit
    open(unit=iunit, file=filename, status='unknown')
    write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
    iunit = iunit + 1
 
    if (lgcl_sst) then
        filename            = 'sst_var.wrk'
        io_var_work_unit(3) = iunit
        open(unit=iunit, file=filename, status='unknown')
        write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
        iunit = iunit + 1
 
        ! open files for sst data
        filename              = 'sst.ob'
        obs_read              = obs_read + 1
        obs_type   (obs_read) = 'sst'
        io_obs_unit(obs_read) = iunit
        obs_intt   (obs_read) = 10
        open(unit=iunit, file=filename, status='old')
        write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
        iunit = iunit + 1
 
        ! caonst
        if (lgcl_nst) then
            filename              = 'nst.ob'
            obs_read              = obs_read + 1
            obs_type   (obs_read) = 'nst'
            io_obs_unit(obs_read) = iunit
            obs_intt   (obs_read) = 10
            open(unit=iunit, file=filename, status='old')
            write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
            iunit = iunit + 1
        endif
 
        ! filename                   = 'sst_obs.wrk'
        ! io_obs_work_unit(obs_read) = iunit
        ! open(unit=iunit, file=filename, status='unknown')
        ! write(io_list, '(5x,a,i4,4x,a,a)') 'open unit',iunit, 'named ', filename
        ! iunit = iunit + 1
    endif
 
    ! open units for xbt data
    filename              = 'xbt.ob'
    obs_read              = obs_read + 1
    obs_type   (obs_read) = 'xbt'
    io_obs_unit(obs_read) = iunit
    obs_intt   (obs_read) = 45  ! 10, lgchen: previously 45, to slow when too many obs data
    open(unit=iunit, file=filename, status='old')
    write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
    iunit = iunit + 1
 
    filename                   = 'xbt_obs.wrk'
    io_obs_work_unit(obs_read) = iunit
    open(unit=iunit, file=filename, status='unknown')
    write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
    iunit = iunit + 1
 
    ! open units for ctd data
    if (lgcl_ctd) then
        filename              = 'ctd.ob'
        obs_read              = obs_read + 1
        obs_type   (obs_read) = 'ctd'
        io_obs_unit(obs_read) = iunit 
        obs_intt   (obs_read) = 45  ! 10, lgchen: previously 45, to slow when too many obs data
        open(unit=iunit, file=filename, status='old')
        write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
        iunit = iunit + 1
 
        filename                   = 'ctd_obs.wrk'
        io_obs_work_unit(obs_read) = iunit
        open(unit=iunit, file=filename, status='unknown')
        write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
        iunit = iunit + 1
    endif
 
    if (lgcl_altm) then
        lrec_altm = imt*jmt*4
 
        filename = 'alpha.dat'
        ioalpha  = iunit
        open(unit=iunit, file=filename, status='old', form='unformatted', access='direct', recl=imt*jmt*ksh*4)
        write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
        iunit = iunit + 1
 
      ! ! ccaoalpha
      ! filename = 'dtdzm.dat'
      ! iodtdzm = iunit
      ! open(unit=iunit, file=filename, status='old', form='unformatted')
      ! write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
      ! iunit = iunit + 1
 
      ! filename = 'dsdzm.dat'
      ! iodsdzm = iunit
      ! open(unit=iunit, file=filename, status='old', form='unformatted')
      ! write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
      ! iunit = iunit + 1
 
        if (lgcl_geosat) then
            filename    = 'geo_msfh.dat'
            io_geo_msfh = iunit
            open(unit=iunit, file=filename, status='old', form='unformatted')
            write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
            iunit = iunit + 1

            ! caonew
            filename    = 'geosat_date.dat'
            io_geo_date = iunit
            open(unit=iunit, file=filename, status='old', form='unformatted')
            write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
            iunit = iunit + 1

            filename              = 'geosat.ob'
            obs_read              = obs_read + 1
            obs_type   (obs_read) = 'geo'
            io_obs_unit(obs_read) = iunit
            obs_intt   (obs_read) = 5
            open(unit=iunit, file=filename, status='old', access='direct', form='unformatted', recl=lrec_altm)
            write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
            iunit = iunit + 1
 
            filename                   = 'geo_obs.wrk'
            io_obs_work_unit(obs_read) = iunit
            open(unit=iunit, file=filename, status='unknown')
            write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
            iunit = iunit + 1
 
            filename                   = 's_geo_obs.wrk'
            obs_read                   = obs_read + 1
            obs_type        (obs_read) = 'sgo'
            obs_intt        (obs_read) = 5
            io_obs_work_unit(obs_read) = iunit
            open(unit=iunit,file=filename,status='unknown')
            write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
            iunit = iunit + 1
        endif
 
        if (lgcl_ers1_c .or. lgcl_ers1_g) then
            if (lgcl_ers1_c) filename = 'e1c_msfh.dat'
            if (lgcl_ers1_g) filename = 'e1g_msfh.dat'
            io_es1_msfh = iunit
            open(unit=iunit, file=filename, status='old', form='unformatted')
            write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
            iunit = iunit + 1

            !caonew
            if (lgcl_ers1_c) filename = 'ers1_c_date.dat'
            if (lgcl_ers1_g) filename = 'ers1_g_date.dat'
            io_es1_date = iunit
            open(unit=iunit, file=filename, status='old', form='unformatted')
            write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
            iunit = iunit + 1


            if (lgcl_ers1_c) filename = 'ers1_c.ob'
            if (lgcl_ers1_g) filename = 'ers1_g.ob'
            obs_read              = obs_read + 1
            obs_type   (obs_read) = 'es1'
            io_obs_unit(obs_read) = iunit
            obs_intt   (obs_read) = 5
            open(unit=iunit, file=filename, status='old', access='direct', form='unformatted', recl=lrec_altm)
            write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
            iunit = iunit + 1
 
            filename                   = 'es1_obs.wrk'
            io_obs_work_unit(obs_read) = iunit
            open(unit=iunit, file=filename, status='unknown')
            write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
            iunit = iunit + 1
 
            filename                   = 's_es1_obs.wrk'
            obs_read                   = obs_read + 1
            obs_type(obs_read)         = 'se1'
            obs_intt(obs_read)         = 5
            io_obs_work_unit(obs_read) = iunit
            open(unit=iunit, file=filename, status='unknown')
            write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
            iunit = iunit + 1
        endif
 
        if (lgcl_ers2) then
            filename    = 'es2_msfh.dat'
            io_es2_msfh = iunit
            open(unit=iunit, file=filename, status='old', form='unformatted')
            write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
            iunit = iunit + 1

            ! caonew
            filename    = 'ers2_date.dat'
            io_es2_date = iunit
            open(unit=iunit, file=filename, status='old', form='unformatted')
            write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
            iunit = iunit + 1

            filename              = 'ers2.ob'
            obs_read              = obs_read + 1
            obs_type   (obs_read) = 'es2'
            io_obs_unit(obs_read) = iunit
            obs_intt   (obs_read) = 5
            open(unit=iunit, file=filename, status='old', access='direct', form='unformatted', recl=lrec_altm)
            write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
            iunit = iunit + 1
 
            filename                   = 'es2_obs.wrk'
            io_obs_work_unit(obs_read) = iunit
            open(unit=iunit, file=filename, status='unknown')
            write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
            iunit = iunit + 1
 
            filename                   = 's_es2_obs.wrk'
            obs_read                   = obs_read + 1
            obs_type        (obs_read) = 'se2'
            obs_intt        (obs_read) = 5
            io_obs_work_unit(obs_read) = iunit
            open(unit=iunit, file=filename, status='unknown')
            write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
            iunit = iunit + 1
        endif
 
        if (lgcl_jason1) then
            filename    = 'js1_msfh.dat'
            io_js1_msfh = iunit
            open(unit=iunit, file=filename, status='old', form='unformatted')
            write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
            iunit = iunit + 1

            ! caonew
            filename    = 'jason1_date.dat'
            io_js1_date = iunit
            open(unit=iunit, file=filename, status='old', form='unformatted')
            write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
            iunit = iunit + 1

            filename              = 'jason1.ob'
            obs_read              = obs_read + 1
            obs_type   (obs_read) = 'js1'
            io_obs_unit(obs_read) = iunit
            obs_intt   (obs_read) = 5
            open(unit=iunit, file=filename, status='old', access='direct', form='unformatted', recl=lrec_altm)
            write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
            iunit = iunit + 1
 
            filename                   = 'js1_obs.wrk'
            io_obs_work_unit(obs_read) = iunit
            open(unit=iunit, file=filename, status='unknown')
            write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
            iunit = iunit + 1
 
            filename                   = 's_js1_obs.wrk'
            obs_read                   = obs_read + 1
            obs_type        (obs_read) = 'sj1'
            obs_intt        (obs_read) = 5
            io_obs_work_unit(obs_read) = iunit
            open(unit=iunit, file=filename, status='unknown')
            write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
            iunit = iunit + 1
        endif
 
        if (lgcl_topex) then
            filename    = 'tpx_msfh.dat'
            io_tpx_msfh = iunit
            open(unit=iunit,file=filename,status='old', form='unformatted')
            write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
            iunit = iunit + 1

            ! caonew
            filename    = 'topex_date.dat'
            io_tpx_date = iunit
            open(unit=iunit, file=filename, status='old', form='unformatted')
            write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
            iunit = iunit + 1

            filename              = 'topex.ob'
            obs_read              = obs_read + 1
            obs_type   (obs_read) = 'tpx'
            io_obs_unit(obs_read) = iunit
            obs_intt   (obs_read) = 5
            open(unit=iunit, file=filename, status='old', access='direct', form='unformatted', recl=lrec_altm)
            write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
            iunit = iunit + 1
 
            filename                   = 'tpx_obs.wrk'
            io_obs_work_unit(obs_read) = iunit
            open(unit=iunit, file=filename, status='unknown')
            write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
            iunit = iunit + 1
 
            filename                   = 's_tpx_obs.wrk'
            obs_read                   = obs_read + 1
            obs_type        (obs_read) = 'stp'
            obs_intt        (obs_read) = 5
            io_obs_work_unit(obs_read) = iunit
            open(unit=iunit, file=filename, status='unknown')
            write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
            iunit = iunit + 1
        endif
    endif
 
    ! open files for bias
    if (lgcl_bias) then
        filename     = 'temp_bias.dat'
        io_temp_bias = iunit
        open(unit=iunit, file=filename, status='unknown', access='direct', form='unformatted', recl=lrec_sgl)
        write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
        iunit = iunit + 1
 
        filename     = 'salt_bias.dat'
        io_salt_bias = iunit
        open(unit=iunit, file=filename, status='unknown', access='direct', form='unformatted', recl=lrec_sgl)
        write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
        iunit = iunit + 1

        filename = 'G5_t.dat'
        io_g5_t  = iunit
        open(unit=iunit, file=filename, status='old', access='direct', form='unformatted', recl=nx*ny*4)
        write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
        iunit = iunit + 1

      ! filename = 'G5_s.dat'
      ! io_g5_s = iunit
      ! open(unit=iunit, file=filename, status='old', access='direct', form='unformatted', recl=nx*ny*4)
      ! write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
      ! iunit = iunit + 1

        filename = 'G1_t.dat'
        io_g1_t  = iunit
        open(unit=iunit, file=filename, status='old', access='direct', form='unformatted', recl=lrec_sgl)
        write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
        iunit = iunit + 1

      ! filename = 'G1_s.dat'
      ! io_g1_s = iunit
      ! open(unit=iunit, file=filename, status='old', access='direct', form='unformatted', recl=lrec_sgl)
      ! write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
      ! iunit = iunit + 1
    endif
 
    ! caoalpha
    ! open file for bias & altimetry assim
    ! if (lgcl_bias .or. lgcl_altm .or. lgkl_ctd) then
        filename = 'dsdt.dat'
        io_dsdt  = iunit
        open(unit=iunit, file=filename, status='old', access='direct', form='unformatted', convert='little_endian', recl=lrec_sgl)
        write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
        iunit = iunit + 1
    ! endif

    ! open file for vertical correlation table
    if (lgcl_vert) then
       filename = 'pop40_ver_tbl.txt'
       io_vert  = iunit
       open(unit=iunit, file=filename, status='old')
       write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
       iunit = iunit + 1
    endif
 
    ! open files for levitus temp and salt
    filename        = 'levitus_temp.dat'
    io_levitus_temp = iunit
   ! open(unit=iunit, file=filename, status='old', access='direct', form='unformatted', recl=lrec_sgl)
    open(unit=iunit, file=filename, status='old', access='direct', form='unformatted', convert='little_endian', recl=lrec_sgl)
    write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
    iunit = iunit + 1
 
    filename        = 'levitus_salt.dat'
    io_levitus_salt = iunit
   ! open(unit=iunit, file=filename, status='old', access='direct', form='unformatted', recl=lrec_sgl)
    open(unit=iunit, file=filename, status='old', access='direct', form='unformatted', convert='little_endian', recl=lrec_sgl)
    write(io_list, '(5x,a,i4,4x,a,a)') 'open unit', iunit, 'named ', filename
    iunit = iunit + 1
 
    ! close files.dat 
    close(unit=io_list)

    ! caonst
    if (myrank == 0) then
        print *, ' total obs data ', obs_read
        print *, '  '
    end if

    return
end subroutine open_files
 
 
SUBROUTINE CLOSE_FILES
    use io_units
    use mpi
 
    io_status = 0
    read(io_list, '(a)') filename
 
    do while (io_status .eq. 0)
        read(io_list,'(a)',iostat=io_status) filename
      ! close(filename)
    end do
 
    if (myrank == 0) then
        write (stdout, '(/,12x,a,/)') ' closed files for assimilation'
    end if
 
    return
end subroutine close_files
 
 
subroutine find_obs_index(obs)
    use io_units
    use obs_types

    use mpi
 
    character(len=3) :: obs

    do kr = 1, obs_read
        if (obs .eq. obs_type(kr)) then
            obs_index = kr

            if (myrank == 0) then
                write(6, *) 'find_obs_index: obs_index=', kr  &
                    , ', obs_read=', obs_read                &
                    , ', obs=', obs                          &
                    , ', obs_type=', obs_type(kr)
            end if

            return
        end if
    end do
 
    write(stdout, '(2x,'' observation type not found'')')
    stop 
 
    return
end subroutine find_obs_index
