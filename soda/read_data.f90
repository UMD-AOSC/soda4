subroutine read_model_direct
    use io_units
    use model_grid
    
    use mpi
 
    real, dimension(imt, jmt) :: buffer
    real                      :: sum
    integer                   :: irec, iprec, nsum
 
    irec  = 0
    iprec = 0
 
    ! read surface height field
 
    irec = irec + 1
  ! print *, ' reading surface height field, irec = ', irec 
    read(io_model, rec=irec)  buffer
 
    ! -- cm -> m
 
!$omp parallel do
    do i = 1, imt
        do j = 1, jmt
          ! sfh_first_guess(i, j) = buffer(i, j)/100.0
            sfh_first_guess(i, j) = buffer(i, j)
        end do
    end do
!$omp end parallel do
 
    ! calculate rms
    nsum   = 0
    sum    = 0.0
    sum1   = 0.0
    rmssfh = 0.0
 
!$omp parallel do reduction(+:sum,nsum)
    do i = 1, imt
        do j = 1, jmt
            if (sfh_first_guess(i, j) .ne. 0.0) then
              nsum = nsum + 1
              sum  = sum + sfh_first_guess(i, j)
            endif
        enddo
    enddo
!$omp end parallel do
    sum = sum/real(nsum)
 
!$omp parallel do reduction(+:rmssfh)
    do i = 1, imt
        do j = 1, jmt
            if (sfh_first_guess(i, j) .ne. 0) then
                rmssfh = rmssfh + (sfh_first_guess(i, j) - sum) * (sfh_first_guess(i, j) - sum)
            end if
        enddo
    enddo
!$omp end parallel do

    rmssfh = sqrt(rmssfh/real(nsum))

    if (myrank == 0) then 
        print *, ' model_sfh & rms', sfh_first_guess(181, 65), rmssfh
    end if
 
    ! now read temp first guess
 
!$omp parallel do private(iprec,buffer)
    do k = 1, km
        iprec = irec + k
      ! print *, 'iprec = ', iprec
        read(io_model, rec=iprec)  buffer
        do i = 1, imt
            do j = 1, jmt
                temp_first_guess(i, j, k) = buffer(i, j)
            end do
        end do
    end do
!$omp end parallel do

    irec = irec + km
  ! print *, 'done reading temp first guess, irec = ', irec
 
    ! and salt

  ! print *, 'start reading salt first guess, irec = ', irec

!$omp parallel do private(iprec,buffer)
    do k = 1, km
        iprec = irec + k
        read(io_model, rec=iprec)  buffer
        do j = 1, jmt
            do i = 1, imt
              ! salt_first_guess(i, j, k) = buffer(i, j)*1000.  ! for POP2
                salt_first_guess(i, j, k) = buffer(i, j)*1.     ! for MOM5
            end do
        end do
    end do
!$omp end parallel do

    irec = irec + km
  ! print *, 'done reading salt first guess, irec = ', irec

    if (myrank == 0) then 
        print *, '  '
        print *, ' model first guess  '
      ! do k = 1, km
        do k = 1, 3
            write(6, *) k, temp_first_guess(180, 65, k), salt_first_guess(180, 65, k)
        end do
    end if
 
    return 
end subroutine read_model_direct

 
 
 
subroutine read_xbt(nid, kr)
    ! Read in xbt or ctd data and put it in file 

    use params
    use io_units
    use obs_types
    use one_d
 
    integer :: nid
 
  ! write(6, *) kr, io_obs_unit(kr), io_obs_work_unit(kr)
 
    do icnt = 1, 500000
        td(:) = 0.
 
        read(io_obs_unit(kr), 831, end=111) id, xla, xlo, n, (td(i), i=1, n)
 
        ! -- data date filter 
 
        ! xbt data are assumed to be sequential...
        !  so if the date is greater than the window...
        !   return
 
        if (id .gt. nid + obs_intt(kr)) return
 
        if (id .ge. nid - obs_intt(kr)) then
            ! -- exclude equator and Greenwich meridian
            if (xlo .ge. 0.0 .and. xlo .le. 360.0) then
                write(io_obs_work_unit(kr), 831) id, xla, xlo, n, (td(i), i=1, n)
            end if
        end if
    end do
 
! 831 format(i6, 2f9.3, i3, 22f8.3)  ! 22 is for old POP model
831   format(i6, 2f9.3, i3, 34f8.3)  ! lgchen: should be 34 
 
111 return
end subroutine read_xbt
 

 
 
subroutine read_levitus
    use model_grid
    use io_units

    use mpi
 
    real, dimension(imt, jmt) :: buffer
 
    irec_offset = (month - 1)*km
    if (myrank == 0) then
        write(stdout, '(/2x,''read levitus with offset irec = '',i4/)') irec_offset
    end if
 
    do k = 1, km
        irec = k + irec_offset
        read(io_levitus_temp, rec=irec) buffer
        levitus_temp(:, :, k) = buffer(:, :)
    end do
  ! if (myrank == 0) then
  !     write(*, *) "levitus_temp(:, 40, 40)", levitus_temp(:, 40, 40)
  ! end if
 
    do k = 1, km
        irec = k + irec_offset
        read(io_levitus_salt, rec=irec) buffer
        levitus_salt(:, :, k) = buffer(:, :)
    end do
  ! if (myrank == 0) then
  !     write(*, *) "levitus_salt(:, 40, 40)", levitus_salt(:, 40, 40)
  ! end if
 
    return
end subroutine read_levitus

 
 
 
subroutine read_cvrt
    use model_grid
    use io_units
 
    do k = 1, kmb
        read(io_vert, 111) (cvrt(l, k), l=1, kmb)
    enddo
 
! 111 format(22(f7.5, 1x))
111   format(34(f7.5, 1x))  ! lgchen: should be 34, 22 is for old POP model
 
    return
end subroutine read_cvrt
 
 

 
! ccaoalpha
subroutine read_dsdt
    use model_grid
    use io_units

    use mpi

    real, dimension(imt, jmt) :: a2
    
    irec_offset = (month - 1)*km
    if (myrank == 0) then
        write(stdout, '(/2x,''read dS/dT with offset irec = '',i4/)') irec_offset
    end if

    do k = 1, km
        irec = k + irec_offset
        read(io_dsdt, rec=irec) a2
        dsdt(:, :, k) = a2(:, :)
    enddo

    return
end subroutine read_dsdt
 
 
 

! caoggbias
subroutine read_bias
    use bias

    ! caosbias
    real, dimension(nx , ny ) :: a5_t
    real, dimension(imt, jmt) :: a1_t
  ! real, dimension(nx , ny ) :: a5_s
  ! real, dimension(imt, jmt) :: a1_s

    irec = 1

    do i_p = 1, np
        do i_z = 1, kmc
            read(io_g5_t, rec=irec) a5_t
          ! read(io_g5_s, rec=irec) a5_s

            do i_y = 1, ny
                do i_x = 1, nx
                    G5_t(i_x, i_y, i_z, i_p) = -a5_t(i_x, i_y)
                  ! G5_s(i_x, i_y, i_z, i_p) = -a5_s(i_x, i_y)/1000.0
                enddo
            enddo
   
            read(io_g1_t, rec=irec) a1_t
          ! read(io_g1_s, rec=irec) a1_s

            do i_y = 1, jmt
                do i_x = 1, imt
                    G1_t(i_x, i_y, i_z, i_p) = -a1_t(i_x, i_y)
                  ! G1_s(i_x, i_y, i_z, i_p) = -a1_s(i_x, i_y)/1000.0
                enddo
            enddo

            irec = irec + 1
        enddo
    enddo

    return
end subroutine read_bias          
