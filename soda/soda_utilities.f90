subroutine make_kmt
    use model_grid
 
    do j = 1, jmt; do i = 1, imt;
        do k = 1, km 
            if (temp_first_guess(i, j, k) .eq. 0.) exit
        end do
        kmt(i, j) = k - 1
    end do; end do

    return
end subroutine make_kmt

 
 
 
subroutine smth5pt(temp, lev)
!=======================================================================
!
! Xianhe Cao  4/18/96
!
! this subroutine does 5 point smoothing for a 2_D temperature field
!
!=======================================================================
 
    use model_grid
 
    real, dimension(imt, jmt) :: temp, tmps
 
    tmps = 0.0
 
    do i = 2, imt - 1; do j = 2, jmt - 1;
        sum = 0.0
        npt = 0
 
        if (lev .le. kmt(i, j)) then
            npt = npt + 1
            sum = sum + temp(i, j)
 
            if (lev .le. kmt(i-1, j)) then
                npt = npt + 1
                sum = sum + temp(i-1, j)
            endif
            if (lev .le. kmt(i+1, j)) then
                npt = npt + 1
                sum = sum + temp(i+1, j)
            endif
            if (lev .le. kmt(i, j-1)) then
                npt = npt + 1
                sum = sum + temp(i, j-1)
            endif
            if (lev .le. kmt(i, j+1)) then
                npt = npt + 1
                sum = sum + temp(i, j+1)
            endif
        end if
 
        ! write(6, *) i, j, npts
        if (npt .eq. 0) then
            tmps(i, j) = 0.
        else
            tmps(i, j) = sum/float(npt)
        end if
    end do; end do;
 
    do j = 2, jmt - 1
        tmps(1  , j) = tmps(imt-1, j)
        tmps(imt, j) = tmps(2    , j)
    end do
 
    tmps(:, 1  ) = tmps(:, 2    )
    tmps(:, jmt) = tmps(:, jmt-1)
 
    temp(:, :) = tmps(:, :)
 
    return
end subroutine smth5pt

 
 
 
function jday(mon, iday, iyr)
!=======================================================================
!
!     compute the julian day corresponding to the
!     day on the gregorian calender
!
!=======================================================================
 
    dimension dpm(12)
    data dpm /31.0, 28.0, 31.0, 30.0, 31.0, 30.0, 31.0, 31.0, 30.0, 31.0, 30.0, 31.0/

    dpm(2) = 28.0
    if (mod(real(iyr), 4.) .eq. 0.) dpm(2) = 29.0
 
    ! first calculate days without leap years
 
    ! cao??11/8/95  for 1950 XXX
    ! cao    iyrs = iyr-1950
    ! cao    days = 3282.0+real(iyrs)*365.0
 
    iyrs = iyr - 1970
    days = 587.0 + real(iyrs)*365.0
    ! cao?? need to think about iyrs+?? e.x.: iyr-1950, should iyrs+1
    ! caojday num_leap = int(real(iyrs+1)/4.)
    num_leap = floor(real(iyrs+1)/4.)
    days     = days + real(num_leap)
 
    ! now sum up for days this year
    sum = 0.
    if (mon .gt. 1)then
        do l = 1, mon - 1
            sum = sum + dpm(l)
        end do
        days = days + sum
    end if
 
    jday = int(days) + iday
 
    return
end function jday

 
 
 
subroutine atlpac(xlat, xlon, lxap)
!=================================================================
!
! X. Cao 12/9/99
!
!   to make a mark to the location of a point in Caribbean area
! (xlat.gt.-2..and.xlat.lt.32..and.xlon.gt.245..and.xlon.lt.295.)
! to indicate if it is in Atlantic ocean (1) or in Pacific ocean (2)
! or out of Caribbean area (0)
!
!=================================================================
 
    lxap = 0

    ! -- Atlantic ? Pacific?
    if (xlat .gt. -2. .and. xlat .le. 8.5) then
        if (xlon .lt. 285.) then
            lxap = 2
        else
            lxap = 1
        endif
    endif
 
    if (xlat .gt. 8.5 .and. xlat .le. 15.5) then
        if (xlon .lt. 276.) then
            lxap = 2
        else
            lxap = 1
        endif
    endif
 
    if (xlat .gt. 15.5 .and. xlat .le. 19.5) then
        if (xlon .lt. 270.) then
            lxap = 2
        else
            lxap = 1
        endif
    endif
 
    if (xlat .gt. 19.5 .and. xlat .le. 32.0) then
        if (xlon .lt. 258.) then
            lxap = 2
        else
            lxap = 1
        endif
    endif
 
    return
end subroutine atlpac

 
 
 
function indp(value, array, ia)
!=======================================================================
!
!     indp = index of nearest data point within "array" corresponding to
!            "value".
!
!     inputs:
!
!     value  = arbitrary data...same units as elements in "array"
!     array  = array of data points  (must be monotonically increasing)
!     ia     = dimension of "array"
!
!     output:
!
!     indp =  index of nearest data point to "value"
!             if "value" is outside the domain of "array" then indp = 1
!             or "ia" depending on whether array(1) or array(ia) is
!             closest to "value"
!
!             note: if "array" is dimensioned array(0:ia) in the calling
!                   program, then the returned index should be reduced
!                   by one to account for the zero base.
!
!     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
!
!     example:
!
!     let model depths be defined by the following:
!     parameter (km=5)
!     dimension z(km)
!     data z /5.0, 10.0, 50.0, 100.0, 250.0/
!
!     k1 = indp (12.5, z, km)
!     k2 = indp (0.0, z, km)
!
!     k1 would be set to 2, & k2 would be set to 1 so that
!     z(k1) would be the nearest data point to 12.5 and z(k2) would
!     be the nearest data point to 0.0
!
!=======================================================================
 
    use io_units
    dimension array(ia)
 
    do i = 2, ia
        if (array(i) .lt. array(i-1)) then
            write (stdout,*)  &
                ' => Error: array must be monotonically increasing in "indp"',  &
                '           when searching for nearest element to value=', value

            write (stdout,*) '           array(i) < array(i-1) for i=', i
            write (stdout,*) '           array(i) for i=1..ia follows:'
            do ii = 1, ia
                write(stdout, *) 'i=', ii, ' array(i)=', array(ii)
            enddo
            stop '=>indp'
        endif
    enddo

    if (value .lt. array(1) .or. value .gt. array(ia)) then
        if (value .lt. array(1 )) indp = 1
        if (value .gt. array(ia)) indp = ia
        return
    else
        do i = 2, ia
            if (value .le. array(i)) then
                indp = i
                if (array(i)-value .gt. value-array(i-1)) indp = i-1
                exit
            endif
        enddo
    endif

    return
end function indp

 
 

function great_circle(rlon1, rlat1, rlon2, rlat2)
    use constants
 
    dlon = rlon2 - rlon1
    dlat = rlat2 - rlat1
    a    = (sin(dlat/2))**2 + cos(rlat1)*cos(rlat2)*(sin(dlon/2))**2
  ! c = 2. * arcsin(min(1.,sqrt(a)))
    c = 1.
    great_circle = radius*c
 
    return
end function great_circle
