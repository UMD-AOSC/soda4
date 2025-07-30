subroutine psd_slt(i, j, tmp, slt, inversion)
!===================================================================
!
! Xianhe Cao 11/25/97
!
! from temp to salt by climatology t-s relationship
!
!  tmp = temperature, input
!  slt = salinity, output
!  inversion = temp inversion (1:yes, o:no; if yes, no salt calculated)
!
!===================================================================
 
    use model_grid
 
    ! - check for temp inversion (if yes: no salt updating)
 
    do k = 1, kmb
        if (levitus_temp(i, j, k+1) .gt. levitus_temp(i, j, k)) then 
            inversion = 1
            return
        endif
    enddo
 
    call saltclim(i, j, tmp, slt)
 
    return
 
end subroutine psd_slt

 
 
 
subroutine saltclim(ig, jg, t, s)
!===================================================================
!
!  Gennady Chepurin  6/4/97
!
!  calculate salinity from observed (or analysed) temperature 
!  and climate mean T-S relationship at ig,jg model grid point
!
!   input:
!         t      =  observed temperature
!         ig,jg  =  model's grid coordinates
!
!   output:
!         s      =  salinity
!
!===================================================================
 
    use model_grid
 
    integer depth, depths
    real dt0, dt10, ds10
    real tcl(km), scl(km)
 
    ! monotonize climat temperature profile
    tcl(1) = levitus_temp(ig, jg, 1)
    scl(1) = levitus_salt(ig, jg, 1)
    n1     = 1

    do k = 2, km
       if (levitus_temp(ig, jg, k) .gt. tcl(n1)) cycle 

       n1      = n1 + 1
       tcl(n1) = levitus_temp(ig, jg, k)
       scl(n1) = levitus_salt(ig, jg, k)
    end do
 
    ! -- find maximum depth
    depths = n1
    do k = 1, n1
        if (tcl(k) .le. 0. .or. scl(k) .le. 0.) depths = k-1
    end do
 
    if (depths .le. 0) then
        stop 'pseudo depths < 0'
    else if (depths .eq. 1) then
        s = scl(1)
        return
    else if (depths .eq. 2) then
        s = scl(1) + (t - tcl(1)) * (scl(2) - scl(1))/(tcl(2) - tcl(1))
        return
    end if
 
    ! -- check if the temp out of climate temperature's limits
    !    When t > tcl salinity is set to surface salinity
 
    if (t .gt. tcl(1)) then 
        s = scl(1)
        return
    end if
 
    if (t .lt. tcl(depths)) then
        dt10 = tcl(depths-1) - tcl(depths)
        ds10 = scl(depths-1) - scl(depths)
        dt0  = t - tcl(depths)

        if (dt10 .ne. 0.0) then
            s = scl(depths) + dt0*ds10/dt10
        else
            s = scl(depths)
        endif

        return
    end if
 
    ! -- select model levels
    depth = 0
    do k = 1, n1
        if (t .ge. tcl(k+1) .and. t .le. tcl(k)) then
            depth = k
            exit
        end if
    end do
 
    ! -- interpolation for the intermediate temperature
    ds10 = scl(depth+1) - scl(depth)
    dt10 = tcl(depth+1) - tcl(depth)
    dt0  = t - tcl(depth)
    if (dt10 .ne. 0.0) then
        s = scl(depth) + dt0*ds10/dt10
    else
        s = scl(depth)
    endif
 
    return
end subroutine saltclim
