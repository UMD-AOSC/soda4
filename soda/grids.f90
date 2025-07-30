subroutine read_grid_z
    use model_grid
    use io_units
 
    real, dimension(km) :: delta_z
 
    do k = 1, km
        read(io_grid_z, *) delta_z(k)
    end do
 
    sum = 0.
    do k = 1, km
        zt(k) = sum + delta_z(k)/2.
        sum   = sum + delta_z(k)
    end do
 
    return
end subroutine read_grid_z

 
 
 
subroutine generate_grid_xy
    use model_grid
 
    rlon = 0. 
    do i = 1, imt
        xt(i) = rlon
        rlon  = rlon + 1.
    end do
 
    rlat = -89.5
    do j = 1, jmt
        yt(j) = rlat
        rlat  = rlat + 1.
    end do
 
    return
end subroutine generate_grid_xy
 

 
 
subroutine generate_umd_grid
    use model_grid
 
    rlon = 0.5
    do i = 1, imt
        xt(i) = rlon
        rlon  = rlon + 1.
    end do
 
    return
end subroutine generate_umd_grid 

 
 
 
subroutine getbatch                      
 
!=======================================================================
!
!     to get batches for analysis
!    
!     nbatch  = number of batches
!     xl,xr,yb,yo = boundary array of batches (degree)
!
!=======================================================================
 
    use model_grid
    use batches
    use constants
    
    use mpi
 
    ! implicit none
 
    batch       = 0
    buffer_zone = bzone
 
    ! we should probably do this the right way
    ! ALSO this only works for odd number of grid points !!!!!
 
    ! km_per_deg = 111.
 
    ! xdeg = 855./111.
    ! ydeg = 633./111.
 
    do j = 1, jmt
    do i = 1, imt
        batch = batch + 1 
 
        !----------------------------------------------------------------------
        !     each batch is lgd*lgd grids + 300km each side
        !----------------------------------------------------------------------
 
        ! First...a batch without a buffer zone
        x_center(batch, 1) = xt(i)
        y_center(batch, 1) = yt(j)
        x_west  (batch, 1) = xt(i)
        x_east  (batch, 1) = xt(i)
        y_south (batch, 1) = yt(j) 
        y_north (batch, 1) = yt(j) 
 
        i_model_west  (batch, 1) = i
        i_model_east  (batch, 1) = i
        j_model_south (batch, 1) = j
        j_model_north (batch, 1) = j
        i_model_center(batch, 1) = i
        j_model_center(batch, 1) = j
 
        ! Now do for the batch with a buffer zone
        x_center(batch, 2) = xt(i)
        y_center(batch, 2) = yt(j)
        x_west  (batch, 2) = xt(i) - buffer_zone/(km_per_deg*cos(yt(j)*deg_to_radians))
        x_east  (batch, 2) = xt(i) + buffer_zone/(km_per_deg*cos(yt(j)*deg_to_radians))
        y_south (batch, 2) = yt(j) - buffer_zone/km_per_deg 
        y_north (batch, 2) = yt(j) + buffer_zone/km_per_deg 
 
        if (x_west(batch, 2) .lt. 0.   ) x_west(batch, 2) = 360. + x_west(batch, 2)
        if (x_east(batch, 2) .gt. 360. ) x_east(batch, 2) = x_east(batch, 2) - 360.
 
        i_model_west  (batch, 2) = indp(x_west  (batch, 2), xt, imt)
        i_model_east  (batch, 2) = indp(x_east  (batch, 2), xt, imt)
        j_model_south (batch, 2) = indp(y_south (batch, 2), yt, jmt)
        j_model_north (batch, 2) = indp(y_north (batch, 2), yt, jmt)
        i_model_center(batch, 2) = indp(x_center(batch, 2), xt, imt)
        j_model_center(batch, 2) = indp(y_center(batch, 2), yt, jmt)
	
        if (myrank == 0) then 
	  write(stdout,'(3i6,10f8.2)') batch, i, j, &
	              x_west(batch,2), x_center(batch, 2), x_east(batch, 2), &
		      x_center(batch, 2)-x_east(batch, 2), yt(j), &
		      buffer_zone/(km_per_deg*cos(yt(j)*deg_to_radians))
        end if
 	
    end do
    end do


    if (myrank == 0) then 
        print *, ' '
        print *, ' total number of batches is  ', batch
        print *, ' '
    end if
    if (batch .ne. nbatch) stop 'calculated number patches /= input'
 
    return
end subroutine getbatch
