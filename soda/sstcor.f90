subroutine sstcor(nid)    
!=======================================================================
!
!     replace the model mixed layer temp & salt with climate sst & sss
!     the mixed layer depth is estimated each analysis time to be
!     the depth at which the temperature is more than 0.18 deg
!     different from the surface temperature.
!
!     the input value of month is the number of months of sst
!     data to skipped.
!
!     the value of month is changed immediately following input
!     to be the month from the beggining of sst data at which
!     the data assimilation will be carried out.
!
!
!     inputs:
!
!     j = latitude
!     t = first guess temperature
!     s = first guess salinity
!
!     outputs:
!
!     t2 = corrected temp: derived from a weighted average
!              of the model and climate sst's
!     s2 = corrected salt
!
!     others:
!
!     kd     = level below which temp is more than 1 deg less than
!              surface temp.
!     month  = the month since jan 83.
!     alpha  = the weight given to model sst
!     ssta   = sst
!     sssa   = sss
!     store  = the three dimensional temperature field which
!              includes sst corrections from climate.
!
!=======================================================================
 
    use model_grid
    use mpi
    
    parameter(alpha=0.25, diff=0.18)
 
    double precision t1, s1, p1, tpot, rou0
    real, dimension(km) :: rou
    real, dimension(imt, jmt) :: dummy
 
    do j = 1, jmt; do i = 1, imt
    
      if (sst_obs(i, j).ne.temp_first_guess(i, j, 1)) then
      
        if (kmt(i, j) .gt. 1) then
            do k = 1, kmt(i, j)
                t1 = temp_first_guess(i, j, k)
                s1 = salt_first_guess(i, j, k)
                p1 = zt(k)
              ! s1 = s1*1000. + 35.
 
                if (t1 .gt. -2.) then
                    call potem(t1, s1, p1, tpot)
                    call unesco(tpot, s1, p1, rou0)
 
                    rou0   = rou0 - 1.0d3 + 2.5d-2
                    rou(k) = rou0
                endif
            end do
 
            kd = 0
            do k = 2, kmt(i, j)
                delta_dens = rou(k) - rou(1)
                if (delta_dens .gt. diff) then 
                    kd = k - 1
                    exit
                endif
            end do
        
            if (kd .ne. 0) then
                do k = 1, kd
                    ctemp_first_guess(i, j, k) = sst_obs(i, j)*(1.0 - alpha) + temp_first_guess(i, j, k)*alpha

                  ! csalt_first_guess          =   ! previously this should be a big bug!
                    csalt_first_guess(i, j, k) = sss_obs(i, j)*(1.0 - alpha) + salt_first_guess(i, j, k)*alpha 
                  ! csalt_first_guess(i, j, k) = salt_first_guess(i, j, k)
                end do
            end if
        end if  ! end of "if (kmt(i, j) .gt. 1) then"
	
      endif
      
    end do; end do;
 
    if (lgcl_sst) ctemp_first_guess(:, :, 1) = sst_obs(:, :)
 
    return
end subroutine sstcor
