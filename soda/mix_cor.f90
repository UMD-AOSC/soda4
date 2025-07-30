subroutine mixlayer_depth(nid)
!
!=======================================================================
!
!    This subroutine calculate mixed layer depth from temperature and
!    salinity first guess fields based on the ensity criteria and
!    save it in the mldk array to be used by the mixlayer_cor subroutine to
!    correct temperature and salinity in the mixed layer
!
!=======================================================================
! 
    use model_grid
    use mpi
    
    parameter(diff=0.05)
 
    double precision t1, s1, p1, tpot, rou0
    real, dimension(km) :: rou
 
    mldk = 0
    do j = 1, jmt; do i = 1, imt
      if (kmt(i, j) .gt. 1) then
        do k = 1, kmt(i, j)
          t1 = temp_first_guess(i, j, k)
          s1 = salt_first_guess(i, j, k)
          p1 = zt(k)/10.0
          if (t1 .gt. -2.) then
            call potem(t1, s1, p1, tpot)
            call unesco(tpot, s1, p1, rou0)
            rou0   = rou0 - 1.0d3 + 2.5d-2
            rou(k) = rou0
          endif
        enddo
        kd = 0
        do k = 2, kmt(i, j)
          delta_dens = rou(k) - rou(1)
          if (delta_dens .gt. diff) then 
            kd = k - 1
            exit
          endif
        enddo	  
        mldk(i,j) = kd
      endif
    enddo; enddo
    
    return
end subroutine mixlayer_depth

	
subroutine mixlayer_cor(nid)
!
!=======================================================================
!
!    This subroutine correct temperature and salinity in the mixed
!    layer calculated by mixlayer_depth subroutine and saved in the
!    mldk array
!
!=======================================================================
! 
    use model_grid
    use mpi
    
    double precision sss
    parameter(alpha=0.0)
	  
    do j = 1, jmt; do i = 1, imt
      if (sst_obs(i, j).ne.temp_first_guess(i, j, 1)) then
        if(mldk(i,j).ne.0) then
!	   write(*,*) 'i=',i,' j=',j,' mldk=',mldk(i,j)
	  sss = 0.0
	  do kk = 1, mldk(i,j)
	    sss = sss + csalt_first_guess(i, j, kk)
	  enddo
	  sss = sss/mldk(i,j)
	  do k = 1, mldk(i,j)
            ctemp_first_guess(i, j, k) = sst_obs(i, j)*(1.0 - alpha) + temp_first_guess(i, j, k)*alpha
            csalt_first_guess(i, j, k) = sss*(1.0 - alpha) + salt_first_guess(i, j, k)*alpha 
!            csalt_first_guess(i, j, k) = salt_first_guess(i, j, k) 
	  enddo
	endif
      endif
    enddo; enddo
    
    return
end subroutine mixlayer_cor


