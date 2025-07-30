subroutine dchco(a, lda, n, ipvt, rcond, z)
 
!=======================================================================
!
!     this routine is from netlib, by 12/2/1991
!
!     dgeco factors a double precision matrix by gaussian elimination
!     and estimates the condition of the matrix.
!
!     if  rcond  is not needed, dgefa is slightly faster.
!     to solve  a*x = b , follow dgeco by dgesl.
!     to compute  inverse(a)*c , follow dgeco by dgesl.
!     to compute  determinant(a) , follow dgeco by dgedi.
!     to compute  inverse(a) , follow dgeco by dgedi.
!
!     inputs:   
!
!        a   = double precision(lda, n)
!              the matrix to be factored.
!
!        lda = integer
!              the leading dimension of the array  a .
!
!        n   = integer
!              the order of the matrix  a .
!
!     outputs:   
!
!        a     = an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                the factorization can be written  a = l*u  where
!                l  is a product of permutation and unit lower
!                triangular matrices and  u  is upper triangular.
!
!        ipvt  = integer(n)
!                an integer vector of pivot indices.
!
!        rcond = double precision
!                an estimate of the reciprocal condition of  a .
!                for the system  a*x = b , relative perturbations
!                in  a  and  b  of size  epsilon  may cause
!                relative perturbations in  x  of size  epsilon/rcond .
!                if  rcond  is so small that the logical expression
!                           1.0 + rcond .eq. 1.0
!                is true, then  a  may be singular to working
!                precision.  in particular,  rcond  is zero  if
!                exact singularity is detected or the estimate
!                underflows.
!
!        z     = double precision(n)
!                a work vector whose contents are usually unimportant.
!                if  a  is close to a singular matrix, then  z  is
!                an approximate null vector in the sense that
!                norm(a*z) = rcond*norm(a)*norm(z) .
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     linpack dgefa
!     blas daxpy,ddot,dscal,dasum
!     fortran dabs,dmax1,dsign
!
!=======================================================================
 
    integer lda, n, ipvt(1)
    double precision a(lda, 1), z(1)
    double precision rcond
 
    ! internal variables
 
    double precision ddot, ek, t, wk, wkm
    double precision anorm, s, dasum, sm, ynorm
    integer info, j, k, kb, kp1, l
 
#ifdef timing
    call tic('toa', 'dchco')
#endif
 
    !----------------------------------------------------------------------
    !     compute 1-norm of a
    !----------------------------------------------------------------------
 
    anorm = 0.0d0
    do j = 1, n
        anorm = dmax1(anorm, dasum(n, a(1, j), 1))
    end do
 
    !----------------------------------------------------------------------
    !     factor
    !----------------------------------------------------------------------
 
    call dgefa(a, lda, n, ipvt, info)
  
    !----------------------------------------------------------------------
    !     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
    !     estimate = norm(z)/norm(y) where  a*z = y  and  trans(a)*y = e .
    !     trans(a)  is the transpose of a .  the components of  e  are
    !     chosen to cause maximum local growth in the elements of w  where
    !     trans(u)*w = e .  the vectors are frequently rescaled to avoid
    !     overflow.
    !
    !     solve trans(u)*w = e
    !----------------------------------------------------------------------
 
    ek = 1.0d0
    do j = 1, n
         z(j) = 0.0d0
    end do

    do k = 1, n
        if (z(k) .ne. 0.0d0) ek = dsign(ek, -z(k))

        if (.not. (dabs(ek - z(k)) .le. dabs(a(k, k)))) then
            s = dabs(a(k, k))/dabs(ek - z(k))
            call dscal(n, s, z, 1)
            ek = s*ek
        end if

        wk  = ek - z(k)
        wkm = -ek - z(k)
        s   = dabs(wk )
        sm  = dabs(wkm)

        if (a(k, k) .eq. 0.0d0) then
            wk  = 1.0d0
            wkm = 1.0d0
        else
            wk  = wk  / a(k, k)
            wkm = wkm / a(k, k)
        end if

        kp1 = k + 1
        if (.not. (kp1 .gt. n)) then  ! go to 90
            do j = kp1, n
                sm   = sm + dabs(z(j) + wkm*a(k, j))
                z(j) = z(j) + wk*a(k, j)
                s    = s + dabs(z(j))
            end do

            if (.not. (s .ge. sm)) then
                t  = wkm - wk
                wk = wkm
                do j = kp1, n
                    z(j) = z(j) + t*a(k, j)
                end do
            end if
        end if

        z(k) = wk
    end do  ! end of "do k = 1, n"


    s = 1.0d0 / dasum(n, z, 1)
    call dscal(n, s, z, 1)
 
    !----------------------------------------------------------------------
    !     solve trans(l)*y = w
    !----------------------------------------------------------------------
 
    do kb = 1, n
        k = n + 1 - kb
        if (k .lt. n) z(k) = z(k) + ddot(n-k, a(k+1, k), 1, z(k+1), 1)
        if (.not. (dabs(z(k)) .le. 1.0d0)) then
            s = 1.0d0 / dabs(z(k))
            call dscal(n, s, z, 1)
        end if

        l    = ipvt(k)
        t    = z(l)
        z(l) = z(k)
        z(k) = t
    end do

    s = 1.0d0 / dasum(n, z, 1)
    call dscal(n, s, z, 1)
 
    ynorm = 1.0d0
 
    !----------------------------------------------------------------------
    !     solve l*v = y
    !----------------------------------------------------------------------
 
    do k = 1, n
        l    = ipvt(k)
        t    = z(l)
        z(l) = z(k)
        z(k) = t
        if (k .lt. n) call daxpy(n-k, t, a(k+1, k), 1, z(k+1), 1)

        if (dabs(z(k)) .le. 1.0d0) cycle

        s = 1.0d0 / dabs(z(k))
        call dscal(n, s, z, 1)
        ynorm = s*ynorm
    end do

    s = 1.0d0 / dasum(n, z, 1)
    call dscal(n, s, z, 1)
    ynorm = s*ynorm
 
    !----------------------------------------------------------------------
    !     solve  u*z = v
    !----------------------------------------------------------------------
 
    do kb = 1, n
        k = n + 1 - kb
        if (.not. (dabs(z(k)) .le. dabs(a(k, k)))) then
            s = dabs(a(k, k))/dabs(z(k))
            call dscal(n, s, z, 1)
            ynorm = s*ynorm
        end if

        if (a(k, k) .ne. 0.0d0) z(k) = z(k) / a(k, k)
        if (a(k, k) .eq. 0.0d0) z(k) = 1.0d0
        t = -z(k)
        call daxpy(k-1, t, a(1, k), 1, z(1), 1)
    end do
 
    !----------------------------------------------------------------------
    !     make znorm = 1.0
    !----------------------------------------------------------------------
 
    s = 1.0d0 / dasum(n, z, 1)
    call dscal(n, s, z, 1)
    ynorm = s*ynorm
 
    if (anorm .ne. 0.0d0) rcond = ynorm / anorm
    if (anorm .eq. 0.0d0) rcond = 0.0d0
 
#ifdef timing
    call toc ('toa', 'dchco')
#endif
 
    return
end subroutine dchco  ! dchco(a, lda, n, ipvt, rcond, z)




subroutine dgefa(a, lda, n, ipvt, info)
 
!=======================================================================
!
!     dgefa factors a double precision matrix by gaussian elimination.
!
!     dgefa is usually called by dgeco, but it can be called
!     directly with a saving in time if  rcond  is not needed.
!     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
!
!     inputs:  
!
!        a   = double precision(lda, n)
!              the matrix to be factored.
!
!        lda = integer
!              the leading dimension of the array  a .
!
!        n   = integer
!              the order of the matrix  a .
!
!     outputs:  
!
!        a    = an upper triangular matrix and the multipliers
!               which were used to obtain it.
!               the factorization can be written  a = l*u  where
!               l  is a product of permutation and unit lower
!               triangular matrices and  u  is upper triangular.
!
!        ipvt = integer(n)
!               an integer vector of pivot indices.
!
!        info = integer
!                = 0  normal value.
!                = k  if  u(k,k) .eq. 0.0 .  this is not an error
!                     condition for this subroutine, but it does
!                     indicate that dgesl or dgedi will divide by zero
!                     if called.  use  rcond  in dgeco for a reliable
!                     indication of singularity.
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     blas daxpy,dscal,idamax
!
!=======================================================================
 
    integer lda, n, ipvt(1), info
    double precision a(lda, 1)
 
    ! internal variables
 
    double precision t
    integer idamax, j, k, kp1, l, nm1
 
    !----------------------------------------------------------------------
    !     gaussian elimination with partial pivoting
    !----------------------------------------------------------------------
 
    info = 0
    nm1  = n - 1
    if (.not. (nm1 .lt. 1)) then
        do k = 1, nm1
            kp1 = k + 1
 
            !----------------------------------------------------------------------
            !     find l = pivot index
            !----------------------------------------------------------------------
 
            l       = idamax(n-k+1, a(k, k), 1) + k - 1
            ipvt(k) = l
 
            !----------------------------------------------------------------------
            !     zero pivot implies this column already triangularized
            !----------------------------------------------------------------------
 
            if (.not. (a(l, k) .eq. 0.0d0)) then
                ! interchange if necessary
                if (.not. (l .eq. k)) then
                    t       = a(l, k)
                    a(l, k) = a(k, k)
                    a(k, k) = t
                end if
 
                ! compute multipliers
 
                t = -1.0d0 / a(k, k)
                call dscal(n-k, t, a(k+1, k), 1)
 
                ! row elimination with column indexing
                do j = kp1, n
                    t = a(l, j)
                    if (.not. (l .eq. k)) then
                        a(l, j) = a(k, j)
                        a(k, j) = t
                    end if
                    call daxpy(n-k, t, a(k+1, k), 1, a(k+1, j), 1)
                end do

                cycle
            end if

            info = k
        end do  ! end of "do k = 1, nm1"
    end if

    ipvt(n) = n
    if (a(n, n) .eq. 0.0d0) info = n

    return
end subroutine dgefa  ! dgefa(a, lda, n, ipvt, info)




integer function idamax(n, dx, incx)
!=======================================================================
!     finds the index of element having max. absolute value.
!     jack dongarra, linpack, 3/11/78.
!=======================================================================
 
    double precision dx(1), dmax
    integer i, incx, ix, n
 
    idamax = 0
    if(n .lt. 1) return

    idamax = 1
    if(n .eq. 1) return

    if(.not. (incx .eq. 1)) then
        ! code for increment not equal to 1
        ix   = 1
        dmax = dabs(dx(1))
        ix   = ix + incx
        do i = 2, n
            if(.not. (dabs(dx(ix)) .le. dmax)) then
                idamax = i
                dmax   = dabs(dx(ix))
            end if
            ix = ix + incx
        end do

        return
    end if

    ! code for increment equal to 1
    dmax = dabs(dx(1))

    do i = 2, n
        if(dabs(dx(i)) .le. dmax) cycle
        idamax = i
        dmax   = dabs(dx(i))
    end do

    return
end function idamax  ! idamax(n, dx, incx)




double precision function dasum(n, dx, incx)
!=======================================================================
!     takes the sum of the absolute values.
!     jack dongarra, linpack, 3/11/78.
!=======================================================================
 
    double precision dx(1), dtemp
    integer i, incx, m, mp1, n, nincx
 
    dasum = 0.0d0
    dtemp = 0.0d0

    if(n .le. 0) return

    if(.not. (incx .eq. 1)) then
        ! code for increment not equal to 1
        nincx = n*incx
        do i = 1, nincx, incx
            dtemp = dtemp + dabs(dx(i))
        end do

        dasum = dtemp

        return
    end if

    ! code for increment equal to 1
    ! clean-up loop
    m = mod(n, 6)

    if (.not. (m .eq. 0)) then
        do i = 1, m
            dtemp = dtemp + dabs(dx(i))
        end do

        if (n .lt. 6) go to 60
    end if

    mp1 = m + 1

    do i = mp1, n, 6
        dtemp = dtemp + dabs(dx(i)) + dabs(dx(i+1)) + dabs(dx(i+2)) + dabs(dx(i+3)) + dabs(dx(i+4)) + dabs(dx(i+5))
    end do

60  dasum = dtemp

    return
end function dasum  ! dasum(n, dx, incx)



    
subroutine daxpy(n, da, dx, incx, dy, incy)
!=======================================================================
! constant times a vector plus a vector.
! uses unrolled loops for increments equal to one.
! jack dongarra, linpack, 3/11/78.
!=======================================================================
 
    double precision dx(1), dy(1), da
    integer i, incx, incy, ix, iy, m, mp1, n
 
    if (n .le. 0) return
    if (da .eq. 0.0d0) return

    if (.not. (incx .eq. 1 .and. incy .eq. 1)) then
        ! code for unequal increments or equal increments not equal to 1
        ix = 1
        iy = 1
        if (incx .lt. 0) ix = (-n+1)*incx + 1
        if (incy .lt. 0) iy = (-n+1)*incy + 1

        do i = 1, n
            dy(iy) = dy(iy) + da*dx(ix)
            ix     = ix + incx
            iy     = iy + incy
        end do

        return
    end if
 
    ! code for both increments equal to 1
    ! clean-up loop
 
    m = mod(n, 4)

    if (.not. (m .eq. 0)) then
        do i = 1, m
            dy(i) = dy(i) + da*dx(i)
        end do

        if (n .lt. 4) return
    end if 

    mp1 = m + 1

    do i = mp1, n, 4
        dy(i)   = dy(i  ) + da*dx(i  )
        dy(i+1) = dy(i+1) + da*dx(i+1)
        dy(i+2) = dy(i+2) + da*dx(i+2)
        dy(i+3) = dy(i+3) + da*dx(i+3)
    end do

      return
end subroutine daxpy  ! daxpy(n, da, dx, incx, dy, incy)




double precision function ddot(n, dx, incx, dy, incy)
!=======================================================================
! forms the dot product of two vectors.
! uses unrolled loops for increments equal to one.
! jack dongarra, linpack, 3/11/78.
!=======================================================================
 
    double precision dx(1), dy(1), dtemp
    integer i, incx, incy, ix, iy, m, mp1, n
 
    ddot  = 0.0d0
    dtemp = 0.0d0

    if (n .le. 0) return

    if (.not. (incx .eq. 1 .and. incy .eq. 1)) then
        ! code for unequal increments or equal increments
        ! not equal to 1
        ix = 1
        iy = 1
        if (incx .lt. 0) ix = (-n+1)*incx + 1
        if (incy .lt. 0) iy = (-n+1)*incy + 1

        do i = 1, n
            dtemp = dtemp + dx(ix)*dy(iy)
            ix    = ix + incx
            iy    = iy + incy
        end do

        ddot = dtemp

        return
    end if

    !----------------------------------------------------------------------
    ! code for both increments equal to 1
    ! clean-up loop
    !----------------------------------------------------------------------
 
    m = mod(n, 5)

    if (.not. (m .eq. 0)) then
        do i = 1, m
            dtemp = dtemp + dx(i)*dy(i)
        end do

        if (n .lt. 5) go to 60
    end if

    mp1 = m + 1
    do i = mp1, n, 5
        dtemp = dtemp + dx(i)*dy(i) + dx(i+1)*dy(i+1) + dx(i+2)*dy(i+2) + dx(i+3)*dy(i+3) + dx(i+4)*dy(i+4)
    end do

60  ddot = dtemp

    return
end function ddot  ! ddot(n, dx, incx, dy, incy)




subroutine dchdi(a, lda, n, ipvt, det, work, job)
!=======================================================================
!
!     this routine is from netlib, by 12/2/1991
!
!     dgedi computes the determinant and inverse of a matrix
!     using the factors computed by dgeco or dgefa.
!
!     inputs:  
!
!        a    = double precision(lda, n)
!               the output from dgeco or dgefa.
!
!        lda  = integer
!               the leading dimension of the array  a .
!
!        n    = integer
!               the order of the matrix  a .
!
!        ipvt = integer(n)
!               the pivot vector from dgeco or dgefa.
!
!        work = double precision(n)
!               work vector.  contents destroyed.
!
!        job  = integer
!                = 11   both determinant and inverse.
!                = 01   inverse only.
!                = 10   determinant only.
!
!     outputs:  
!
!        a   = inverse of original matrix if requested.
!              otherwise unchanged.
!
!        det = double precision(2)
!              determinant of original matrix if requested.
!              otherwise not referenced.
!              determinant = det(1) * 10.0**det(2)
!              with  1.0 .le. dabs(det(1)) .lt. 10.0
!              or  det(1) .eq. 0.0 .
!
!     error condition:
!
!        a division by zero will occur if the input factor contains
!        a zero on the diagonal and the inverse is requested.
!        it will not occur if the subroutines are called correctly
!        and if dgeco has set rcond .gt. 0.0 or dgefa has set
!        info .eq. 0 .
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     blas daxpy,dscal,dswap
!     fortran dabs,mod
!
!=======================================================================
 
    integer lda, n, ipvt(1), job
    double precision a(lda, 1), det(2), work(1)
 
    ! internal variables
    double precision t
    double precision ten
    integer i, j, k, kb, kp1, l, nm1
 
#ifdef timing
    call tic('toa', 'dchdi')
#endif
 
    !----------------------------------------------------------------------
    ! compute determinant
    !----------------------------------------------------------------------
 
    if (.not. (job/10 .eq. 0)) then
        det(1) = 1.0d0
        det(2) = 0.0d0
        ten    = 10.0d0

        do i = 1, n
            if (ipvt(i) .ne. i) det(1) = -det(1)
            det(1) = a(i, i)*det(1)
          ! ...exit
            if (det(1) .eq. 0.0d0) exit 

            do while (.not. (dabs(det(1)) .ge. 1.0d0))
                det(1) = ten*det(1)
                det(2) = det(2) - 1.0d0
            end do

            do while (.not. (dabs(det(1)) .lt. ten))
                det(1) = det(1)/ten
                det(2) = det(2) + 1.0d0
            end do 
        end do
    end if
 
    !----------------------------------------------------------------------
    ! compute inverse(u)
    !----------------------------------------------------------------------
 
    if (.not. (mod(job, 10) .eq. 0)) then
        do k = 1, n
            a(k, k) = 1.0d0/a(k, k)
            t       = -a(k, k)
            call dscal(k-1, t, a(1, k), 1)
            kp1 = k + 1
            if (.not. (n .lt. kp1)) then
                do j = kp1, n
                    t      = a(k, j)
                    a(k,j) = 0.0d0
                    call daxpy(k, t, a(1, k), 1, a(1, j), 1)
                end do
            end if
        end do
 
        !----------------------------------------------------------------------
        ! form inverse(u)*inverse(l)
        !----------------------------------------------------------------------
        nm1 = n - 1
        if (.not. (nm1 .lt. 1)) then
            do kb = 1, nm1
                k   = n - kb
                kp1 = k + 1

                do i = kp1, n
                   work(i) = a(i, k)
                   a(i, k) = 0.0d0
                end do

                do j = kp1, n
                   t = work(j)
                   call daxpy(n, t, a(1, j), 1, a(1, k), 1)
                end do

                l = ipvt(k)
                if (l .ne. k) call dswap(n, a(1, k), 1, a(1, l), 1)
            end do
        end if
    end if
 
#ifdef timing
    call toc ('toa', 'dchdi')
#endif
 
    return
end subroutine dchdi  ! dchdi(a, lda, n, ipvt, det, work, job)




subroutine dscal(n, da, dx, incx)
!=======================================================================
! scales a vector by a constant.
! uses unrolled loops for increment equal to one.
! jack dongarra, linpack, 3/11/78.
!=======================================================================
 
    double precision da, dx(1)
    integer i, incx, m, mp1, n, nincx
 
    if (n .le. 0) return

    if (.not. (incx .eq. 1)) then
        ! code for increment not equal to 1
        nincx = n*incx

        do i = 1, nincx, incx
            dx(i) = da*dx(i)
        end do

        return
    end if 

    ! code for increment equal to 1
    ! clean-up loop
    m = mod(n, 5)

    if (.not. (m .eq. 0)) then
        do i = 1, m
            dx(i) = da*dx(i)
        end do

        if (n .lt. 5) return
    end if

    mp1 = m + 1

    do i = mp1, n, 5
        dx(i  ) = da*dx(i  )
        dx(i+1) = da*dx(i+1)
        dx(i+2) = da*dx(i+2)
        dx(i+3) = da*dx(i+3)
        dx(i+4) = da*dx(i+4)
    end do

    return
end subroutine dscal  ! dscal(n, da, dx, incx)




subroutine dswap(n, dx, incx, dy, incy)
!=======================================================================
! interchanges two vectors.
! uses unrolled loops for increments equal one.
! jack dongarra, linpack, 3/11/78.
!=======================================================================
 
    double precision dx(1), dy(1), dtemp
    integer i, incx, incy, ix, iy, m, mp1, n
 
    if (n .le. 0) return

    if (.not. (incx .eq. 1 .and. incy .eq. 1)) then
        ! code for unequal increments or equal increments not equal to 1
        ix = 1
        iy = 1
        if (incx .lt. 0) ix = (-n+1)*incx + 1
        if (incy .lt. 0) iy = (-n+1)*incy + 1

        do i = 1, n
            dtemp  = dx(ix)
            dx(ix) = dy(iy)
            dy(iy) = dtemp
            ix     = ix + incx
            iy     = iy + incy
        end do

        return
    end if

    ! code for both increments equal to 1
    ! clean-up loop
    m = mod(n, 3)

    if (.not. (m .eq. 0)) then
        do i = 1, m
            dtemp = dx(i)
            dx(i) = dy(i)
            dy(i) = dtemp
        end do

        if (n .lt. 3) return
    end if

    mp1 = m + 1

    do i = mp1, n, 3
        dtemp   = dx(i)
        dx(i)   = dy(i)
        dy(i)   = dtemp

        dtemp   = dx(i+1)
        dx(i+1) = dy(i+1)
        dy(i+1) = dtemp

        dtemp   = dx(i+2)
        dx(i+2) = dy(i+2)
        dy(i+2) = dtemp
    end do

    return
end subroutine dswap
