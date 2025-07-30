module mpi
!==============================================================================
! PURPOSE: General MPI procedures
!
! HISTORY: 
!     10/01/2015, Ligang Chen, Created
!==============================================================================
    use params

    implicit none

    include 'mpif.h'
  ! public
  ! public initialize_mpi, finalize_mpi

    integer, save :: nprocs, nprocs_x, nprocs_y
    integer, save :: myrank, myrank_x, myrank_y
    integer, save :: is, ie, js, je

    integer, parameter :: mpibufsize=1000 !600 !(this worked as 'safe' with 480 procs on Gaea) 
                                          !200 !1000  !STEVE: this fixes the problem of bad output when using over 6 nodes default=1000,mom2(mpich2)=200
    integer, save :: nij1     ! STEVE: this is the number of gridpoints to run on this (myrank) processor
    integer, save :: nij1max  ! STEVE: the largest number of gridpoints on any 1 processor
    integer, allocatable, save :: nij1node(:)
    integer, save :: bs, be

    integer, save :: ierr
    
    integer, save :: MPI_r_size

contains
subroutine initialize_mpi
  ! use mpi

    implicit none
    integer :: i, n

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
  ! write(6, '(A,I3.3,A,I3.3)') 'Hello from MYRANK ', myrank, '/', nprocs - 1

!   if (nprocs == 512) then
!       nprocs_y = 16
!       nprocs_x = 32
!   else if (nprocs == 768) then
!       nprocs_y = 32
!       nprocs_x = 48
!   else if (nprocs == 1024) then
!       nprocs_y = 24
!       nprocs_x = 32       
!   end if

!   myrank_x = mod(myrank, nprocs_x)
!   myrank_y = myrank/nprocs_xa

!   i_mod = mod(360, nprocs_x)
!   nx_max = (360 - i_mod) / nprocs_x + 1
!   if (myrank_x < i_mod) then
!       is = myrank_x*nprocs_x + 1
!       ie = is + nx_max - 1
!   else
!       is = i_mod*nx_max + (myrank_x-i_mod)*(nx_max-1) + 1
!       ie = is + (nx_max - 1) - 1
!   end if 

!   j_mod = mod(180, nprocs_y)
!   ny_max = (180 - j_mod) / nprocs_y + 1
!   if (myrank_y < j_mod) then
!       js = myrank_y*nprocs_y + 1
!       je = js + ny_max - 1
!   else
!       js = j_mod*ny_max + (myrank_y-j_mod)*(ny_max-1) + 1
!       je = js + (ny_max - 1) - 1
!   end if 
!  
!   write(6, *) 'nprocs = ', nprocs, ', nprocs_y = ', nprocs_y, ', nprocs_x = ', nprocs_x  &
!       , ', myrank_y = ', myrank_y, ', myrank_x = ', myrank_x  &
!       , ', is = ', is, ', ie = ', ie, ', js = ', js, ', je = ', je


    i       = mod(imt*jmt, nprocs)
    nij1max = (imt*jmt - i)/nprocs + 1
  ! if (mpibufsize > nij1max) then
  !     write(*, *) "mpibufsize > nij1max :: ", mpibufsize, nij1max
  !     write(6, *) "using scatter/gather grd_mpi_fast."
  ! else
  !     write(*, *) "mpibufsize <= nij1max :: ", mpibufsize, nij1max
  !     write(6, *) "using scatter/gather grd_mpi_safe."
  ! end if

    if (myrank < i) then
        nij1 = nij1max
    else
        nij1 = nij1max - 1
    end if
  ! write(6, '(A,I3.3,A,I6)') 'MYRANK ',myrank,' number of grid points: nij1= ',nij1

    allocate(nij1node(nprocs))
    do n = 1, nprocs
        if (n-1 < i) then
            nij1node(n) = nij1max
        else
            nij1node(n) = nij1max - 1
        end if
    end do

    if (myrank < i) then
        bs = myrank    *nij1max + 1
        be = (myrank+1)*nij1max
    else
        bs = i*nij1max + (myrank-i)*(nij1max-1) + 1
        be = bs + nij1max - 1 - 1
    end if 
  ! write(6, '(A,I3.3,A,I6)') 'MYRANK = ', myrank, ', bs = ', bs, ', be = ', be
    if (myrank == 0 .or. myrank == 511) then
        write(6, '(A,I3.3,A,I6,A,I6)') 'MYRANK = ', myrank, ', bs = ', bs, ', be = ', be
    end if


  ! if (r_size == r_dble) then
  !     MPI_r_size = MPI_DOUBLE_PRECISION
  ! else if(r_size == r_sngl) then
  !     MPI_r_size = MPI_REAL
  ! end if

    return
end subroutine initialize_mpi


subroutine finalize_mpi
  ! use mpi

    implicit none

  ! integer :: ierr

    call MPI_FINALIZE(ierr)

    return
end subroutine finalize_mpi

end module mpi
