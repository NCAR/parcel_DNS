   subroutine calmean(var,rm)
  !! calculate mean of fieldfield and subtract it from all elements
   implicit none
   include 'param.inc'
   include 'mpif.h' 
   real, dimension(n1d,n3d,n2dp) :: var
   integer :: ii,ix,iy,iz
   real :: rm,tmp
   integer :: status(MPI_STATUS_SIZE)
   integer :: mype,ierror
   common/mpi/mype

   rm=0.d0
   do 100 iy=1,n2dp
      do 100 iz=1,n3d
         do 100 ix=1,n1d
            rm=rm+var(ix,iz,iy)
100 continue
    
!   if ( MPI == 1 ) then
      call mpi_allreduce(rm,tmp,1,mpi_real,mpi_sum,mpi_comm_world,ierror)
      rm=tmp
!   endif
 !! calculate mean
      rm=rm/(n2dp*n3d*n1d)
 !! subtract the mean
   do 200 iy=1,n2dp
      do 200 iz=1,n3d
         do 200 ix=1,n1d
            var(ix,iz,iy)=var(ix,iz,iy) -rm
200 continue  
  return
  end
