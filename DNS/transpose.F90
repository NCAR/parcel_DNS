subroutine mpitranspose(x,n1,n2,n3p,xt,n3,n2p,npe,xsend,xrecv)

  implicit none
 
  include 'mpif.h'

  integer :: status(MPI_STATUS_SIZE),npe
  integer :: n1,n2,n3,n2p,n3p
  integer :: i,j,k,ip,ito,ifrom,jblock,kblock,koff,n123p
  complex, dimension(n1,n2,n3p) :: x
  complex, dimension(n1,n3,n2p) :: xt
  complex, dimension(n1,n2p,n3p,npe) :: xsend,xrecv

  integer :: mype,ierror
  common/mpi/mype

  n123p=n1*n2p*n3p

  do ip=1,npe
     do k=1,n3p
        do j=1,n2p
           jblock=j+(ip-1)*n2p
           do i=1,n1
              xsend(i,j,k,ip) = x(i,jblock,k)
           enddo
        enddo
     enddo
  enddo
  
  call mpi_alltoall(xsend,n123p,MPI_COMPLEX,xrecv,n123p,MPI_COMPLEX,MPI_COMM_WORLD,ierror)
 
  do ip=1,npe
     do j=1,n2p
        do k=1,n3p
           do i=1,n1
              xt(i,k+(ip-1)*n3p,j) = xrecv(i,j,k,ip)
           enddo
        enddo
     enddo
  enddo
  
end subroutine mpitranspose


subroutine serialtranspose(x,xt,n1,n2,n3)

  implicit none

  integer :: n1,n2,n3
  integer :: i,j,k
  complex, dimension(n1,n2,n3) :: x
  complex, dimension(n1,n3,n2) :: xt

  do j=1,n2
     do k=1,n3
        do i=1,n1
           xt(i,k,j) = x(i,j,k)
        enddo
     enddo
  enddo

end subroutine serialtranspose
