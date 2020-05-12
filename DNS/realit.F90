
subroutine realit(zk)

! Enforces the reality condition on the plane kx=0
! by writing on modes with L(ikx,iky,ikz) = 0.

  implicit none
  include 'param.inc'
  include 'mpif.h'

  integer :: n2h,n3h
  integer :: ikx,iky,ikz,kz,inkz,ky,inky 
  complex :: zk(iktx,ikty,iktzp)

  integer :: ikza
  real :: kyy,kzz
  integer, parameter :: iktyh = ikty/2
  integer :: mype,ierror
  common/mpi/mype

!if ( MPI == 1 ) then
  integer :: status(MPI_STATUS_SIZE)
  complex :: buf1(iktyh,iktzp-1),buf2(iktyh)
  integer :: nto,nfrom,nbuf1,nbuf2,nph
!endif

!if ( MPI == 1 ) then
  nph   = npe/2
  nbuf1 = iktyh*(iktzp-1)
  nbuf2 = iktyh
!endif

  n2h = n2/2
  n3h = n3/2

! First, negative ky axis; no communication required
! Set zk(0,-ky,0) = conjg(zk(0,ky,0)) 
  if (mype.eq.0) then
     ikx = 1
     ikz = 1
     do ky=1,n2h-1
        iky  = n2h+1-ky
        inky = n2h+1+ky
        zk(ikx,inky,ikz) = conjg( zk(ikx,iky,ikz) )           
     enddo
  endif


!if ( MPI == 0 ) then
! Negative kz axis. Set Set zk(0,0,-kz) = conjg(zk(0,0,kz))  
!  ikx = 1
!  iky = 1
!  do kz=1,n3h-1
!     ikz  = n3h+1-kz
!     inkz = n3h+1+kz
!     zk(ikx,iky,inkz) = conjg( zk(ikx,iky,ikz) )           
!  enddo

! Do rest of kx=0 plane
!  ikx = 1
!  do kz=1,n3h-1
!     ikz  = n3h+1-kz
!     inkz = n3h+1+kz
!     do ky=1,n2h-1
!        iky  = n2h+1-ky 
!        inky = n2h+1+ky 
!        zk(ikx,inky,inkz) = conjg( zk(ikx,iky,ikz) )
!        zk(ikx,inky,ikz)  = conjg( zk(ikx,iky,inkz) )
!     enddo
!  enddo
!else !/* MPI */
! Next, send kz>0,ky>=0 to kz<0,ky<=0
! Write zk(0,-ky,-kz) = conjg(zk(0,ky, kz)) 
! Have to do it in 2 parts

  if (mype.le.nph-1) then  
     nto  = npe-mype-1
     buf1 = zk(1,1:ikty/2,2:iktzp)
     call mpi_send(buf1,nbuf1,MPI_COMPLEX,nto,1,MPI_COMM_WORLD,ierror)
  else
     nfrom = npe-mype-1
     call mpi_recv(buf1,nbuf1,MPI_COMPLEX,nfrom,1,MPI_COMM_WORLD,status,ierror)
     ikx = 1
     do kz=2,iktzp
        ikz  = kz-1
        inkz = iktzp+2-kz
        iky = 1
        zk(ikx,iky,inkz) = conjg( buf1(iky,ikz) )  ! do ky=0
        do ky=1,n2h-1
           iky  = n2h+1-ky 
           inky = n2h+1+ky 
           zk(ikx,inky,inkz) = conjg( buf1(iky,ikz) )
        enddo
     enddo
  endif

  if ((mype.gt.0).and.(mype.le.nph-1)) then  
     nto  = npe-mype
     buf2 = zk(1,1:ikty/2,1)
     call mpi_send(buf2,nbuf2,MPI_COMPLEX,nto,1,MPI_COMM_WORLD,ierror)
  endif
  if (mype.gt.nph) then
     nfrom = npe-mype
     call mpi_recv(buf2,nbuf2,MPI_COMPLEX,nfrom,1,MPI_COMM_WORLD,status,ierror)
     ikx  = 1
     inkz = 1
     iky  = 1
     zk(ikx,iky,inkz) = conjg( buf2(iky) )  ! do ky=0
     do ky=1,n2h-1
        iky  = n2h+1-ky 
        inky = n2h+1+ky 
        zk(ikx,inky,inkz) = conjg( buf2(iky) )
     enddo
  endif

! Finally, send kz<0,ky>=0 to kz>0,ky<=0
! Write zk(0,-ky,kz) = conjg(zk(0,ky,-kz)) 
! Have to do it in 2 parts

  if (mype.ge.nph) then  
     nto  = npe-mype-1
     buf1 = zk(1,1:ikty/2,2:iktzp)
     call mpi_send(buf1,nbuf1,MPI_COMPLEX,nto,1,MPI_COMM_WORLD,ierror)
  else
     nfrom = npe-mype-1
     call mpi_recv(buf1,nbuf1,MPI_COMPLEX,nfrom,1,MPI_COMM_WORLD,status,ierror)
     ikx = 1
     do kz=2,iktzp
        ikz  = kz-1
        inkz = iktzp+2-kz
        iky = 1
        do ky=1,n2h-1
           iky  = n2h+1-ky 
           inky = n2h+1+ky 
           zk(ikx,inky,inkz) = conjg( buf1(iky,ikz) )
        enddo
     enddo
  endif

  if (mype.gt.nph) then
     nto  = npe-mype
     buf2 = zk(1,1:ikty/2,1)
     call mpi_send(buf2,nbuf2,MPI_COMPLEX,nto,1,MPI_COMM_WORLD,ierror)
  endif
  if ((mype.gt.0).and.(mype.le.nph-1)) then
     nfrom = npe-mype
     call mpi_recv(buf2,nbuf2,MPI_COMPLEX,nfrom,1,MPI_COMM_WORLD,status,ierror)
     ikx  = 1
     inkz = 1
     iky  = 1
     do ky=1,n2h-1
        iky  = n2h+1-ky 
        inky = n2h+1+ky 
        zk(ikx,inky,inkz) = conjg( buf2(iky) )
     enddo
  endif
!endif !MPI

  return
end subroutine realit
