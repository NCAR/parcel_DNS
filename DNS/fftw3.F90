!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c FFTS
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine fftwrk(zr,zk)

  use scratch_mod
  implicit none
  include 'param.inc'
!  include 'fftw3.f'

  complex, dimension(iktx,ikty,iktzp) :: zk
  real,    dimension(n1d,n3d,n2dp)    :: zr
  real :: norm
  integer :: i,j,k,ikstart
  integer, parameter :: n1n3=n1d*n3d
  integer, parameter :: iktxz=iktx*iktz,iktxy=iktx*ikty,iktyp=ikty/npe
  integer, parameter :: iktxyzp=iktx*ikty*iktzp

  integer*8 :: plan2_rk,plan2_kr,plan1_rk,plan1_kr
  common/FTW2D/plan2_rk,plan2_kr
  common/FTW1D/plan1_rk,plan1_kr

  integer :: mype
  common/mpi/mype

!  if( .not. scratch_alloue ) then
!        allocate( zkt(iktx,iktz,iktyp), zk1(iktx,ikty,iktzp), zk2(iktx,ikty,iktzp) )
!  endif

  norm = float(n1*n2*n3) !#grids

  ! Do 2D (x,z) transforms at n2 levels all at once 
  ! Note, output data has dimensions iktx,iktz,iktyp 
  ! zk1 is scratch
  ! real to k-space
  call rfftwnd_f77_real_to_complex(plan2_rk,n2pe,zr,1,n1n3,zk1,1,iktxz)

  ! 2D transforms are in-place so output is in zk
  ! but has dimensions iktx,iktz,iktyp
  ! copy to zkt in prep for transpose
!!  zkt=zk ! this doesn't work on sharcnet!  Do this instead:  (sharcnet: shared hierarchical academic research computing network)
   do i=1,iktxyzp
     zkt(i,1,1)=zk(i,1,1)
   enddo

  ! Transpose zkt(iktx,iktz,iktyp) -> zk(iktx,ikty,iktzp)
  ! Transpose output from previous step to have dimensions iktx,ikty,iktzp
!if (MPI == 1 ) then
  call mpitranspose(zkt,iktx,iktz,iktyp,zk,ikty,iktzp,npe,zk1,zk2) ! k-space to real
!else
!  call serialtranspose(zkt,zk,iktx,iktz,ikty)
!endif

  ! Do remaining 1D (y) transforms at iktx*iktz rows
  do k=1,iktzp
     ikstart=1+(k-1)*iktxy
     call fftw_f77(plan1_rk,iktx,zk(ikstart,1,1),iktx,1,zk1,iktx,1)
  enddo
  
  ! Normalize
  zk=zk/norm
end subroutine fftwrk


subroutine fftwkr(zr,zk)

  use scratch_mod
  implicit none
  include 'param.inc'
  include 'mpif.h'

  complex, dimension(iktx,ikty,iktzp) :: zk
  real,    dimension(n1d,n3d,n2dp)    :: zr
  integer :: i,k,ikstart
  integer, parameter :: n1n3=n1d*n3d
  integer, parameter :: iktxz=iktx*iktz,iktxy=iktx*ikty,iktyp=ikty/npe
  integer, parameter :: iktxyzp=iktx*ikty*iktzp

  integer*8 :: plan2_rk,plan2_kr,plan1_rk,plan1_kr
  common/FTW2D/plan2_rk,plan2_kr
  common/FTW1D/plan1_rk,plan1_kr

  integer :: mype,ierror
!  integer :: mype
  common/mpi/mype


!  call mpi_barrier(MPI_COMM_WORLD,ierror)
!  if (mype .eq. 0) print*, 'In fftwkr before allocate'
!  call mpi_barrier(MPI_COMM_WORLD,ierror)


!  if( .not. scratch_alloue ) then
!        allocate( zkt(iktx,iktz,iktyp), zk1(iktx,ikty,iktzp), zk2(iktx,ikty,iktzp) )
!  endif

!  call mpi_barrier(MPI_COMM_WORLD,ierror)
!  if (mype .eq. 0) print*, 'In fftwkr after allocate' 
!  call mpi_barrier(MPI_COMM_WORLD,ierror)


  call realit(zk)

  ! Do 1D (y) transforms at iktx*iktz rows
  do k=1,iktzp
     ikstart=1+(k-1)*iktxy
     call fftw_f77(plan1_kr,iktx,zk(ikstart,1,1),iktx,1,zk1,iktx,1)
  enddo

  ! Transpose zk(iktx,ikty,iktzp) -> zkt(iktx,iktz,iktyp)
  ! Note, output of transpose has dimensions iktx,iktz,iktyp but store it in zk
!if (MPI == 1) then
  call mpitranspose(zk,iktx,ikty,iktzp,zkt,iktz,iktyp,npe,zk1,zk2)
!else
!  call serialtranspose(zk,zkt,iktx,ikty,iktz)
!endif

  ! Copy transposed array back into zk in preparation for in-place 2D transforms
!! zk=zkt ! this doesnt work on sharcnet
  do i=1,iktxyzp
     zk(i,1,1)=zkt(i,1,1)
  enddo

  ! Finally do 2D (x,z) transforms at n2 levels all at once
  ! Note, input data in zk has dimensions iktx,iktz,iktyp
  ! zk1 is scrtatch
  call rfftwnd_f77_complex_to_real(plan2_kr,n2pe,zk,1,iktxz,zk1,1,n1n3)

end subroutine fftwkr
