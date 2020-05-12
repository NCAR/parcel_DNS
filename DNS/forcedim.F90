
  subroutine forset(KFMIN,KFMAX,ikxf,ikyf,ikzf,ikzfp,nf,L,KXA,KYA,KZA,DIMK) 


! Use this once to initialize arrays containing the forced-mode
! coordinates ik(x,y,z)f(nfmax). Nfmax is the maximum allowable
! number of forced modes.
  implicit none
  include 'param.inc'

  integer :: ikx,iky,ikz,ikza,L(IKTX,IKTY,IKTZP),ikxf(nfmax),ikyf(nfmax),ikzf(nfmax),ikzfp(nfmax),nf
  real :: kfmin,kfmax,k,kx,ky,kz,DIMK,KXA(IKTX),KYA(IKTY),KZA(IKTZ),tmp

  integer :: mype
  common/mpi/mype

  nf = 0

  DO 340 IKZ=1,IKTZP
     IKZA = MYPE*IKTZP+IKZ
     KZ = KZA(IKZA)
     DO 340 IKY=1,IKTY
        KY = KYA(IKY)
        DO 340 IKX=1,IKTX
           IF (L(IKX,IKY,IKZ).NE.1)        go to 340
           KX = KXA(IKX)
           K  = SQRT( KX*KX+KY*KY+KZ*KZ )
           if (k/DIMK.gt.kfmax .or. k/DIMK.lt.kfmin) go to 340
              nf = nf + 1
              if (nf.gt.nfmax) then
                 print*,'You need to increase NFMAX.'
                 stop
              endif
           ikxf(nf) = ikx
           ikyf(nf) = iky
           ikzf(nf) = ikza
           ikzfp(nf) = ikz
340 continue

  return
  end



  subroutine force(delt,zxnk,zynk,zznk,ikxf,ikyf,ikzf,ikzfp,nf,L,KXA,KYA,KZA,nt,nout,DIMK) 

! Use this after call convol. It's followed up with a 
! call proj(nzxk,nzyk,nzzk) to remove divergence of the forcing.

!  use fred_mod
  use dyn_mod
  implicit none
  include 'param.inc'
  include 'mpif.h'

  integer :: ikxf(nfmax),ikyf(nfmax),ikzf(nfmax),ikzfp(nfmax),nf,i,L(iktx,ikty,iktzp),nt,nout
  complex, dimension(iktx,ikty,iktzp) :: zxnk,zynk,zznk
  complex :: zi
  real :: kxa(iktx),kya(ikty),kza(iktz),mult,ke,phase,kx,ky,kz,wk,vz,dimk,tmp
  real :: efornew
  real(8) :: delt
  integer :: mype,ierror
  common/mpi/mype
  
  zi    = CMPLX(0.,1.)

  KE = 0.0
  DO 40 I=1,NF  !!! indexes of these forcing modes in globe
     KX = KXA(IKXF(i))
     KY = KYA(IKYF(i))
     KZ = KZA(IKZF(i))
     WK  =  KX*KX+KY*KY+KZ*KZ

     VZ = real( zxnk(ikxf(i),ikyf(i),ikzfp(i))*CONJG(zxnk(ikxf(i),ikyf(i),ikzfp(i))) )
     KE = KE + VZ/WK
     VZ = real( zynk(ikxf(i),ikyf(i),ikzfp(i))*CONJG(zynk(ikxf(i),ikyf(i),ikzfp(i))) )
     KE = KE + VZ/WK
     VZ = real( zznk(ikxf(i),ikyf(i),ikzfp(i))*CONJG(zznk(ikxf(i),ikyf(i),ikzfp(i))) )
     KE = KE + VZ/WK

40 CONTINUE

!if ( MPI == 1 ) then
    call mpi_allreduce(ke,tmp,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierror)
    ke=tmp
!endif


 EFORnew = KE + 2*edr*delt !new forcing

  MULT = SQRT(EFORnew/KE)

  IF ( mod(nt,nout).eq. 0 .and. mype .eq. 0) print*,'MULTIPLICATION FACTOR = ',MULT,EFORnew,KE

  PHASE = 0.       ! Must not modify phase of forced modes

  do 100 i=1,nf
     zxnk(ikxf(i),ikyf(i),ikzfp(i)) = MULT*CEXP(ZI*PHASE)*zxnk(ikxf(i),ikyf(i),ikzfp(i))
     zynk(ikxf(i),ikyf(i),ikzfp(i)) = MULT*CEXP(ZI*PHASE)*zynk(ikxf(i),ikyf(i),ikzfp(i))
     zznk(ikxf(i),ikyf(i),ikzfp(i)) = MULT*CEXP(ZI*PHASE)*zznk(ikxf(i),ikyf(i),ikzfp(i))
100 continue

  return
  end
