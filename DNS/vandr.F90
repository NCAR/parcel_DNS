  SUBROUTINE rspace (zxk,zyk,zzk,ttk,zxr,zyr,zzr,ttr,nxk,nyk,nzk,ntk,nxr,nyr,nzr,ntr,    &
                      uk,vk,wk,ur,vr,wr,l,zi,kxa,kya,kza,time)

! CALCULATES REAL-SPACE FIELDS. calc k-space velo from vorticity field & real space fields
!  use fred_mod   
  implicit none
  include 'param.inc'
  include 'mpif.h'
     
  INTEGER :: L(iktx,ikty,iktzp),IKX,IKY,IKZ,ikza

  real,    dimension(n1d,n3d,n2dp)    :: ur,vr,wr,nxr,nyr,nzr,zxr,zyr,zzr,ttr,ntr
  REAL :: TIME,KX,KY,KZ,KXA(IKTX),KYA(IKTY),KZA(IKTZ),v,s,tmp,wk2

  complex, dimension(iktx,ikty,iktzp) :: uk,vk,wk,nxk,nyk,nzk,zxk,zyk,zzk,ttk,ntk
  complex :: zi

  integer :: mype,ierror
  common/mpi/ mype

  v = 0.
  s = 0.
  CALL VELO (ZXK,ZYK,ZZK,UK,VK,WK,L,KXA,KYA,KZA)
  if (mype .eq. 0) then
    PRINT*,'                '
    PRINT*,'VELOCITY FIELD DIVERGENCE:'
  endif 
!  call mpi_barrier(MPI_COMM_WORLD,ierror)
  CALL DIVERG (UK,VK,WK,L,KXA,KYA,KZA)

! Calculate longitudinal derivatives
  DO 5 ikz = 1, iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     DO 5 IKY = 1, IKTY
        KY = KYA(IKY)
        DO 5 IKX = 1, IKTX
           KX = KXA(IKX)
           UK(IKX,IKY,IKZ)  = ZI*KX* UK(IKX,IKY,IKZ)
           VK(IKX,IKY,IKZ)  = ZI*KY* VK(IKX,IKY,IKZ)
           WK(IKX,IKY,IKZ)  = ZI*KZ* WK(IKX,IKY,IKZ)
 5 CONTINUE

! Go to real space
  call fftwkr(ur,uk)
  call fftwkr(vr,vk)
  call fftwkr(wr,wk)

  do 15 ikz = 1, n3d
     do 15 iky = 1, n2dp
        do 15 ikx = 1, n1d
           v = v + ur(ikx,ikz,iky)**2
           s = s + ur(ikx,ikz,iky)**3
           v = v + vr(ikx,ikz,iky)**2
           s = s + vr(ikx,ikz,iky)**3
           v = v + wr(ikx,ikz,iky)**2
           s = s + wr(ikx,ikz,iky)**3
 15 continue

!if ( MPI == 1 ) then
    call mpi_reduce(v,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    v=tmp
    call mpi_reduce(s,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    s=tmp
!endif  

  if (mype .eq. 0) then
    v = v/(3.*n1*n2*n3)
    s = -s/(3.*n1*n2*n3)
    s = s/v**1.5
    print*,'Skewness of velocity derivatives = ',s
  endif
  nxk = cmplx(0.,0.)  
  nyk = cmplx(0.,0.)  
  nzk = cmplx(0.,0.)  
  uk = cmplx(0.,0.)  
  vk = cmplx(0.,0.)  
  wk = cmplx(0.,0.)  
  
  return
  end


  SUBROUTINE VELO (ZX,ZY,ZZ,U,V,W,L,KXA,KYA,KZA)

! CALCULATES K-SPACE VELOCITY FROM K-SPACE VORTICITY.
! CURL (VORTICITY) = - LAPLACIAN (VELOCITY) IF VELOCITY
! IS SOLENOIDAL.

  implicit none
  include 'param.inc'
    
  integer :: L(iktx,ikty,iktzp),ikx,iky,ikz,ikza
  real :: KX,KY,KZ,K2,KXA(IKTX),KYA(IKTY),KZA(IKTZ)
  complex, dimension(iktx,ikty,iktzp) :: zx,zy,zz,u,v,w
  complex :: zi,c1,c2,c3

  integer :: mype
  common/mpi/ mype
  
  zi = CMPLX(0.,1.)
  u = cmplx(0.,0.)  
  v = cmplx(0.,0.)  
  w = cmplx(0.,0.)  

  DO 70 ikz = 1, iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     DO 70 IKY=1,IKTY
        KY = KYA(IKY)
        DO 70 IKX=1,IKTX
           IF (L(IKX,IKY,IKZ).NE.1) GO TO 70
           KX = KXA(IKX)
           K2 = KX*KX+KY*KY+KZ*KZ
           C1 = +  KY*ZZ(IKX,IKY,IKZ) -  KZ*ZY(IKX,IKY,IKZ) !eps_ijk(dk/dxj)
           C2 = +  KZ*ZX(IKX,IKY,IKZ) -  KX*ZZ(IKX,IKY,IKZ)
           C3 = +  KX*ZY(IKX,IKY,IKZ) -  KY*ZX(IKX,IKY,IKZ)
           U(IKX,IKY,IKZ) = zi * C1 / K2 ! = vorticity/k**2
           V(IKX,IKY,IKZ) = zi * C2 / K2
           W(IKX,IKY,IKZ) = zi * C3 / K2

 70 CONTINUE

  RETURN
  END
