
  SUBROUTINE convol (zxk,zyk,zzk,ttk,qvk,nxk,nyk,nzk,ntk,nrk,l,uk,vk,wk,ur,vr,wr,wxk,wyk,wzk,wtk,wrk,         &
                     zxr,zyr,zzr,ttr,qvr,urt,vrt,wrt,nxr,nyr,nzr,ntr,nqv,zi,kxa,kya,kza,ak,bk,ck,ar,br,cr,    &
                     time0,nt,ndumpd,ndropreal,nddone)

! CALCULATES CONVOLUTION SUMS, CALLS FFT'S, micro, ETC.! advection
  use dyn_mod
  use mic_mod
  implicit none
  include 'param.inc'
  include 'mpif.h'

  ! --- argument ----
  integer :: nt,ikx,iky,ikz,ikza,L(iktx,ikty,iktzp),i,j,k
  integer :: nbufr,nbufs,nto,nfrom,ntofrom(npe+2)
  real :: kxa(iktx),kya(ikty),kza(iktz),kx,ky,kz,tmp

  ! --- local ---
  real,     dimension(n1+1,n3+1,n2pe+1) :: urt,vrt,wrt       
  real,     dimension(n1d,n3d,n2dp)       :: ur,vr,wr,zxr,zyr,zzr,ttr,qvr,nxr,nyr,nzr,ntr,nqv,ar,br,cr
  complex, dimension(iktx,ikty,iktzp) :: uk,vk,wk,zxk,zyk,zzk,ttk,qvk,nxk,nyk,nzk,ntk,nrk,wxk,wyk,wzk,wtk,wrk,ak,bk,ck
  complex :: c1,c2,c3,zi

  integer :: status(MPI_STATUS_SIZE)
  integer :: mype,ierror
  common/mpi/mype

! DECLARATIONS FOR MICROPHYSICS
  real :: time0
  integer :: ndumpd
  integer :: ndropreal,nddone

!  character*128 :: formatvel
  real :: maxu,maxuvw,absur,absvr,abswr,maxur,maxvr,maxwr,maxii,maxjj,maxkk


  nxk = cmplx(0.,0.)
  nyk = cmplx(0.,0.)
  nzk = cmplx(0.,0.) 
  wxk = zxk
  wyk = zyk
  wzk = zzk
  
  if (thermo .eq. 1) then
     ntk = cmplx(0.,0.)
     nrk = cmplx(0.,0.)
     wtk = ttk
     wrk = qvk 
     call fftwkr(ttr,ttk)
     call fftwkr(qvr,qvk)
  endif      !thermo
if (iturb .eq. 1) then
  CALL VELO (ZXK,ZYK,ZZK,UK,VK,WK,L,KXA,KYA,KZA)

  CALL FFTWKR ( UR,UK)
  CALL FFTWKR ( VR,VK)
  CALL FFTWKR ( WR,WK)
  CALL FFTWKR (ZXR,ZXK)
  CALL FFTWKR (ZYR,ZYK)
  CALL FFTWKR (ZZR,ZZK)
endif !iturb
!-----------------------------------------------------------
! MICROPHYSICAL ADJUSTMENT -- COMPLETING TIMESTEP (NT - 1)
  if ( gomic .gt. 0) then
     call micro(urt,vrt,wrt,ttr,qvr,time0,ndumpd,ndropreal,nt,nddone)
! OUTPUT PDF OF SCALAR FIELDS
!     IF(MOD(NT-1,NOUT).EQ.0) THEN
!        CALL CALCPDFSCAL(TTR,1.,N1D,N2D,N3D,-.005,.01,0)
!        CALL CALCPDFSCAL(QVR,1.,N1D,N2D,N3D,-.0000025,.000005,0)
!     ENDIF

! SAVE VELOCITY FIELD FOR NEXT STEP IN MICRO AND DIMENSIONALIZE
! SUBTRACT AVERAGES OFF SCALAR FIELDS
  
     urt = 0.
     vrt = 0.
     wrt = 0.
  if (iturb .eq. 1) then
!    if ( MPI == 1 ) then

! Shift arrays one index so they can be padded with extra elements (2 for first and third direction, but 1 only for second direction)       
        urt(1:n1,1:n3,1:n2pe) = ur(1:n1,1:n3,1:n2pe) 
        vrt(1:n1,1:n3,1:n2pe) = vr(1:n1,1:n3,1:n2pe) 
        wrt(1:n1,1:n3,1:n2pe) = wr(1:n1,1:n3,1:n2pe) 

        nbufr = n1*n3*3 
        nbufs = n1*n3*3 
        uvecs(1:n1,1:n3,1) = ur(1:n1,1:n3,1) 
        uvecs(1:n1,1:n3,2) = vr(1:n1,1:n3,1) 
        uvecs(1:n1,1:n3,3) = wr(1:n1,1:n3,1)

        do i=1,npe
           ntofrom(i+1) = i-1
        enddo
        ntofrom(1) = npe -1
        ntofrom(npe+2) = 0
        nto = ntofrom(mype+1)
        nfrom = ntofrom(mype+3)

! Send left edge of slab to processor to left for droplet trajectory interpolation. Also, receive left edge of slab to the right.
        call mpi_sendrecv(uvecs,nbufs,MPI_REAL,nto,nto,uvecr,nbufr,MPI_REAL,nfrom,mype,MPI_COMM_WORLD,status,ierror)

! Added left edge of right-neighbour to urt-arrays for droplet interpolation.
        urt(1:n1,1:n3,n2pe+1) = uvecr(:,:,1) 
        vrt(1:n1,1:n3,n2pe+1) = uvecr(:,:,2) 
        wrt(1:n1,1:n3,n2pe+1) = uvecr(:,:,3) 

!    else!MPI

!        urt(1:n1,1:n3,1:n2pe) = ur(1:n1,1:n3,1:n2pe) 
!        urt(1:n1,1:n3,n2pe+1) = ur(1:n1,1:n3,1) 
        
!        vrt(1:n1,1:n3,1:n2pe) = vr(1:n1,1:n3,1:n2pe) 
!        vrt(1:n1,1:n3,n2pe+1) = vr(1:n1,1:n3,1) 

!        wrt(1:n1,1:n3,1:n2pe) = wr(1:n1,1:n3,1:n2pe)
!        wrt(1:n1,1:n3,n2pe+1) = wr(1:n1,1:n3,1) 

!    endif ! MPI

! Impose periodicity in remaining two directions

     urt(n1+1,1:n3,:) = urt(1,1:n3,:)
     vrt(n1+1,1:n3,:) = vrt(1,1:n3,:)
     wrt(n1+1,1:n3,:) = wrt(1,1:n3,:)

     urt(:,n3+1,:) = urt(:,1,:)
     vrt(:,n3+1,:) = vrt(:,1,:)
     wrt(:,n3+1,:) = wrt(:,1,:)

! Test to see if periodicity is done correctly using MPI. Values should be the same when npe = 8
!     if (mype .eq. 0) print*, 'ur(1,1,1), wr(1,1,1) =', mype, ur(1,1,1), wr(1,1,1) 
!     if (mype .eq. npe-1) print*, 'urt(n1+1,n3+1,n2pe+1), wrt(n1+1,n3+1,n2pe+1) =', mype, urt(n1+1,n3+1,n2pe+1), wrt(n1+1,n3+1,n2pe+1) 

  endif !iturb

! Subtract averages of scalar fields
     if (thermo .eq. 1) then
        ttr = ttr - rmth
        qvr = qvr - rmqv
     endif

  ENDIF ! GOMIC
!-----------------------------------------------

  if (iturb .eq. 1) then
! 	OUTPUT PDF OF CHOSEN FIELDS
   	 if (mod(nt-1,nout).eq.0 ) then
!	 MODIFY THESE LATER
!  	   call calcpdfvel(ur,n1d,n2d,n3d,-.5,1.,0)
!  	   call calcpdfvel(vr,n1d,n2d,n3d,-.5,1.,0)
!   	   call calcpdfvel(wr,n1d,n2d,n3d,-.5,1.,0)
         
!  	   call calcpdfvort(zxr,n1d,n2d,n3d,-250.,500.,0)
!  	   call calcpdfvort(zzr,n1d,n2d,n3d,-250.,500.,0)

   	  	maxu  = - 10.
   	   DO 109 K=1,n3d
     	     DO 109 J=1,n2dp
     	       DO 109 I=1,n1d
       		    maxu = max(maxu,abs(ur(i,k,j)),abs(vr(i,k,j)),abs(wr(i,k,j)))
109  CONTINUE
!     	   if ( MPI == 1 ) then
          	call mpi_reduce(maxu,tmp,1,MPI_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierror)
          	maxu=tmp
!     	   endif  
     	   if (mype .eq. 0) then
        	print*,'  '
        	print*,'***MAX COURANT NUMBER : ', maxu*delt/h,maxu
        	print*,'  '
!        	write(97,*) '***MAX COURANT NUMBER : ', maxu*delt/h,maxu
     	   endif
    	 endif ! nt 
   endif!iturb


! PROCEED WITH CALCULATION OF ADVECTIVE TERMS
     if(iturb .eq. 1) then
  	nxr = wr*zyr - vr*zzr
  	nyr = ur*zzr - wr*zxr
  	nzr = vr*zxr - ur*zyr

  	CALL FFTWRK (NXR,NXK)
  	CALL FFTWRK (NYR,NYK)
  	CALL FFTWRK (NZR,NZK)
  	do ikz=1,iktzp
  	   ikza = mype*iktzp+ikz
   	   kz = kza(ikza)
   	   do iky=1,ikty
        	ky = kya(iky)
        	do ikx=1,iktx
        		kx = kxa(ikx)
           		c1 = ky*nzk(ikx,iky,ikz) - kz*nyk(ikx,iky,ikz)
           		c2 = kz*nxk(ikx,iky,ikz) - kx*nzk(ikx,iky,ikz)
           		c3 = kx*nyk(ikx,iky,ikz) - ky*nxk(ikx,iky,ikz)
           		nxk(ikx,iky,ikz) = - zi*c1*L(ikx,iky,ikz)
           		nyk(ikx,iky,ikz) = - zi*c2*L(ikx,iky,ikz)
           		nzk(ikx,iky,ikz) = - zi*c3*L(ikx,iky,ikz)

           		zxk(ikx,iky,ikz) = wxk(ikx,iky,ikz)  
           		zyk(ikx,iky,ikz) = wyk(ikx,iky,ikz)
           		zzk(ikx,iky,ikz) = wzk(ikx,iky,ikz)

        	enddo
     	   enddo
  	enddo
     endif
     if (thermo .eq. 1) then
     	ar = ur*qvr
     	br = vr*qvr
     	cr = wr*qvr
 
     	ur = ur*ttr
     	vr = vr*ttr
     	wr = wr*ttr

     	call fftwrk(ar,ak)
     	call fftwrk(br,bk)
    	call fftwrk(cr,ck)
    	call fftwrk(ur,uk)
    	call fftwrk(vr,vk)
    	call fftwrk(wr,wk)

        do ikz=1,iktzp
           ikza = mype*iktzp+ikz
           kz = kza(ikza)
           do iky=1,ikty
              ky = kya(iky)
              do ikx=1,iktx
                 kx = kxa(ikx)
                 ntk(ikx,iky,ikz) = - zi * ( kx*uk(ikx,iky,ikz) + ky*vk(ikx,iky,ikz) + kz*wk(ikx,iky,ikz) )*L(ikx,iky,ikz)
                 nrk(ikx,iky,ikz) = - zi * ( kx*ak(ikx,iky,ikz) + ky*bk(ikx,iky,ikz) + kz*ck(ikx,iky,ikz) )*L(ikx,iky,ikz)
              enddo
           enddo
        enddo
        if (gomic .gt. 0) then
           call fftwrk(ttr,ttk)
           call fftwrk(qvr,qvk)
    	else
 	   ttk=wtk
	   qvk=wrk
	endif!gomic
  endif!thermo 

999 continue
  return
  END 
