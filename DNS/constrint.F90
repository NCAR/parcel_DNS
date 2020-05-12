
    subroutine constr (zxk,zyk,zzk,ttk,rrk,nxk,nyk,nzk,ntk,nrk,l,uk,vk,wk,ur,vr,wr,zxr,zyr,zzr,ttr,qvr,   &
                    nxr,nyr,nzr,ntr,nrr,zi,kxa,kya,kza,ak,bk,ck,ar,br,cr,urt,vrt,wrt,ivelinit)

! calculates convolution sums (got the nonlinear terms and transform back to k-space), calls fft's, etc.
    use dyn_mod
    implicit none
    include 'param.inc'
    include 'mpif.h'

    ! ---- argument -----
    integer :: ikx,iky,ikz,ikza,l(iktx,ikty,iktzp),i,j,k,ivelinit
    integer :: nbufr,nbufs,nto,nfrom,ntofrom(npe+2)
    real :: kxa(iktx),kya(ikty),kza(iktz),kx,ky,kz
    real,     dimension(n1+1,n3+1,n2pe+1) :: urt,vrt,wrt       
    real,     dimension(n1d,n3d,n2dp)   :: ur,vr,wr,zxr,zyr,zzr,ttr,qvr,nxr,nyr,nzr,ntr,nrr,ar,br,cr
    complex, dimension(iktx,ikty,iktzp) :: uk,vk,wk,zxk,zyk,zzk,ttk,rrk,nxk,nyk,nzk,ntk,nrk,ak,bk,ck

    ! ------ local -----
    complex :: c1,c2,c3,zi

    integer :: status(mpi_status_size)
    integer :: mype,ierr
    common/mpi/mype

    nxr = 0.
    nyr = 0.
    nzr = 0.
    ntr = 0.
    nrr = 0.
    call velo (zxk,zyk,zzk,uk,vk,wk,l,kxa,kya,kza)
    if (mype .eq. 0) print*,'calling first fft'
    call fftwkr (ur,uk)
    call fftwkr (vr,vk)
    call fftwkr (wr,wk)

    call fftwkr (zxr,zxk)
    call fftwkr (zyr,zyk)
    call fftwkr (zzr,zzk)


    if (ivelinit .eq. 1) then

       urt = 0.
       vrt = 0.
       wrt = 0.

!       if ( mpi == 1 ) then

	! shift arrays one index so they can be padded with extra elements (2 for first and third direction, but 1 only for second direction)       
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

	! send left edge of slab to processor to left for droplet trajectory interpolation. also, receive left edge of slab to the right.
          call mpi_sendrecv(uvecs,nbufs,mpi_real,nto,nto,uvecr,nbufr,mpi_real,nfrom,mype,mpi_comm_world,status,ierr)

	!!! interpolation
	! added left edge of right-neighbour to urt-arrays for droplet interpolation.
          urt(1:n1,1:n3,n2pe+1) = uvecr(:,:,1) 
          vrt(1:n1,1:n3,n2pe+1) = uvecr(:,:,2)  
          wrt(1:n1,1:n3,n2pe+1) = uvecr(:,:,3) 

!       else

!          urt(1:n1,1:n3,1:n2pe) = ur(1:n1,1:n3,1:n2pe)
!          urt(1:n1,1:n3,n2pe+1) = ur(1:n1,1:n3,1)  
        
!          vrt(1:n1,1:n3,1:n2pe) = vr(1:n1,1:n3,1:n2pe)  
!          vrt(1:n1,1:n3,n2pe+1) = vr(1:n1,1:n3,1) 

!          wrt(1:n1,1:n3,1:n2pe) = wr(1:n1,1:n3,1:n2pe)  
!          wrt(1:n1,1:n3,n2pe+1) = wr(1:n1,1:n3,1)  

!       endif ! mpi

	!!!periodicity
	! impose periodicity in remaining two directions

       urt(n1+1  ,1:n3,1:n2pe+1) = urt(1     ,1:n3,1:n2pe+1)
       vrt(n1+1  ,1:n3,1:n2pe+1) = vrt(1     ,1:n3,1:n2pe+1)
       wrt(n1+1  ,1:n3,1:n2pe+1) = wrt(1     ,1:n3,1:n2pe+1)

       urt(1:n1+1,n3+1,1:n2pe+1) = urt(1:n1+1,1   ,1:n2pe+1)
       vrt(1:n1+1,n3+1,1:n2pe+1) = vrt(1:n1+1,1   ,1:n2pe+1)
       wrt(1:n1+1,n3+1,1:n2pe+1) = wrt(1:n1+1,1   ,1:n2pe+1)
	!!! for test only
	! test to see if periodicity is done correctly using mpi. values should be the same when npe = 8
	!     if (mype .eq. 0) print*, 'ur(1,1,1), wr(1,1,1) =', mype, ur(1,1,1), wr(1,1,1) 
	!     if (mype .eq. npe-1) print*, 'urt(n1+1,n3+1,n2pe+1), wrt(n1+1,n3+1,n2pe+1) =', mype, urt(n1+1,n3+1,n2pe+1), wrt(n1+1,n3+1,n2pe+1) 

    endif ! ivelinit 

    nxr = wr*zyr - vr*zzr
    nyr = ur*zzr - wr*zxr
    nzr = vr*zxr - ur*zyr

    ! print*,'first real-complex fft'
    call fftwrk (nxr,nxk)
    call fftwrk (nyr,nyk)
    call fftwrk (nzr,nzk)

    call fftwrk (zxr,zxk)
    call fftwrk (zyr,zyk)
    call fftwrk (zzr,zzk)

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
             nxk(ikx,iky,ikz) = - zi * c1 * l(ikx,iky,ikz)
             nyk(ikx,iky,ikz) = - zi * c2 * l(ikx,iky,ikz)
             nzk(ikx,iky,ikz) = - zi * c3 * l(ikx,iky,ikz)
          enddo
       enddo
    enddo


    if (thermo .eq. 1) then
       call fftwkr (ttr,ttk)
       call fftwkr (qvr,rrk)
       ar = ur*qvr
       br = vr*qvr
       cr = wr*qvr

       ur = ur*ttr
       vr = vr*ttr
       wr = wr*ttr     

       call fftwrk (ttr,ttk)
       call fftwrk (qvr,rrk)

       call fftwrk ( ur,uk)
       call fftwrk ( vr,vk)
       call fftwrk ( wr,wk)
       call fftwrk ( ar,ak)
       call fftwrk ( br,bk)
       call fftwrk ( cr,ck)

       do ikz=1,iktzp
          ikza = mype*iktzp+ikz
          kz = kza(ikza)
          do iky=1,ikty
             ky = kya(iky)
             do ikx=1,iktx
                kx = kxa(ikx)
                ntk(ikx,iky,ikz) = - zi * ( kx*uk(ikx,iky,ikz)  &
                                        + ky*vk(ikx,iky,ikz)  &
                                        + kz*wk(ikx,iky,ikz)  ) * l(ikx,iky,ikz)
                nrk(ikx,iky,ikz) = - zi * ( kx*ak(ikx,iky,ikz)  &
                                        + ky*bk(ikx,iky,ikz)  &
                                        + kz*ck(ikx,iky,ikz)  ) * l(ikx,iky,ikz)
             enddo
          enddo
       enddo
    endif

999 continue
    return
  end 
