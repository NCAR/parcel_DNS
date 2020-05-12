
    SUBROUTINE initturb (seed,zx,zy,zz,tt,qv,ampk,ampp,l,irest,kxa,kya,kza,ktrunc,ts,dimk)

! INITIALISES STUFF LIKE WAVENUMBERS, INDICES, SPECTRA, PHASES, ETC.
    use dyn_mod
!    use fred_mod
    implicit none
    include 'param.inc'
    include 'mpif.h'
 

    complex, dimension(iktx,ikty,iktzp) :: zx,zy,zz,tt,qv
    INTEGER :: IKX,IKY,IKZ,ikza,L(iktx,ikty,iktzp),IREST,JJ,idum
    REAL :: KXA(IKTX),KYA(IKTY),KZA(IKTZ),TS,RAN1,twopi,DIMK
 !   real, dimension(n1d,n3,n2dp) :: ttr,qvr

    integer :: ktp,nyy,nzz,nmx,nmy,nmz,kx1,ky1,kz1,ntime,seed

    COMPLEX :: ZI

    real :: kx,ky,kz,wk,ek,ampk,phase,kinen,r1,r2,vk,ktrunc,tmp
    real :: ampp,poten,tk

!    real, allocatable, dimension(:,:,:)    :: zxr,zyr,zzr,ur,vr,wr!ttr,qvr,

    integer :: mype,ierror
    common/mpi/ mype

!    allocate( zxr(n1d,n3d,n2dp), zyr(n1d,n3d,n2dp), zzr(n1d,n3d,n2dp))!, ttr(n1d,n3d,n2dp), qvr(n1d,n3d,n2dp))
!    allocate( ur(n1d,n3d,n2dp), vr(n1d,n3d,n2dp), wr(n1d,n3d,n2dp) )

    zi = cmplx(0.,1.)
    twopi = 4.*asin(1.)

! Initialize wavenumber arrays
    do  ikx = 1,iktx !iktx=n/2+1
       kxa(ikx) = dimk*float(ikx-1)
    enddo
    do iky = 1,ikty           !(1,n)
       jj = iky - 1           !(0,n-1)
       if (iky.gt.kty)   jj = jj - 2*kty   !(iky > n/2)
       if (iky.eq.kty+1) jj = 0            !(iky = n/2+1 )
       if (iky.gt.2*kty) jj = 0            !(iky > n)
       kya(iky) = dimk*float(jj)           !(-n/2,n/2)
    enddo
    do ikz = 1,iktz
       jj = ikz - 1
       if (ikz.gt.ktz)   jj = jj - 2*ktz
       if (ikz.eq.ktz+1) jj = 0
       if (ikz.gt.2*ktz) jj = 0
       kza(ikz) = dimk*float(jj)
    enddo
    L  = 1
    zx = cmplx(0.,0.)
    zy = cmplx(0.,0.)
    zz = cmplx(0.,0.)
    tt = cmplx(0.,0.)
    qv = cmplx(0.,0.)
    KINEN = 0.
    POTEN = 0.
! L(ikx,iky,ikz) is unity for retained modes and zero beyond the truncation, at k=0, and for half the modes on the plane kx=0.

    do 5 ikz = 1,iktzp ! iktzp = iktz/npe . trunc. w# for each processes
       ikza = mype*iktzp+ikz
       kz = kza(ikza)
       do 5 iky = 1,ikty
          ky = kya(iky)
          do 5 ikx = 1,iktx
             kx = kxa(ikx)
             wk = kx*kx + ky*ky + kz*kz
             wk = sqrt( wk )
	! Modes N/2(aussi impose dans fftg)
             if (iky.eq.kty+1)                          L(ikx,iky,ikz) = 0
             if (ikza.eq.ktz+1)                         L(ikx,iky,ikz) = 0
	! Symetrie (on kx=0 plane)
             if (kx.lt.0)                               L(ikx,iky,ikz) = 0
             if (kx.eq.0 .and. ky.lt.0)                 L(ikx,iky,ikz) = 0
             if (kx.eq.0 .and. ky.eq.0 .and. kz.lt.0)   L(ikx,iky,ikz) = 0
	! Mode 0 et Troncature
             if (wk.eq.0 .or. wk.gt.ifix(dimk*ktrunc+0.5)-0.5)       L(ikx,iky,ikz) = 0  
	! Mode a l'exterierur du domaine(aussi impose dans fftg)
             if (iky.gt.2*kty)                          L(ikx,iky,ikz) = 0
             if (ikza.gt.2*ktz)                         L(ikx,iky,ikz) = 0
5   continue

! Restart from output file? 
    IF (IREST.ne.0) THEN  ! Begin Restart mode
       if (mype .eq. 0)    print*,'RESTART FROM OUTPUT DATA. '
       if (netcdf == 1 ) then
          call ncreadrst(zx,zy,zz,tt,qv,ts)
       else
          print*, 'Sorry, do not know how to restart without NetCDF'
          stop
       endif
       if (mype .eq. 0) then
          print*,'RESTART #',IREST
          print*,'Just read time= ',ts
       endif

    ELSE ! irest = 0, Initialize all the arrays:
       idum = seed
       IF(iturb .eq. 1) then
          DO 30 ikz = 1,iktzp
             ikza = mype*iktzp+ikz
             kz = kza(ikza)

             DO 30 IKY=1,IKTY
                KY = KYA(IKY)
                DO 30 IKX=1,IKTX
                   KX = KXA(IKX)

                   IF (L(ikx,iky,ikz) .eq. 1) THEN       
                      WK = SQRT(KX*KX+KY*KY+KZ*KZ)
                      !EK = INITIAL ENERGY SPECTRUM.
                      EK = AMPK*EXP(-1.5*WK/DIMK)
                      !EK = AMPK*EXP(-WK/DIMK)
                      !VK = INITIAL VORTICITY SPECTRUM.
                      VK = WK**2 * EK
                  
                      PHASE = ran1(idum)*twopi
                      ZX(IKX,IKY,IKZ) = SQRT(VK)*DIMK/WK * CEXP(ZI*PHASE)
                      PHASE = ran1(idum)*twopi
                      ZY(IKX,IKY,IKZ) = SQRT(VK)*DIMK/WK * CEXP(ZI*PHASE)
                      PHASE = ran1(idum)*twopi
                      ZZ(IKX,IKY,IKZ) = SQRT(VK)*DIMK/WK * CEXP(ZI*PHASE)

                      if (kx .eq. 0.) then
                         ZX(IKX,IKY,IKZ) = CMPLX(real(ZX(IKX,IKY,IKZ)),0.)
                         ZY(IKX,IKY,IKZ) = CMPLX(real(ZY(IKX,IKY,IKZ)),0.)
                         ZZ(IKX,IKY,IKZ) = CMPLX(real(ZZ(IKX,IKY,IKZ)),0.)
                      endif
                   ENDIF ! L
30   CONTINUE
       ENDIF !iturb
       if (thermo .eq. 1) then
          DO 39 ikz = 1,iktzp
             ikza = mype*iktzp+ikz
             kz = kza(ikza)
             DO 39 IKY=1,IKTY
                KY = KYA(IKY)
                DO 39 IKX=1,IKTX
                   KX = KXA(IKX)
                   IF (L(ikx,iky,ikz) .eq. 1) THEN
                      WK = SQRT(KX*KX+KY*KY+KZ*KZ)
                      TK = AMPK*EXP(-1.5*WK/DIMK)!TK = INITIAL TEMPERATURE SPECTRUM.
                      !TK = AMPK*EXP(-WK/DIMK)
                    
                      PHASE = ran1(idum)*twopi
                      TT(IKX,IKY,IKZ) = SQRT(TK)*DIMK/WK * CEXP(ZI*PHASE)
                      PHASE = ran1(idum)*twopi
                      QV(IKX,IKY,IKZ) = SQRT(TK)*DIMK/WK * CEXP(ZI*PHASE)
                      if (kx.eq.0. .and. IREST.eq.0) then
                         TT(IKX,IKY,IKZ) = CMPLX(real(TT(IKX,IKY,IKZ)),0.)
                         QV(IKX,IKY,IKZ) = CMPLX(real(QV(IKX,IKY,IKZ)),0.)
                      endif
                   ENDIF ! L
39   CONTINUE
       endif ! thermo
       if (iturb .eq. 1) then
          CALL PROJ (ZX,ZY,ZZ,L,KXA,KYA,KZA)
          if (mype .eq. 0) print*,'Divergence after projection.'
          CALL DIVERG (ZX,ZY,ZZ,L,KXA,KYA,KZA)

! NORMALISE AMPLITUDES.
          DO 31 ikz = 1,iktzp
             ikza = mype*iktzp+ikz
             kz = kza(ikza)
             DO 31 IKY=1,IKTY
                KY = KYA(IKY)
                DO 31 IKX=1,IKTX
                   KX = KXA(IKX)
                   WK = SQRT( KX*KX+KY*KY+KZ*KZ )
                   IF (L(ikx,iky,ikz).ne.1) GO TO 31
                   R1     =  real( ZX(IKX,IKY,IKZ) )**2
                   R2     = IMAG( ZX(IKX,IKY,IKZ) )**2
                   KINEN  = KINEN + (R1 + R2)/WK**2
                   R1     =  real( ZY(IKX,IKY,IKZ) )**2
                   R2     = IMAG( ZY(IKX,IKY,IKZ) )**2
                   KINEN  = KINEN + (R1 + R2)/WK**2
                   R1     =  real( ZZ(IKX,IKY,IKZ) )**2
                   R2     = IMAG( ZZ(IKX,IKY,IKZ) )**2
                   KINEN  = KINEN + (R1 + R2)/WK**2

31    CONTINUE
       endif !iturb
       if (thermo .eq. 1) then
          DO 49 ikz = 1,iktzp
             ikza = mype*iktzp+ikz
             kz = kza(ikza)
             DO 49 IKY=1,IKTY
                KY = KYA(IKY)
                DO 49 IKX=1,IKTX
                   KX = KXA(IKX)
                   WK = SQRT( KX*KX+KY*KY+KZ*KZ )
                   R1     = real( TT(IKX,IKY,IKZ) )**2
                   R2     = IMAG( TT(IKX,IKY,IKZ) )**2
                   POTEN  = POTEN + (R1 + R2)
49    CONTINUE
       endif ! thermo
       IF (expsiv.NE.0. .AND. temgrad.NE.0.) THEN
          POTEN = expsiv*POTEN/temgrad
       ENDIF

!       if ( MPI == 1 ) then
          call mpi_allreduce(kinen,tmp,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierror)
          kinen=tmp
          call mpi_allreduce(poten,tmp,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierror)
          poten=tmp
!       endif!mpi
! Normalise amplitudes now that kinen and poten are known (new method for Fortran 90)
       if (iturb .eq. 1 ) then
          zx = zx * sqrt(ampk/(kinen))
          zy = zy * sqrt(ampk/(kinen))
          zz = zz * sqrt(ampk/(kinen))
       endif !iturb

       if (thermo .eq. 1 .and. poten .ne. 0) then
          tt = tt * sqrt(ampp/(poten))
          qv = qv * sqrt(ampp/(poten))
       !subtract the mean
       endif !!thermo&poten
    endif ! irest

!    deallocate( zxr, zyr, zzr,  ur, vr, wr )!ttr, qvr,

    END



    SUBROUTINE DIVERG (ZX,ZY,ZZ,L,KXA,KYA,KZA)

! CALCULATES DIVERGENCE OF A VECTOR ZX,Y,Z
!    use fred_mod
    implicit none
    include 'param.inc'
    include 'mpif.h'
  
    INTEGER :: ikx,iky,ikz,ikza,L(iktx,ikty,iktzp)
    REAL :: KX,KY,KZ,KXA(IKTX),KYA(IKTY),KZA(IKTZ),div,divtot!tmp
    complex, dimension(iktx,ikty,iktzp) :: zx,zy,zz
    complex :: c1,zi

    integer :: mype,ierror
    common/mpi/ mype
  
    zi = CMPLX(0.,1.)

    DIV = 0. !divergence of the vorticity

    DO 70 ikz = 1, iktzp
        ikza = mype*iktzp+ikz
        kz = kza(ikza)
        DO 70 IKY=1,IKTY
            KY = KYA(IKY)
            DO 70 IKX=1,IKTX
           KX = KXA(IKX)
           C1 = KX*ZX(IKX,IKY,IKZ) + KY*ZY(IKX,IKY,IKZ) + KZ*ZZ(IKX,IKY,IKZ)
           C1 = C1 * L(IKX,IKY,IKZ)
           DIV = DIV + 2.*IMAG(zi*C1)
 70 CONTINUE

!    if ( MPI == 1 ) then
    call mpi_reduce(div,divtot,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierror)
!    divtot=tmp
!    endif!mpi

    if (mype .eq. 0) PRINT*,'DIVERGENCE = ',DIVTOT

    RETURN
    END



  subroutine proj(zx,zy,zz,L,kxa,kya,kza)

! Fourier-space determination of the solenoidal part of a vector zx,y,z.

  implicit none
  include 'param.inc'

  integer :: ikx,iky,ikz,ikza,L(iktx,ikty,iktzp)
  real :: kx,ky,kz,k2,kxa(iktx),kya(ikty),kza(iktz)
  complex, dimension(iktx,ikty,iktzp) :: zx,zy,zz
  complex :: c1,c2,c3

  integer :: mype
  common/mpi/mype

  do ikz = 1,iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     do iky = 1,ikty
        ky = kya(iky)
        do ikx = 1,iktx
	   kx = kxa(ikx)     
           if (L(ikx,iky,ikz).eq.1) then
              k2 = max(kx*kx + ky*ky + kz*kz, 1.e-15)
              c1 =  (k2-kx*kx)*zx(ikx,iky,ikz) - kx*ky*zy(ikx,iky,ikz) - kx*kz*zz(ikx,iky,ikz)
              c2 = -ky*kx*zx(ikx,iky,ikz) + (k2-ky*ky)*zy(ikx,iky,ikz) - ky*kz*zz(ikx,iky,ikz)
              c3 = -kz*kx*zx(ikx,iky,ikz) - kz*ky*zy(ikx,iky,ikz) + (k2-kz*kz)*zz(ikx,iky,ikz)
              zx(ikx,iky,ikz) = c1 / k2
              zy(ikx,iky,ikz) = c2 / k2
              zz(ikx,iky,ikz) = c3 / k2
           endif
        enddo
     enddo
  enddo

  return
  end subroutine proj



  FUNCTION ran1(idum)

  implicit none
  integer :: idum
  real ::ran1
  integer, parameter  :: ia=16807,im=2147483647,iq=127773,ir=2836,ntab=32,ndiv=1+(im-1)/ntab
  real, parameter  :: am=1./im,eps=1.2e-7,rnmx=1.-eps
  integer :: j,k,iv(ntab),iy
  save iv,iy
  data iv /ntab*0/, iy /0/
  if (idum.le.0.or.iy.eq.0) then
      idum=max(-idum,1)
      do j=ntab+8,1,-1
         k=idum/iq
         idum=ia*(idum-k*iq)-ir*k
         if (idum.lt.0) idum=idum+im
            if (j.le.ntab) iv(j)=idum
      enddo
      iy=iv(1)
  endif
  k=idum/iq
  idum=ia*(idum-k*iq)-ir*k
  if (idum.lt.0) idum=idum+im
  j=1+iy/ndiv
  iy=iv(j)
  iv(j)=idum
  ran1=min(am*iy,rnmx)
  return
  end
