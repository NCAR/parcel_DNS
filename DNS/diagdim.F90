
  SUBROUTINE spec (zx,zy,zz,tt,qv,nm,spz,l,time,kxa,kya,kza,vvk,jmax,dimk)

! CALCULATES SPECTRA.
  use dyn_mod
  implicit none
  include 'param.inc'
  include 'mpif.h'
     
  ! ---- argument ---
  integer :: ikx,iky,ikz,ikza,i,j,L(iktx,ikty,iktzp),nspz,nspz4
  complex, dimension(iktx,ikty,iktzp) :: zx,zy,zz,tt,qv
  real :: kx,ky,kz,wk,vz,e,z,t,q,avz,avv,avt,avq,time,kxa(iktx),kya(ikty),kza(iktz)
  real, dimension(ktx,4)  :: spz
  integer, dimension(ktx) :: nm
  
  ! ---- local ----
  real :: VMAX,VMAXT,JMAX,VMAXTO,TAYMIC,REYTAY,INTE,VVK,REYINT,EDDYT,DIMK  
  real :: rmsvel,dissrate,taymic2,reytay2,kmaxeta,tmp
    

  real, allocatable, dimension(:,:)  :: spztot
  integer, allocatable, dimension(:) :: ntot

  integer :: mype,ierror
  common/mpi/mype 
  

  allocate( spztot(ktx,4) )
  allocate( ntot(ktx) )

  nspz  = ktx
  nspz4 = ktx*4

  avz = 0.0
  avv = 0.0
  avt = 0.0
  avq = 0.0

  spz = 0.
  nm = 0. 

  vmaxto = 0.
  DO 20 ikz = 1, iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     DO 20 IKY=1,IKTY
        KY = KYA(IKY)
        DO 20 IKX=1,IKTX
           IF (L(IKX,IKY,IKZ).NE.1)    GO TO 20
           KX = KXA(IKX)
           WK = SQRT(KX*KX+KY*KY+KZ*KZ)
           
           J = INT(WK/DIMK+0.5)
           IF (J.LE.0 .OR. J.GT.KTX) THEN
              PRINT*,'SPEC: SCREW-UP.   k= ',j
           ENDIF
           if (iturb .eq. 1) then
             VZ = real( ZX(IKX,IKY,IKZ)*CONJG(ZX(IKX,IKY,IKZ)) ) ! in x dir
             SPZ(J,1) = SPZ(J,1) + VZ/WK**2	! KE = U^2
             SPZ(J,2) = SPZ(J,2) + VZ		! enstrophy Vort^2
             VMAXTO = VMAXTO + VZ
             VZ = real( ZY(IKX,IKY,IKZ)*CONJG(ZY(IKX,IKY,IKZ)) ) ! in y dir
             SPZ(J,1) = SPZ(J,1) + VZ/WK**2
             SPZ(J,2) = SPZ(J,2) + VZ
             VMAXTO = VMAXTO + VZ
             VZ = real( ZZ(IKX,IKY,IKZ)*CONJG(ZZ(IKX,IKY,IKZ)) ) ! in z dir
             SPZ(J,1) = SPZ(J,1) + VZ/WK**2
             SPZ(J,2) = SPZ(J,2) + VZ
             VMAXTO = VMAXTO + VZ
           endif !iturb

           if (thermo .eq. 1) then
              VZ = real( TT(IKX,IKY,IKZ)*CONJG(TT(IKX,IKY,IKZ)) )
              SPZ(J,3) = SPZ(J,3) + VZ          ! TK = T'^2
              VZ = real( qv(IKX,IKY,IKZ)*CONJG(qv(IKX,IKY,IKZ)) )
              SPZ(J,4) = SPZ(J,4) + VZ          ! QK = qv'^2
           endif 
  
            NM(J)   = NM(J) + 2 
20 CONTINUE

    if (mype .eq. 0) then
       PRINT *,'KINETIC ENERGY, ENSTROPHY, POTENTIAL ENERGY & EQUIVALENT FOR QV AT T = ',TIME
       WRITE(80,*) TIME,'  = TIME'
    endif

    VMAX = -100.0
    INTE = 0.0

!    if(iturb .eq. 1) then
      do j=1,ktx-1
        avz=avz+spz(j,1)
	avv=avv+spz(j,2)
	avt=avt+spz(j,3)
        avq=avq+spz(j,4)
!       avz=sum(spz(1:ktx-1,1))
!       avv=sum(spz(1:ktx-1,2))
      enddo
!    endif
!    if(thermo .eq. 1) then
!      do j=1,ktx-1
!	avt=avt+spz(j,3)
!	avq=avq+spz(j,4)
!      enddo
!	avt=sum(spz(1:ktx-1,3))
!	avq=sum(spz(1:ktx-1,4))
 !   endif
!  DO  J=1,KTX-1
!     if (iturb .eq. 1) then
!        E = SPZ(J,1)            ! Kinetic energy, here <U**2>=E
!        Z = SPZ(J,2)            ! Enstrophy, here <V**2>=Z
!        AVZ = AVZ + E
!        AVV = AVV + Z
!     endif
!     if (thermo .eq. 1) then
!        T = SPZ(J,3)
!        Q = SPZ(J,4)
!        avt = avt + T
!        avq = avq + Q
!     endif
!  ENDDO

!if ( MPI == 1 ) then
  if (iturb .eq. 1 .or. thermo .eq. 1) then
    call mpi_reduce(spz,spztot,nspz4,MPI_REAL,   MPI_SUM,0,MPI_COMM_WORLD,ierror)

    call mpi_reduce( nm,  ntot, nspz,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierror)

  endif!iturb&thermo

  if (iturb .eq. 1) then
    call mpi_reduce(avz,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    avz=tmp
    call mpi_reduce(avv,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    avv=tmp
    call mpi_reduce(VMAXTO,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    VMAXTO=tmp
  endif!iturb
  if (thermo .eq. 1) then
     call mpi_reduce(avt,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierror)
     avt=tmp
     call mpi_reduce(avq,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierror)
     avq=tmp
  endif!thermo
!else
!  spztot=spz
!  ntot=nm
!endif!mpi

  if (mype.eq.0) then
!	print*,'ntot',ntot !mark
    if (iturb .ne. 1 .and. thermo .eq. 1) then
     do j=1,ktx-1 
        if (ntot(j) .ne. 0) then
           write(80,5000) float(j),(spztot(j,i),i=1,4),ntot(j)
           write(6,5000) float(j),(spztot(j,i),i=1,4),ntot(j)
        endif
     enddo
    elseif(iturb .eq. 1  ) then
     do j=1,ktx-1 			! Calculate Kolmogorov wavenumber (1/L)
        e = spztot(j,1)				
        vmaxt = (dimk*j)**2*E ! (2pi*h*J/L)**2 * E
        if (vmaxt .gt. vmax) then
           vmax = vmaxt
           jmax = j
!           print*, 'e,vmax,jmax=',e,vmax,jmax
        endif  
        inte = inte + e/(dimk*j) 	! Integral scale

        if (ntot(j) .ne. 0) then
           write(80,5000) float(j),(spztot(j,i),i=1,4),ntot(j)
           write(6,5000) float(j),(spztot(j,i),i=1,4),ntot(j)
        endif
     enddo
     write(80,*) '           '
     call flush(80)
! Calculate integral scale and Taylor microscale
     INTE = 4.*ASIN(1.)*INTE/AVZ     !! DEFINED AS IN VINCENT AND MENEGUZZI(1991) *2*PI (integral(k-1*KE)/interal(KE))
     EDDYT = INTE/SQRT(2*AVZ)
     REYINT = INTE*SQRT(2*AVZ)/visc
     TAYMIC = SQRT(AVZ/VMAXTO)          !! DEFINED AS IN VINCENT AND MENEGUZZI(1991) (integral(KE)/integral(Enstrophy)
      
     REYTAY = TAYMIC*SQRT(2*AVZ)/visc
     VVK = spztot(NINT(JMAX),1)*DIMK*JMAX

     WRITE(6,5010) AVZ,AVV,AVT,AVQ
     WRITE(80,*) '           '
     WRITE(6,*) '           '
     WRITE(6,*) 'APPROX. KOLMOGOROV WAVENUMBER   = ',JMAX,4.*asin(1.)/(DIMK*JMAX)
!     WRITE(6,*) 'ENERGY AT KOLMOGOROV WAVENUMBER = ',spztot(NINT(JMAX),1)
     WRITE(6,*) 'TAYLOR MICROSCALE               = ',TAYMIC
     WRITE(6,*) 'ASSOCIATED REYNOLDS NUMBER      = ',REYTAY
!     WRITE(6,*) 'INTEGRAL SCALE                  = ',INTE
!     WRITE(6,*) 'ASSOCIATED REYNOLDS NUMBER      = ',REYINT
!     WRITE(6,*) 'EDDY TURNOVER TIME              = ',EDDYT

!     WRITE(97,*) 'APPROX. KOLMOGOROV WAVENUMBER   = ',JMAX,4.*asin(1.)/(DIMK*JMAX)
!     WRITE(97,*) 'ENERGY AT KOLMOGOROV WAVENUMBER = ',spztot(NINT(JMAX),1)
!     WRITE(97,*) 'TAYLOR MICROSCALE               = ',TAYMIC
!     WRITE(97,*) 'ASSOCIATED REYNOLDS NUMBER      = ',REYTAY
!     WRITE(97,*) 'INTEGRAL SCALE                  = ',INTE
!     WRITE(97,*) 'ASSOCIATED REYNOLDS NUMBER      = ',REYINT
!     WRITE(97,*) 'EDDY TURNOVER TIME              = ',EDDYT
      
     rmsvel = sqrt(2.*AVZ/3.)
     dissrate = visc*2.*AVV ! 2 niu AVV = eps. AVV = 1/2 * eps/niu
     taymic2 = rmsvel*sqrt(15*visc/dissrate)
     Reytay2 = taymic2*rmsvel/visc
     kmaxeta = dimk*IKTY/2*(visc**2/(2.*AVV))**.25         ! IKTY = N 
     WRITE(98,5020) TIME,rmsvel,dissrate,taymic2,reytay2,eddyt,kmaxeta,VMAXTO,AVV, TAYMIC, REYTAY, sqrt(2*AVZ) 
    endif !iturb
  endif ! mype

  deallocate( spztot, ntot )
  RETURN

5000  FORMAT(1X,F4.0,4X,4(E15.8,1x),10X,I8)
5010  FORMAT(1X, 'TOTALS   ', 4(E14.8,2x), '     ', E14.8)
5020  FORMAT(1X,12(E15.8,1x))

  END 


  SUBROUTINE OUT (ZX,ZY,ZZ,TT,QV,UX,UY,UZ,time,nstop,ns,spz,l,kxa,kya,kza,vvk,kkol,dimk)

! PRINTS STUFF OUT.

!  use fred_mod
  use dyn_mod
  implicit none
  include 'param.inc'
  include 'mpif.h'

  CHARACTER*90 :: AER
  CHARACTER*10 :: BER
  INTEGER :: IKX,IKY,IKZ,ikza,NS(KTX)
  INTEGER :: L(IKTX,IKTY,IKTZP)
  complex, dimension(iktx,ikty,iktzp) :: zx,zy,zz,ux,uy,uz,tt,qv
  REAL :: hh,KE,PE,TS, whatever,FR
  REAL :: E,V,KX,KY,KZ,TIME,WK,SPZ(KTX,4),vz,vh
  REAL :: KXA(IKTX),KYA(IKTY),KZA(IKTZ),VVK,KKOL,DIMK,tmp
  INTEGER :: I,NSTOP!,IO,iout

  integer :: mype,ierror
  common/mpi/mype 


  !IO = IO + 1
  BER = '__________'
  AER = BER//BER//BER//BER//BER//BER//BER//BER//BER


 
   FR = 0.
   whatever=0.
   KE = 0.
   PE = 0.
   hh = 0.
   V  = 0.
   SPZ= 0.
   if ( iturb .eq. 1) then
  CALL PROJ (ZX,ZY,ZZ,L,KXA,KYA,KZA)
  CALL VELO (ZX,ZY,ZZ,UX,UY,UZ,L,KXA,KYA,KZA)
   endif

   DO 40 IKZ=1,IKTZP
      ikza=mype*iktzp+ikz
      KZ = KZA(IKZA)
      DO 40 IKY=1,IKTY
         KY = KYA(IKY)
         DO 40 IKX=1,IKTX
            KX = KXA(IKX)
            WK = KX*KX + KY*KY + KZ*KZ

            IF (L(IKX,IKY,IKZ).NE.1) GO TO 40
            if (iturb .eq. 1) then
             VZ = real( ZX(IKX,IKY,IKZ)*CONJG(ZX(IKX,IKY,IKZ)) )
             KE = KE + VZ/WK
             V  = V + VZ
             FR = FR + VZ
             VZ = real( ZY(IKX,IKY,IKZ)*CONJG(ZY(IKX,IKY,IKZ)) )
             KE = KE + VZ/WK
             V  = V + VZ
             FR = FR + VZ
             VZ = real( ZZ(IKX,IKY,IKZ)*CONJG(ZZ(IKX,IKY,IKZ)) )
             KE = KE + VZ/WK
             V  = V + VZ

             VH = real( ZX(IKX,IKY,IKZ)*CONJG(UX(IKX,IKY,IKZ)) )
             hh  = hh + vh
             VH = real( ZY(IKX,IKY,IKZ)*CONJG(UY(IKX,IKY,IKZ)) )
             hh  = hh + vh
             VH = real( ZZ(IKX,IKY,IKZ)*CONJG(UZ(IKX,IKY,IKZ)) )
             hh  = hh + vh
            endif
            if (thermo .eq. 1) then
               VH = real( TT(IKX,IKY,IKZ)*CONJG(TT(IKX,IKY,IKZ)) )
               PE = PE  + VH
               VH = real( qv(IKX,IKY,IKZ)*CONJG(qv(IKX,IKY,IKZ)) )
               whatever = whatever  + VH
            endif

40 CONTINUE
! V = 2.*V
! hh = 2.*hh

!if ( MPI == 1 ) then
  if (iturb .eq. 1) then
    call mpi_reduce(ke,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    ke=tmp
    call mpi_reduce(v,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    v=tmp
    call mpi_reduce(hh,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    hh=tmp
  endif
  if (thermo .eq. 1) then
    call mpi_reduce(pe,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    pe=tmp
    call mpi_reduce(whatever,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    whatever=tmp
    call mpi_reduce(fr,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    fr=tmp
  endif
!endif!mpi

  if (mype .eq. 0) then  
     IF (expsiv.NE.0. .AND. temgrad.NE.0.) THEN
        PE = expsiv*PE/(temgrad)
        E  = (PE + KE)
        FR = SQRT( FR/(expsiv*temgrad) )
     ELSE
        E = 99999999.
        FR= 99999999.
     ENDIF

!     E1(IO) = TIME
!     E2(IO) = KE
!     E3(IO) = PE    ! serves here as scalar energy for temperature
!     E4(IO) = VVK   ! kolmogorov wavenumber X energy at this wavenumber
!     E5(IO) = V     ! enstrophy
!     E6(IO) = vh!hh    ! helicity
!     E7(IO) = whatever  ! serves here as scalar energy for qv
!     E8(IO) = KKOL  ! kolmogorov wave#
     e1 = time
     e2 = ke
     e3 = pe
     e4 = vvk
     e5 = v
     e6 = vh
     e7 = whatever
     e8 = kkol 
     PRINT*,'    ' 
     PRINT*,'    ' 
     PRINT*,'    ' 
     PRINT*,'________________________________________________________'
     WRITE( 6,5042) AER
     WRITE( 6,5043) 
!     WRITE(6,5044) E1(IO),E2(IO),E3(IO),E4(IO),E5(IO),E6(IO),E7(IO),E8(IO)
     WRITE(6, 5044) e1,e2,e3,e4,e5,e6,e7,e8!TIME, KE, PE, VVK, V, vh, whatever, KKOL
     WRITE(81,5044) e1,e2,e3,e4,e5,e6,e7,e8!TIME, KE, PE, VVK, V, vh, whatever, KKOL
!     IF (IOUT.EQ.1) THEN
!        WRITE( 6,5043)
!        DO 17 I= 1,IO
!17        WRITE( 6,5044) E1(I),E2(I),E3(I),E4(I),E5(I),E6(I),E7(I),E8(I)
!        DO 18 I= 1,IO
!18         WRITE(81,5044) E1(I),E2(I),E3(I),E4(I),E5(I),E6(I),E7(I),E8(I)
!     ENDIF
     Print*,'Froude = ',FR

     PRINT *,'                             '
     PRINT *,'                             '
     PRINT *,'                             '

  endif !mype

  RETURN

5042  FORMAT(1X,A90)
5043  FORMAT(7X,'T',8x,'KE',10x,'PE',10X,'VVK',10x,'V',10x,'H',10x,'S2',6x,'KKOL')
5044 format(1x,(F8.3),2x,3(e10.3,2x),(e10.3),2x,f8.3,2x,e10.3,2x,f4.1)

  END
