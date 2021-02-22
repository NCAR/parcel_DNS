      program parcel_model !parcel model main program
        !aerosol activation for parcel of 1m^3 in volume
      implicit none
      real*8 :: tau,source,qvpp,qvs,PP,esat,ks,temp,dp,time,temp0,upp,time_prep
      real*8 :: rhoa,rhoa0,sp,sp2,vol,avgconc,cql,h,thetapp,cdp,cd,seq,seq1,seq2
      real*8 :: exner,racp,p0,p1,sumrp,edr
      real*8   :: lwc , rm,rm1,rm0,curv,solu,taup
      real*8 :: aa11,aa22,vtemp,vtemp1,delt
      real*8 :: thetap,qvp,deltaqp,dtheta,gamma,gamma1
      character*1 :: name
      real*8, allocatable, dimension(:) :: rad,rad1,dsd,mass,dr3,rad_wet,udrop,fv
      integer :: nbins,nbinsout,iinit,ifinal
      integer :: disp,GCCN
      integer :: iter,ntmic,ntot,i
      real*8 :: ndrop
!      integer, parameter :: NN =64
      real*8 :: diffvnd1,diffvnd2,ka1,ka2,Red
      real, parameter :: diffvnd = 2.55d-5              ! Coefficient of diffusion of water vapour in air [m**2/s]
      real, parameter :: ka = 2.48d-2                   ! Thermal conductivity of air [J/msK]
      real,parameter :: pi=3.1415
      real,parameter ::  KK = 8.54d-11
      real,parameter :: grav=9.8
      real,parameter :: visc = 1.78d-5!0.16d-4
      real,parameter :: lat = 2.5d6!2.477d6
      real,parameter :: up=2.d0
      real,parameter :: ra=287.0
      real,parameter :: cp=1004!1005.0
      real,parameter :: rv= 467!461.5
      real,parameter :: m_w=18.d-3 !molecular weight of water
      real,parameter :: m_s=132.14d-3!molecular weight of solute; ammonium sulfate=132.14d-3; NaCl = 58.44d-3
      real,parameter :: rhow=1000.0
      real,parameter :: eps=ra/rv
      real,parameter :: vh=2. !van hoff factor
      real,parameter :: latovercp=lat/cp
      real,parameter :: rho_ccn = 1726.d0 !kg/m**3 for ammonium sulfate =1726.d0 for NaCl=2160.d0
      real(8), parameter :: Nsc = visc/diffvnd
      real(8),parameter :: kwt = 1.233d8
196   format(1x,6(f16.8,2x))
145   format(1x,100(e16.8,2x))

      time_prep=0.d0
      delt= 1.d-04
      ntot=600./delt
      !setup output files
      OPEN(UNIT=16,FILE='new.out',ACCESS='APPEND')!parcel mean variables
      OPEN(UNIT=50,FILE='new.dsd',ACCESS='APPEND')!number concentration of each bin
      OPEN(UNIT=51,FILE='new.rad',ACCESS='APPEND')!droplet size of each bin
      OPEN(UNIT=60,FILE='new.test',ACCESS='APPEND')!output the tested variables
!--------------initialize aerosols---------
      disp = 30
      GCCN = 0 !=1 insert Giant CCN =2 monotonic seeding r=1micron;!=3 add 3 mode lognormal seeding distribution by Cooper et al. 1997
      !disp=35 Xue10 urban, 30 Xue10 marine, 31 JN17 polluted, 32 NJ17 pristine
      nbins = 100
      allocate (dsd(nbins),rad(nbins),rad1(nbins),mass(nbins),dr3(nbins),rad_wet(nbins))
      allocate (udrop(nbins),fv(nbins))
      rm =0.d0!1.d-5
      rad=rm
      dsd=0.d0
      dr3=0.d0
!     get the initial aerosol size spectrum
      call iaerosol(disp,rad,mass,dsd,nbins,ndrop,iinit,ifinal,rm,nbinsout,GCCN)
      print*,'nbinsout,GCCN',nbinsout,GCCN
      write(50,145) 0.,(dsd(i), i=1,nbinsout)
      write(51,145) 0.,(rad(i), i=1,nbinsout)
      print*,'total droplet number',ndrop,sum(dsd(1:nbinsout)),'dispersion parameter',disp
!-------------initialize variables------------
      rad1=rad
      h = 1.d-2 !.01m=1cm
      vol = h**3
      cd=sum(rad(1:nbinsout)**3*dsd(1:nbinsout))*4.d0/3.d0*pi*rhow/vol!condensation
      temp=284.3d0!281.5d0!284.3!293.15
      p1=93850.0d0!90650.d0!93850.0!95000.0
      p0=1.d5
      sp = -14.39d-2!0.d0!-14.39d-2
      lwc=0.d0
      pp=p1
      sumrp=0.d0!mark new cdp sumrp=sum(r**2*dr)
      racp= ra/cp
      exner=(pp/p0)**racp
      thetapp=temp/exner
      cql=4.0d0*pi*rhow/vol !!mark new cql(rhoa*vol)
      esat = 2.53d11*exp(-5.42d3/temp)
!      ks=1.0/(lat**2*eps*rhow/(Ka*Ra*Temp**2)+Ra*Temp*rhow/(eps*diffvnd*esat))
      ks = 1.d0/(rhow*Rv*temp/(esat*diffvnd)+rhow*Lat/(Ka*Temp)*(Lat/(Rv*Temp)-1))
      qvs = eps*esat/(PP-esat)
      qvpp= (sp+1.d0)*qvs !kg/m^3
      rhoa=pp/(Ra*(1+18.d0/29.d0*qvpp)*temp)
      rhoa0=rhoa
      write(16,*) 0.d0-time_prep,up*delt*ntmic,Sp,ndrop,pp,temp,& !8
                thetapp,qvpp,qvs,lwc,rhoa,deltaqp
!--------------first guess the radius of wet aerosol--------!
!--------------start with r_wet=1.5*r_d-----------------------!
      rad_wet=rad1*1.5d0
      diffvnd1=1.d-5*(0.015*Temp-1.9)
      ka1=1.5d-11*temp**3-4.8d-8*temp**2+1.d-4*temp-3.9d-4
      do i=1,nbinsout
         seq=sp+.01d0
	 seq1=seq+.01d0
	 seq2=seq+.01d0
      do while(abs(seq-sp) .gt. 1.d-7 .and. seq2 .ne. seq) 
	 seq2=seq1
	 seq1=seq
         diffvnd2=diffvnd1*1.d0/(rad_wet(i)/(rad_wet(i)+0.104d-6)+diffvnd1/(rad_wet(i)*0.036)*sqrt(2.d0*pi/(Ra*temp)))
         ka2=ka1*1.d0/(rad_wet(i)/(rad_wet(i)+.216d-6)+ka1/(rad_wet(i)*.7*rhoa*cp)*sqrt(2.d0*pi/(Ra*temp)))
         ks = 1.d0/(rhow*Rv*temp/(esat*diffvnd2)+rhow*Lat/(Ka2*Temp)*(Lat/(Rv*Temp)-1))
         curv=2.d0*7.61d-2/(Rv*rhow*temp*rad_wet(i))  
         solu=vh*m_w/m_s*rho_ccn*rad1(i)**3/rhow!4.3d0*2.d0*mass(i)/132.14d0*1.d-6 !solute effect coefficient ms=132.14 for ammonium sulfate !m^3
         seq=exp(curv-solu/(rad_wet(i)**3-rad1(i)**3))-1.d0
         if(seq .gt. sp) then
	   rad_wet(i)=rad_wet(i)-rad1(i)*1.d-6
	 elseif(seq .lt. sp) then 
	   rad_wet(i)=rad_wet(i)+rad1(i)*1.d-6
	 endif
      enddo
      enddo
      rad=rad_wet
!--------------spin-up------------------
      upp=0.d0
      do 200 ntmic=1,int(time_prep/delt*2.d0)
         time = ntmic*delt/2.d0-time_prep
         pp=rhoa*Ra*(1.d0+18.d0/29.d0*qvpp)*temp
         exner = (PP/P0)**RACP
         sumrp=sum(dr3(1:nbinsout)*dsd(1:nbinsout))/3.d0
         deltaqp=cql*sumrp/rhoa!new cdp condensation
         cd=sum(rad(1:nbinsout)**3*dsd(1:nbinsout))*4.d0/3.d0*pi*rhow/vol/rhoa!!mark
         temp0=temp
         temp=temp0-grav/cp*delt/2.0d0*upp+latovercp*deltaqp!mark new
         esat=2.53d11*exp(-5.42d3/temp)
         qvs = eps*esat/(PP-esat)
         rhoa= rhoa*(-grav*upp/(Ra*temp)*delt/2.d0-(temp-temp0)/temp)+rhoa
          thetapp=thetapp+latovercp*deltaqp/exner
          qvpp    = qvpp - deltaqp
          sp = qvpp/qvs-1.d0 !mark new sp
          rm=0.d0
          diffvnd1=1.d-5*(0.015*Temp-1.9)
          ka1=1.5d-11*temp**3-4.8d-8*temp**2+1.d-4*temp-3.9d-4
          do i = 1,nbinsout
             curv=2.d0*7.61d-2/(Rv*rhow*temp*rad(i)) !curvature effect coefficient !in m
             solu=vh*m_w/m_s*rho_ccn*rad1(i)**3/rhow!4.3d0*2.d0*mass(i)/132.14d0*1.d-6 !solute effect coefficient ms=132.14 for ammonium sulfate !m^3
             diffvnd2=diffvnd1*1.d0/(rad(i)/(rad(i)+0.104d-6)+diffvnd1/(rad(i)*0.036)*sqrt(2.d0*pi/(Ra*temp)))
             ka2=ka1*1.d0/(rad(i)/(rad(i)+.216d-6)+ka1/(rad(i)*.7*rhoa*cp)*sqrt(2.d0*pi/(Ra*temp)))
             ks =1.d0/(rhow*Rv*temp/(esat*diffvnd2)+rhow*Lat/(Ka2*Temp)*(Lat/(Rv*Temp)-1))
!!caculate equillibrium supersat.
             seq=exp(curv-solu/(rad(i)**3-rad1(i)**3))-1.d0
             rm0=rad(i)*3.0d0*delt/2.d0*ks*(sp-seq)+rad(i)**3 !!r^3 scheme
             if(rm0 .gt. rad1(i)**3) then
                dr3(i)=rm0-rad(i)**3
                rad(i)=rm0**(1.d0/3.d0)
             else
                dr3(i)=0.d0
                rad(i)=rad1(i)
             endif
          enddo
          rm=(sum(rad(1:nbinsout)**3*dsd(1:nbinsout))/ndrop)**(1.d0/3.0d0)
          if(mod(ntmic,int(time_prep/delt*2.d0)/1000) .eq. 0 .or. int(time_prep/delt*2.d0) .le. 1000) then
             write(16,*) time,0,Sp,ndrop,pp,temp,thetapp,qvpp,qvs,lwc,rhoa,deltaqp
             write(50,145) time,(dsd(i), i=1,nbinsout)
             write(51,145) time,(rad(i), i=1,nbinsout)
	  endif
  200 enddo
!--------------evolution------------------
      do 100 ntmic=1,ntot
         time = ntmic*delt
!         dp = rhoa*grav*Up*delt
!         PP = PP-DP
         pp=rhoa*Ra*(1+m_w/29.d-3*qvpp)*temp
         exner = (PP/P0)**RACP         
!         sumrp = sum(rad(1:nbinsout)*dsd(1:nbinsout))!mark old cdp
         sumrp=sum(dr3(1:nbinsout)*dsd(1:nbinsout))/3.d0
         do i = 1,nbinsout
             taup=kwt*rad(i)**2/grav
             udrop(i)=taup*grav
             Red=2.d0*rad(i)*udrop(i)/visc
             if(Red**(1.d0/2.d0)*Nsc**(1.d0/3.d0) .lt. 1.4) then
                 fv(i) = 1.0+0.108*(Nsc**(1.d0/3.d0)*sqrt(Red))**2
             else
                 fv(i) = .78 + .308*Nsc**(1.d0/3.d0)*sqrt(Red)
             endif !ventilation effect coefficient
         enddo
         if(ntmic.eq.ntot) print*,fv
         fv(i)=1.d0
	 deltaqp=cql*sumrp!new cdp condensation
!         print*,'cdold,new',(sum(rad(1:nbinsout)**3*dsd(1:nbinsout))*4.d0/3.d0*pi*rhow)/vol-cd,deltaqp
         cd=sum(rad(1:nbinsout)**3*dsd(1:nbinsout))*4.d0/3.d0*pi*rhow/vol!!mark
!         cdp=delt*cql*ks*sumrp!mark old sp
         temp0=temp
         temp=temp0-grav/cp*delt*up+latovercp*deltaqp!mark new
         esat=2.53d11*exp(-5.42d3/temp)
         qvs = eps*esat/(PP-esat)
!          rhoa=pp/(Ra*temp)
         rhoa= rhoa*(-grav*up/(Ra*temp)*delt-(temp-temp0)/temp)+rhoa
          thetapp=thetapp+latovercp*deltaqp/exner
          qvpp    = qvpp - deltaqp
	  sp = qvpp/qvs-1.d0 !mark new sp
          rm=0.d0
          diffvnd1=1.d-5*(0.015*Temp-1.9)
          ka1=1.5d-11*temp**3-4.8d-8*temp**2+1.d-4*temp-3.9d-4
          do i = 1,nbinsout
!             curv=3.3d-7/(temp*rad(i))!
             curv=2.d0*7.61d-2/(Rv*rhow*temp*rad(i)) !curvature effect coefficient !in m
             solu=vh*m_w/m_s*rho_ccn*rad1(i)**3/rhow !solute effect coefficient ms=132.14 for ammonium sulfate !m^3
             diffvnd2=diffvnd1*1.d0/(rad(i)/(rad(i)+0.104d-6)+diffvnd1/(rad(i)*0.036)*sqrt(2.d0*pi/(Ra*temp)))
             ka2=ka1*1.d0/(rad(i)/(rad(i)+.216d-6)+ka1/(rad(i)*.7*rhoa*cp)*sqrt(2.d0*pi/(Ra*temp)))
             ks =1.d0/(rhow*Rv*temp/(esat*diffvnd2)+rhow*Lat/(Ka2*Temp)*(Lat/(Rv*Temp)-1))
!!--------------caculate equillibrium supersat.

             seq=exp(curv-solu/(rad(i)**3-rad1(i)**3))-1.d0
             rm0=rad(i)*3.0d0*delt*ks*(sp-seq)*fv(i)+rad(i)**3 !!r^3 scheme    

	      if(rm0 .gt. rad1(i)**3) then
	        dr3(i)=rm0-rad(i)**3
	        rad(i)=rm0**(1.d0/3.d0)
	      else
		dr3(i)=0.d0
                rad(i)=rad1(i)
	      endif
          enddo 
            lwc=sum(rad(1:nbinsout)**3*dsd(1:nbinsout))*4.d0/3.d0*pi*rhow/vol/rhoa
            rm=(sum(rad(1:nbinsout)**3*dsd(1:nbinsout))/ndrop)**(1.d0/3.0d0)
            if(mod(ntmic,1000) .eq. 0) then
              write(16,*) time,up*delt*ntmic,Sp,ndrop,pp,temp,&
                thetapp,qvpp,qvs,lwc,rhoa,deltaqp
              write(50,145) time,(dsd(i), i=1,nbinsout)
              write(51,145) time,(rad(i), i=1,nbinsout)
            endif
  100 enddo

      close(unit=16)
      deallocate (dsd, rad, rad1,mass,dr3,rad_wet)
      deallocate (udrop,fv)
      close(unit=50)


      end program parcel_model

    SUBROUTINE IAEROSOL(disp,rad,mass,nrad,nbins,ndrop,iinit,ifinal,rm,nbinsout,GCCN)
!---- This subroutine determines the initial size of aerosols at each size bin.
  implicit none

  ! --- argument ---
  integer :: nbins,nbinsout,GCCN
  real(8), dimension(nbins) :: rad,nrad,mass
  integer :: disp,iinit,ifinal
  real(8) :: rm,ndrop

  ! --- local --

  integer :: i
  real(8), allocatable, dimension(:) :: wid,dNdlogr,dNdr
  real(8) :: r1,n1,logsig,rmin,rmax,rho_ccn
  real(8) :: logrmin,logrmax,rad_power,bin_factor
  real,parameter :: pi=3.14159265
 145   format(1x,100(e16.8,2x))
 111 format(a20,i3)

! Set everything to zero.
rad = 0.0
ndrop=0.0 
rm=0.d0
nrad=0.
nbinsout=nbins
allocate(wid(nbins),dNdlogr(nbins),dNdr(nbins))
dNdlogr = 0.0
dNdr = 0.0d0
wid =0.d0
mass=0.d0
rho_ccn=1726 !kg/m**3 for Ammonium sulfate
!size dispersion
  if (disp .eq. 30) then !lulin maritime case
     nbinsout=39
     rmin = 6.d-9
     rad(1)=rmin
     mass(1)=4.d0/3.d0*pi*rmin**3*rho_ccn 
     bin_factor=2.d0
     wid(1)=rad(1)*(bin_factor**(1.d0/3.d0)-1.d0)
     do i=2,nbinsout
        mass(i)=mass(1)*bin_factor**i
        rad_power=real(i)/3.d0
        rad(i)=rad(1)*bin_factor**rad_power
        wid(i)=rad(i)-rad(i-1)
     enddo
     do i=1,nbinsout
        n1=133.d0
        r1=0.0039d-6
        logsig=.657d0
        logsig=log(10.d0)*logsig
        dNdlogr(i) = n1/(sqrt(2.0d0*pi) *logsig) * exp(-((log(rad(i))-log(r1))/(sqrt(2.0d0)*logsig))**2)
        n1=66.6d0
        r1=.133d-6
        logsig=.21d0
        logsig=log(10.d0)*logsig
        dNdlogr(i)= dNdlogr(i)+n1/(sqrt(2.0d0*pi) *logsig) * exp(-((log(rad(i))-log(r1))/(sqrt(2.0d0)*logsig))**2)
        n1=3.06d0
        r1=.29d-6
        logsig=.396d0
        logsig=log(10.d0)*logsig
        dNdlogr(i)= dNdlogr(i)+n1/(sqrt(2.0d0*pi) *logsig) * exp(-((log(rad(i))-log(r1))/(sqrt(2.0d0)*logsig))**2)
	if (GCCN .eq. 3) then !add 3 mode lognormal seeding distribution by Cooper et al. 1997
	   !mode 1
	   n1=100.d0
	   r1=.15d-6
	   logsig=.2d0
	   logsig=log(10.d0)*logsig
           dNdlogr(i)= dNdlogr(i)+n1/(sqrt(2.0d0*pi) *logsig) * exp(-((log(rad(i))-log(r1))/(sqrt(2.0d0)*logsig))**2)
	   !mode 2
           n1=100.d0*1.7d-4
           r1=.5d-6
           logsig=.4d0
           logsig=log(10.d0)*logsig
           dNdlogr(i)= dNdlogr(i)+n1/(sqrt(2.0d0*pi) *logsig) * exp(-((log(rad(i))-log(r1))/(sqrt(2.0d0)*logsig))**2)
	   !mode 3
           n1=100.d0*3.d-7
           r1=5.d-6
           logsig=.6d0
           logsig=log(10.d0)*logsig
           dNdlogr(i)= dNdlogr(i)+n1/(sqrt(2.0d0*pi) *logsig) * exp(-((log(rad(i))-log(r1))/(sqrt(2.0d0)*logsig))**2)
	endif!GCCN=3
     enddo
     do i=1,nbinsout
        rm=rm+rad(i)**3*dNdlogr(i)*wid(i)/rad(i)
        ndrop=ndrop+dNdlogr(i)*wid(i)/rad(i)      
        dNdr(i)=dNdlogr(i)/rad(i)
        nrad(i)=dNdr(i)*wid(i)
     enddo
  elseif (disp .eq. 31) then !Jensen&Nugent 2017 modified polluted case
     nbinsout=30
     rmin = 2.0514d-9
     rad(1)=rmin
     rho_ccn=1726.d0!ammonium sulfate
     mass(1)=4.d0/3.d0*pi*rmin**3*rho_ccn
     bin_factor=2.0d0 !mass increment
     wid(1)=log(rad(1))-log(rad(1)/bin_factor**(1.d0/3.d0))
     do i=2,nbinsout
        mass(i)=mass(1)*bin_factor**i
        rad_power=real(i)/3.d0
        rad(i)=rad(1)*bin_factor**rad_power
        wid(i)=log(rad(i))-log(rad(i-1))
        n1=48.d0!160!48.d0
        r1=.029d-6
        logsig=1.36d0
        logsig=log(logsig)
        dNdlogr(i)= dNdlogr(i)+n1/(sqrt(2.0d0*pi) *logsig) *exp(-((log(rad(i))-log(r1))/(sqrt(2.0d0)*logsig))**2)
        n1=125.d0!380!125.d0
        r1=.071d-6
        logsig=1.57d0
        logsig=log(logsig)
        dNdlogr(i)= dNdlogr(i)+n1/(sqrt(2.0d0*pi) *logsig)*exp(-((log(rad(i))-log(r1))/(sqrt(2.0d0)*logsig))**2)
     enddo
     dNdr(1:nbinsout)=dNdlogr(1:nbinsout)/rad(1:nbinsout)
     nrad(1:nbinsout)=dNdlogr(1:nbinsout)*wid(1:nbinsout)
     ndrop=sum(nrad(1:nbinsout))
     rm = sum(rad(1:nbinsout)**3*nrad(1:nbinsout))


  elseif (disp .eq. 32) then !Jensen&Nugent 2017 pristine case
     nbinsout=30
     rmin = 2.1d-9
     rad(1)=rmin
     rho_ccn=1726.d0!ammonium sulfate
     mass(1)=4.d0/3.d0*pi*rmin**3*rho_ccn
     bin_factor=2.0d0 !mass increment
     wid(1)=log(rad(1))-log(rad(1)/bin_factor**(1.d0/3.d0))
     do i=2,nbinsout
        mass(i)=mass(1)*bin_factor**i
        rad_power=real(i)/3.d0
        rad(i)=rad(1)*bin_factor**rad_power
        wid(i)=log(rad(i))-log(rad(i-1))
        n1=125.d0
        r1=.011d-6
        logsig=1.2d0
        logsig=log(logsig)
        dNdlogr(i)= dNdlogr(i)+n1/(sqrt(2.0d0*pi) *logsig)*exp(-((log(rad(i))-log(r1))/(sqrt(2.0d0)*logsig))**2)
        n1=65.d0
        r1=.06d-6
        logsig=1.7d0
        logsig=log(logsig)
        dNdlogr(i)= dNdlogr(i)+n1/(sqrt(2.0d0*pi) *logsig)*exp(-((log(rad(i))-log(r1))/(sqrt(2.0d0)*logsig))**2)
     enddo
        dNdr(1:nbinsout) = dNdlogr(1:nbinsout)/rad(1:nbinsout)
        nrad(1:nbinsout) = dNdlogr(1:nbinsout)*wid(1:nbinsout)
        ndrop=sum(nrad(1:nbinsout))
        rm = sum(rad(1:nbinsout)**3*nrad(1:nbinsout))
  elseif (disp .eq. 35) then !Lulin 2010 rural
     rmin = 6.d-9
     rad(1)=rmin
     rho_ccn=1726.d0
     mass(1)=4.d0/3.d0*pi*rmin**3*rho_ccn
     bin_factor=2.d0 !mass increment
     wid(1)=rad(1)*(bin_factor**(1.d0/3.0d0)-1)
     do i=2,35
        mass(i)=mass(1)*bin_factor**i
        rad_power=real(i)/3.0
        rad(i)=rad(1)*bin_factor**rad_power
        wid(i)=rad(i)-rad(i-1)
        n1=6650.d0
        r1= 0.00739d-6
        logsig= .225d0
        logsig=log(10.d0)*logsig
        dNdlogr(i) = n1/(sqrt(2.0d0*pi) *logsig) * exp(-((log(rad(i))-log(r1))/(sqrt(2.0d0)*logsig))**2)
        n1= 147.d0    
        r1= .0269d-6
        logsig= .557
        logsig=log(10.d0)*logsig
        dNdlogr(i) = dNdlogr(i)+n1/(sqrt(2.0d0*pi) *logsig) * exp(-((log(rad(i))-log(r1))/(sqrt(2.0d0)*logsig))**2)
        n1= 1990.d0
        r1=.0419d-6
        logsig= .266
        logsig=log(10.d0)*logsig
        dNdlogr(i) = dNdlogr(i)+n1/(sqrt(2.0d0*pi) *logsig) * exp(-((log(rad(i))-log(r1))/(sqrt(2.0d0)*logsig))**2)
     enddo
     do i=1,35
        rm=rm+rad(i)**3*dNdlogr(i)*wid(i)/rad(i)
        ndrop=ndrop+dNdlogr(i)*wid(i)/rad(i)
        dNdr(i)=dNdlogr(i)/rad(i)
        nrad(i)=dNdr(i)*wid(i)
     enddo
 
  elseif(disp .eq. 39 ) then !Jaenicke1988
     rmin = 6.d-9
     rmax = 5.d-6
     logrmin=10.d0**floor(log10(rmin))
     logrmax=10.d0**floor(log10(rmax))
     iinit=Nint(rmin/logrmin)
     ifinal=(floor(log10(rmax))-floor(log10(rmin)))*9+floor(rmax/logrmax)
     r1=0.
     n1=0.
     rm=0.
     logsig=0.
     do i=iinit,ifinal
          rad(i) = real(mod(i,9))*10.d0**(-7+i/9)
          if (mod(i,9) .eq. 0) then
             rad(i) = 9.0d0*10.d0**(-7+i/9-1)
          endif
          n1=133.d0
          r1=0.0039d-6
          logsig=.657d0
          dNdlogr(i) = n1/(sqrt(2.0d0*pi) *logsig) * exp(-((log10(rad(i))-log10(r1))/(sqrt(2.0d0)*logsig))**2)
          n1=66.6d0
          r1=.133d-6
          logsig=.21d0
          dNdlogr(i)= dNdlogr(i)+n1/(sqrt(2.0d0*pi) *logsig) * exp(-((log10(rad(i))-log10(r1))/(sqrt(2.0d0)*logsig))**2)
          n1=3.06d0
          r1=.29d-6
          logsig=.396d0
          dNdlogr(i)= dNdlogr(i)+n1/(sqrt(2.0d0*pi) *logsig) * exp(-((log10(rad(i))-log10(r1))/(sqrt(2.0d0)*logsig))**2)
          rm=rm+rad(i)**3*dNdlogr(i)*10.0d0**(floor(log10(rad(i))))/rad(i)/log(10.d0)
          ndrop=ndrop+dNdlogr(i)*10.0d0**(floor(log10(rad(i))))/rad(i)/log(10.d0)
     enddo

  elseif (disp .eq. 37) then !Jaenicke1988 maritime case !6d-7~5d-4cm
     !r<1d-5
     rmin = 6.d-7
     rmax = 5.d-4
     logrmin=10.d0**floor(log10(rmin))
     logrmax=10.d0**floor(log10(rmax))
     iinit=Nint(rmin/logrmin)
     ifinal=(floor(log10(rmax))-floor(log10(rmin)))*9+floor(rmax/logrmax)
     r1=0.
     n1=0.
     rm=0.
     logsig=0.
     do i=iinit,ifinal
          rad(i) = real(mod(i,9))*10.d0**(-7+i/9)
          if (mod(i,9) .eq. 0) then
             rad(i) = 9.0d0*10.d0**(-7+i/9-1)
          endif
          if(rad(i) .lt. 1.d-5) then
             n1 = 1.33d2
             r1 = 3.9d-7
             logsig = .657d0
             dNdlogr(i) = n1/(sqrt(2.0d0*pi) *logsig) * exp(-((log10(rad(i))-log10(r1))/(sqrt(2.0d0)*logsig))**2)
          elseif(rad(i) .lt. 1.d-4) then
             n1 = 6.66d1
             r1 = 1.33d-5
             logsig = .21d0
             dNdlogr(i) = n1/(sqrt(2.0*pi) *logsig) * exp(-((log10(rad(i))-log10(r1))/(sqrt(2.0d0)*logsig))**2)
             r1 = 2.9d-5
             logsig = .396d0
             dNdlogr(i) = dNdlogr(i) +n1/(sqrt(2.0*pi) *logsig) * exp(-((log10(rad(i))-log10(r1))/(sqrt(2.0d0)*logsig))**2)
          else !r>1.d-4
             n1=2.84d-1
             r1=2.0d-4
             logsig = .5d0
             dNdlogr(i) = n1/(sqrt(2.0*pi) *logsig) * exp(-((log10(rad(i))-log10(r1))/(sqrt(2.0d0)*logsig))**2)
          endif
          rm=rm+rad(i)**3*dNdlogr(i)*10.0d0**(floor(log10(rad(i))))/rad(i)/log(10.d0)
          ndrop=ndrop+dNdlogr(i)*10.0d0**(floor(log10(rad(i))))/rad(i)/log(10.d0)

     enddo

  endif ! disp
  if(GCCN .eq. 1) then
     do i=1,42
        rad(nbinsout+i)=.8d-6+0.2d-6*real(i-1)
        mass(nbinsout+i)=4.d0/3.d0*pi*rad(nbinsout+i)**3*rho_ccn
     enddo
        nrad(nbinsout+1)=0.1118d0
        nrad(nbinsout+2)=.06849d0!1micron
        nrad(nbinsout+3)=.0384d0
        nrad(nbinsout+4)=.02182d0
        nrad(nbinsout+5)=.0133d0
        nrad(nbinsout+6)=.8496d-2
        nrad(nbinsout+7)=.5486d-2!2micron
        nrad(nbinsout+8)=.3805d-2
        nrad(nbinsout+9)=.2593d-2
        nrad(nbinsout+10)=.1919d-2
        nrad(nbinsout+11)=.1278d-2
        nrad(nbinsout+12)=.9884d-3!3micron
        nrad(nbinsout+13)=.7779d-3
        nrad(nbinsout+14)=.5195d-3
        nrad(nbinsout+15)=.4005d-3
        nrad(nbinsout+16)=.3769d-3
        nrad(nbinsout+17)=.2653d-3!4micron
        nrad(nbinsout+18)=.2124d-3
        nrad(nbinsout+19)=.1378d-3
        nrad(nbinsout+20)=.1214d-3
        nrad(nbinsout+21)=.1009d-3
        nrad(nbinsout+22)=.1222d-3!5micron
        nrad(nbinsout+23)=.5064d-4
        nrad(nbinsout+24)=.383d-4
        nrad(nbinsout+25)=.5547d-4
        nrad(nbinsout+26)=.2145d-4
        nrad(nbinsout+27)=.1295d-4!6micron
        nrad(nbinsout+28)=.4323d-4
        nrad(nbinsout+29)=.2626d-4
        nrad(nbinsout+30)=.305d-4
        nrad(nbinsout+31)=.4385d-5
        nrad(nbinsout+32)=.4372d-5!7micron
        nrad(nbinsout+33)=.4465d-5
        nrad(nbinsout+34)=.4395d-5
        nrad(nbinsout+35)=.4427d-5
        nrad(nbinsout+36)=.4411d-5
        nrad(nbinsout+37)=.0d0    !8micron
        nrad(nbinsout+38)=.0d0
        nrad(nbinsout+39)=.0d0
        nrad(nbinsout+40)=.4522d-5
        nrad(nbinsout+41)=.0d0
        nrad(nbinsout+42)=.4542d-5!9micron
        nbinsout=nbinsout+42
  elseif (GCCN==2) then
        rad(nbinsout+1)=1.d-6
	mass=4.d0/3.d0*pi*rad(nbinsout+1)**3*rho_ccn
        nrad(nbinsout+1)=10
        nbinsout=nbinsout+1
  endif!GCCN
  ndrop=sum(nrad(1:nbinsout))
  write(60,145) 0.,(dNdlogr(i), i=1,nbinsout)
  rm = (sum(rad(1:nbinsout)**3*nrad(1:nbinsout))/ndrop)**(1.d0/3.d0)

  print*,"rm", rm,ndrop,10.d0**(floor(log10(rad(7))))
deallocate(wid,dNdlogr,dNdr)

 299   format(1x,(e12.6),3(f8.6))
 199   format(1x,3(f8.6))


  end SUBROUTINE IAEROSOL

real function log2(x)
  implicit none
  real(8), intent(in) :: x

  log2 = log(x) / log(2.d0)
end function
