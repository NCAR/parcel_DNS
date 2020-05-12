  subroutine micro(uflow,vflow,wflow,ttr,qvr,time0,ndumpd,ndropreal,nt,nddone)

! Main program for condensational growth of a spectrum of 
! droplets randomly distributed in 3D space and sedimenting. 
! All variables ending with a p refer to parcel variables 
! Variables with a double pp are parcel variables at a later timestep 
! Grid variables as theta or qv are perturbation variables 
! AS OF MARCH 19, 1996: THETA IS TEMPERATURE PERTURBATION(NOT POT. TEMP.) 
! PARCEL MODEL STILL WRITTEN WITH POT. TEMP. 
! TO BE USED WITH MAIN THAT CONTAINS FORCING DUE TO W PERTURBATIONS
!      use fred_mod
      use dyn_mod
      use thermo_mod
      use hydro_mod
      use mic_mod
      implicit none

    ! cell structure
      interface
         integer function ICELL(jx,jy,jz,m,my)
            implicit none 
!	    intent: variable in a function that can be passed in & out. intent(in) can be entered not changed.
            integer, intent(in) :: jx,jy,jz,m,my
         end function ICELL

	 real function rangen(idum)
	    integer, intent(in) ::idum
	 end function rangen
      end interface

      include 'param.inc'
      include 'mpif.h'

      ! argument
      real, dimension(n1+1,n3+1,n2pe+1) :: uflow,vflow,wflow
      real, dimension(n1d,n3d,n2dp) :: ttr,qvr
      real, dimension(n1,n3,n2pe)    :: s,sumrordeltaq
      real :: time0,time
      integer :: ndropreal, nt, nddone!, nrstd

      !local
      integer :: iter
      integer :: seed
      integer :: ymaxid,yminid,tmpi
      integer :: ndumpd,ntdump,id
      integer :: nsendl,nsendr,nrecvl,nrecvr,nto,nfrom,ntofrom(npe+2),naway
      integer :: nbufsendr,nbufsendl,nbufrecvr,nbufrecvl,ndroptot
      integer, allocatable, dimension(:) :: nawayid
      integer :: ikx,iky,ikz
      integer, parameter :: nawaymax = 2*ndropp

      real :: red,dist,dxp,dyp,dzp,tmpv,rtmp,ltmp
      real(8) :: cql,vtemp,tempd,solu,curv,lwctmp,tmp
      real :: taup,wr,Vs,Vbox,yrelmax,yrelmin
      real, allocatable, dimension(:,:) :: infosendl,infosendr,inforecvl,inforecvr

      real :: onepmin,onepmax,twopmin,twopmax
      real :: fourpmin, fourpmax, tenpmin, tenpmax, twenpmin, twenpmax
      real :: vmax, maxvel 

      integer :: status(MPI_STATUS_SIZE)
      integer :: mype,ierror
      common/mpi/mype

    ! For plotting droplet's trajectory

      integer :: ifol,sendfol,idpfoli,ndumpdrop
      integer, allocatable, dimension(:) :: idpfols
      real, allocatable, dimension(:) :: xfols,yfols,zfols
      real(8), allocatable, dimension(:) :: rfols,rfolr !,ufols,vfols,wfols
      integer, allocatable, dimension(:) :: idpfolr
      real, allocatable, dimension(:) :: xfolr,yfolr,zfolr !,ufolr,vfolr,wfolr

    ! for the linked lists and collision detection

      integer, parameter :: maxcol = 10
      integer, parameter :: ncele = 23
      integer, parameter :: ngrele = 6
      real, allocatable, dimension(:,:) :: ncdat
      real, allocatable, dimension(:,:) :: ngrdat
      integer :: nc, nctot(npe), ofs, numgr, numgrtot(npe)
      integer :: imap,jx,jy,jz,i,tcell,k,l,jcell0,jcell,j,nabor
      integer :: xcell,ycell,zcell
      integer :: colpars,colparso,colparr,ndropcol
      real :: celli,cell,a,b,c,q,p,t,root1,root2,coll,tc,rr,vv,f,d
      real :: rxi,ryi,rzi,vxij,vyij,vzij,rxij,ryij,rzij,vxi,vyi,vzi,ncc
      integer, allocatable, dimension(:) :: list
      integer :: map(mapsiz)

    ! For double-collision check and collisional growth

      integer :: numcoll, numngr, kcol, zz, ndc, istatus
      integer :: ndrems, ndremr, cols, nrtest, ii, tempcol
      real :: curtc, xcoll, ycoll, zcoll, pos
      real(8):: frac,r1cube,r2cube,r3cube
      real, allocatable, dimension(:,:) :: nctemp, ngrtemp, ncall
      real, allocatable, dimension(:) :: sorttemp
      real(8), allocatable, dimension(:) :: rcol, r_ccncol
      real, allocatable, dimension(:) :: xcol, ycol, zcol, ucol, vcol, wcol
      integer, allocatable, dimension(:) :: idpc, dblcol, icolcomms, icolcommr, icol, idpcol
      real, allocatable, dimension(:) :: xc, yc, yrelc, zc, udropc, vdropc, wdropc, uudropc, vvdropc, wwdropc, rc,r_ccnc

    ! For generating Droplet Size Distribution

      integer :: entr


 66   format(//)
 99   format(1x,2(e12.6),3(f8.6))
 142  format(1x,f16.6,2(i10,1x),15(e15.8,1x))
 143  format(1x,6(i10,1x))
 144  format(1x,f16.6,2(i10,1x),3(e15.8,1x))
 145  format(1x,f16.6,2x,100(i10,1x))
 188  format(1x,3(i3,2x),e12.6)
 199  format(1x,3(f8.6))
 196  format(1x,(f16.6,2x),9(e16.10,2x))
 197  format(1x,3(f8.6,2x))
 198  format(1x,3(e12.6,2x))
 233  format(1x,4(e13.6,2x))
 237  format(i10)
 247  format(f8.6)
 257  format(1x,2(f8.6,2x),2(i10,1x))
 274  format(1x,9(i10,1x))
 277  format(1x,2(f9.6))
 298  format(1x,10(i12))
 299  format(1x,(e12.6),3(f8.6))
 347  format(f8.6,1x,i8)
 497  format(1x,4(d20.14,2x))
 498  format(1x,f8.6,2x,f8.6,2x,i8,2x,i8)
 500  format(1x,i10,1X,4(e15.8,1x))
 501  format(1x,i6,1x,6(e15.8,1x))
 502  format(1x,i6,12(e15.8,1x))
 557  format(1x,3(e14.8,2x),i10)
 567  format(1x,5(f11.8,2x),i10)
 597  format(1x,4(f8.6,2x),i10)
 598  format(1x,f8.6,2x,f8.6,2x,f8.6,2x,i8,2x,i8,2x,i8)
 600  format(1x,10(f8.4,2x)) !for velocity output



      time = time0 + (nt-1)*delt

!! move droplets
      uudrop = udrop !drop velo from last time step
      vvdrop = vdrop
      wwdrop = wdrop
      nddone = 0

      ixp(1:ndropreal)=int(x(1:ndropreal)*oneoverh) + 1
      iyp(1:ndropreal)=int(y(1:ndropreal)*oneoverh) - slabw + 1
      izp(1:ndropreal)=int(z(1:ndropreal)*oneoverh) + 1
      do id = 1,ndropreal !mark
          if(ixp(id) .eq. 0) then
                  ixp(id)=1
          elseif(ixp(id) .eq. N+1) then
                  ixp(id)=N
          endif

          if(izp(id) .eq. 0) then
                  izp(id)=1
          elseif(izp(id) .eq. N+1) then
                  izp(id)=N
          endif

          if(iyp(id) .eq. 0)then
                  iyp(id)=1
          elseif(iyp(id) .eq. 3) then
                  iyp(id)=2
          endif
       if(min(ixp(id),izp(id),iyp(id)) .le. 0 .or. max(ixp(id),izp(id)) .gt. N .or. iyp(id) .gt. 2) then
               print*,'x,y,z index of drop location, id, r,dr3'
               print*,ixp(id),iyp(id),izp(id),x(id),y(id)-mype*n2pe*h,z(id), id, r(id),dr3(id)!markd
       endif
       enddo
   if(thermo .eq. 1 .and. gomic .eq. 2) then
      pp=rhoa*Ra*(1+m_w/29.d-3*qvpp)*temp
      exner = (PP/P0)**RACP
      fv = 1.0
      ! beginning of macroscopic calculations
      !compute sum of radii
      sumrp = 0.0 !parcel sum of radius
      do id = 1,ndropreal
             Red=sqrt(udrop(id)**2+vdrop(id)**2+wdrop(id)**2)*2.d0*r(id)/visc
             if(Red**(1.d0/2.d0)*Nsc**(1.d0/3.d0) .lt. 1.4) then
                 fv(id) = 1.0+0.108*(Nsc**(1.d0/3.d0)*sqrt(Red))**2
             else
                 fv(id) = .78 + .308*Nsc**(1.d0/3.d0)*sqrt(Red)
             endif !ventilation effect coefficient
!!         sumrp = sumrp+r(id)*fv(id) 
      enddo
!      sumrp = sum(r(1:ndropreal)*fv(1:ndropreal))
      sumrp=sum(dr3(1:ndropreal))/3.d0
      call mpi_allreduce(sumrp,tmp,1,mpi_real8,mpi_sum,mpi_comm_world,ierror)
      sumrp = tmp
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!use new P to solve set of conservation eq for T&Qv using Hall's method!!!
      !!/!for the super. (see Grabowski,1989 JAS) uncentered in time formulation!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! Cd = Cdp*S
      cql = 4.0d0*pi*rhow/(rhoa*vol)
      !Cdp=delt*cql*ks*sumrp/real(N1*N2*N3)
      deltaqp = cql*sumrp/real(n1*n2*n3)
      temp0=temp !temperature at previous timestep
      temp=temp0-grav/cp*delt*up+latovercp*deltaqp
      esat=2.53d11*exp(-5.42d3/temp)
      qvs = eps*esat/(PP-Esat)
      rhoa= rhoa*(-grav*up/(Ra*temp)*delt-(temp-temp0)/temp)+rhoa
      thetap = thetapp !old thetapp
      !pp=rhoa*Ra*(1+m_w/29.d-3*qvpp)*temp
      !exner = (pp/P0)**RACP
      thetapp = thetap+latovercp*deltaqp/exner
      qvp    = qvpp !old qvp
      qvpp=qvp-deltaqp
      sp = qvpp/qvs-1.d0
      diffvnd1=1.d-5*(0.015*Temp-1.9)
      ka1=1.5d-11*temp**3-4.8d-8*temp**2+1.d-4*temp-3.9d-4
      ks =1.d0/(rhow*Rv*temp/(esat*diffvnd1)+rhow*Lat/(Ka1*Temp)*(Lat/(Rv*Temp)-1))
     !!sumdr3*cql/3.d0 at each gridpoint
       sumrordeltaq = 0.d0 
       do 120 id = 1,ndropreal  
     !      if(min(ixp(id),izp(id),iyp(id)) .le. 0 .or. max(ixp(id),izp(id),iyp(id)) .gt. N) then
      !             print*,'x,y,z index of drop location, id, r,dr3'
       !            print*,ixp(id),iyp(id),izp(id), id, r(id),dr3(id)!markd
        !   endif
           sumrordeltaq(ixp(id),izp(id),iyp(id)) = sumrordeltaq(ixp(id),izp(id),iyp(id)) + cql*dr3(id)/3.d0!mark
120     continue
     !!calulate S and Cd at each gridpoint
     !! Add averages to scalar fields     
      temp0 = thetap*exner !temp without condensation
      !temp0=temp0-grav/cp*delt*up !temperature after updraft cooling without latent heat release 
      esat0 = 2.53d11*exp(-5.42d3/temp0)
      factor = 5.42d3*esat0/temp0**2
      ttr = ttr + rmth
      qvr = qvr + rmqv
      do 130 iky= 1,n2pe
      do 130 ikz= 1,n3
      do 130 ikx= 1,n1
            ttr(ikx,ikz,iky)=ttr(ikx,ikz,iky) &
                 + latovercp*(sumrordeltaq(ikx,ikz,iky)-deltaqp)
            tempd = temp0+ttr(ikx,ikz,iky)+LatoverCp*deltaqp
            esat = 2.53d11*exp(-5.42d3/tempd)
            qvr(ikx,ikz,iky)=qvr(ikx,ikz,iky)-(sumrordeltaq(ikx,ikz,iky)-deltaqp)
            qvs  = eps*esat/(pp-esat)
            s(ikx,ikz,iky)=(qvr(ikx,ikz,iky)+qvpp)/qvs-1.d0
 130  continue
      rmth = 0.d0
      rmqv = 0.d0
      rmth = sum(ttr(1:n1,1:n3,1:n2pe))
      rmqv = sum(qvr(1:n1,1:n3,1:n2pe))
      call mpi_allreduce(rmth,tmpv,1,mpi_real,mpi_sum,mpi_comm_world,ierror)
      rmth = tmpv
      call mpi_allreduce(rmqv,tmpv,1,mpi_real,mpi_sum,mpi_comm_world,ierror)
      rmqv = tmpv
      rmth = rmth/real(n1*n3*n2)
      rmqv = rmqv/real(n1*n3*n2)
      !temp = temp + rmth
      thetapp = temp/exner !add back to mean theta value.
      !qvpp  = qvpp + rmqv
      lwc = 0.0d0
    !!adjust size of droplets & including ventilation coeffient fv
      do 140 id=1,ndropreal
         tempd=temp0+ttr(ixp(id),izp(id),iyp(id))+LatoverCp*deltaqp!+LatoverCp*sumrordeltaq(ixp(id),izp(id),iyp(id)) !temperature of that grid
            !total temperature at the droplet location
         curv=2.d0*7.61d-2/(Rv*rhow*tempd*r(id))
	 if(isolu .ne. 0) then
         solu=2.0d0*m_w/m_s*rho_ccn*r_ccn(id)**3/rhow
         seq=exp(curv-solu/(r(id)**3-r_ccn(id)**3))-1.d0!equillibrium supersat.
	 else
	 seq = curv
	 endif
         if(abs(seq) .ge. 1) print*, 'big seq=',seq,'tempd,r(id)',tempd,r(id)
         vtemp=r(id)*3.0d0*delt*ks*(s(ixp(id),izp(id),iyp(id))-seq)*fv(id)+r(id)**3
         if(vtemp .gt. r_ccn(id)**3 .and. r(id) .lt. rtrunc) then
            dr3(id) = vtemp-r(id)**3
            r(id) = vtemp**(1.d0/3.0d0)
         else
            dr3(id)=0.d0
         endif
         lwc = lwc+r(id)**3
 140   continue
         call mpi_reduce(lwc,lwctmp,1,mpi_real8,mpi_sum,0,mpi_comm_world,ierror)
      if(mype.eq.0) lwc = lwctmp*(cql/3.d0)/real(n1*n3*n2)
      
    !!adjust t and qv field
!    if(mype .eq. 0 ) print*,'gegeda',r(10) !mark
  if (mype.eq.0 .and. mod(nt-1,micdtout) .eq. 0 )  then
    write(82,196) time,Sp,pp,temp,thetapp,qvpp,lwc,qvs,deltaqp,rhoa
  endif !mype
   endif !thermo==1 condensate

      do 160 id = 1,ndropreal
         if (r(id) .lt. rtrunc_min) then
	   taup=delt
	 elseif (r(id) .lt. 4.d-5) then
           taup=kwt*r(id)**2/grav
         else
           taup=(kwt2*r(id)-.1227)/grav
         endif !terminal

      dtovertaup(id) = delt/taup
      dxp = oneoverh*x(id) - real(ixp(id)) + 1
      dyp = oneoverh*y(id) - real(iyp(id)) - slabw + 1
      dzp = oneoverh*z(id) - real(izp(id)) + 1

!!! flow velocity at droplet position using linear interpolation
      flowu(id) = &
        (1.-dxp)*(1.-dyp)*(1.-dzp)*uflow(ixp(id),izp(id),iyp(id)) &
        + dxp*(1.-dyp)     *(1.-dzp)*uflow(ixp(id)+1,izp(id),iyp(id)) &
        + dxp*dyp          *(1.-dzp)*uflow(ixp(id)+1,izp(id),iyp(id)+1) &
        + (1.-dxp)*dyp     *(1.-dzp)*uflow(ixp(id),izp(id),iyp(id)+1) &
        + (1.-dxp)*(1.-dyp)*dzp     *uflow(ixp(id),izp(id)+1,iyp(id)) &
        + dxp*(1.-dyp)     *dzp     *uflow(ixp(id)+1,izp(id)+1,iyp(id)) &
        + dxp*dyp          *dzp     *uflow(ixp(id)+1,izp(id)+1,iyp(id)+1) &
        + (1.-dxp)*dyp     *dzp     *uflow(ixp(id),izp(id)+1,iyp(id)+1)

      UDROP(id)=dtovertaup(id)*flowu(id)+UDROP(id)*(1.-dtovertaup(id))

      flowv(id) = &
        (1.-dxp)*(1.-dyp)*(1.-dzp)*vflow(ixp(id),izp(id),iyp(id)) &
        + dxp*(1.-dyp)     *(1.-dzp)*vflow(ixp(id)+1,izp(id),iyp(id)) &
        + dxp*dyp          *(1.-dzp)*vflow(ixp(id)+1,izp(id),iyp(id)+1) &
        + (1.-dxp)*dyp     *(1.-dzp)*vflow(ixp(id),izp(id),iyp(id)+1) &
        + (1.-dxp)*(1.-dyp)*dzp     *vflow(ixp(id),izp(id)+1,iyp(id)) &
        + dxp*(1.-dyp)     *dzp     *vflow(ixp(id)+1,izp(id)+1,iyp(id)) &
        + dxp*dyp          *dzp     *vflow(ixp(id)+1,izp(id)+1,iyp(id)+1) &
        + (1.-dxp)*dyp     *dzp     *vflow(ixp(id),izp(id)+1,iyp(id)+1)

      VDROP(id)=dtovertaup(id)*flowv(id)+VDROP(id)*(1.-dtovertaup(id))

      floww(id) = &
        (1.-dxp)*(1.-dyp)*(1.-dzp)*wflow(ixp(id),izp(id),iyp(id)) &
        + dxp*(1.-dyp)     *(1.-dzp)*wflow(ixp(id)+1,izp(id),iyp(id)) &
        + dxp*dyp          *(1.-dzp)*wflow(ixp(id)+1,izp(id),iyp(id)+1) &
        + (1.-dxp)*dyp     *(1.-dzp)*wflow(ixp(id),izp(id),iyp(id)+1) &
        + (1.-dxp)*(1.-dyp)*dzp     *wflow(ixp(id),izp(id)+1,iyp(id)) &
        + dxp*(1.-dyp)     *dzp     *wflow(ixp(id)+1,izp(id)+1,iyp(id)) &
        + dxp*dyp          *dzp     *wflow(ixp(id)+1,izp(id)+1,iyp(id)+1) &
        + (1.-dxp)*dyp     *dzp     *wflow(ixp(id),izp(id)+1,iyp(id)+1)
      WDROP(id)=dtovertaup(id)*floww(id)+WDROP(id)*(1.-dtovertaup(id))-grav*delt
160   continue
    ! Possibility to save droplet trajectories
          if (everynd .ge. 1) then
           ndumpdrop = 10              ! Save positions of droplets every ndumpdrop timesteps
           if (mod(nt,ndumpdrop) .eq. 0) then

              sendfol = ndroppe*4/everynd

              allocate( idpfols(sendfol) )
              allocate( xfols(sendfol), yfols(sendfol), zfols(sendfol) , rfols(sendfol))!, ufols(sendfol), vfols(sendfol), wfols(sendfol) )
              allocate( idpfolr(sendfol*npe) )
              allocate( xfolr(sendfol*npe), yfolr(sendfol*npe), zfolr(sendfol*npe), rfolr(sendfol*npe))!, ufolr(sendfol*npe), vfolr(sendfol*npe), wfolr(sendfol*npe) )

              ifol = 0
              idpfols = 0
              xfols = 0.0
              yfols = 0.0
              zfols = 0.0
!              ufols = 0.0
!              vfols = 0.0
!              wfols = 0.0
              rfols = 0.0
              do 170 id = 1, ndropreal
                 if (mod(idp(id),everynd) .eq. 0) then
                   ifol = ifol + 1
                   idpfols(ifol) = idp(id)
                   xfols(ifol) = x(id)
                   yfols(ifol) = y(id)
                   zfols(ifol) = z(id)
!                   ufols(ifol) = udrop(id)
!                   vfols(ifol) = vdrop(id)
!                   wfols(ifol) = wdrop(id)
                   rfols(ifol) = r(id)
                 endif
170            continue ! id
              call mpi_gather(idpfols,sendfol,MPI_INTEGER,idpfolr,sendfol,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
              call mpi_gather(xfols,sendfol,MPI_REAL,xfolr,sendfol,MPI_REAL,0,MPI_COMM_WORLD,ierror)
              call mpi_gather(yfols,sendfol,MPI_REAL,yfolr,sendfol,MPI_REAL,0,MPI_COMM_WORLD,ierror)
              call mpi_gather(zfols,sendfol,MPI_REAL,zfolr,sendfol,MPI_REAL,0,MPI_COMM_WORLD,ierror)
!              call mpi_gather(ufols,sendfol,MPI_REAL,ufolr,sendfol,MPI_REAL,0,MPI_COMM_WORLD,ierror)
!              call mpi_gather(vfols,sendfol,MPI_REAL,vfolr,sendfol,MPI_REAL,0,MPI_COMM_WORLD,ierror)
!              call mpi_gather(wfols,sendfol,MPI_REAL,wfolr,sendfol,MPI_REAL,0,MPI_COMM_WORLD,ierror)
              call mpi_gather(rfols,sendfol,MPI_REAL8,rfolr,sendfol,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
              if (mype .eq. 0) then
                 idpfol = 0
                 xfol = 0.0
                 yfol = 0.0
                 zfol = 0.0
!                 ufol = 0.0
!                 vfol = 0.0
!                 wfol = 0.0
                 rfol = 0.0
                 do 180 i = 1,ndroppe/everynd*npe*4
                   idpfoli = idpfolr(i)
                   if (idpfoli .ne. 0) then
                       idpfol(idpfoli/everynd) = idpfoli
                       xfol(idpfoli/everynd) = xfolr(i)
                       yfol(idpfoli/everynd) = yfolr(i)
                       zfol(idpfoli/everynd) = zfolr(i)
!                       ufol(idpfoli/everynd) = ufolr(i)
!                       vfol(idpfoli/everynd) = vfolr(i)
!                       wfol(idpfoli/everynd) = wfolr(i)
                       rfol(idpfoli/everynd) = rfolr(i)
                   endif
180              continue ! i
                 do i = 1,ndrop/everynd
                   write(41,500) idpfol(i),xfol(i),yfol(i),zfol(i),rfol(i) !,ufol(i),vfol(i),wfol(i)
                 enddo
              endif ! mype .eq. 0
              deallocate( idpfols )
              deallocate( xfols, yfols, zfols,rfols) !, ufols, vfols, wfols, rfols )
              deallocate( xfolr, yfolr, zfolr,rfolr) !, ufolr, vfolr, wfolr, rfolr )

           endif ! mod

      endif ! tracking droplets

    ! Collision detection scheme &  collision calculation(enable =1)
      if (colcal .eq. 1) then
           colpars = 0
           infocols = 0.0
           if(npe .ne. m) then
           do 190 id=1,ndropreal
              yrel(id) = y(id) - slabw*h
              ! check if droplets on the left boundary of each slab can be a potential collisional partner of its left neighbour
              if (yrel(id) .le. cellw) then
              colpars = colpars + 1
              infocols(colpars,1) = real(idp(id))
              infocols(colpars,2) = x(id)
              infocols(colpars,3) = y(id)
              infocols(colpars,4) = z(id)
              infocols(colpars,5) = udrop(id)
              infocols(colpars,6) = vdrop(id)
              infocols(colpars,7) = wdrop(id)
              infocols(colpars,8) = flowu(id)
              infocols(colpars,9) = flowv(id)
              infocols(colpars,10) = floww(id)
              infocols(colpars,11) = r(id)
              infocols(colpars,12) = dtovertaup(id)
              infocols(colpars,13) = uudrop(id)
              infocols(colpars,14) = vvdrop(id)
              infocols(colpars,15) = wwdrop(id)
              infocols(colpars,16) = real(id)
              infocols(colpars,17) = r_ccn(id)
              endif!yrel 
190	   continue
           else!npe=m
              colpars=ndropreal
              infocols(1:ndropreal,1) = real(idp(1:ndropreal))
              infocols(1:ndropreal,2) = x(1:ndropreal)
              infocols(1:ndropreal,3) = y(1:ndropreal)
              infocols(1:ndropreal,4) = z(1:ndropreal)
              infocols(1:ndropreal,5) = udrop(1:ndropreal)
              infocols(1:ndropreal,6) = vdrop(1:ndropreal)
              infocols(1:ndropreal,7) = wdrop(1:ndropreal)
              infocols(1:ndropreal,8) = flowu(1:ndropreal)
              infocols(1:ndropreal,9) = flowv(1:ndropreal)
              infocols(1:ndropreal,10) = floww(1:ndropreal)
              infocols(1:ndropreal,11) = r(1:ndropreal)
              infocols(1:ndropreal,12) = dtovertaup(1:ndropreal)
              infocols(1:ndropreal,13) = uudrop(1:ndropreal)
              infocols(1:ndropreal,14) = vvdrop(1:ndropreal)
              infocols(1:ndropreal,15) = wwdrop(1:ndropreal)
              do id=1,ndropreal
              infocols(id,16) = real(id) 
              enddo
              infocols(1:ndropreal,17)= r_ccn(1:ndropreal)
           endif!npe=m
! exchange droplet infos on the left boundary of each slab with its left neighbours. only do left boundary to avoid repeated collision counts.

           do i=1,npe
              ntofrom(i+1) = i-1
           enddo
           ntofrom(1) = npe -1
           ntofrom(npe+2) = 0

           nto = ntofrom(mype+1)
           nfrom = ntofrom(mype+3)

           call mpi_sendrecv(colpars,1,MPI_INTEGER,nto,nto,colparr,1,MPI_INTEGER,nfrom,mype,MPI_COMM_WORLD,status,ierror)
           nbufsendl = colpars*17
           nbufrecvr = colparr*17

           call mpi_sendrecv(infocols(1:colpars,:),nbufsendl,MPI_REAL,nto,nto,infocolr(1:colparr,:),nbufrecvr,MPI_REAL,nfrom,mype,MPI_COMM_WORLD,status,ierror)

! Finally, add droplets from right neighbour to list of potential collision partners.
           ndropcol = colparr + ndropreal

           allocate( list(ndropcol) )
           allocate( idpc(ndropcol), xc(ndropcol), yc(ndropcol), yrelc(ndropcol), zc(ndropcol) )
           allocate( udropc(ndropcol), vdropc(ndropcol),wdropc(ndropcol), rc(ndropcol),r_ccnc(ndropcol) )
           allocate( uudropc(ndropcol), vvdropc(ndropcol), wwdropc(ndropcol) )
           allocate( dtovertaupc(ndropcol), flowuc(ndropcol), flowvc(ndropcol), flowwc(ndropcol) )

           idpc(1:ndropreal) = idp(1:ndropreal)

           xc(1:ndropreal) = x(1:ndropreal)
           yc(1:ndropreal) = y(1:ndropreal)
           zc(1:ndropreal) = z(1:ndropreal)

           udropc(1:ndropreal) = udrop(1:ndropreal)
           vdropc(1:ndropreal) = vdrop(1:ndropreal)
           wdropc(1:ndropreal) = wdrop(1:ndropreal)
           uudropc(1:ndropreal) = uudrop(1:ndropreal)
           vvdropc(1:ndropreal) = vvdrop(1:ndropreal)
           wwdropc(1:ndropreal) = wwdrop(1:ndropreal)
           flowuc(1:ndropreal) = flowu(1:ndropreal)
           flowvc(1:ndropreal) = flowv(1:ndropreal)
           flowwc(1:ndropreal) = floww(1:ndropreal)

           rc(1:ndropreal) = r(1:ndropreal)
           r_ccnc(1:ndropreal) = r_ccn(1:ndropreal)
           dtovertaupc(1:ndropreal) = dtovertaup(1:ndropreal)

           idpc(ndropreal+1:ndropcol) = int(infocolr(1:colparr,1))
           xc(ndropreal+1:ndropcol) = infocolr(1:colparr,2)
           if (mype .eq. npe-1) then
              yc(ndropreal+1:ndropcol) = infocolr(1:colparr,3) + real(n2)*h
           else
              yc(ndropreal+1:ndropcol) = infocolr(1:colparr,3)
           endif!mype
           zc(ndropreal+1:ndropcol) = infocolr(1:colparr,4)

           udropc(ndropreal+1:ndropcol) = infocolr(1:colparr,5)
           vdropc(ndropreal+1:ndropcol) = infocolr(1:colparr,6)
           wdropc(ndropreal+1:ndropcol) = infocolr(1:colparr,7)
           flowuc(ndropreal+1:ndropcol) = infocolr(1:colparr,8)
           flowvc(ndropreal+1:ndropcol) = infocolr(1:colparr,9)
           flowwc(ndropreal+1:ndropcol) = infocolr(1:colparr,10)

           rc(ndropreal+1:ndropcol) = infocolr(1:colparr,11)
           dtovertaupc(ndropreal+1:ndropcol) = infocolr(1:colparr,12)
           uudropc(ndropreal+1:ndropcol) = infocolr(1:colparr,13)
           vvdropc(ndropreal+1:ndropcol) = infocolr(1:colparr,14)
           wwdropc(ndropreal+1:ndropcol) = infocolr(1:colparr,15)

!**************************************************************** 
!      ICELL ( jx, jy, jz) = 1 + MOD ( jx - 1 + M, M )
!     :                          + MOD ( jz - 1 + M, M ) * M
!     :                          + MOD ( jy + MY, MY + 1 ) * M * M

!    ** FIND HALF THE NEAREST NEIGHBOURS OF EACH CELL **^M

           DO 200 jz = 1, M
              DO 200 jy = 1, MY
                 DO 200 jx = 1, M

                   IMAP = ( ICELL ( jx, jy, jz, M, MY ) - 1 ) * 13

                   MAP( IMAP + 1  ) = ICELL( jx + 1, jy    , jz    , M, MY )
                   MAP( IMAP + 2  ) = ICELL( jx + 1, jy + 1, jz    , M, MY )
                   MAP( IMAP + 3  ) = ICELL( jx    , jy + 1, jz    , M, MY )
                   MAP( IMAP + 4  ) = ICELL( jx - 1, jy + 1, jz    , M, MY )
                   MAP( IMAP + 5  ) = ICELL( jx + 1, jy    , jz - 1, M, MY )
                   MAP( IMAP + 6  ) = ICELL( jx + 1, jy + 1, jz - 1, M, MY )
                   MAP( IMAP + 7  ) = ICELL( jx    , jy + 1, jz - 1, M, MY )
                   MAP( IMAP + 8  ) = ICELL( jx - 1, jy + 1, jz - 1, M, MY )
                   MAP( IMAP + 9  ) = ICELL( jx + 1, jy    , jz + 1, M, MY )
                   MAP( IMAP + 10 ) = ICELL( jx + 1, jy + 1, jz + 1, M, MY )
                   MAP( IMAP + 11 ) = ICELL( jx    , jy + 1, jz + 1, M, MY )
                   MAP( IMAP + 12 ) = ICELL( jx - 1, jy + 1, jz + 1, M, MY )
                   MAP( IMAP + 13 ) = ICELL( jx    , jy    , jz + 1, M, MY )

  200      CONTINUE


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  set up the linked list and head of chain arrays
!
!    *******************************************************************^M
!    ** ROUTINE TO SET UP LINKED LIST AND THE HEAD OF CHAIN ARRAYS    **^M
!    **                                                               **^M
!    ** PRINCIPAL VARIABLES:                                          **^M
!    **                                                               **^M
!    ** INTEGER N                  NUMBER OF ATOMS                    **^M
!    ** INTEGER M                  NUMBER OF CELLS IN EACH DIRECTION  **^M
!    ** INTEGER NCELL              TOTAL NUMBER OF CELLS (M**3)       **^M
!    ** INTEGER LIST(N)            LINKED LIST OF ATOMS               **^M
!    ** REAL    RX(N),RY(N),RZ(N)  POSITIONS                          **^M
!    ** REAL    RCUT               THE CUTOFF DISTANCE FOR THE FORCE  **^M
!    **                                                               **^M
!    ** USAGE:                                                        **^M
!    **                                                               **^M
!    ** EACH ATOM IS SORTED INTO ONE OF THE M**3 SMALL CELLS.         **^M
!    ** THE FIRST ATOM IN EACH CELL IS PLACED IN THE HEAD ARRAY.      **^M
!    ** SUBSEQUENT ATOMS ARE PLACED IN THE LINKED LIST ARRAY.         **^M
!    ** ATOM COORDINATES ARE ASSUMED TO BE BETWEEN -0.5 AND +0.5.     **^M
!    ** THE ROUTINE IS CALLED EVERY TIMESTEP BEFORE THE FORCE ROUTINE.**^M
!    *******************************************************************^M
!    ** ZERO HEAD OF CHAIN ARRAY                                      **^M

           head = 0
           list = 0
           celli = real ( m )
           cell=1.0/((N*h)/celli) !1/2delx

!    ** SORT ALL ATOMS **^M

           do 210 i = 1, ndropcol
	      if(rc(i) .ge. rtrunc_min .and. rc(i) .le. rtrunc) then !treat rc(i) beyond the radius range as ghost particles
              yrelc(i) = yc(i) - slabw*h

              xcell = int(xc(i) * cell) !#series of cell (x/2delx)
              if (xcell .gt. m .or. xcell .lt. -1) then 
                 print*, 'mype', mype, 'has droplet', idpc(i),'whose xcell', xcell, 'is out of boundary','udrop',uudropc(i)
              elseif (xcell .eq. m) then
                 xcell = xcell - 1
              elseif (xcell .eq. -1) then
                 xcell = xcell + 1
              endif !xcell

              zcell = int(zc(i) * cell)
              if (zcell .gt. m .or. zcell .lt. -1) then 
                 print*, 'mype', mype, 'has droplet', idpc(i),'whose zcell', zcell, 'is out of boundary','wdrop',wwdropc(i)
              elseif (zcell .eq. m) then
                 zcell = zcell - 1
              elseif (zcell .eq. -1) then
                 zcell = zcell + 1
              endif !zcell

              ycell = int(yrelc(i) * cell)
              if(ycell .gt. my+1 .or. ycell .lt. -1) then
                 print*, 'mype',mype,'has droplet', idpc(i), 'with too big ycell', ycell,'vdrop',vvdropc(i)
              elseif (ycell .eq. my+1) then
                 ycell = ycell - 1
              elseif (ycell .eq. -1) then
                 ycell = ycell + 1
              endif !ycell

              tcell = 1 + xcell + zcell*m + ycell*m*m

              if (tcell .gt. m*m*(my+1)) then
                 print*, '         '
                 print*, 'mype has tcell too big, tcell, i, idp,  x(i), y(i), yrel(i), z(i) =', mype, tcell, i, idpc(i), xc(i), yc(i), yrelc(i), zc(i)
                 print*, 'xcell, ycell),zcell: ', xcell, ycell, zcell
                 print*, '         '
              elseif (tcell .lt. 0) then
                 print*, 'mype has tcell <0 , tcell, i, idp, x(i), y(i), yrel(i), z(i) =', mype, tcell, i, idpc(i), xc(i), yc(i), yrelc(i), zc(i)
              endif !tcell

              list(i)     = head(tcell)
              head(tcell) = i
	    endif !r_trunc_min<r<r_trunc_max
  210       continue

     !write out head and list
!     if(mype.eq.0) then
!           write(57,298) (head(id),id=1,ncell)
!           write(57,66)
!           write(57,298) (list(id),id=1,ndrop)
!     endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   *******************************************************************  !!
!!   *******************************************************************  !!
!!   ** LOOP OVER ALL CELLS TO UPDATS DISTURBANCE FLOW  ****************  !!
!!   **                         by Sisi Chen 2016       ****************  !!
!!   *******************************************************************  !!
!!   *******************************************************************  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           if (ihydro .eq. 1) then
              allocate( udist(ndropcol), vdist(ndropcol), wdist(ndropcol) )
              allocate( uudist(ndropcol), vvdist(ndropcol), wwdist(ndropcol) )
              allocate( distsendr1(colparr,3), distrecvr2(colparr,3) )
              allocate( distrecvl1(colpars,3), distsendl2(colpars,3))
              allocate( ucharact(ndropcol), udist_syn(ndropcol), uudist_syn(ndropcol) )
		  uudist_syn    = 0.
    	      	  uudist        = 0.
	          vvdist        = 0.
	          wwdist        = 0.
	          maxudist      = 1.
	          udist_iter    = 0
	          maxuall       = 1.0
		  maxuall_old   = 0.5
		  maxuall_old2  = 1.5

	      do id = 1,ndropcol
                if(rc(id) .gt. 4.d-5 ) then !terminal velocity
			ucharact(id) = rc(id)*kwt2 - 0.1227
		else
			ucharact(id) = rc(id)**2*kwt
		endif
	      enddo

          do 7100 while (maxuall .gt. 1.e-5  .AND. udist_iter .lt. 200 .AND. (maxuall-maxuall_old)*(maxuall-maxuall_old2) .ne. 0 )
		maxuall_old2 = maxuall_old
	        maxuall_old = maxuall
                udist = 0.
                vdist = 0.
                wdist = 0.
              DO 8500 tcell = 1, ncell
                 i = head(tcell)

!  ** LOOP OVER ALL DROPLETS IN THE CELL **^M
                do 8100 while (i .gt. 0 ) 

                   rxi = xc(i)
                   ryi = yc(i)
                   rzi = zc(i)

                   vxi = uudropc(i)
                   vyi = vvdropc(i)
                   vzi = wwdropc(i) !  ** LOOP OVER ALL DROPLETS BELOW I IN THE CURRENT CELL **^M

                   j = list(i)

               do 8200 while (j .GT. 0 )

                       rxij = rxi - xc(j)
                       if (abs(rxij) .gt. 2*cellw) rxij = -( rxij/abs(rxij) ) * ( N*H - abs(rxij) )
                       ryij = ryi - yc(j)
                       rzij = rzi - zc(j)
                       if (abs(rzij) .gt. 2*cellw) rzij = -( rzij/abs(rzij) ) * ( N*H - abs(rzij) )

                       rxji = -rxij
                       ryji = -ryij
                       rzji = -rzij

                       rr = sqrt(rxij**2+ryij**2+rzij**2)
                  !disturbance flow caused by droblet i on j
                       if (rr .le. 30.*rc(i) .and. rr .gt. (1.0+1.d-4)*(rc(i)+rc(j)) ) then

                           distdxji = uudropc(i)-flowuc(i)-uudist(i) !Vp
                           distdyji = vvdropc(i)-flowvc(i)-vvdist(i)
                           distdzji = wwdropc(i)-flowwc(i)-wwdist(i)

                           distaji = (0.75d0*rc(i)/rr - 0.75d0*(rc(i)/rr)**3) !constant1

                           distbxji = rxji*(distdxji*rxji+distdyji*ryji+distdzji*rzji)/(rr*rr)
                           distbyji = ryji*(distdxji*rxji+distdyji*ryji+distdzji*rzji)/(rr*rr)
                           distbzji = rzji*(distdxji*rxji+distdyji*ryji+distdzji*rzji)/(rr*rr)

                           distcji = (0.75d0*rc(i)/rr + 0.25d0*(rc(i)/rr)**3) !constant2

                           udisttemp = distaji*distbxji+distcji*distdxji
                           vdisttemp = distaji*distbyji+distcji*distdyji
                           wdisttemp = distaji*distbzji+distcji*distdzji

                           udist(j) = udist(j) + udisttemp
                           vdist(j) = vdist(j) + vdisttemp
                           wdist(j) = wdist(j) + wdisttemp


                           distdxij = uudropc(j)-flowuc(j)-uudist(j)
                           distdyij = vvdropc(j)-flowvc(j)-vvdist(j)
                           distdzij = wwdropc(j)-flowwc(j)-wwdist(j)

                           distaij = (0.75d0*rc(j)/rr - 0.75d0*(rc(j)/rr)**3)

                           distbxij = rxij*(distdxij*rxij+distdyij*ryij+distdzij*rzij)/(rr**2)
                           distbyij = ryij*(distdxij*rxij+distdyij*ryij+distdzij*rzij)/(rr**2)
                           distbzij = rzij*(distdxij*rxij+distdyij*ryij+distdzij*rzij)/(rr**2)

                           distcij = (0.75d0*rc(j)/rr + 0.25d0*(rc(j)/rr)**3)

                           udisttemp = distaij*distbxij+distcij*distdxij
                           vdisttemp = distaij*distbyij+distcij*distdyij
                           wdisttemp = distaij*distbzij+distcij*distdzij

                           udist(i) = udist(i) + udisttemp
                           vdist(i) = vdist(i) + vdisttemp
                           wdist(i) = wdist(i) + wdisttemp

                       endif!rr

                       j = list(j)
8200		     continue!j loop
!  ** LOOP OVER NEIGHBOURING CELLS **^M

                   jcell0 = 13 * (tcell - 1)

                   DO 8000 nabor = 1, 13

                       jcell = map ( jcell0 + nabor )

!  ** LOOP OVER ALL DROPLETS IN NEIGHBOURING CELLS **^M

                       j = head(jcell)

                       do 8300 while ( j .NE. 0 )

                           rxij = rxi - xc(j)
                           if (abs(rxij) .gt. 2*cellw) rxij = -( rxij/abs(rxij) ) * ( N*H - abs(rxij) )
                           ryij = ryi - yc(j)
                           rzij = rzi - zc(j)
                           if (abs(rzij) .gt. 2*cellw) rzij = -( rzij/abs(rzij) ) * ( N*H - abs(rzij) )

                           rxji = -rxij
                           ryji = -ryij
                           rzji = -rzij

                           rr = sqrt(rxij*rxij+ryij*ryij+rzij*rzij)

                        !compute disturbance flow
                        !disturbance flow caused by droblet i on j
                           if (rr .le. 30.*rc(i).and. rr .gt. (1.0+1.0d-4)*(rc(i)+rc(j)) ) then

                              distdxji = uudropc(i)-flowuc(i)-uudist(i)
                              distdyji = vvdropc(i)-flowvc(i)-vvdist(i)
                              distdzji = wwdropc(i)-flowwc(i)-wwdist(i)

                              distaji = (0.75d0*rc(i)/rr - 0.75d0*(rc(i)/rr)**3)

                              distbxji = rxji*(distdxji*rxji+distdyji*ryji+distdzji*rzji)/(rr*rr)
                              distbyji = ryji*(distdxji*rxji+distdyji*ryji+distdzji*rzji)/(rr*rr)
                              distbzji = rzji*(distdxji*rxji+distdyji*ryji+distdzji*rzji)/(rr*rr)

                              distcji = (0.75d0*rc(i)/rr + 0.25d0*(rc(i)/rr)**3)

                              udisttemp = distaji*distbxji+distcji*distdxji
                              vdisttemp = distaji*distbyji+distcji*distdyji
                              wdisttemp = distaji*distbzji+distcji*distdzji 

                              udist(j) = udist(j) + udisttemp
                              vdist(j) = vdist(j) + vdisttemp
                              wdist(j) = wdist(j) + wdisttemp

                           !disturbance flow caused by droblet j on i

                              distdxij = uudropc(j)-flowuc(j)-uudist(j)
                              distdyij = vvdropc(j)-flowvc(j)-vvdist(j)
                              distdzij = wwdropc(j)-flowwc(j)-wwdist(j)

                              distaij = (0.75d0*rc(j)/rr - 0.75d0*(rc(j)/rr)**3)

                              distbxij = rxij*(distdxij*rxij+distdyij*ryij+distdzij*rzij)/(rr*rr)
                              distbyij = ryij*(distdxij*rxij+distdyij*ryij+distdzij*rzij)/(rr*rr)
                              distbzij = rzij*(distdxij*rxij+distdyij*ryij+distdzij*rzij)/(rr*rr)

                              distcij = (0.75d0*rc(j)/rr + 0.25d0*(rc(j)/rr)**3)

                              udisttemp = distaij*distbxij+distcij*distdxij
                              vdisttemp = distaij*distbyij+distcij*distdyij
                              wdisttemp = distaij*distbzij+distcij*distdzij 

                              udist(i) = udist(i) + udisttemp
                              vdist(i) = vdist(i) + vdisttemp
                              wdist(i) = wdist(i) + wdisttemp 

                           endif!rr

                           j = list(j)

8300			continue !j loop
8000               continue

!		endif!cellopen
                   i = list(i)
8100          continue !i loop

8500        continue

             !update droplet velocity by adding disturbance flow

              distsendr1(1:colparr,1) = udist(ndropreal+1:ndropcol)
              distsendr1(1:colparr,2) = vdist(ndropreal+1:ndropcol)
              distsendr1(1:colparr,3) = wdist(ndropreal+1:ndropcol)

              nbufright = colparr*3 !from me to right
              nbufleft = colpars*3 !from left to me

             !first send udist data to its right slab to update velo
              nto = ntofrom(mype+3)
              nfrom = ntofrom(mype+1)
              call mpi_sendrecv(distsendr1,nbufright,MPI_REAL,nto,nto,distrecvl1,nbufleft,MPI_REAL,nfrom,mype,MPI_COMM_WORLD,status,ierror)
!!vectorize
              do id = 1,colpars

                 idinfocols= nint(infocols(id,16))

                 udist(idinfocols) = udist(idinfocols)+distrecvl1(id,1)
                 vdist(idinfocols) = vdist(idinfocols)+distrecvl1(id,2)
                 wdist(idinfocols) = wdist(idinfocols)+distrecvl1(id,3)

                 distsendl2(id,1) = udist(idinfocols)
                 distsendl2(id,2) = vdist(idinfocols)
                 distsendl2(id,3) = wdist(idinfocols)

              enddo

             !then send the updated velo back to its left slab

              nto = ntofrom(mype+1)
              nfrom = ntofrom(mype+3)
              call mpi_sendrecv(distsendl2,nbufleft,MPI_REAL,nto,nto,distrecvr2,nbufright,MPI_REAL,nfrom,mype,MPI_COMM_WORLD,status,ierror)
	      
              udist(ndropreal+1:ndropcol) = distrecvr2(1:colparr,1)
              vdist(ndropreal+1:ndropcol) = distrecvr2(1:colparr,2)
              wdist(ndropreal+1:ndropcol) = distrecvr2(1:colparr,3)

		udist_syn=sqrt(udist*udist+vdist*vdist+wdist*wdist)
!mark the max function should change as well!
	         maxudist = maxval(abs((udist_syn(1:ndropcol)-uudist_syn(1:ndropcol))/ucharact(1:ndropcol)))

	      call mpi_allreduce(maxudist,maxuall,1,mpi_real,mpi_max,mpi_comm_world,ierror)
              uudist = udist
              vvdist = vdist
              wwdist = wdist
	      uudist_syn = udist_syn
              udist_iter = udist_iter + 1 
7100        continue
              do id= 1,ndropcol
                        if (abs(udist(id))>1.) then
                                print*,'udist=',udist(id),'id',idpc(id)
                                udist(id)=0.
                        endif
                        if (abs(vdist(id))>1.) then
                                print*,'vdist=',vdist(id),'id',idpc(id)
                                vdist(id)=0.
                        endif
                        if (abs(wdist(id))>1.) then
                                print*,'wdist=',wdist(id),'id',idpc(id)
                                wdist(id)=0.
                        endif
              enddo !id=ndropcol
             !update the velocity
              udropc(1:ndropcol) = udropc(1:ndropcol)+dtovertaupc(1:ndropcol)*udist(1:ndropcol)
              vdropc(1:ndropcol) = vdropc(1:ndropcol)+dtovertaupc(1:ndropcol)*vdist(1:ndropcol)
              wdropc(1:ndropcol) = wdropc(1:ndropcol)+dtovertaupc(1:ndropcol)*wdist(1:ndropcol)

              udrop(1:ndropreal) = udropc(1:ndropreal)
              vdrop(1:ndropreal) = vdropc(1:ndropreal)
              wdrop(1:ndropreal) = wdropc(1:ndropreal) 

              deallocate( distsendr1, distrecvr2 )
              deallocate( distrecvl1, distsendl2 )
              deallocate( udist, vdist, wdist )
		      deallocate( uudist, vvdist, wwdist )
		      deallocate( ucharact, udist_syn, uudist_syn)

           endif!ihydro=1


!!!  !!!!!!!collision detection

           allocate( ncdat(maxcol,ncele) )
           allocate( ngrdat(maxcol,ngrele) )

           ncdat = 0.
           ngrdat = 0.
           nc=0
           wr=0.
           numgr = 0

!  ** LOOP OVER ALL CELLS **^M

           do 5000 tcell = 1, ncell

              i = head(tcell)

!  ** LOOP OVER ALL MOLECULES IN THE CELL **^M

 1000         if ( i .gt. 0 ) then 

                 rxi = xc(i)
                 ryi = yc(i)
                 rzi = zc(i)
                 vxi = udropc(i)
                 vyi = vdropc(i)
                 vzi = wdropc(i)

!  ** LOOP OVER ALL MOLECULES BELOW I IN THE CURRENT CELL **^M

                 j = list(i)

 2000            if (j .gt. 0 ) then

                   rxij = rxi - xc(j) !2x(2xltx)
                   if (abs(rxij) .gt. 2*cellw) rxij = -( rxij/abs(rxij) ) * ( N*H - abs(rxij) ) !periodicity
                   ryij = ryi - yc(j)
                   rzij = rzi - zc(j)
                   if (abs(rzij) .gt. 2*cellw) rzij = -( rzij/abs(rzij) ) * ( N*H - abs(rzij) )
                   rr = sqrt(rxij**2+ryij**2+rzij**2) !distance in 3D

                   vxij = vxi - udropc(j)
                   vyij = vyi - vdropc(j)
                   vzij = vzi - wdropc(j)
                   vv = sqrt(vxij**2+vyij**2+vzij**2)

                   a = vv**2
                   b = 2*(rxij*vxij+ryij*vyij+rzij*vzij)
                   c = rr**2 - (rc(i)+rc(j))**2

                   wr=(b/2.)/rr

                   onepmin = 0.99*(rc(i)+rc(j))
                   onepmax = 1.01*(rc(i)+rc(j))

                   if (rr .lt. onepmax .and. rr .gt. onepmin) then !at contact
                       numgr = numgr + 1
                       ngrdat(numgr,1) = time
                       ngrdat(numgr,2) = real(idpc(i))
                       ngrdat(numgr,3) = real(idpc(j))
                       ngrdat(numgr,4) = rc(i)
                       ngrdat(numgr,5) = rc(j)
                       ngrdat(numgr,6) = wr
                    !  ngrdat(numgr,7) = (vzi-wdrop(j))
                    !  ngrdat(numgr,8) = rzij
                    !  ngrdat(numgr,9) = ryij
                    !  ngrdat(numgr,10)= rxij
                    !  ngrdat(numgr,11)= rr
                    !  ngrdat(numgr,12)= onepmax
                   endif

                   if (a .ne. 0.0 .and. b**2 .ge. 4.0*a*c) then
                       !if (b**2 .ge. 4.0*a*c) then
                           d = sqrt(b**2-4.0*a*c)
                           if (b .gt. 0.0) then
                                root1 = (-b-d)/(a+a)
                           else
                                root1 = (-b+d)/(a+a)
                           endif
                           root2 = (c/a)/root1
			  if(root1 .gt. 0 .and. root2 .gt. 0) then
				tc=min(root1,root2)
!				if (root1 .lt. root2) then
!				  tc = root1
!				else
!				  tc = root2
!				endif
!
                              if (tc .le. delt ) then
                                 nc = nc + 1
                                 ncdat(nc,1) = time
                                 ncdat(nc,2) = real(idpc(i))
                                 ncdat(nc,3) = real(idpc(j))
                                 ncdat(nc,4) = tc
                                 ncdat(nc,5) = xc(i)
                                 ncdat(nc,6) = yc(i)
                                 ncdat(nc,7) = zc(i)
                                 ncdat(nc,8) = udropc(i)
                                 ncdat(nc,9) = vdropc(i)
                                 ncdat(nc,10) = wdropc(i)
                                 ncdat(nc,11) = xc(j)
                                 ncdat(nc,12) = yc(j)
                                 ncdat(nc,13) = zc(j)
                                 ncdat(nc,14) = udropc(j)
                                 ncdat(nc,15) = vdropc(j)
                                 ncdat(nc,16) = wdropc(j)
                                 ncdat(nc,17) = rc(i)
                                 ncdat(nc,18) = rc(j)
                                 ncdat(nc,19) = real(i)
                                 ncdat(nc,20) = real(j)
                                 ncdat(nc,22) = r_ccnc(i)
                                 ncdat(nc,23) = r_ccnc(j)
                          endif!tc<delt

			endif!root1 and root2 >0
 !                      endif ! p*p .ge. 4.0*q*t

                   endif ! a .ne. 0.0

                   j = list(j)
                   GO TO 2000

                 endif ! 2000 - j > 0

!  ** LOOP OVER NEIGHBOURING CELLS **^M

                 jcell0 = 13 * (tcell - 1)

                 do 4000 nabor = 1, 13

                   jcell = map ( jcell0 + nabor )

!  ** LOOP OVER ALL MOLECULES IN NEIGHBOURING CELLS **^M

                   j = head(jcell)

 3000              if ( j .ne. 0 ) then

                       rxij = rxi - xc(j)

                       if (abs(rxij) .gt. 2*cellw) rxij = -( rxij/abs(rxij) ) * ( N*H - abs(rxij) )
                       ryij = ryi - yc(j)
                       rzij = rzi - zc(j)
                       if (abs(rzij) .gt. 2*cellw) rzij = -( rzij/abs(rzij) ) * ( N*H - abs(rzij) )
                       rr = sqrt(rxij**2+ryij**2+rzij**2)

                       vxij = vxi - udropc(j)
                       vyij = vyi - vdropc(j)
                       vzij = vzi - wdropc(j)
                       vv = sqrt(vxij**2+vyij**2+vzij**2)

                       a = vv**2
                       b = 2*(rxij*vxij+ryij*vyij+rzij*vzij)
                       c = rr**2 - (rc(i)+rc(j))**2

                       wr=(b/2.)/rr

                       onepmin = 0.99*(rc(i)+rc(j))
                       onepmax = 1.01*(rc(i)+rc(j))

                       if (rr .lt. onepmax .and. rr .gt. onepmin) then
                          numgr = numgr + 1
                          ngrdat(numgr,1) = time
                          ngrdat(numgr,2) = real(idpc(i))
                          ngrdat(numgr,3) = real(idpc(j))
                          ngrdat(numgr,4) = rc(i)
                          ngrdat(numgr,5) = rc(j)
                          ngrdat(numgr,6) = wr
                       !  ngrdat(numgr,7) = (vzi-wdrop(j))
                       !  ngrdat(numgr,8) = rzij
                       !  ngrdat(numgr,9) = ryij
                       !  ngrdat(numgr,10)= rxij
                       !  ngrdat(numgr,11)= rr
                       !  ngrdat(numgr,12)= onepmax
                       endif

                       if (a .ne. 0.0 .and. b**2 .ge. 4.0*a*c) then
			   !if (b**2 .ge. 4.0*a*c) then
				d = sqrt(b**2-4.0*a*c)
				  if (b .gt. 0.0) then
					root1 = (-b-d)/(a+a)
				  else
					root1 = (-b+d)/(a+a)
				  endif
				  root2 = (c/a)/root1
                              if (root1 .gt. 0.0 .and. root2 .gt. 0.0 ) then
				 tc=min(root1,root2)
!                                 if (root1 .lt. root2) then
!                                     tc = root1
!                                 else
!                                     tc = root2
!                                 endif

				 if (tc .lt. delt) then
                                      nc = nc + 1
                                      ncdat(nc,1) = time
          			      ncdat(nc,2) = real(idpc(i))
                                      ncdat(nc,3) = real(idpc(j))
                                      ncdat(nc,4) = tc
                                      ncdat(nc,5) = xc(i)
                                      ncdat(nc,6) = yc(i)
                                      ncdat(nc,7) = zc(i)
                                      ncdat(nc,8) = udropc(i)
                                      ncdat(nc,9) = vdropc(i)
                                      ncdat(nc,10) = wdropc(i)
                                      ncdat(nc,11) = xc(j)
                                      ncdat(nc,12) = yc(j)
                                      ncdat(nc,13) = zc(j)
                                      ncdat(nc,14) = udropc(j)
                                      ncdat(nc,15) = vdropc(j)
                                      ncdat(nc,16) = wdropc(j)
                                      ncdat(nc,17) = rc(i)
                                      ncdat(nc,18) = rc(j)
                                      ncdat(nc,19) = real(i)
                                      ncdat(nc,20) = real(j)
                                      ncdat(nc,22) = r_ccnc(i)
                                      ncdat(nc,23) = r_ccnc(j)
                              endif!tc<delt

                           endif!root1 and root2 >0

 !                          endif ! p*p .ge. 4.0*q*t

                       endif ! a .ne. 0.0

                       j = list(j)

                       GO TO 3000

                   endif ! 3000 - j .ne. 0

 4000              CONTINUE

                   i = list(i)

                   GO TO 1000

              endif ! 1000 - i > 0

 5000      CONTINUE

           deallocate(list)
!            if ( MPI == 1 ) then
              call mpi_Allgather(nc,1,MPI_INTEGER,nctot,1,MPI_INTEGER,MPI_COMM_WORLD,ierror)
              call mpi_Allgather(numgr,1,MPI_INTEGER,numgrtot,1,MPI_INTEGER,MPI_COMM_WORLD,ierror)
!           endif

           numcoll = sum(nctot)
           if (numcoll > 0) then

              if (colgr .ge. 1) then
                 kcol = 0
                   dcids = 0
                   ndc = 0

                 if (mype == 0) then
                   allocate( ncall(numcoll,ncele) )
                   allocate( nctemp(maxcol,ncele) )
                   allocate( dblcol(numcoll) )
                   allocate( sorttemp(ncele) )
                   dblcol = 0
                   if (nctot(1) .gt. 0) then
                       if (nc .ne. nctot(1) ) print*, 'nc not equal to nctot(1) on mype 0', nc, nctot(1)
                       kcol = kcol + nc
                       ncall(1:nc,:) = ncdat(1:nc,:)
                   endif
                 endif ! mype == 0

                 do i = 1,npe-1
                   if (nctot(i+1) .gt. 0) then
                       if (mype == 0) then
                           call mpi_recv(nctemp(1:nctot(i+1),:),nctot(i+1)*ncele,MPI_REAL,i,i,MPI_COMM_WORLD,status,istatus)
                           ncall(kcol+1:kcol+nctot(i+1),:) = nctemp(1:nctot(i+1),:)
                           kcol = kcol + nctot(i+1)
                       elseif (mype == i) then
                           call mpi_send(ncdat(1:nc,:),nc*ncele,MPI_REAL,0,i,MPI_COMM_WORLD,istatus)
                       endif
                       call mpi_barrier(MPI_COMM_WORLD,ierror)
                   endif
                 enddo ! i

                 if (mype == 0) then
                   if (numcoll > 1) then
                   !check which pair collides first sort out the collision time from small to large.
                       do i = 2, numcoll
                           j = i - 1
                           curtc = ncall(i,4)
                           sorttemp = ncall(i,:)
                           do while (j >= 1 .and. ncall(j,4) .gt. curtc)
                              ncall(j+1,:) = ncall(j,:)
                              j = j - 1
                           enddo
                           ncall(j+1,:) = sorttemp
                       enddo
                   endif !numcoll
                   do i = 1,numcoll-1
                       if ( dblcol(i) == 0) then
                           write(42,142) time, (int(ncall(i,l)),l=2,3), (ncall(i,zz),zz=4,18)
			   
                           do j = i+1,numcoll

                              if (int(ncall(i,2)) == int(ncall(j,2)) .or. int(ncall(i,2)) == int(ncall(j,3))) then
                                 print*, 'double collision with IDs:', int(ncall(i,2)), int(ncall(i,3)), int(ncall(j,2)), int(ncall(j,3))
                                 dblcol(j) = 1
                                     ndc = ndc + 1
                                     dcids(ndc,1) = int(ncall(j,2))
                                     dcids(ndc,2) = int(ncall(j,3))

                              elseif (int(ncall(i,3)) == int(ncall(j,2)) .or. int(ncall(i,3)) == int(ncall(j,3))) then
                                 print*, 'double collision with IDs:', int(ncall(i,2)), int(ncall(i,3)), int(ncall(j,2)), int(ncall(j,3))
                                 dblcol(j) = 1
                                     ndc = ndc + 1
                                     dcids(ndc,1) = int(ncall(j,2))
                                     dcids(ndc,2) = int(ncall(j,3))
                              endif!ncall
                           enddo !j
                       endif !dblcol
                   enddo !i
                   if (dblcol(numcoll) == 0) then ! still have to output data of last collision
                       write(42,142) time,(int(ncall(numcoll,l)),l=2,3), (ncall(numcoll,zz),zz=4,18)
                   endif

                   deallocate ( ncall, nctemp, dblcol, sorttemp )
                 endif ! mype == 0

              elseif (colgr == 0) then !colgr

                 if (mype == 0) then
                   if (nctot(1) .gt. 0) then
                       if (nc .ne. nctot(1) ) print*, 'nc not equal to nctot(1) on mype 0', nc, nctot(1)
                       do i = 1, nctot(1)
                           write(42,142) time,(int(ncdat(i,l)),l=2,3),(ncdat(i,zz),zz=4,18)
                       enddo
                   endif !nctot
                 endif !mype

                 do i = 1,npe-1
                   if (nctot(i+1) .gt. 0) then

                       if (mype == 0) then
                           call mpi_recv(ncdat(1:nctot(i+1),:),nctot(i+1)*ncele,MPI_REAL,i,i,MPI_COMM_WORLD,status,istatus)
                           do j = 1, nctot(i+1)
                              write(42,142) time,(int(ncdat(j,l)),l=2,3), (ncdat(j,zz),zz=4,18)
                           enddo
                       elseif (mype == i) then
                           call mpi_send(ncdat(1:nc,:),nc*ncele,MPI_REAL,0,i,MPI_COMM_WORLD,istatus)
                       endif !mype

                       call mpi_barrier(MPI_COMM_WORLD,ierror)
                   endif !nctot
                 enddo ! i

              endif ! colgr .ge. 1 or .eq. 0

              if (colgr .ge. 1) then
                 call mpi_Bcast(ndc,1,MPI_INTEGER,0,MPI_COMM_WORLD,istatus)
                 if (ndc > 0) then
                   call mpi_Bcast(dcids,10,MPI_INTEGER,0,MPI_COMM_WORLD,istatus)
                   if (nc > 0) then
                       do 220 j = 1, ndc
                           do 220 i = 1, nc
                              if (dcids(j,1) .eq. int(ncdat(i,2)) .and. dcids(j,2) .eq. int(ncdat(i,3))) then
                                 print*, 'located double collision on mype', mype, 'with ids:', dcids(j,1), dcids(j,2)
                                 ncdat(i,21) = 1
                              endif
 220                continue
                   endif ! nc > 0
                 endif ! ndc > 0
                 if (colgr == 1 .and. ireloc == 1) then
                 seed  = seedmic*(mype+1)+nt
                   do i = 1,nc !reallocate half of the collided droplets to random spots
                       if (ncdat(i,21) == 0) then
                           id = nint(ncdat(i,19))
                           x(id) = rangen(seed)*h*N1*1.0d0
                           y(id) = real(mype*h*n2pe) + rangen(seed)*h*n2pe*1.0d0
			               yrel(id) = y(id) - slabw*h
                           if (yrel(id) .le. 0.0) then ! HAS TOO SMALL Of A YREL at ini
                            yrel(id) = 0.0 - yrel(id) + h/100.0
                            y(id) = yrel(id) + slabw*h
                           elseif (yrel(id) .ge. real(n2pe)*h) then ! HAS TOO BIG OF A YREL at ini
                            yrel(id) =  2.0d0*h*real(n2pe) - yrel(id) - h/100.0
                            y(id) = yrel(id) + slabw*h
                           endif
                           z(id) = rangen(seed)*h*N3*1.0d0
			   udrop(id)=0.
			   vdrop(id)=0.
                           if (r(id) .lt. rtrunc_min) then !non-inertia particles
                             wdrop(id)=floww(id)
		           elseif (r(id) .lt. 4.d-5 ) then
		             wdrop(id) = -kwt*r(id)**2
		           else
		             wdrop(id) = -kwt2*r(id) + 0.1227
			   endif
                       endif!ncdat
                   enddo
                 endif !colgr == 1
         if (colgr == 2 ) then
                 ndrems = 0
                 cols = 0
                 allocate(icolcomms(nc+1), icolcommr(nctot(ntofrom(mype+3))+1), icol(nc+nctot(ntofrom(mype+3))+1))
                 icolcomms = 0
                 icolcommr = 0
                 icol = 0
                 if (nc > 0) then
                   allocate(rcol(nc), r_ccncol(nc))
                   allocate(xcol(nc), ycol(nc), zcol(nc), ucol(nc), vcol(nc), wcol(nc), idpcol(nc))
                   do i = 1,nc
                       if (ncdat(i,21) == 0) then
                           nddone = nddone + 1
                           cols = cols + 1
                           icol(cols) = int(ncdat(i,19))
                           idpcol(nddone) = int(ncdat(i,2))
                           r1cube = ncdat(i,17)**3
                           r2cube = ncdat(i,18)**3
                           r3cube = r1cube+r2cube
                           rcol(nddone) = r3cube**(1.d0/3.d0)
                           frac = r2cube/r3cube !ncdat(i,17)/(ncdat(i,17)+ncdat(i,18))
                           xcoll = ncdat(i,5) + ncdat(i,4)*ncdat(i,8)
                           pos = ncdat(i,11) + ncdat(i,4)*ncdat(i,14)
                           xcol(nddone) = xcoll + frac*(pos - xcoll)
                           ycoll = ncdat(i,6) + ncdat(i,4)*ncdat(i,9)
                           pos = ncdat(i,12) + ncdat(i,4)*ncdat(i,15)
                           ycol(nddone) = ycoll + frac*(pos - ycoll)
                           zcoll = ncdat(i,7) + ncdat(i,4)*ncdat(i,10)
                           pos = ncdat(i,13) + ncdat(i,4)*ncdat(i,16)
                           zcol(nddone) = zcoll + frac*(pos - zcoll)
                           ucol(nddone) = (ncdat(i,8)*r1cube + ncdat(i,14)*r2cube)/r3cube
                           vcol(nddone) = (ncdat(i,9)*r1cube + ncdat(i,15)*r2cube)/r3cube
                           wcol(nddone) = (ncdat(i,10)*r1cube + ncdat(i,16)*r2cube)/r3cube
                           !for merging the CCNs
                           r1cube = ncdat(i,22)**3
                           r2cube = ncdat(i,23)**3
                           r3cube = r1cube+r2cube
                           r_ccncol(nddone) = r3cube**(1.d0/3.d0)
                           
                           !then make the merged droplet go backwards by time of tc in order to go forwards for a full timestep together with other droplets in following step
                           xcol(nddone) = xcol(nddone)-ucol(nddone)*ncdat(i,4)
                           ycol(nddone) = ycol(nddone)-vcol(nddone)*ncdat(i,4)
                           zcol(nddone) = zcol(nddone)-wcol(nddone)*ncdat(i,4)
                           if ((ncdat(i,12) - slabw*h) > real(n2pe*h) ) then
                              ndrems = ndrems + 1 !# of droplets need to be removed in the right slab
                              icolcomms(ndrems) = int(ncdat(i,3))!id
                           else
                              cols = cols + 1 !# of droplets need to be removed in this slab
                              icol(cols) = int(ncdat(i,20))
                           endif
                       endif ! ncdat = 0
                   enddo ! i
                 endif ! nc > 0
                 nto = ntofrom(mype+3)
                 nfrom = ntofrom(mype+1)
                 call mpi_sendrecv(ndrems,1,MPI_INTEGER,nto,nto,ndremr,1,MPI_INTEGER,nfrom,mype,MPI_COMM_WORLD,status,ierror)
                 call mpi_sendrecv(icolcomms(1:ndrems),ndrems,MPI_INTEGER,nto,nto,icolcommr(1:ndremr),ndremr,MPI_INTEGER,nfrom,mype,MPI_COMM_WORLD,status,ierror) 
                 if (ndremr > 0) then
                   nrtest = 0
                   outer: do ii = 1,ndremr
                       inner: do id = 1,ndropreal
                           if (idp(id) .eq. icolcommr(ii)) then
                              nrtest = nrtest + 1
                              cols = cols + 1
                              icol(cols) = id
                              exit inner
                           endif
                       enddo inner ! id
                   enddo outer ! ii
                   if (nrtest .ne. ndremr) then
                       print*, 'nrtest .ne. ndremr', nrtest, ndremr
                       stop
                   endif
                 endif ! ndremr > 0
                 if (cols > 1) then !reorder large id follows small id
                   do i = 2, cols
                       j = i - 1
                       tempcol = icol(i)
                       do while (j >= 1 .and. icol(j) > tempcol)
                           icol(j+1) = icol(j)
                           j = j - 1
                       enddo
                       icol(j+1) = tempcol
                   enddo
                 endif ! cols
                 if (cols > 0) then
                   do i = cols,1,-1 !remove all merged droplets
                       if (icol(i) > ndropreal) then
                           print*, 'something odd, since icol(i) > ndropreal', icol(i), ndropreal
                       endif
                       idp(icol(i)) = idp(ndropreal)
                       x(icol(i)) = x(ndropreal)
                       y(icol(i)) = y(ndropreal)
                       z(icol(i)) = z(ndropreal)
                       udrop(icol(i)) = udrop(ndropreal)
                       vdrop(icol(i)) = vdrop(ndropreal)
                       wdrop(icol(i)) = wdrop(ndropreal)
                       r(icol(i)) = r(ndropreal)
                       r_ccn(icol(i)) = r_ccn(ndropreal)
                       ndropreal = ndropreal - 1
                   enddo
                   if (nc > 0) then !put back half of the merged droplets and put them at the end
                       idp(ndropreal+1:ndropreal+nddone)=idpcol(1:nddone)
                       x(ndropreal+1:ndropreal+nddone)=xcol(1:nddone)
                       y(ndropreal+1:ndropreal+nddone)=ycol(1:nddone)
                       z(ndropreal+1:ndropreal+nddone)=zcol(1:nddone)
                       udrop(ndropreal+1:ndropreal+nddone)=ucol(1:nddone)
                       vdrop(ndropreal+1:ndropreal+nddone)=vcol(1:nddone)
                       wdrop(ndropreal+1:ndropreal+nddone)=wcol(1:nddone)
                       r(ndropreal+1:ndropreal+nddone)=rcol(1:nddone)
                       r_ccn(ndropreal+1:ndropreal+nddone)=r_ccncol(1:nddone)
                       ndropreal = ndropreal + nddone
                   endif ! nc > 0
                 endif ! cols > 0
                 if (nc > 0) then
                   deallocate(rcol, r_ccncol)
                   deallocate(xcol, ycol, zcol, ucol, vcol, wcol, idpcol)
                 endif
                 deallocate (icolcomms, icolcommr, icol)

                 endif ! colgr == 2
              endif ! colgr .ge. 1
           endif ! numcoll

           numngr = sum(numgrtot)
           if (numngr > 0) then
              if (mype == 0) then
                 if (numgrtot(1) .gt. 0) then
                   if (numgr .ne. numgrtot(1) ) print*, 'numgr not equal to numgrtot(1) on mype 0', numgr, numgrtot(1)
                   do i = 1, numgrtot(1)
                       write(44,144) time,(int(ngrdat(i,l)),l=2,3),(ngrdat(i,zz),zz=4,ngrele)
                   enddo
                 endif
              endif ! mype == 0
              do i = 1,npe-1
                 if (numgrtot(i+1) .gt. 0) then
                   if (mype == 0) then
                       call mpi_recv(ngrdat(1:numgrtot(i+1),:),numgrtot(i+1)*ngrele,MPI_REAL,i,i,MPI_COMM_WORLD,status,istatus)
                       do j = 1, numgrtot(i+1)
                           write(44,144) time,(int(ngrdat(j,l)),l=2,3), (ngrdat(j,zz),zz=4,ngrele)
                       enddo
                   elseif (mype == i) then
                       call mpi_send(ngrdat(1:numgr,:),numgr*ngrele,MPI_REAL,0,i,MPI_COMM_WORLD,istatus)
                   endif
                   call mpi_barrier(MPI_COMM_WORLD,ierror)
                 endif
              enddo ! i
           endif ! numngr

           deallocate(ncdat )
           deallocate(ngrdat )
           deallocate(idpc, xc, yc, yrelc, zc, udropc, vdropc, wdropc, rc, r_ccnc)
           deallocate(uudropc,vvdropc,wwdropc)
           deallocate(dtovertaupc)
           deallocate(flowuc,flowvc,flowwc)

      endif ! Close loop of if (colcal .eq. 1) then ! update droplet positions

!!DSD
           if (colgr == 2 .or.  (thermo .eq. 1 .and. gomic .ge. 2)) then !dsd output
                if (mod(nt-1,dsdout) .eq. 0) then
                  if(mype .eq. 0) print*,'begin output dsd'
                  dsd = 0
!                  dsd_log2 = 0
                  do id = 1,ndropreal
                    entr = Nint(r(id)/1.d-6)
                    if (entr > nbins) then
                        print*, 'entr > nbins, so have to fix it', entr, nbins
                        entr = nbins
                    endif
                    !  if (entr > nbins) print*, 'entr > nbins', entr, nbins
		    if (entr .eq. 0) then
			dsd(1) = dsd(1) + 1
		    else
                        dsd(entr) = dsd(entr) + 1
		    endif
		  !   entr = Nint(log2(r(id)/1.e-6))
		  !  if (entr .lt. 0) then
		!	dsd_log2(1) = dsd_log2(1) + 1
		!    else
		!        dsd_log2(entr+1)=dsd_log2(entr+1)+1
		!    endif

                  enddo

                  call mpi_reduce(dsd,dsdtot,nbins,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierror)
!		  call mpi_reduce(dsd_log2,dsdtot_log2,nbins,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierror)

                  if (mype == 0) then
                     write(45,145) time, (dsdtot(i), i=1,nbins)
!		     write(46,145) time, (dsdtot_log2(i), i=1,nbins)
                  endif
                endif ! mod
           endif ! colgr ==2 thermo gomic 

!           deallocate( flowu, flowv, floww, dtovertaup, uudrop, vvdrop, wwdrop )
      ltmp = h*real(N1)
      do id = 1,ndropreal
           X(id) = X(id) + UDROP(id)*delt
           !dist = X(id) - dmid
           !if(abs(dist).ge.dmid) X(id) = X(id) - (abs(dist)/dist)*ltmp
           if(X(id) .ge. ltmp) x(id)=x(id)-ltmp
           if(x(id) .lt. 0.0 ) x(id)=x(id)+ltmp
           Y(id) = Y(id) + VDROP(id)*delt
           Z(id) = Z(id) + WDROP(id)*delt
           !dist = Z(id) - dmid
           !if(abs(dist).ge.dmid) Z(id) = Z(id) - (abs(dist)/dist)*ltmp
           if(z(id) .ge. ltmp) z(id)=z(id)-ltmp
           if(z(id) .lt. 0.0) z(id)=z(id)+ltmp
      enddo
      ! cross boundary movement
      allocate( infosendl(nawaymax,9), infosendr(nawaymax,9), inforecvl(nawaymax,9), inforecvr(nawaymax,9) )
      allocate( nawayid(nawaymax) )
      ymaxid = 0
      yminid = 0
      yrelmax = 0.0
      yrelmin = 0.0
      nsendl = 0
      nsendr = 0
      naway = 0
      nawayid = 0
      infosendl = 0.0
      infosendr = 0.0

      do 9000 id = 1,ndropreal
           yrel(id) = y(id) - slabw*h
           IF (yrel(id) .gt. yrelmax) then
              yrelmax = yrel(id)
              ymaxid = idp(id)
              if (yrelmax .gt. h*real(n2pe+2)) then
                 print*, 'MYPE', mype,'HAS TOO BIG A YREL! idp=',ymaxid,'yrel=',yrelmax,'over ',(yrelmax-(N*h/real(npe)+2*h))/h,'h,','udrop(id)=',udrop(id),vdrop(id),wdrop(id),'y',y(id),'r',r(id)
              endif
           ENDIF
           if (yrel(id) .lt. yrelmin) then
              yrelmin = yrel(id)
              yminid = idp(id)
              if (yrelmin .lt. -2*h) then
                 print*, 'MYPE', mype,'HAS TOO SMALL A YREL! idp=',yminid,'yrel=',yrelmin,'over ',yrelmin/h,'h,', 'udrop(id)=',udrop(id),vdrop(id),wdrop(id),'y',y(id),'r',r(id)
              endif
           endif
           !!!!left
           if (yrel(id) .lt. 0.d0) then
              if (mype .eq. 0) y(id) = y(id) + ltmp
              nsendl = nsendl + 1
              infosendl(nsendl,1) = real(idp(id))
              infosendl(nsendl,2) = x(id)
              infosendl(nsendl,3) = y(id)
              infosendl(nsendl,4) = z(id)
              infosendl(nsendl,5) = udrop(id)
              infosendl(nsendl,6) = vdrop(id)
              infosendl(nsendl,7) = wdrop(id)
              infosendl(nsendl,8) = r(id)
              infosendl(nsendl,9) = r_ccn(id)
              naway = naway + 1
              nawayid(naway) = id
              !!!!right
           elseif (yrel(id) .ge. ltmp/npe) then
              if (mype .eq. npe-1) y(id) = y(id) - ltmp
              nsendr = nsendr + 1
              infosendr(nsendr,1) = real(idp(id))
              infosendr(nsendr,2) = x(id)
              infosendr(nsendr,3) = y(id)
              infosendr(nsendr,4) = z(id)
              infosendr(nsendr,5) = udrop(id)
              infosendr(nsendr,6) = vdrop(id)
              infosendr(nsendr,7) = wdrop(id)
              infosendr(nsendr,8) = r(id)
              infosendr(nsendr,9) = r_ccn(id)
              naway = naway + 1
              nawayid(naway) = id
           endif ! yrel < 0.0 or yrel > slabw*h 
9000  continue
      if (naway .gt. nawaymax) then
           print*, 'mype', mype, 'has too large naway = ', naway, nawaymax
           stop
      endif

      if (naway .ne. nsendl+nsendr) print*,'naway not equal to nsendl+nsendr by mype',mype

! Transferring droplets that cross slab edges boundaries 
! Communicate with left neighbour to let him know how many droplets are leaving me to go to him and also communicate 
! with right neighbour to know how many droplets are coming from him to me. This call is then followed by the actual 
! transfer of the droplets. Leftward motion
      do i=1,npe
         ntofrom(i+1) = i-1
      enddo
      ntofrom(1) = npe -1
      ntofrom(npe+2) = 0
      nto = ntofrom(mype+1)
      nfrom = ntofrom(mype+3)
      call mpi_sendrecv(nsendl,1,MPI_INTEGER,nto,nto,nrecvr,1,MPI_INTEGER,nfrom,mype,MPI_COMM_WORLD,status,ierror)

      nbufsendl = nsendl*9
      nbufrecvr = nrecvr*9

      call mpi_sendrecv(infosendl(1:nsendl,:),nbufsendl,MPI_REAL,nto,nto,inforecvr(1:nrecvr,:),nbufrecvr,MPI_REAL,nfrom,mype,MPI_COMM_WORLD,status,ierror)

      if (nrecvr .ge. 1) then
           idp(ndropreal+1:ndropreal+nrecvr) = int(inforecvr(1:nrecvr,1))
           x(ndropreal+1:ndropreal+nrecvr) = inforecvr(1:nrecvr,2)
           y(ndropreal+1:ndropreal+nrecvr) = inforecvr(1:nrecvr,3)
           z(ndropreal+1:ndropreal+nrecvr) = inforecvr(1:nrecvr,4)
           udrop(ndropreal+1:ndropreal+nrecvr) = inforecvr(1:nrecvr,5)
           vdrop(ndropreal+1:ndropreal+nrecvr) = inforecvr(1:nrecvr,6)
           wdrop(ndropreal+1:ndropreal+nrecvr) = inforecvr(1:nrecvr,7)
           r(ndropreal+1:ndropreal+nrecvr) = inforecvr(1:nrecvr,8)
           r_ccn(ndropreal+1:ndropreal+nrecvr) = inforecvr(1:nrecvr,9)
      endif
      ndropreal = ndropreal + nrecvr
! Communicate with right neighbour to let him know how many droplets are leaving me to go to him and also communicate 
! with left neighbour to know how many droplets are coming from him to me. This call is the followed by the actual 
! transfor of the droplets. Rightward motion
      nto = ntofrom(mype+3)
      nfrom = ntofrom(mype+1)
      call mpi_sendrecv(nsendr,1,MPI_INTEGER,nto,nto,nrecvl,1,MPI_INTEGER,nfrom,mype,MPI_COMM_WORLD,status,ierror)

      nbufsendr = nsendr*9
      nbufrecvl = nrecvl*9

      call mpi_sendrecv(infosendr(1:nsendr,:),nbufsendr,MPI_REAL,nto,nto,inforecvl(1:nrecvl,:),nbufrecvl,MPI_REAL,nfrom,mype,MPI_COMM_WORLD,status,ierror)

      if (nrecvl .ge. 1) then
           idp(ndropreal+1:ndropreal+nrecvl) = int(inforecvl(1:nrecvl,1))
           x(ndropreal+1:ndropreal+nrecvl) = inforecvl(1:nrecvl,2)
           y(ndropreal+1:ndropreal+nrecvl) = inforecvl(1:nrecvl,3)
           z(ndropreal+1:ndropreal+nrecvl) = inforecvl(1:nrecvl,4)
           udrop(ndropreal+1:ndropreal+nrecvl) = inforecvl(1:nrecvl,5)
           vdrop(ndropreal+1:ndropreal+nrecvl) = inforecvl(1:nrecvl,6)
           wdrop(ndropreal+1:ndropreal+nrecvl) = inforecvl(1:nrecvl,7)
           r(ndropreal+1:ndropreal+nrecvl) = inforecvl(1:nrecvl,8)
           r_ccn(ndropreal+1:ndropreal+nrecvl) = inforecvl(1:nrecvl,9)
      endif

      ndropreal = ndropreal + nrecvl

      do i = naway,1,-1
           idp(nawayid(i)) = idp(ndropreal)
           x(nawayid(i)) = x(ndropreal)
           y(nawayid(i)) = y(ndropreal)
           z(nawayid(i)) = z(ndropreal)
           udrop(nawayid(i)) = udrop(ndropreal)
           vdrop(nawayid(i)) = vdrop(ndropreal)
           wdrop(nawayid(i)) = wdrop(ndropreal)
           r(nawayid(i)) = r(ndropreal)
           r_ccn(nawayid(i))=r_ccn(ndropreal)
           ndropreal = ndropreal - 1
      enddo

!mark test if any index is out of bound
      ixp(1:ndropreal)=int(x(1:ndropreal)*oneoverh) + 1
      iyp(1:ndropreal)=int(y(1:ndropreal)*oneoverh) - slabw + 1
      izp(1:ndropreal)=int(z(1:ndropreal)*oneoverh) + 1

      do id = 1,ndropreal !mark

          if(ixp(id) .eq. 0) then
                  ixp(id)=1
          elseif(ixp(id) .eq. N+1) then
                  ixp(id)=N
          endif

          if(izp(id) .eq. 0) then
                  izp(id)=1
          elseif(izp(id) .eq. N+1) then
                  izp(id)=N
          endif

          if(iyp(id) .eq. 0)then
                  iyp(id)=1
          elseif(iyp(id) .eq. 3) then
                  iyp(id)=2
          endif

!          if(min(ixp(id),izp(id),iyp(id)) .le. 0 .or. max(ixp(id),izp(id)) .gt. N .or. iyp(id) .gt. 2)  then
!                  print*,'after adjustment x,y,z index of drop location, id, r,dr3'
!                  print*, ixp(id),iyp(id),izp(id), x(id), y(id)-mype*n2pe*h , mype, z(id), id, r(id), dr3(id)!markd
      !    endif
      enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Calculate max-CFL number for droplets based on maximum droplet velocity.
      if (mod(nt-1,micdtout) .eq. 0) then
           maxvel = 0.0
           maxvel = max(maxvel, maxval(abs(udrop(1:ndropreal))), &
                          maxval(abs(vdrop(1:ndropreal))), &
                          maxval(abs(wdrop(1:ndropreal))))
!           if ( MPI == 1 ) then
              call mpi_reduce(maxvel,tmpv,1,MPI_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierror)
              maxvel=tmpv
    !  call mpi_reduce(vmax,tmpv,1,MPI_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierror)
    !  vmax=tmpv

!           endif
           if (mype .eq. 0) then
              print*,' '
              print*,'***MAX DROPLET COURANT NUMBER : ', maxvel*delt/h,maxvel
              print*,' '
!              write(97,*) '***MAX DROPLET COURANT NUMBER : ', maxvel*delt/h,maxvel
           endif

      endif 
      if (ndropreal .ge. int(1.25*real(ndropp))) then
           print*, 'mype', mype, 'has ndropreal > 1.25*ndropp'
      endif

      if (ndropreal .ge. int(1.5*real(ndropp))) then
           print*, 'mype', mype, 'has ndropreal > 1.5*ndropp'
      endif

      if (ndropreal .ge. int(1.9*real(ndropp))) then
           print*, 'mype', mype, 'has ndropreal > 1.9*ndropp'
      endif

      call mpi_reduce(ndropreal,ndroptot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierror)

      if (mype .eq. 0 .and. ndroptot .ne. ndrop .and. gomic .le. 1) then
           print*, '------------------------------'
           print*, ' ndroptot not equal to ndrop'
           print*, ' ndroptot = ', ndroptot, 'ndrop = ', ndrop
           print*, '------------------------------'
           stop
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!---------------------Adaptive time step--------------------------!!!!!!!!!
!      if(gomic .ge. 2 .and. mod(nt,micdtout) .eq. 0) then 
!         if(ndropreal .gt. 0) rmin=minval(r(1:ndropreal))
!         if(MPI .eq. 1) then
!           call mpi_allreduce(rmin,rtmp,1,mpi_real,mpi_min,mpi_comm_world,ierror)
!           rmin=rtmp
!         endif
!         if(delt*1.05 .lt.  min(0.15d0*kwt*rmin**2/real(grav),0.01d0*(EDR*real(N)*h)**(-1./3.)*h)) delt=1.05*delt
!         if(mype .eq. 0) print*, 'delt', delt
!       endif !gomic
!------------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! OUTPUT microphysical restart file
      if (mod((nt-1),ndumpd) .eq. 0) then
           ntdump=(nt-1)/ndumpd
           if (netcdf .eq. 1) then
              call ncdumpdrop(ndropreal,ntdump,time,nddone )
              if (mype .eq. 0) print*, 'time of ncdumpdrop = ',time
!		if (mype .eq. 0) print*, 'PP,rhoa,thetapp',PP,rhoa,thetapp
           else
              print*, 'Sorry, do not know how to dump droplet restart files without NetCDF'
           endif
      endif ! ndump


      deallocate( nawayid )
!      deallocate( ixp, iyp, izp )
      deallocate( infosendl, infosendr, inforecvl, inforecvr )

      RETURN

      END

! Integer function ICELL(jx,jy,jz,m,my) 
! Calculate cell-number of cell with indices jx,jy,jz of a block that has m-cells 
! in x- and z-direction and my-cells in y-direction

      integer function ICELL(jx,jy,jz,m,my)

           implicit none
           integer, intent(in) :: jx,jy,jz,m,my
           ICELL = 1 + MOD ( jx - 1 + M, M )  &
                              + MOD ( jz - 1 + M, M ) * M &
                              + MOD ( jy + MY, MY + 1 ) * M * M

      end function ICELL


