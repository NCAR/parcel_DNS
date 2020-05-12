! triply periodic boussinesq with centred difference and crank-nicholson damping
! original code by peter bartello
! subsequent modification by paul vaillancourt and charmaine franklin
! further modified by kevin zwijsen 
! f90, mpi: fall 2010

! set preprocessor flags:
! netcdf: 1 to use netcdf i/o.  note: at the moment, this is required.
! mpi: 1 to compile with mpi; 0 for serial/openmp

!     notes: 
!     kx,ky,kz are wavenumbers. 
!     ikx,y,z  are corresponding array indices.
!     ktx,y,z  are truncation wavenumbers.  
!     iktx,y,z are corresponding indices.
!     zx,y,z is vorticity. u,v,w is velocity. 
!     nzx,y,z is nonlinear term in vorticity equation.  
!     complex stuff is in k-space.
!     suffix o/n means old/new for leap-frog     


program main 

  use dyn_mod
!  use fred_mod
  use scratch_mod
  use mic_mod
  use netcdf_mod
  implicit none 

  include 'param.inc'
!  include 'fftw3.f'

  complex, target, allocatable, dimension(:,:,:) :: zxok,zyok,zzok,ttok,qvok,zxnk,zynk,zznk,ttnk,qvnk
  complex, target, allocatable, dimension(:,:,:) :: nzxk,nzyk,nzzk,nttk,nqvk,uk,vk,wk
  complex, target, allocatable, dimension(:,:,:) :: rhzx,rhzy,rhzz,rhtt,rhqv,ak,bk,ck
  real,    target, allocatable, dimension(:,:,:) :: mdiff,tdiff,qdiff

  complex :: avzx,avzy,avzz,avtt,avqv
  complex :: termzx,termzy,termzz,termtt,termqv,tzx,tzy,tzz,ttt,tqv
  complex :: u,v,w,c1,c2,c3,zi

  real, allocatable, dimension(:,:,:) :: urt,vrt,wrt
  real, allocatable, dimension(:,:) :: spz
  real, allocatable, dimension(:) :: kxa,kya,kza
  real :: efor,bf2,ke,pe,etot,robert,r1,r2
  real :: ts=0.,time
  real(8) :: d2, delt_drop, delt_drop10,delt_drop15,delt_drop20,delt_dyn
  real(8) :: cql
  real :: time1(2),time2(2),time3,etime
  real :: k,kx,ky,kz,wk2,kfmin,kfmax,ktrunc
  real :: tmp,vz,vvk,kkol

  integer :: ikx,iky,ikz,ikza
  integer, allocatable, dimension(:,:,:) :: l
  integer, allocatable, dimension(:) :: ns !?
  integer, allocatable, dimension(:) :: ikxf,ikyf,ikzf,ikzfp

  integer :: nproj
  integer :: ndump,ndumpd,nbig,rsflag,nrsp,ntdump,ntdumpd!,nrst,nrstd
  integer :: jw,irest,seed,ampfor
  integer :: nf,nftot
  integer :: i,linear,ilap2,void
  integer :: nt=0!,io=0

  character*11 :: nameexp 
 
! declarations for microphysical model
  real :: time0,dimk,dxp,dyp,dzp	
!  integer, allocatable, dimension(:) :: ixp,iyp,izp
  integer :: id,ndropreal,nddone

! fftw stuff
  integer, parameter :: fftw_forward         = -1, fftw_backward                  = 1
  integer, parameter :: fftw_real_to_complex = -1, fftw_complex_to_real           = 1
  integer, parameter :: fftw_estimate        =  0, fftw_measure=1,fftw_use_wisdom = 16
  integer, parameter :: fftw_in_place        =  8
  integer*8 :: plan2_rk,plan2_kr,plan1_rk,plan1_kr
  common/ftw2d/plan2_rk,plan2_kr
  common/ftw1d/plan1_rk,plan1_kr

  integer, parameter :: iktyp=ikty/npe


!if (netcdf == 1) then
  include 'netcdf.inc'


! mpi stuff
  include 'mpif.h'
  integer mype,size,ierror
  common/mpi/ mype 

  allocate( zxok(iktx,ikty,iktzp), zyok(iktx,ikty,iktzp), zzok(iktx,ikty,iktzp), &
          ttok(iktx,ikty,iktzp), qvok(iktx,ikty,iktzp), &
          zxnk(iktx,ikty,iktzp), zynk(iktx,ikty,iktzp), zznk(iktx,ikty,iktzp), &
          ttnk(iktx,ikty,iktzp), qvnk(iktx,ikty,iktzp) )
  allocate( nzxk(iktx,ikty,iktzp), nzyk(iktx,ikty,iktzp), nzzk(iktx,ikty,iktzp), &
          nttk(iktx,ikty,iktzp), nqvk(iktx,ikty,iktzp), &
          uk(iktx,ikty,iktzp), vk(iktx,ikty,iktzp), wk(iktx,ikty,iktzp) )
  allocate( mdiff(iktx,iktz,iktzp), tdiff(iktx,iktz,iktzp), qdiff(iktx,iktz,iktzp) )
  allocate( rhzx(iktx,iktz,iktzp), rhzy(iktx,iktz,iktzp), rhzz(iktx,iktz,iktzp), rhtt(iktx,iktz,iktzp), &
          rhqv(iktx,iktz,iktzp), ak(iktx,iktz,iktzp), bk(iktx,iktz,iktzp), ck(iktx,iktz,iktzp) )

  allocate( urt(n1+2,n2+2,n3+2),vrt(n1+2,n2+2,n3+2),wrt(n1+2,n2+2,n3+2) )

  allocate( x(ndroppe), y(ndroppe), z(ndroppe), udrop(ndroppe), vdrop(ndroppe), wdrop(ndroppe) )
  allocate( r(ndroppe), r_ccn(ndroppe))

!  allocate( e1(5000), e2(5000), e3(5000), e4(5000), e5(5000), e6(5000), e7(5000), e8(5000) )

  allocate( spz(ktx,4) )
  allocate( kxa(iktx),kya(iktz),kza(iktz) )
  allocate( l(iktx,ikty,iktzp) )
  allocate( ns(ktx) )
  allocate( ikxf(nfmax),ikyf(nfmax),ikzf(nfmax),ikzfp(nfmax) )

  allocate( zkt(iktx,iktz,iktyp), zk1(iktx,ikty,iktzp), zk2(iktx,ikty,iktzp) )


! ----------------------------------------------------------------------------------------
! initialize mpi
! ----------------------------------------------------------------------------------------
  call mpi_init(ierror)
  call mpi_comm_size(mpi_comm_world,size,ierror)
  call mpi_comm_rank(mpi_comm_world,mype,ierror)
  if (size .ne. npe) then
    print*, 'size',size,'npe',npe,'wrong number of processors!'
    stop
  endif

! ----------------------------------------------------------------------------------------
! parameters
! ----------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------
  nameexp=  'RUN06262018'				! name of output files
!  gomic= 1   !microphysics flag (0 = no mic.phy., 1 = no mic. restart file, >=2 = read mic. restart-file)
!  iturb = 1  ! 0 = turbulence off; 1 = on !!only allow 0 or 1 because of the constraints in outputting.F90
  colcal =   int(gomic/2)           ! collision calculation flag. if = 0 ignore collision calculations.
  colgr  =   2*colcal               ! Collisional growth flag. 0 = do nothing, 1 = double-collision check, 2 = collisional growth 
  ihydro   =  int(gomic/2)          ! =1 hydro =0 nonhydro
  isolu    = 1                      ! =1 solute effect on =0 off
  edr      = 0.05                   ! expected dissipation rate
  nproj    =  1500                  ! call proj every nproj timesteps.
  nout     =  1500                  ! frequency of output from turbulence model, such as spectra. output every nout timesteps.
  micdtout =  1500                  ! frequency of output from selected droplets.
  dsdout   =  3000                  ! frequency of dsd output
!  nstop    =  int(real(N)*h/(kwt*rad2**2*delt)*1.1) !avoid periodic collision in gravitational case
  nstop    = 320000                ! number of timesteps in current restart block.
  rsflag   = 0                      ! flag for outputting real space fields
! ----------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------
  if (iturb .eq. 0) then 
          irest    = gomic-1  				! current or initial restart block. if = 0 then turbulence does not read restart file
  else
          irest = gomic
  endif
  h        = PI/1.3d0*(visc**3./edr)**.25d0 !dx
  delt_drop=    0.1d0*kwt*rtrunc_min**2/real(grav)
  delt_dyn = 0.06d0*(edr*real(N)*h)**(-1./3.)*h
  if(gomic .eq. 0 ) then
     delt = delt_dyn
  else
     delt = min(delt_drop,delt_dyn) ! actual timestep of simulation is minimum of cfl-based timesteps for droplets and dynamics.
  endif
  if(mype .eq. 0) print*,'delt',delt,'delt_drop',delt_drop,'delt_dyn',delt_dyn
  ndump    = nstop/nrst             ! dump restart file every ndump timesteps
  ndumpd   = nstop/nrstd            ! dump droplet restart file every ndumpd timesteps
  nbig     = nstop/nrst             ! dump real space fields every nbig timesteps
  nrsp     = rsflag*(nstop/nbig+1)  ! number of real space fields to dump
  etot     = edr*delt*20.d0
  pe       =   0.5d0*etot           ! potential energy
  pe       = pe*real(iturb)
  ke       =   0.5d0*etot           ! kinetic energy
  ke       = ke*real(iturb)
  ampfor   = iturb                  ! forcing flag. if = 1 means forcing is applied at every timestep to maintain energy.
  linear   =   0                    ! linear flag. if = 1 ignore nonlinear terms in vorticity equation.
  kfmin    =   0.0                  ! minimum wavenumber over which forcing is applied.
  kfmax    =   1.50                 ! maximum wavenumber over which forcing is applied.
  seed     =   seedturb*(mype+1)    ! seeding of random number generator.
  ktrunc   =   float(n)/2.          ! truncation wavenumber. when = float(n)/3 then aliasing is removed.
  dimk     =   2.d0*pi/(n*h)        ! dimensions of wavenumber k
  !dimu     =   1.		    !dimensions of velocity u
  !dimt     =   1.		    ! dimensions of temperature t
  !diml     =   1.		    ! dimensions of length l
  ilap2    =   2                    ! order of laplacian dissipation/hyperviscosity (2 = regular viscosity)
  robert   =   0.001                ! robert filter to ensure coupling between subsequent timesteps in leap-frog
!  expsiv   =   grav/283.15d0* 0.d0  ! = g/T0. thermal expansivity. put to zero for no buoyancy.
!  temgrad  =   expsiv * 0.d0        ! vertical temperature gradient
  zi       =   cmplx(0.,1.)         ! complex parameter i
  d2       =   2.d0*delt            ! 
! -----------------------------------------------------------------------------------  
! open files for diagnostics and microphysics
! -----------------------------------------------------------------------------------
  if (mype .eq. 0) then
    jw = nstop/nout + 1
    print*, 'files will contain',jw,' output times.'
    open(unit=14,file='exp'//nameexp,form='formatted',position='append')
    open (41,file=nameexp//'.track',form='formatted',position='append')
    open (42,file=nameexp//'.coldat',form='formatted',position='append')
!    open (43,file=nameexp//'.ngrdat',form='formatted',position='append')
    open (44,file=nameexp//'.dispdat',form='formatted',position='append')
    open (45,file=nameexp//'.dsd',form='formatted',position='append')
    open (46,file=nameexp//'.dsdlog',form='formatted',position='append')
    open (57,file=nameexp//'.list',form='formatted',position='append')
    open (80,file=nameexp//'.spc',form='formatted',position='append')
    open (81,file=nameexp//'.eng',form='formatted',position='append')
    open (82,file=nameexp//'.nc',form='formatted',position='append')
!    open (93,file=nameexp//'.pdfvort',position='append')
!    open (95,file=nameexp//'.pdfvel',position='append')
!    open (97,file=nameexp//'.info',position='append')
    open (98,file=nameexp//'.para',position='append')
!    open (94,file=nameexp//'.velorvort',position='append')
!    open (61,file=nameexp//'.np',form='formatted',position='append')
    print*, 'done opening files'
  endif ! mype 

! ----------------------------------------------------------------------------------------------
! test for some requirements
! ----------------------------------------------------------------------------------------------

  if (mype .eq. 0) then
     if (mod(n,npe) .ne. 0) then
        print*, 'number of gridpoints n must be an exact multiple of number of processors'
        stop
     endif
     if (mod(ndrop,npe) .ne. 0) then
        print*, 'must have equal amount of droplets on each processor at the beginning'
        stop 
     endif
     if (mod(n,2*npe) .ne. 0) then
        print*, 'must have even number of N for npe'
        stop
     endif
  endif 

! ----------------------------------------------------------------------------------------------
! set up a few last things
! ----------------------------------------------------------------------------------------------

! initialize netcdf
if (netcdf ==1) then
  call ncprep(nrsp)
endif
  if (mype.eq.0) print*,'initializing fftw'

  call rfftw2d_f77_create_plan(plan2_rk,n1,n3,fftw_real_to_complex,fftw_measure+fftw_in_place)
  call rfftw2d_f77_create_plan(plan2_kr,n1,n3,fftw_complex_to_real,fftw_measure+fftw_in_place)
  call    fftw_f77_create_plan(plan1_rk,n2,   fftw_forward,        fftw_measure+fftw_in_place)
  call    fftw_f77_create_plan(plan1_kr,n2,   fftw_backward,       fftw_measure+fftw_in_place)

! ----------------------------------------------------------------------------------------------
! initialize fields, l, etc.
! ----------------------------------------------------------------------------------------------
  if (mype .eq. 0) write(*,*) "going in initturb &thermo"
!  div = 0.
  allocate( uvecs(n1,n3,3),uvecr(n1,n3,3) )
  call initturb(seed,zxok,zyok,zzok,ttok,qvok,ke,pe,l,irest,kxa,kya,kza,    &
                 ktrunc,ts,dimk)
  if (mype .eq. 0) write(*,*) "done with initturb &thermo"

! ----------------------------------------------------------------------------------------------
! initialize microphysics.
! ----------------------------------------------------------------------------------------------
  if (gomic .gt. 0) then
     allocate( idp(ndroppe), ixp(ndroppe), iyp(ndroppe), izp(ndroppe) )
     allocate( flowu(ndroppe), flowv(ndroppe), floww(ndroppe),dtovertaup(ndroppe), uudrop(ndroppe), vvdrop(ndroppe), wwdrop(ndroppe))
     allocate(fv(ndroppe))
     allocate(dr3(ndroppe))
     allocate(yrel(ndroppe))
     allocate(idpfol(ndrop/everynd))
     allocate(xfol(ndrop/everynd),yfol(ndrop/everynd),zfol(ndrop/everynd))
     allocate(rfol(ndrop/everynd))
     allocate(infocols(ndroppe,17),infocolr(ndroppe,17))
     allocate(head(ncell+m*m))
     allocate(maxutemp(npe))
     allocate(dcids(5,2))
     allocate(dsd(nbins),dsd_log2(nbins))
     cellw    =   N*H/real(m) !2delx
     slabw    =   real(mype*n2pe) !starting slab of each mype
     vol      =   h**3
     if(mype .eq. 0 ) allocate(dsdtot(nbins),dsdtot_log2(nbins))
     call initmicro(time0,ndropreal,nddone)
  endif ! gomic


! ------------------------------------------------------------------------------
! print out all parameters and write some information about simulation to file
! ------------------------------------------------------------------------------

  if (mype .eq. 0) then 
         !   open(unit=14,file='exp'//nameexp,position='append')
            write(14,*)   'dim model flag on or off(0)'
            write(14,*)  'sedimentation=', kwt
            write(14,*)  'turbulence is', iturb
            write(14,*)  'microphysical startup', gomic
            write(14,*)  'hydrodynamic effect=', ihydro
            write(14,*) 'disperse system=', disp
            write(14,*) 'condensation =', thermo
            write(14,*) 'cooling term =', gammand
            write(14,*) 'buoyancy term =', expsiv
            write(14,*) 'collision mech colgr=', colgr
               !  write(14,*)  'colgr = 0 : just collision statistics'
                ! write(14,*)  'colgr = 1 : double collision check'
                ! write(14,*)  'colgr = 2 : double collision check + collisional growth'
                    write(14,*) 'irest = ',irest
                    write(14,*) 'gomic = ',gomic
                    write(14,*) '----------------------------------------'
                    write(14,*) '                '
                    write(14,*) '                '
! check if criteria for crank nicolson is satisfied
!       r1 = delt*visc*(dimk*ktrunc)**ilap2
!       if(gomic.ne.0) r1 = kappand*delt*(dimk*ktrunc)**ilap2
!       if(r1.ge.1.and.ilap2.eq.2) then
!          print*,'*****problem: must decrease time step : ',r1
!          stop
!        endif

  time = ts + nt*delt
  bf2 = expsiv*temgrad

    print*,'                '
    print*,'                '
    print*,'parameters: ----------------------------------------'
    print*,'boussinesq'
    print*,'                      n = ',n
    print*,'                     kt = ',ktrunc
    print*,'              viscosity = ',visc
    if (linear .eq. 1 .and. iturb .ne. 0) print*,'nonlinear terms switched off'
    print*,'               timestep = ',delt
    print*,'     integration length = ',nstop*delt, ' = ',nstop,' delt.'
    print*,'       output frequency = ',nout*delt, ' = ',nout,' delt.'
    if (irest.eq.0 .and. iturb .ne. 0) then
      print*,' starting from ics '
      print*,' initial kinetic energy = ',ke
      print*,'   "    potential  "    = ', pe
      print*,'   "      total    "    = ', ke + pe
      print*,'        random i.c.s '
      print*,'        random no. seed = ',seed
    elseif (irest .ne. 0) then
      print*,'  starting from restart file'
      print*,'  restart from time     = ', ts
    endif !irest iturb
     print*,'    thermal expansivity = ', expsiv
     print*,'    vertical t gradient = ', temgrad
     print*,'    brunt-vaisala freq. = ', sqrt(bf2)
     print*,'robert filter parameter = ', robert

     print*,'Number of droplets :',ndrop
     write(14,*)   'From Microetturb'
     write(14,111) 'Initial time =      ',time0
     write(14,111) 'Final time =        ',time0+delt*nstop
     write(14,112) 'Timestep =          ',delt
     write(14,111) 'Vertical velocity = ',UP
     write(14,111) 'Vertical distance = ',delt*nstop*UP
     write(14,122) 'Domain size (m) =   ',n1*h,n2*h,n3*h
     write(14,177) '# of gridpoints =   ',n1,n2,n3
     write(14,133) 'Number of droplets = ',ndrop
     write(14,*)   'Seed =               ',seed
     write(14,*)   'kernel constants   =',(real(n1)*h)**3/((real(ndrop)/2)**2*(time0+delt*nstop))

 111  format(1x,a20,f12.4)
 112  format(1x,a20,e10.3)
 122  format(1x,a20,f5.3,' by ',f5.3,' by ',f5.3)
 133  format(1x,a21,i10)
 177  format(1x,a20,i3,' by ',i3,' by ',i3)

  endif!mype

  if (ampfor .eq. 1.) then
       if (mype .eq. 0) print*,' forcing active over ',kfmin,' <= k =< ',kfmax       
       call forset(kfmin,kfmax,ikxf,ikyf,ikzf,ikzfp,nf,l,kxa,kya,kza,dimk)
       efor = 0.0
       do  i=1,nf
         kx = kxa(ikxf(i))
         ky = kya(ikyf(i))
         kz = kza(ikzf(i))
         wk2  =  kx*kx+ky*ky+kz*kz

         vz = real( zxok(ikxf(i),ikyf(i),ikzfp(i))*conjg(zxok(ikxf(i),ikyf(i),ikzfp(i))) )
         efor = efor + vz/wk2
         vz = real( zyok(ikxf(i),ikyf(i),ikzfp(i))*conjg(zyok(ikxf(i),ikyf(i),ikzfp(i))) )
         efor = efor + vz/wk2
         vz = real( zzok(ikxf(i),ikyf(i),ikzfp(i))*conjg(zzok(ikxf(i),ikyf(i),ikzfp(i))) )
         efor = efor + vz/wk2
       enddo

!         if ( mpi == 1 ) then
          call mpi_allreduce(efor,tmp,1,mpi_real,mpi_sum,mpi_comm_world,ierror)
          efor=tmp
          call mpi_reduce(nf,nftot,1,mpi_integer,mpi_sum,0,mpi_comm_world,ierror)
!         else
!          nftot=nf
!         endif
         if (mype .eq. 0) then
          print*,'        initial energy in forced band = ', efor
          print*,'        there are ',nftot,' forced modes.'
         endif
      else
        if (mype .eq. 0) print*,'        no forcing.'
  endif ! ampfor

      if (mype .eq. 0) then
        print*,'                '
        print*,'                '
        print*,'                '
      endif

  if (mype == 0) print*, 'ndroppe = ', ndroppe

! ----------------------------------------------------------------------------------------------
! diagnostics on ics
! ----------------------------------------------------------------------------------------------

  if(thermo .eq. 1 .or. iturb .eq. 1) then
     if (mype .eq. 0) print*, 'diagnostics in initial conditions'
     call spec(zxok,zyok,zzok,ttok,qvok,ns,spz,l,time,kxa,kya,kza,vvk,kkol,dimk)
     call out (zxok,zyok,zzok,ttok,qvok,uk,vk,wk,time,nstop,ns,spz,l,kxa,kya,kza,vvk,kkol,dimk)
     if (mype .eq. 0)     print*,'  '
     if(iturb .eq. 1)      call rspace (zxok,zyok,zzok,ttok,zxok,zyok,zzok,ttok,nzxk,nzyk,nzzk,nttk,nzxk,nzyk,nzzk,nttk, &
                  uk,vk,wk,uk,vk,wk,l,zi,kxa,kya,kza,time)
     if (rsflag .eq. 1) then
        if (netcdf .eq. 1) then 
           call dumpreal(zxok,zyok,zzok,ttok,qvok,zxok,zyok,zzok,ttok,qvok,uk,vk,wk, &
                          uk,vk,wk,nzxk,l,kxa,kya,kza,1,time,nrsp)
        else
           print*, 'sorry, do not know how to write real space fields without netcdf'
        endif  
     endif

     if (mype .eq. 0 .and. iturb .eq. 1) print*,'vorticity field divergence:'
     if (iturb .eq. 1) call diverg (zxok,zyok,zzok,l,kxa,kya,kza)


!                        first timestep

     if (iturb .eq. 0) then
        nzxk = cmplx(0.,0.)
        nzyk = cmplx(0.,0.)
        nzzk = cmplx(0.,0.)
        zxok = cmplx(0.,0.)
        zyok = cmplx(0.,0.)
        zzok = cmplx(0.,0.)
     endif
     if (thermo .eq. 0) then
        nttk = cmplx(0.,0.)
        nqvk = cmplx(0.,0.)
        ttok = cmplx(0.,0.)
        qvok = cmplx(0.,0.)
     endif

    if (linear.ne.1 ) then !got nonlinear terms, if no turbulence then advection=0

!---------------------------------------------------------
      call constr (zxok,zyok,zzok,ttok,qvok,nzxk,nzyk,nzzk,nttk,nqvk,l,uk,vk,wk,uk,vk,wk,zxok,zyok,zzok,ttok,qvok,  &
                   nzxk,nzyk,nzzk,nttk,nqvk,zi,kxa,kya,kza,ak,bk,ck,ak,bk,ck,urt,vrt,wrt,1)
!---------------------------------------------------------
    else
        nzxk = cmplx(0.,0.)
        nzyk = cmplx(0.,0.)
        nzzk = cmplx(0.,0.)
        nttk = cmplx(0.,0.)
        nqvk = cmplx(0.,0.)

    endif

  endif!iturb&thermo

if (iturb .eq. 1) then

  do 20 ikz=1,iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     do 20 iky=1,ikty
        ky = kya(iky)
        do 20 ikx=1,iktx
           kx = kxa(ikx)
!           if ( l(ikx,iky,ikz).ne.1 ) then
!              zxnk(ikx,iky,ikz) = cmplx(0.,0.)
!              zynk(ikx,iky,ikz) = cmplx(0.,0.)
!              zznk(ikx,iky,ikz) = cmplx(0.,0.)
!           else
              wk2= kx*kx+ky*ky+kz*kz
              k  = sqrt( wk2 )
              r1 = visc*k**ilap2 !viscosity derivative term (without vorticity)

              termzx = nzxk(ikx,iky,ikz)-zxok(ikx,iky,ikz)*r1
              termzy = nzyk(ikx,iky,ikz)-zyok(ikx,iky,ikz)*r1
!              termzx = nzxk(ikx,iky,ikz)+expsiv*zi*ky*ttok(ikx,iky,ikz)-zxok(ikx,iky,ikz)*r1
!              termzy = nzyk(ikx,iky,ikz)-expsiv*zi*kx*ttok(ikx,iky,ikz)-zyok(ikx,iky,ikz)*r1
              termzz = nzzk(ikx,iky,ikz)-zzok(ikx,iky,ikz)*r1

              rhzx(ikx,iky,ikz) = termzx
              zxnk(ikx,iky,ikz) = zxok(ikx,iky,ikz) + delt*termzx
              zxnk(ikx,iky,ikz) = zxnk(ikx,iky,ikz) * l(ikx,iky,ikz)

              rhzy(ikx,iky,ikz) = termzy
              zynk(ikx,iky,ikz) = zyok(ikx,iky,ikz) + delt*termzy
              zynk(ikx,iky,ikz) = zynk(ikx,iky,ikz) * l(ikx,iky,ikz)

              rhzz(ikx,iky,ikz) = termzz
              zznk(ikx,iky,ikz) = zzok(ikx,iky,ikz) + delt*termzz
              zznk(ikx,iky,ikz) = zznk(ikx,iky,ikz) * l(ikx,iky,ikz)

!           endif
20 continue
endif ! iturb

  if (thermo .eq. 1) then
     do 21 ikz=1,iktzp
        ikza=mype*iktzp+ikz
        kz = kza(ikza)
        do 21 iky=1,ikty  
           ky = kya(iky)
           do 21 ikx=1,iktx
!              if ( l(ikx,iky,ikz).ne.1 ) then
!                 ttnk(ikx,iky,ikz) = cmplx(0.,0.)
!                 qvnk(ikx,iky,ikz) = cmplx(0.,0.)
!              else
                 kx = kxa(ikx)
                 wk2= kx*kx+ky*ky+kz*kz
                 k  = sqrt( wk2 )

                 r1 = kappand*k**ilap2 !kappand=1.38*visc
                 termtt = nttk(ikx,iky,ikz)-ttok(ikx,iky,ikz)*r1+ gammand*wk(ikx,iky,ikz)
                 rhtt(ikx,iky,ikz) = termtt
                 ttnk(ikx,iky,ikz) = ttok(ikx,iky,ikz) + delt*termtt
                 ttnk(ikx,iky,ikz) = ttnk(ikx,iky,ikz) * l(ikx,iky,ikz)!mark get rid of if

                 r2 = diffvnd*k**ilap2
                 termqv = nqvk(ikx,iky,ikz) - qvok(ikx,iky,ikz)*r2
                 rhqv(ikx,iky,ikz) = termqv
                 qvnk(ikx,iky,ikz) = qvok(ikx,iky,ikz) + delt*termqv
                 qvnk(ikx,iky,ikz) = qvnk(ikx,iky,ikz) * l(ikx,iky,ikz)!mark
!              endif
21   continue
  endif!thermo

  nt = 1 ! timestep in current restart block
  time = ts + nt*delt

!---------------------------------------------------------
   if (iturb .eq. 0) then
      nzxk = cmplx(0.,0.)
      nzyk = cmplx(0.,0.)
      nzzk = cmplx(0.,0.)
   endif
   if (thermo .eq. 0) then
      nttk = cmplx(0.,0.)
      nqvk = cmplx(0.,0.)
   endif

  if (linear.ne.1 ) then


  call constr (zxnk,zynk,zznk,ttnk,qvnk,nzxk,nzyk,nzzk,nttk,nqvk,l,uk,vk,wk,uk,vk,wk,zxnk,zynk,zznk,ttnk,qvnk,   &
                nzxk,nzyk,nzzk,nttk,nqvk,zi,kxa,kya,kza,ak,bk,ck,ak,bk,ck,urt,vrt,wrt,0)

!---------------------------------------------------------
  else
        nzxk = cmplx(0.,0.)
        nzyk = cmplx(0.,0.)
        nzzk = cmplx(0.,0.)
        nttk = cmplx(0.,0.)
        nqvk = cmplx(0.,0.)
  endif !linear

 if (iturb .eq. 1) then

  do 30 ikz=1,iktzp
     ikza=mype*iktzp+ikz
     kz = kza(ikza)
     do 30 iky=1,ikty
        ky = kya(iky)
        do 30 ikx=1,iktx
!           if ( l(ikx,iky,ikz).ne.1 ) then
!              zxnk(ikx,iky,ikz) = cmplx(0.,0.)
!              zynk(ikx,iky,ikz) = cmplx(0.,0.)
!              zznk(ikx,iky,ikz) = cmplx(0.,0.)
!           else
              kx = kxa(ikx)
              wk2= kx*kx+ky*ky+kz*kz
              k  = sqrt( wk2 )

              r1 = visc*k**ilap2

              termzx = nzxk(ikx,iky,ikz) - zxnk(ikx,iky,ikz)*r1
              termzy = nzyk(ikx,iky,ikz) - zynk(ikx,iky,ikz)*r1
!              termzx = nzxk(ikx,iky,ikz) + expsiv*zi*ky*ttnk(ikx,iky,ikz) - zxnk(ikx,iky,ikz)*r1
!              termzy = nzyk(ikx,iky,ikz) - expsiv*zi*kx*ttnk(ikx,iky,ikz) - zynk(ikx,iky,ikz)*r1
              termzz = nzzk(ikx,iky,ikz)-zznk(ikx,iky,ikz)*r1
 
              zxnk(ikx,iky,ikz) = zxok(ikx,iky,ikz) +  delt*(termzx+rhzx(ikx,iky,ikz))/2.
              zynk(ikx,iky,ikz) = zyok(ikx,iky,ikz) +  delt*(termzy+rhzy(ikx,iky,ikz))/2.
              zznk(ikx,iky,ikz) = zzok(ikx,iky,ikz) +  delt*(termzz+rhzz(ikx,iky,ikz))/2.
              zxnk(ikx,iky,ikz) = zxnk(ikx,iky,ikz) *  l(ikx,iky,ikz)
              zynk(ikx,iky,ikz) = zynk(ikx,iky,ikz) *  l(ikx,iky,ikz)
              zznk(ikx,iky,ikz) = zznk(ikx,iky,ikz) *  l(ikx,iky,ikz)
           !endif
30 continue

! prepare table for analytical diffusion
  do 35 ikz=1,iktzp
     ikza=mype*iktzp+ikz
     kz = kza(ikza)
     do 35 iky=1,ikty
        ky = kya(iky)
        do 35 ikx=1,iktx
           kx = kxa(ikx)
           wk2= kx*kx+ky*ky+kz*kz
           k  =  sqrt( wk2 )
           mdiff(ikx,iky,ikz) = exp(-visc*delt*k**ilap2) * l(ikx,iky,ikz)
35 continue
 else !iturb=0
    zxnk = cmplx(0.,0.)
    zynk = cmplx(0.,0.)
    zznk = cmplx(0.,0.)
if(thermo .eq.0 )       l = 0
    urt = 0.
    vrt = 0.
    wrt = 0.
 endif !iturb
!vectorize mark
  if (thermo .eq. 1) then
    do 31 ikz=1,iktzp
       ikza=mype*iktzp+ikz
       kz = kza(ikza)
       do 31 iky=1,ikty
          ky = kya(iky)
          do 31 ikx=1,iktx
!             if ( l(ikx,iky,ikz).ne.1 ) then!l=0
!                ttnk(ikx,iky,ikz) = cmplx(0.,0.)
!                qvnk(ikx,iky,ikz) = cmplx(0.,0.)
!             else
                kx = kxa(ikx)
                wk2= kx*kx+ky*ky+kz*kz
                k  = sqrt( wk2 )

!                r1 = visc*k**ilap2

                r1 = kappand*k**ilap2

                termtt = nttk(ikx,iky,ikz) - ttnk(ikx,iky,ikz)*r1 + gammand*wk(ikx,iky,ikz)
                ttnk(ikx,iky,ikz) = ttok(ikx,iky,ikz) + delt*(termtt+rhtt(ikx,iky,ikz))/2.
	        ttnk(ikx,iky,ikz) = ttnk(ikx,iky,ikz) * L(ikx,iky,ikz)!mark remove if

                r2 = diffvnd*k**ilap2

                termqv = nqvk(ikx,iky,ikz) - qvnk(ikx,iky,ikz)*r2
                qvnk(ikx,iky,ikz) = qvok(ikx,iky,ikz) + delt*0.5*(termqv+rhqv(ikx,iky,ikz))
		qvnk(ikx,iky,ikz) = qvnk(ikx,iky,ikz) * L(ikx,iky,ikz) !mark remove if
!             endif
31 continue

! prepare table for analytical diffusion
     do 36 ikz=1,iktzp
        ikza=mype*iktzp+ikz
        kz = kza(ikza)
        do 36 iky=1,ikty
           ky = kya(iky)
           do 36 ikx=1,iktx
              kx = kxa(ikx)
              wk2= kx*kx+ky*ky+kz*kz
              k  =  sqrt( wk2 )
              tdiff(ikx,iky,ikz) = exp(-kappand*delt*k**ilap2) * l(ikx,iky,ikz)
              qdiff(ikx,iky,ikz) = exp(-diffvnd*delt*k**ilap2) * l(ikx,iky,ikz)
36 continue
  else
    qvnk = cmplx(0.,0.)
    ttnk = cmplx(0.,0.)
  endif!thermo




  if(gomic.eq.1 ) then

! initialize droplet velocities
! arrays of indices of droplet positions for interpolation purposes

!     do id = 1,ndropp
        ixp(1:ndropp)=int(x(1:ndropp)*oneoverh) + 1
        iyp(1:ndropp)=int(y(1:ndropp)*oneoverh) - mype*n2pe + 1
        izp(1:ndropp)=int(z(1:ndropp)*oneoverh) + 1
!     enddo

   if (iturb .eq. 1) then
     do 200 id=1,ndropp
        dxp = oneoverh*x(id) - real(ixp(id)) + 1
        dyp = oneoverh*y(id) - real(iyp(id) + mype*n2pe) + 1
        dzp = oneoverh*z(id) - real(izp(id)) + 1
        udrop(id) = & 
         (1.-dxp)*(1.-dyp)*(1.-dzp)*urt(ixp(id),izp(id),iyp(id)) &
         + dxp*(1.-dyp)     *(1.-dzp)*urt(ixp(id)+1,izp(id),iyp(id)) &
         + dxp*dyp          *(1.-dzp)*urt(ixp(id)+1,izp(id),iyp(id)+1) &
         + (1.-dxp)*dyp     *(1.-dzp)*urt(ixp(id),izp(id),iyp(id)+1) &
         + (1.-dxp)*(1.-dyp)*dzp     *urt(ixp(id),izp(id)+1,iyp(id)) &
         + dxp*(1.-dyp)     *dzp     *urt(ixp(id)+1,izp(id)+1,iyp(id)) &
         + dxp*dyp          *dzp     *urt(ixp(id)+1,izp(id)+1,iyp(id)+1) &
         + (1.-dxp)*dyp     *dzp     *urt(ixp(id),izp(id)+1,iyp(id)+1) 

        vdrop(id) = &
         (1.-dxp)*(1.-dyp)*(1.-dzp)*vrt(ixp(id),izp(id),iyp(id)) &
         + dxp*(1.-dyp)     *(1.-dzp)*vrt(ixp(id)+1,izp(id),iyp(id)) &
         + dxp*dyp          *(1.-dzp)*vrt(ixp(id)+1,izp(id),iyp(id)+1) &
         + (1.-dxp)*dyp     *(1.-dzp)*vrt(ixp(id),izp(id),iyp(id)+1) &
         + (1.-dxp)*(1.-dyp)*dzp     *vrt(ixp(id),izp(id)+1,iyp(id)) &
         + dxp*(1.-dyp)     *dzp     *vrt(ixp(id)+1,izp(id)+1,iyp(id)) &
         + dxp*dyp          *dzp     *vrt(ixp(id)+1,izp(id)+1,iyp(id)+1) &
         + (1.-dxp)*dyp     *dzp     *vrt(ixp(id),izp(id)+1,iyp(id)+1) 

        wdrop(id) = & 
         (1.-dxp)*(1.-dyp)*(1.-dzp)*wrt(ixp(id),izp(id),iyp(id)) &
         + dxp*(1.-dyp)     *(1.-dzp)*wrt(ixp(id)+1,izp(id),iyp(id)) &
         + dxp*dyp          *(1.-dzp)*wrt(ixp(id)+1,izp(id),iyp(id)+1) &
         + (1.-dxp)*dyp     *(1.-dzp)*wrt(ixp(id),izp(id),iyp(id)+1) &
         + (1.-dxp)*(1.-dyp)*dzp     *wrt(ixp(id),izp(id)+1,iyp(id)) &
         + dxp*(1.-dyp)     *dzp     *wrt(ixp(id)+1,izp(id)+1,iyp(id)) &
         + dxp*dyp          *dzp     *wrt(ixp(id)+1,izp(id)+1,iyp(id)+1) &
         + (1.-dxp)*dyp     *dzp     *wrt(ixp(id),izp(id)+1,iyp(id)+1) 
        if (r(id) .lt. 4.d-5) then
           wdrop(id) = wdrop(id) - kwt*r(id)**2
        else
           wdrop(id) = wdrop(id) - (kwt2*r(id)-0.1227)
        endif!terminal velo

 200  continue
   else!iturb=0
      udrop = 0.
      vdrop = 0.
      do id = 1,ndropp
         if (r(id) .lt. 4.d-5) then
            wdrop(id) = -kwt*r(id)**2
         else
            wdrop(id) = -kwt2*r(id)+0.1227
         endif!terminal velo
      enddo
   endif!iturb
  endif


! ---------------------------------------------------------------------------------------
! subsequent timesteps
! ---------------------------------------------------------------------------------------

  time3=etime(time1)

  do 50 nt=2,nstop+1
! nt goes to nstop+1 to complete microphysics
     time = ts + nt*delt
     ! max function is necessary for case when irest = 0

     if (mype .eq. 0) print*, 'nt =', nt

     if (linear.ne.1) then
!---------------------------------------------------------
! convol returns the nonlinear term, nzxk,nzyk,nzzk in k-space. and microphysics
     call convol (zxnk,zynk,zznk,ttnk,qvnk,nzxk,nzyk,nzzk,nttk,nqvk,l,uk,vk,wk,uk,vk,wk,rhzx,rhzy,rhzz,rhtt,rhqv,   &
                   zxnk,zynk,zznk,ttnk,qvnk,urt,vrt,wrt,nzxk,nzyk,nzzk,nttk,nqvk,zi,kxa,kya,kza,ak,bk,ck,ak,bk,ck,  &
                   time0,nt,ndumpd,ndropreal,nddone)
!---------------------------------------------------------
     else
        nzxk = cmplx(0.,0.)   
        nzyk = cmplx(0.,0.)   
        nzzk = cmplx(0.,0.)   
        nttk = cmplx(0.,0.)   
        nqvk = cmplx(0.,0.)   
     endif

     if(nt.eq.nstop+1) goto 50      !!!! end after micro. adjustment

     if (iturb .eq. 1) then
     do 40 ikz=1,iktzp
        ikza=mype*iktzp+ikz
        kz = kza(ikza)
        do 40 iky=1,ikty
           ky = kya(iky)
           do 40 ikx=1,iktx
              kx = kxa(ikx)

              termzx = nzxk(ikx,iky,ikz)
              termzy = nzyk(ikx,iky,ikz)
!              termzx = nzxk(ikx,iky,ikz)+expsiv*zi*ky*ttnk(ikx,iky,ikz)
!              termzy = nzyk(ikx,iky,ikz)-expsiv*zi*kx*ttnk(ikx,iky,ikz)
              termzz = nzzk(ikx,iky,ikz)
!
              tzx =  zxok(ikx,iky,ikz)*mdiff(ikx,iky,ikz)**2 + d2*mdiff(ikx,iky,ikz)*termzx
              tzy =  zyok(ikx,iky,ikz)*mdiff(ikx,iky,ikz)**2 + d2*mdiff(ikx,iky,ikz)*termzy
              tzz =  zzok(ikx,iky,ikz)*mdiff(ikx,iky,ikz)**2 + d2*mdiff(ikx,iky,ikz)*termzz
!             the robert filter. (see asselin,r.a.,1972,mon.wea.rev.,100,p.487-490)

              avzx = zxnk(ikx,iky,ikz) + robert* (tzx - 2.*zxnk(ikx,iky,ikz) + zxok(ikx,iky,ikz))
              avzy = zynk(ikx,iky,ikz) + robert* (tzy - 2.*zynk(ikx,iky,ikz) + zyok(ikx,iky,ikz))
              avzz = zznk(ikx,iky,ikz) + robert* (tzz - 2.*zznk(ikx,iky,ikz) + zzok(ikx,iky,ikz))

              zxok(ikx,iky,ikz) = avzx*l(ikx,iky,ikz)
              zxnk(ikx,iky,ikz) =  tzx*l(ikx,iky,ikz)
              zyok(ikx,iky,ikz) = avzy*l(ikx,iky,ikz)
              zynk(ikx,iky,ikz) =  tzy*l(ikx,iky,ikz)
              zzok(ikx,iky,ikz) = avzz*l(ikx,iky,ikz)
              zznk(ikx,iky,ikz) =  tzz*l(ikx,iky,ikz)
40 continue
     endif !iturb

     if (thermo .eq. 1) then
        do 41 ikz=1,iktzp
           ikza=mype*iktzp+ikz
           kz = kza(ikza)
           do 41 iky=1,ikty
              ky = kya(iky)
              do 41 ikx=1,iktx
                 kx = kxa(ikx)

                 termtt = nttk(ikx,iky,ikz) + gammand*wk(ikx,iky,ikz)
                 termqv = nqvk(ikx,iky,ikz)

                 ttt =  ttok(ikx,iky,ikz)*tdiff(ikx,iky,ikz)**2 + d2*tdiff(ikx,iky,ikz)*termtt
                 tqv =  qvok(ikx,iky,ikz)*qdiff(ikx,iky,ikz)**2 + d2*qdiff(ikx,iky,ikz)*termqv


!                 the robert filter. (see asselin,r.a.,1972,mon.wea.rev.,100,p.487-490)

                 avtt = ttnk(ikx,iky,ikz) + robert* (ttt - 2.*ttnk(ikx,iky,ikz) + ttok(ikx,iky,ikz))
                 avqv = qvnk(ikx,iky,ikz) + robert* (tqv - 2.*qvnk(ikx,iky,ikz) + qvok(ikx,iky,ikz))

                 ttok(ikx,iky,ikz) = avtt*l(ikx,iky,ikz)
                 ttnk(ikx,iky,ikz) =  ttt*l(ikx,iky,ikz)
                 qvok(ikx,iky,ikz) = avqv*l(ikx,iky,ikz)
                 qvnk(ikx,iky,ikz) =  tqv*l(ikx,iky,ikz)
41      continue
     endif !thermo

     if (ampfor.eq.1) then
        call force(delt,zxnk,zynk,zznk,ikxf,ikyf,ikzf,ikzfp,nf,l,kxa,kya,kza,nt,nout,dimk)
     endif !ampfor
! write to restart file
     if (mod(nt,ndump) .eq. 0) then
        ntdump=nt/ndump
        if (netcdf .eq. 1) then
           if(iturb .eq. 1) call proj(zxnk,zynk,zznk,l,kxa,kya,kza)
           call ncdumprst(zxnk,zynk,zznk,ttnk,qvnk,uk,ntdump,time)
        else
           print*, 'sorry, do not know how to dump restart files without netcdf'
        endif 
     endif ! ndump
 
! write real space fields
     if (rsflag .eq. 1 .and. mod(nt,nbig) .eq. 0 ) then
        ntdumpd=nt/nbig+1
        if (netcdf .eq. 1) then
           call dumpreal(zxnk,zynk,zznk,ttnk,qvnk,zxnk,zynk,zznk,ttnk,qvnk,uk,vk,wk, &
                          uk,vk,wk,nzxk,l,kxa,kya,kza,ntdumpd,time,nrsp)
        else
           print*, 'sorry, do not know how to write real space fields without netcdf'
        endif 
     endif   !rsflag


     if ( mod(nt,nproj+1) .eq. 0 ) then 
        if (iturb .eq. 1) call proj (zxok,zyok,zzok,l,kxa,kya,kza)
     endif !nt

     if ( mod(nt-1,nout) .ne. 0 ) go to 50
     call spec(zxok,zyok,zzok,ttok,qvok,ns,spz,l,time,kxa,kya,kza,vvk,kkol,dimk)

     call out (zxok,zyok,zzok,ttok,qvok,uk,vk,wk,time,nstop,ns,spz,l,kxa,kya,kza,vvk,kkol,dimk)
     if (mype .eq. 0 .and. iturb .eq. 1) then
         print*,'  '
         print*,'kmax*eta = ',dimk*n/2.*(visc**2/(2.*e5))**.25!dimk*n/2.*(visc**2/(2.*e5(io)))**.25
!         write(97,*) 'kmax*eta = ',dimk*n/2.*(visc**2/(2.*e5))**.25!dimk*n/2.*(visc**2/(2.*e5(io)))**.25
     endif

     if (iturb .eq. 1) call rspace (zxnk,zynk,zznk,ttnk,zxnk,zynk,zznk,ttnk,nzxk,nzyk,nzzk,nttk,nzxk,nzyk,nzzk,nttk, &
                        uk,vk,wk,uk,vk,wk,l,zi,kxa,kya,kza,time)

     if (mype .eq.  0 .and. iturb .eq. 1) then
        print*,'  '
        print*,'vorticity field divergence:'
     endif
     if(iturb .eq. 1) call diverg (zxnk,zynk,zznk,l,kxa,kya,kza)
50 continue


  time = time - delt

! compute diagnostics 
  call spec(zxok,zyok,zzok,ttok,qvok,ns,spz,l,time,kxa,kya,kza,vvk,kkol,dimk)
  call out (zxok,zyok,zzok,ttok,qvok,uk,vk,wk,time,nstop,ns,spz,l,kxa,kya,kza,vvk,kkol,dimk)
  if (iturb .eq. 1)    call rspace (zxnk,zynk,zznk,ttnk,zxnk,zynk,zznk,ttnk,nzxk,nzyk,nzzk,nttk,nzxk,nzyk,nzzk,nttk, &
                uk,vk,wk,uk,vk,wk,l,zi,kxa,kya,kza,time)
  print*,'      '

  if (mype .eq. 0) then 
     close(14)
     close(41)
     close(42)
     close(43)
     close(44)
     close(45)
     close(57)
     close(61)
     close(80)
     close(81)
     close(82)
     close(93)
     close(94)
     close(95)
     close(97)

     print*,'              '
     print*,'              '
     print*,'              '
     print*,'              '
     print*,'______________________________________________________ '
  endif ! mype

333  format(1x,a10)

   time3 = etime(time2)
   time3 = time2(1) - time1(1)


  if (mype .eq. 0) then
     print*,'     '
     write(6,5000) time3, time3/60., time3/3600., time3/86400.
  endif

!if (mpi == 1) then 
  call mpi_finalize(ierror)
!endif

5000 format(1x,'cpu time required for main loop = ',f7.0,' s = ',f7.1,' m = ',f7.2,' h = ',f7.3,' d.')

  end 

