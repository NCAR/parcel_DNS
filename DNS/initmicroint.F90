 
  SUBROUTINE INITMICRO(time0,ndropreal,nddone)
! This subroutine initializes variables for parcel calculations
! and fields for detailed model
!!! judge and tell us whether it has restarted block if yes print info out
  use thermo_mod
  use mic_mod
  implicit none
  include 'param.inc'
  include 'mpif.h'

  integer :: seed,id
  integer :: ndropreal,nddone
  real(8) :: rm
  real(8) :: cql,tmp
  character*128 :: form155
  real    :: time0

  integer :: mype,ierror
  common/mpi/ mype
      
  time0=0.0d0
  oneoverh = 1.d0/h
 
  if(gomic .eq. 1 ) then
! ----------------------------------------------------------
! Initialize size and position of droplets
     seed  = seedmic*(mype+1) 			! Seeding of random number generator.
     call idrops(seed,rm)
     udrop = 0.0  
     vdrop = 0.0  
     wdrop = 0.0  
     ndropreal = ndropp
     nddone = 0
     rmth = 0.d0
     rmqv = 0.d0
!     deallocate(integsup,integsup2)
     if (thermo .eq. 1) then
! Initialize parcel variables
       TEMP = T0
       pp = p1
       exner = (pp/p0)**RaCp
       thetapp= TEMP/exner
   !!  assume initial supersaturation is the quasi-steady S. 
   !!  See Korolev & Mazin (2003) appendix C or Rogers&Yau(1989)
    !!  assume qvpp=qvs in first estimate of source and tau
       !ks = 1.0d0/(rhow*Rv*temp/(esat*diffvnd)+rhow*Lat/(Ka*Temp)*(Lat/(Rv*Temp)-1))
       esat = 2.53d11*exp(-5.42d3/temp) !saturation vapor pressure
       qvs = eps*esat/(PP-esat)
       sp = sp0 !initial Supersaturation set to 0
       qvpp=(sp+1.0d0)*qvs
       rhoa = pp/(Ra*(1+18.d0/29.d0*qvpp)*temp)
       cql = 4.0d0*pi*rhow/(rhoa*vol)
     endif !thermo

  ELSE ! GOMIC >1; Begin microphysics restart mode 
     if (mype .eq. 0) then 
        print*,'RESTART WITH MYCROPHYSICAL FILE. '
     endif
     if (netcdf == 1) then
        call ncreaddrop(ndropreal,time0,nddone)
        if(mype .eq. 0) print*,'sp initial=', sp !
     else
        print*, 'Sorry, do not know how to restart without NetCDF'
        stop
     endif
 if(mype.eq.0) print*,'qvpp,thetapp,rmth,rmqv',qvpp,thetapp,rmth,rmqv !mark
     exner = (PP/P0)**RaCp
     temp = thetapp*exner
     esat = 2.53d11*exp(-5.42d3/temp) !saturation vapor pressure
     cql = 4.0d0*pi*rhow/(rhoa*vol)
     !ks = 1.0d0/(rhow*Rv*temp/(esat*diffvnd)+rhow*Lat/(Ka*Temp)*(Lat/(Rv*Temp)-1))
     if (mype .eq. 0)   print*,'Just read time = ',time0, 'for droplet restart'

  ENDIF ! GOMIC .EQ. 1
  !calculate initial LWC
  lwc = 0.0
  do id = 1,ndropreal
     lwc = lwc + r(id)**3
  enddo
  lwc=lwc*4.d0/3.d0*pi*rhow
  call mpi_reduce(lwc,tmp,1,mpi_real8,mpi_sum,0,mpi_comm_world,ierror)
  if (mype .eq. 0) then
      lwc=tmp/(n1*n2*n3*vol*rhoa)
      write(14,*) 'Initial LWC = ', LWC
  endif
     

  return
  end



