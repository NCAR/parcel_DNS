! -----------------------------------------------------------------
! Set model resolution

      integer, parameter :: N = 128
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! Set number of processors (must be even)

      integer, parameter :: MPI = 1
      integer, parameter :: NPE = 64
! -----------------------------------------------------------------

      integer, parameter :: N1 = N, N2 = N, N3 = N 			! Resolution in the three dimensions
      integer, parameter :: N1D = N1+2, N2D = N2, N3D = N3   
      integer, parameter :: N2DP = N2D/NPE, N2PE = N2/NPE	! Real-space points per processor in 2nd direction  
      integer, parameter :: KTX = N1/2, KTY = N2/2, KTZ = N3/2		! Truncation wavenumbers
      integer, parameter :: IKTX = KTX + 1, IKTY = N2D, IKTZ = N3D 	! Array indices of truncation wavenumbers
      integer, parameter :: IKTZP = IKTZ/NPE				! Same but per processor for the third direction	
      integer, parameter :: NFMAX = 64					! Maximum number of forced wavenumbers 

! -----------------------------------------------------------------

! -----------------------------------------------------------------
! Use NetCDF for generating restart files and real space fields

      integer, parameter :: netcdf = 1 

!------------------------------------------------------------------------------------------------
! Parameters for microphysics. All units are in metre, kilogram and/or second
!------------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------------
! General parameters

      integer, parameter :: disp = 6094
      integer, parameter :: ndrop = 390016		! Number of droplets ndrop=390016 for disp=6094nature
     							! ndrop = 780032 for disp=12188 double nature; 
							! ndrop=434944 for disp=6796 seeding&6800GCCN
							! ndrop = 479872 for disp=7498 double seeding
      integer, parameter :: ndropp = ndrop/npe		! Number of droplets per processor at initialization
      integer, parameter :: ndroppe = ndropp * 3	! Possible number of droplets in a slab during simulation (estimate)
      integer, parameter :: everynd = -1000		! Every 'everynd'-th droplet is followed
      integer, parameter :: seedturb = -200
      integer, parameter :: seedmic = -200
      integer, parameter :: nbins = 100
      real(8), parameter :: up = 2.d0			! Vertical velocity of air parcel
      real(8), parameter :: pi = 3.1415926535d0		! Pi
      real(8), parameter :: grav = 9.8d0			! Gravitational constant	
      real(8), parameter :: visc   =   1.6d-5				! kinematic viscosity of turbulent flow field [m**2/s]    
      real, parameter :: rtrunc = 50.0d-6		!truncated radius above which droplet stops grow and has no interaction with flow
      real, parameter :: rtrunc_min = 5.d-6		!truncated minimum radius below which droplets won't interact with other droplets and regarded as passive tracer.

!------------------------------------------------------------------------------------------------
! Parameters used for collision detection		

      integer, parameter :: m = N/2 			
      integer, parameter :: my = m/npe			!!  m should be multiple of npe (i.e. N is multiple of twice npe)
      integer, parameter :: ncell = m * m * my		!! number of cells
      integer, parameter :: mapsiz = 13 * ncell		!! used for linked list cell mapping

!------------------------------------------------------------------------------------------------
! Parameters that are function of temperature and/or pressure, values at 283.15K and 90kPa	

      real(8), parameter :: Cw = 4187.0                    ! Specific heat for water at 10degreeC [J/kgK]
      real(8), parameter :: Cp = 1004.0			! Specific heat for dry air at constant pressure [J/kgK]
      real(8), parameter :: Cv = 718.0			! Specific heat for dry air at constant volume [J/kgK]
      real(8), parameter :: diffvnd = 2.55d-5		! Coefficient of diffusion of water vapour in air [m**2/s]
      real(8), parameter :: kappand = 2.22d-5		! Coefficient of thermal diffusivity
      real(8), parameter :: Ra = 287.0			! Individual gas constant for dry air [J/kgK]
      real(8), parameter :: Rv = 467			! Individual gas constant for water vapour [J/kgK]
      real(8), parameter :: eps = Ra/Rv			! Epsilon = R'/Rv = Ra/Rv
      real(8), parameter :: Racp = Ra/Cp			! = R'/cp : for adiabatic processes: (T/T0) = (P/P0)**k, k = R'/cp
      real(8), parameter :: m_w = 18.d-3			! molecular weight for water 
      real(8), parameter :: m_s = 132.14d-3		! molecular weight of solute; ammonium sulfate=132.14d-3; NaCl = 58.44d-3
      real(8), parameter :: rho_ccn = 1726.d0 !kg/m**3 for ammonium sulfate =1726.d0 for NaCl=2160.d0
      real(8), parameter :: rhow = 1000.0			! Density of water [kg/m**3] 
      real(8), parameter :: sp0 = 1.59d-2			! initial supersaturation	
      real(8), parameter :: T0 = 281.2d0			! Reference temperature [K] ( = 10 degrees Celsius)
      real(8), parameter :: P0 = 1.0d5			! Reference pressure [Pa] ( = 100kPa = 1000mbar)
      real(8), parameter :: P1 = 90220.d0                 ! Initial pressure inside the cloud (=900hpa)
      real(8), parameter :: lat = 2.477d6 			! Latent heat of condensation of water [J/kg] L(T)=(2500.8-2.36T+.0016T^2-0.00006T^3)*10^3
      real(8), parameter :: LatoverRv= lat/Rv      	! constant B to calculate saturated water vapor
      real(8), parameter :: LatoverCp = lat/cp		! = L/cp (see Yau page 20-24. eqn (2.29))
      real(8), parameter :: ka = 2.48d-2			! Thermal conductivity of air [J/msK]
      real(8), parameter :: KK = 8.54d-11			! ?? (value at T = 280.65K)
      real, parameter :: Gammand = -grav/Cp*1.d0                    ! for cooling term. (=dry adiabatic) put to zero to exclude cooling term.
      real(8), parameter :: Nsc = visc/diffvnd		!Schmidt number
      real(8), parameter :: isedi   = 1.d0          	! put to zero for no sedimentation in microphysical model.
      real(8), parameter :: kwt     = 1.233d8	        !for r < 40micron
      real(8), parameter :: kwt2    = 8.d3           	! for r ~ [40micron, 0.6mm]
      real, parameter :: expsiv   =   grav/283.15d0* 0.d0  ! = g/T0. thermal expansivity. put to zero for no buoyancy.
      real, parameter :: temgrad  =   expsiv * 0.d0        ! vertical temperature gradient
!------------------------------------------------------------------------------------------------
! Additional type declarations

      integer, parameter :: thermo=1  ! thermo flag. if = 0 ignore temperature and qv equations. !!only allow 0 or 1 because of the constraints in outputting.F90
      integer, parameter :: gomic= 2   !microphysics flag (0 = no mic.phy., 1 = no mic. restart file, >=2 = read mic. restart-file)
      integer, parameter :: iturb = 1  ! 0 = turbulence off; 1 = on !!only allow 0 or 1 because of the constraints in outputting.F90

