# aerosol activation code

#parcel model
The model was developed based on Jensen & Nugent 2017 JAS's equations
Collision-coalescence is excluded in the parcel model

#DNS model

This is the version of DNS based on Chen et al. (2020) in ACP: https://www.atmos-chem-phys-discuss.net/acp-2019-886/


The model equations & algorithms can be seen in Chen et al. (2016) http://journals.ametsoc.org/doi/abs/10.1175/JAS-D-15-0203.1 
and Chen et al. (2018) http://journals.ametsoc.org/doi/abs/10.1175/JAS-D-17-0123.1

Modifications were also made in Chen et al. (2018) ACP: https://www.atmos-chem-phys.net/18/7251/2018/acp-18-7251-2018.html


Updates:
May 2020:
Changed the precision of droplet-related variables to double precision. The single precision has caused the round-off truncation error which prevent the droplet from growing.

May 2019:
Debug finished in terms of solute effect. The model now can handle the hygroscopic growth of multiple chemical species.

Feb. 2019: 
Physics: Added aerosol activation code to DNS. Particle can grow from its dry size. 
Numerics: modified the condensational growth equation (use hygroscopicity parameter from Petters&Kreidenweis 2007) in microhydroall.F90. Modified the idrops.F90 to include the calculation of the initial equilibrium droplet size. (iteration)

Early versions:

The first version of this DNS was developed in late 90s using Fortran 77 by Paul Vaillancourt (see his PhD thesis and publications 2001&2002).The code was in serial at that time.


The second version given by Franklin et al. (2005&2007) was parallelized using OpenMP. The MPI parallelization on the turbulence code was finally realized by Mike Waite. And the code was then partially transformed to F90


The model is in git_aerosol/

The test data is in example/

##--------modules loaded---------##

1) ncarenv/1.2   2) intel/17.0.1   3) ncarcompilers/0.4.1   4) mpt/2.19   5) netcdf/4.6.1
#library installed locally
FFTW2.1.5  (location: /glade/work/sisichen/fftw2.1.5intel_openmpi )

##-----how to run the code at Cheyenne----##

1 compile the code by "make" 
2 submit the job by "qsub run_cheyenne"
or simply "./compile_and_run"

##------model config & simulations (need to run 2 spin-ups before the real simulation)-----##

1 turbulence spin-up run: need to get fully-developed turbulence before insert droplets

param.inc: 
    iturb = 1 #turn on turbulence
    gomic = 0 #turn off droplets
    thermo = 1 #turn on thermodynamic calculation (e.g., Temperature & water vapor)
compile and submit
    ./compile_and_run
File output:
    Zk1.out.ncf     #flow field data in Fourier space
    out.ncf         #flow field data in real space if rsflag = 1 in main.F90
    Run_aerosol.*   #various informations with its file explained in main.F90

2 droplet(microphysics) spin-up run: 

Restart file:
    mv Zk1.out.ncf Zk.in.ncf #rename the output to *in.ncf and used as restart file
param.inc:
    gomic = 1 #initiate droplets
compile and submit
    ./compile_and_run
File output:
    Zk1.out.ncf     #flow field data in Fourier space
    out.ncf         #flow field data in real space if rsflag = 1 in main.F90
    Run_aerosol.*   #various informations with its file explained in main.F90
    drop1.out.ncf   #droplet data

3 real simulation

Restart file:
    mv Zk1.out.ncf Zk.in.ncf
    mv drop1.out.ncf drop.in.ncf
param.inc:
    gomic = 2 #droplets from restart files

compile and submit
    ./compile_and_run


##----test run----------------------------------##

The results are in the folder "example"
cd example
spin1: output for spin1
spin2: output for spin2
simulation: output for real simulation

