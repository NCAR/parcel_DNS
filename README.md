If you use the model/model code
Please cite the following paper. Thank you!
Chen, S., Xue, L., and Yau, M.-K.: Impact of aerosols and turbulence on cloud droplet growth: an in-cloud seeding case study using a parcel–DNS (direct numerical simulation) approach, Atmos. Chem. Phys., 20, 10111–10124, https://doi.org/10.5194/acp-20-10111-2020, 2020.

Contact Sisi Chen (sisichen@ucar.edu) if you have any questions. Thank you!

# aerosol activation code
Both models were written in FORTRAN.

#parcel model
The model was developed based on Jensen & Nugent 2017 JAS's equations
Collision-coalescence is excluded in the parcel model
The code is stored in /parcel 
The folder include:
1) run is the script to compile & run the code. The code can be run at your desktop (if you have FORTRAN compiler installed)
2) activation.F90 is the source code of the parcel model. 
3) a sample of the output files using the current configuration 
	3a) new.dsd is the output of the droplet size distribution (number concentration of the droplet at each size bin)
	3b) new.rad is the radius of the droplet at each size bin
  	3c) new.out is the parcel-mean variables of time, height, supersaturation, total number of droplets, pressure, temperature, potential temperature, water vapor mixing ratio, saturated water vapor mixing ratio, liquid water content, air density, deltaqp (some intermediate quantities that you may skip). 
	Look for the line starting with "write(16,*) time,0,Sp,ndrop,pp,temp,”... in the code, and you can find the variable names corresponding to the above quantities. 


#DNS model

This is the version of DNS based on Chen et al. (2020) in ACP: https://acp.copernicus.org/articles/20/10111/2020/acp-20-10111-2020.html


The model equations & algorithms can be seen in Chen et al. (2016) http://journals.ametsoc.org/doi/abs/10.1175/JAS-D-15-0203.1 
and Chen et al. (2018) http://journals.ametsoc.org/doi/abs/10.1175/JAS-D-17-0123.1

Modifications were also made in Chen et al. (2018) ACP: https://www.atmos-chem-phys.net/18/7251/2018/acp-18-7251-2018.html


Updates:
May 2020:
Changed the precision of droplet-related variables to double precision. The single precision has caused the round-off truncation error which prevented the droplet from growing.

May 2019:
Debug finished in terms of solute effect. The model now can handle the hygroscopic growth of multiple chemical species.

Feb. 2019: 
Physics: Added aerosol activation code to DNS. Particle can grow from its dry size. 
Numerics: modified the condensational growth equation (use hygroscopicity parameter from Petters&Kreidenweis 2007) in microhydroall.F90. Modified the idrops.F90 to include the calculation of the initial equilibrium droplet size. (iteration)

Early versions:

The first version of this DNS was developed in late 90s using Fortran 77 by Paul Vaillancourt (see his PhD thesis and publications 2001&2002).The code was in serial at that time.


The second version given by Franklin et al. (2005&2007) was parallelized using OpenMP. The MPI parallelization on the turbulence code was finally realized by Mike Waite. And the code was then partially transformed to F90


##--------modules loaded at Cheyenne---------##

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
    
Edit run_cheyenne (job submission script to Cheyenne, you may modify it based on the job schedulers in your supercomputer)

compile and submit
    ./compile_and_run
    
File output:
    Zk1.out.ncf     #flow field data in Fourier space (netcdf)
    out.ncf         #flow field data in real space if rsflag = 1 in main.F90 (netcdf)
    Run_aerosol.*   #various informations with its file explained/defined in main.F90
    
2 droplet(microphysics) spin-up run: 

Rename restart file:
    mv Zk1.out.ncf Zk.in.ncf #rename the output to *in.ncf and used as restart file
param.inc:
    gomic = 1 #initiate droplets

Edit run_cheyenne (job submission script to Cheyenne, you may modify it based on the job schedulers in your supercomputer)
compile and submit
    ./compile_and_run
    
File output:
    Zk1.out.ncf     #flow field data in Fourier space (netcdf)
    drop1.out.ncf   #droplet data (netcdf)
    out.ncf         #flow field data in real space if rsflag = 1 in main.F90 (netcdf)
    Run_aerosol.*   #various informations with its file explained/defined in main.F90
    
3 real simulation

Rename restart file:
    mv Zk1.out.ncf Zk.in.ncf
    mv drop1.out.ncf drop.in.ncf
param.inc:
    gomic = 2 #droplets from restart files

Edit run_cheyenne (job submission script to Cheyenne, you may modify it based on the job schedulers in your supercomputer)

compile and submit
    ./compile_and_run




