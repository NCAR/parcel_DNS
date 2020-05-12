
subroutine dumpreal(zxk,zyk,zzk,ttk,qvk,zxr,zyr,zzr,ttr,qvr,uk,vk,wk,ur,vr,wr,junk, &
                    L,kxa,kya,kza,ntdump,time,nrsp)

  implicit none
  include 'param.inc'
  integer :: ikx,iky,ikz,ikza,i,j,k,L(iktx,ikty,iktzp),ntdump,nrsp
  real :: kxa(iktx),kya(ikty),kza(iktz),kx,ky,kz,time
  real,    dimension(n1d,n3d,n2dp)       :: ur,vr,wr,zxr,zyr,zzr,ttr,qvr,junk
  complex, dimension(iktx,ikty,iktzp) :: zxk,zyk,zzk,ttk,qvk,uk,vk,wk

  integer :: mype
  common/mpi/mype

  call velo(zxk,zyk,zzk,uk,vk,wk,L,kxa,kya,kza)

! Move vorticity and velocity to physical space.
  call fftwkr(zxr,zxk)
  call fftwkr(zyr,zyk)
  call fftwkr(zzr,zzk)
  if (thermo .eq. 1) then
     call fftwkr(ttr,ttk)
     call fftwkr(qvr,qvk)
  endif
  call fftwkr( ur, uk)
  call fftwkr( vr, vk)
  call fftwkr( wr, wk)

  call ncdumprsp(ur,vr,wr,zxr,zyr,zzr,ttr,qvr,junk,ntdump,time,nrsp)

! Put everything back to fourier for future timestepping
  call fftwrk(zxr,zxk)
  call fftwrk(zyr,zyk)
  call fftwrk(zzr,zzk)
  if (thermo .eq. 1) then
     call fftwrk(ttr,ttk)
     call fftwrk(qvr,qvk)
  endif
  uk = cmplx(0.,0.)
  vk = cmplx(0.,0.)
  wk = cmplx(0.,0.)

  return
end subroutine dumpreal


subroutine ncprep(nrsp)
  use netcdf_mod 
  implicit none
  include 'param.inc'
  include 'netcdf.inc'
  include 'mpif.h'
  integer :: istatus
!  integer :: istatus, status

!! real space output vars
  character*12 :: filename
  character*14 :: filenamed
  integer :: nrsp,nf
  integer :: idx,idy,idz,idt
  integer, dimension(4) :: ncdims

!! restart ourput vars
  integer :: idkx,idky,idkz,idkri,idkt
  integer, dimension(5) :: ncdimsk

!! droplet restart output vars
  integer :: idndroppe,iddt
  integer, dimension(2) :: ncdimsd

  integer :: mype
  common/mpi/mype


  !!! first, prep physical space output file out.ncf

  if (nrsp.gt.0) then

     if (mype.eq.0) then
 
        print*,'Creating netcdf file for output'
!        istatus = nf_create("out.ncf",0,idnc)
        istatus = nf_create("out.ncf",OR(NF_CLOBBER,NF_NETCDF4),idnc)
!        istatus = nf_create("out.ncf",OR(NF_CLOBBER,NF_64BIT_OFFSET),idnc)
        if (istatus .ne. nf_noerr) then 
           print*,'error in nf_create of out.ncf'
           call handle_err(istatus)
        endif
        istatus = nf_def_dim(idnc,"X",n1,idx)
        if (istatus .ne. nf_noerr) then 
           print*,'error in nf_def_dim x'
           call handle_err(istatus)
        endif
        istatus = nf_def_dim(idnc,"Y",n2,idy)
        if (istatus .ne. nf_noerr) then 
           print*,'error in nf_def_dim y'
           call handle_err(istatus)
        endif
        istatus = nf_def_dim(idnc,"Z",n3,idz)
        if (istatus .ne. nf_noerr) then 
           print*,'error in nf_def_dim z'
           call handle_err(istatus)
        endif
        istatus = nf_def_dim(idnc,"T",nrsp,idt)
        if (istatus .ne. nf_noerr) then 
           print*,'error in nf_def_dim t'
           call handle_err(istatus)
        endif

        ncdims(1) = idx
        ncdims(2) = idz
        ncdims(3) = idy
        ncdims(4) = idt

        istatus = nf_def_var(idnc,"TIMES",NF_FLOAT,1,idt,idtimes)
        if (istatus .ne. nf_noerr) then 
           print*,'error nf_def_var TIMES'
           call handle_err(istatus)
        endif

        if (iturb.eq.1) then
           istatus = nf_def_var(idnc,"U",NF_FLOAT,4,ncdims,idu)
           if (istatus .ne. nf_noerr) then 
              print*,'error nf_def_var U'
              call handle_err(istatus)
           endif
           istatus = nf_def_var(idnc,"V",NF_FLOAT,4,ncdims,idv)
           if (istatus .ne. nf_noerr) then 
              print*,'error nf_def_var V'
              call handle_err(istatus)
           endif
           istatus = nf_def_var(idnc,"W",NF_FLOAT,4,ncdims,idw)
           if (istatus .ne. nf_noerr) then 
              print*,'error nf_def_var W'
              call handle_err(istatus)
           endif
           istatus = nf_def_var(idnc,"ZX",NF_FLOAT,4,ncdims,idzx)
           if (istatus .ne. nf_noerr) then 
              print*,'error nf_def_var ZX'
              call handle_err(istatus)
           endif
           istatus = nf_def_var(idnc,"ZY",NF_FLOAT,4,ncdims,idzy)
           if (istatus .ne. nf_noerr) then 
              print*,'error nf_def_var ZY'
              call handle_err(istatus)
           endif
           istatus = nf_def_var(idnc,"ZZ",NF_FLOAT,4,ncdims,idzz)
           if (istatus .ne. nf_noerr) then 
              print*,'error nf_def_var ZZ'
              call handle_err(istatus)
           endif
        endif!iturb
        if (thermo .eq. 1) then
           istatus = nf_def_var(idnc,"TH",NF_FLOAT,4,ncdims,idth)
           if (istatus .ne. nf_noerr) then 
              print*,'error nf_def_var TH'
              call handle_err(istatus)
           endif
           istatus = nf_def_var(idnc,"QV",NF_FLOAT,4,ncdims,idqv)
           if (istatus .ne. nf_noerr) then 
              print*,'error nf_def_var QV'
              call handle_err(istatus)
           endif
        endif !thermo
        istatus = nf_enddef(idnc)
        if (istatus .ne. nf_noerr) then 
           print*,'error enddef for out.ncf'
           call handle_err(istatus)
        endif

        ncstart(1) = 1
        ncstart(2) = 1
        nccount(1) = n1
        nccount(2) = n3
        nccount(3) = n2pe
        nccount(4) = 1

     endif ! mype

!if ( MPI == 1) then 
  ! broadcast variable IDs to all procs
     if(iturb .eq. 1)  then
        call mpi_bcast(idu,1,MPI_INTEGER,0,MPI_COMM_WORLD,istatus) 
        call mpi_bcast(idv,1,MPI_INTEGER,0,MPI_COMM_WORLD,istatus)
        call mpi_bcast(idw,1,MPI_INTEGER,0,MPI_COMM_WORLD,istatus)
        call mpi_bcast(idzx,1,MPI_INTEGER,0,MPI_COMM_WORLD,istatus)
        call mpi_bcast(idzy,1,MPI_INTEGER,0,MPI_COMM_WORLD,istatus)
        call mpi_bcast(idzz,1,MPI_INTEGER,0,MPI_COMM_WORLD,istatus)
     endif!iturb
     if (thermo .eq. 1) then
        call mpi_bcast(idth,1,MPI_INTEGER,0,MPI_COMM_WORLD,istatus)
        call mpi_bcast(idqv,1,MPI_INTEGER,0,MPI_COMM_WORLD,istatus)
     endif!ithermo
!endif!mpi

  endif ! nrsp

  !!! next, prep restart file

  if (mype.eq.0) then

     print*,'Creating netcdf restart file'
     do 800 nf = 1,nrst
        if(nf<10) then
	   write(filename,'("Zk",I1,".out.ncf")')nf
        else
           write(filename,'("Zk",I2,".out.ncf")')nf
        endif
	!write(filename,format_string)nf
!     istatus = nf_create(filename,0,idnck(nf))
     istatus = nf_create(trim(filename),OR(NF_CLOBBER,NF_NETCDF4),idnck(nf))
!     istatus = nf_create(filename,OR(NF_CLOBBER,NF_64BIT_OFFSET),idnck(nf))
     if (istatus .ne. nf_noerr) then 
        print*,'error in nf_create of Zk.out.ncf nf=',nf
        call handle_err(istatus)
     endif
     istatus = nf_def_dim(idnck(nf),"KX",iktx,idkx)
     if (istatus .ne. nf_noerr) then 
        print*,'error nf_def_dim kx'
        call handle_err(istatus)
     endif
     istatus = nf_def_dim(idnck(nf),"KY",ikty,idky)
     if (istatus .ne. nf_noerr) then 
        print*,'error nf_def_dim ky'
        call handle_err(istatus)
     endif
     istatus = nf_def_dim(idnck(nf),"KZ",iktz,idkz)
     if (istatus .ne. nf_noerr) then 
        print*,'error nf_def_dim kz'
        call handle_err(istatus)
     endif
     istatus = nf_def_dim(idnck(nf),"RI",2,idkri) !real&imag
     if (istatus .ne. nf_noerr) then 
        print*,'error nf_def_dim ri'
        call handle_err(istatus)
     endif
     istatus = nf_def_dim(idnck(nf),"T", 1,idkt)
     if (istatus .ne. nf_noerr) then 
        print*,'error nf_def_dim t'
        call handle_err(istatus)
     endif

     ncdimsk(1) = idkx
     ncdimsk(2) = idky
     ncdimsk(3) = idkz
     ncdimsk(4) = idkri  ! dim 4 is real/imag part
     ncdimsk(5) = idkt

     istatus = nf_def_var(idnck(nf),"TIMES",NF_FLOAT,1,idkt,idtimesk)
     if (istatus .ne. nf_noerr) then 
        print*,'error nf_def_var TIMES'
        call handle_err(istatus)
     endif

     if(iturb .eq. 1) then
        istatus = nf_def_var(idnck(nf),"ZXK",NF_FLOAT,5,ncdimsk,idzxk)
        if (istatus .ne. nf_noerr) then 
          print*,'error nf_def_var ZX'
          call handle_err(istatus)
        endif
        istatus = nf_def_var(idnck(nf),"ZYK",NF_FLOAT,5,ncdimsk,idzyk)
        if (istatus .ne. nf_noerr) then 
           print*, 'error nf_def_var ZY'
           call handle_err(istatus)
        endif
        istatus = nf_def_var(idnck(nf),"ZZK",NF_FLOAT,5,ncdimsk,idzzk)
        if (istatus .ne. nf_noerr) then 
           print*,'error nf_def_var ZZ'
           call handle_err(istatus)
        endif
     endif!iturb
     if (thermo .eq. 1) then
        istatus = nf_def_var(idnck(nf),"TTK",NF_FLOAT,5,ncdimsk,idttk)
        if (istatus .ne. nf_noerr) then 
           print*,'error nf_def_var TT'
           call handle_err(istatus)
        endif
        istatus = nf_def_var(idnck(nf),"QVK",NF_FLOAT,5,ncdimsk,idqvk)
        if (istatus .ne. nf_noerr) then 
           print*,'error nf_def_var QV'
           call handle_err(istatus)
        endif
     endif!thermo
 
     istatus = nf_enddef(idnck(nf))
     if (istatus .ne. nf_noerr) then 
        print*,'enddef error for Zk.out.ncf'
        call handle_err(istatus)
     endif

800	enddo

     ncstartk(1) = 1
     ncstartk(2) = 1

     nccountk(1) = iktx
     nccountk(2) = ikty
!     nccountk(3) = iktzp*npe
     nccountk(3) = iktzp
     nccountk(4) = 1
     nccountk(5) = 1
  endif ! mype

!if ( MPI == 1 ) then
  ! broadcast variable IDs to all procs
  if(iturb .eq. 1) then
     call mpi_bcast(idzxk,1,MPI_INTEGER,0,MPI_COMM_WORLD,istatus)
     call mpi_bcast(idzyk,1,MPI_INTEGER,0,MPI_COMM_WORLD,istatus)
     call mpi_bcast(idzzk,1,MPI_INTEGER,0,MPI_COMM_WORLD,istatus)
  endif!iturb
  if (thermo .eq. 1) then
     call mpi_bcast(idttk,1,MPI_INTEGER,0,MPI_COMM_WORLD,istatus)
     call mpi_bcast(idqvk,1,MPI_INTEGER,0,MPI_COMM_WORLD,istatus)
  endif!thermo
!endif!mpi


  !!! next, prep droplet restart file

  if (gomic .ge. 1) then

     if (mype.eq.0) then

        print*,'Creating netcdf droplet restart file'
	do 801 nf = 1,nrstd
        if(nf<10) then
           write(filenamed,'("drop",I1,".out.ncf")')nf
        else
           write(filenamed,'("drop",I2,".out.ncf")')nf
        endif
	  !write(filenamed,'("drop",I1,".out.ncf")')nf
!        istatus = nf_create(filenamed,0,idncd(nf))
        istatus = nf_create(trim(filenamed),OR(NF_CLOBBER,NF_NETCDF4),idncd(nf))
!        istatus = nf_create(filenamed,OR(NF_CLOBBER,NF_64BIT_OFFSET),idncd(nf))
        if (istatus .ne. nf_noerr) then 
           print*,'error nf_create of drop.out.ncf nf=',nf
           call handle_err(istatus)
        endif
        istatus = nf_def_dim(idncd(nf),"NDROPPE",ndroppe*npe,idndroppe)
        if (istatus .ne. nf_noerr) then 
           print*,'error nf_def_dim ndroppe'
           call handle_err(istatus)
        endif
        istatus = nf_def_dim(idncd(nf),"TD",1,iddt)
        if (istatus .ne. nf_noerr) then 
           print*,'error nf_def_dim t'
           call handle_err(istatus)
        endif

        ncdimsd(1) = idndroppe
        ncdimsd(2) = iddt

        istatus = nf_def_var(idncd(nf),"TIMES",NF_FLOAT,1,iddt,idtimesd)
        if (istatus .ne. nf_noerr) then 
           print*,'error nf_def_var TIMES'
           call handle_err(istatus)
        endif

        if(thermo .eq. 1)then
        istatus = nf_def_var(idncd(nf),"PP",NF_DOUBLE,1,iddt,idPP)
        if (istatus .ne. nf_noerr) then
           print*,'error nf_def_var pp'
           call handle_err(istatus)
        endif
        istatus = nf_def_var(idncd(nf),"rhoa",NF_DOUBLE,1,iddt,idrhoa)
        if (istatus .ne. nf_noerr) then 
           print*,'error nf_def_var rhoa'
           call handle_err(istatus)
        endif
        istatus = nf_def_var(idncd(nf),"thetapp",NF_DOUBLE,1,iddt,idthetapp)
        if (istatus .ne. nf_noerr) then 
           print*,'error nf_def_var thetapp'
           call handle_err(istatus)
        endif
        istatus = nf_def_var(idncd(nf),"sp",NF_DOUBLE,1,iddt,idsp)
        if (istatus .ne. nf_noerr) then 
           print*,'error nf_def_var sp'
           call handle_err(istatus)
        endif
        istatus = nf_def_var(idncd(nf),"qvpp",NF_DOUBLE,1,iddt,idqvpp)
        if (istatus .ne. nf_noerr) then 
           print*,'error nf_def_var qvpp'
           call handle_err(istatus)
        endif
        istatus = nf_def_var(idncd(nf),"rmth",NF_FLOAT,1,iddt,idrmth)
        if (istatus .ne. nf_noerr) then 
           print*,'error nf_def_var rmth'
           call handle_err(istatus)
        endif
        istatus = nf_def_var(idncd(nf),"rmqv",NF_FLOAT,1,iddt,idrmqv)
        if (istatus .ne. nf_noerr) then 
           print*,'error nf_def_var rmqv'
           call handle_err(istatus)
        endif
        endif !thermo

        istatus = nf_def_var(idncd(nf),"IDP",NF_FLOAT,2,ncdimsd,ididp)
        if (istatus .ne. nf_noerr) then 
           print*,'error nf_def_var IDP'
           call handle_err(istatus)
        endif
        istatus = nf_def_var(idncd(nf),"XDROP",NF_FLOAT,2,ncdimsd,idxd)
        if (istatus .ne. nf_noerr) then 
           print*,'error nf_def_var XDROP'
           call handle_err(istatus)
        endif
        istatus = nf_def_var(idncd(nf),"YDROP",NF_FLOAT,2,ncdimsd,idyd)
        if (istatus .ne. nf_noerr) then 
           print*,'error nf_def_var YDROP'
           call handle_err(istatus)
        endif
        istatus = nf_def_var(idncd(nf),"ZDROP",NF_FLOAT,2,ncdimsd,idzd)
        if (istatus .ne. nf_noerr) then 
           print*,'error nf_def_var ZDROP'
           call handle_err(istatus)
        endif
        istatus = nf_def_var(idncd(nf),"R",NF_DOUBLE,2,ncdimsd,idr)
        if (istatus .ne. nf_noerr) then 
           print*,'error nf_def_var R'
           call handle_err(istatus)
        endif
        if(thermo .eq. 1)then
           istatus = nf_def_var(idncd(nf),"DR3",NF_DOUBLE,2,ncdimsd,iddr3)
           if (istatus .ne. nf_noerr) then 
              print*,'error nf_def_var DR3'
              call handle_err(istatus)
           endif
           istatus = nf_def_var(idncd(nf),"R_CCN",NF_DOUBLE,2,ncdimsd,idr_ccn)
           if (istatus .ne. nf_noerr) then
                   print*,'error nf_def_var R_CCN'
                   call handle_err(istatus)
           endif
        endif
        istatus = nf_def_var(idncd(nf),"UDROP",NF_FLOAT,2,ncdimsd,idudrop)
        if (istatus .ne. nf_noerr) then 
           print*,'error nf_def_var UDROP'
           call handle_err(istatus)
        endif
        istatus = nf_def_var(idncd(nf),"VDROP",NF_FLOAT,2,ncdimsd,idvdrop)
        if (istatus .ne. nf_noerr) then 
           print*,'error nf_def_var VDROP'
           call handle_err(istatus)
        endif
        istatus = nf_def_var(idncd(nf),"WDROP",NF_FLOAT,2,ncdimsd,idwdrop)
        if (istatus .ne. nf_noerr) then 
           print*,'error nf_def_var WDROP'
           call handle_err(istatus)
        endif

        istatus = nf_enddef(idncd(nf))
        if (istatus .ne. nf_noerr) then 
           print*,'enddef error for drop.out.ncf'
           call handle_err(istatus)
        endif
801	continue

        nccountd(1) = ndroppe*npe
        nccountd(2) = 1
     endif ! mype

  endif ! gomic

  ! broadcast variable IDs to all procs
  call mpi_bcast(ididp,1,MPI_INTEGER,0,MPI_COMM_WORLD,istatus)
  call mpi_bcast(idxd,1,MPI_INTEGER,0,MPI_COMM_WORLD,istatus)
  call mpi_bcast(idyd,1,MPI_INTEGER,0,MPI_COMM_WORLD,istatus)
  call mpi_bcast(idzd,1,MPI_INTEGER,0,MPI_COMM_WORLD,istatus)
  call mpi_bcast(idr,1,MPI_INTEGER,0,MPI_COMM_WORLD,istatus)
  if(thermo .eq. 1) then
          call mpi_bcast(iddr3,1,MPI_INTEGER,0,MPI_COMM_WORLD,istatus)
          call mpi_bcast(idr_ccn,1,MPI_INTEGER,0,MPI_COMM_WORLD,istatus)
  endif
  call mpi_bcast(idudrop,1,MPI_INTEGER,0,MPI_COMM_WORLD,istatus)
  call mpi_bcast(idvdrop,1,MPI_INTEGER,0,MPI_COMM_WORLD,istatus)
  call mpi_bcast(idwdrop,1,MPI_INTEGER,0,MPI_COMM_WORLD,istatus)
end subroutine ncprep


subroutine ncdumprsp(ur,vr,wr,zxr,zyr,zzr,ttr,qvr,junk,ntdump,time,nrsp)

! dump real space fields into out.ncf
  use netcdf_mod 
  implicit none
  include 'param.inc'
  include 'netcdf.inc'
  include 'mpif.h'

  real :: time
  real, dimension(n1d,n3d,n2dp) :: ur,vr,wr,zxr,zyr,zzr,ttr,qvr,junk
  integer :: iproc,nbuf,nsends,nrecvs,ntdump,nrsp,i,j
  integer :: istatus,ip,idvar

!if ( MPI == 1 ) then
  integer :: status(MPI_STATUS_SIZE)
!endif!mpi
  integer :: mype,ierror
  common/mpi/mype

!if ( MPI == 1 ) then
  nbuf       = n1d*n3d*n2dp
     nsends     = (npe-1)*(iturb*6+thermo*2)
     nrecvs     = iturb*6+thermo*2 
!endif!mpi

  ncstart(4) = ntdump

  if (mype.eq.0 .and. nrecvs .ne. 0) then
     print*,'Writing to out.ncf at time ',time
     istatus    = nf_put_vara_real(idnc,idtimes,ntdump,1,time)
     if (istatus .ne. nf_noerr) then
        print*, 'Error in putting time in out.ncf'
        call handle_err(istatus)
     endif
!  endif
    
!  if (mype .eq. 0) print*, 'after writing time'
     print*, 'after writing time'

!  if (mype.eq.0) then
     ncstart(3) = 1
     if (iturb.eq.1)  then
        istatus = nf_put_vara_real(idnc,idu, ncstart,nccount, ur(1:n1,:,:))
        if (istatus .ne. nf_noerr) then
           print*, 'Error in putting ur by mype 0'
           call handle_err(istatus)
        endif

        istatus = nf_put_vara_real(idnc,idv, ncstart,nccount, vr(1:n1,:,:))
        if (istatus .ne. nf_noerr) then
           print*, 'Error in putting vr by mype 0'
           call handle_err(istatus)
        endif

        istatus = nf_put_vara_real(idnc,idw, ncstart,nccount, wr(1:n1,:,:))
        if (istatus .ne. nf_noerr) then
           print*, 'Error in putting wr by mype 0'
           call handle_err(istatus)
        endif
        istatus = nf_put_vara_real(idnc,idzx,ncstart,nccount,zxr(1:n1,:,:))
        if (istatus .ne. nf_noerr) then
           print*, 'Error in putting zxr by mype 0'
           call handle_err(istatus)
        endif
        istatus = nf_put_vara_real(idnc,idzy,ncstart,nccount,zyr(1:n1,:,:))
        if (istatus .ne. nf_noerr) then
           print*, 'Error in putting zyr by mype 0'
           call handle_err(istatus)
        endif
        istatus = nf_put_vara_real(idnc,idzz,ncstart,nccount,zzr(1:n1,:,:))
        if (istatus .ne. nf_noerr) then
           print*, 'Error in putting zzr by mype 0'
           call handle_err(istatus)
        endif
     endif!iturb

     if (thermo .eq. 1) then
        istatus = nf_put_vara_real(idnc,idth,ncstart,nccount,ttr(1:n1,:,:))
        if (istatus .ne. nf_noerr) then
           print*, 'Error in putting ttr by mype 0'
           call handle_err(istatus)
        endif
        istatus = nf_put_vara_real(idnc,idqv,ncstart,nccount,qvr(1:n1,:,:))
        if (istatus .ne. nf_noerr) then
           print*, 'Error in putting qvr by mype 0'
           call handle_err(istatus)
        endif
     endif

     print*, 'after mype = 0 dumped his data'
  endif ! if mype .eq. 0 


  if (nrecvs .ne. 0) then
  do i = 1,npe-1
     if (mype .eq. 0) then
     ! receive from mype=i
        do j = 1,nrecvs
           call mpi_recv(junk,nbuf,MPI_REAL,i,MPI_ANY_TAG,MPI_COMM_WORLD,status,istatus)
           ncstart(3) = i*n2pe+1
           idvar = status(MPI_TAG)
           istatus = nf_put_vara_real(idnc,idvar,ncstart,nccount,junk(1:n1,:,:))
           if (istatus .ne. nf_noerr) then
              print*,'Error in putting variable with ncstart(3),idvar = ',ncstart(3), idvar
              call handle_err(istatus)
           endif
        enddo
        print*, 'mype = 0 has received all data for out.ncf from mype =', i
     elseif (mype .eq. i) then
     ! send to mype=0; 
        if (iturb.eq.1) then
            call mpi_send(ur, nbuf,MPI_REAL,0,idu, MPI_COMM_WORLD,istatus)
            call mpi_send(vr, nbuf,MPI_REAL,0,idv, MPI_COMM_WORLD,istatus)
            call mpi_send(wr, nbuf,MPI_REAL,0,idw, MPI_COMM_WORLD,istatus)
            call mpi_send(zxr,nbuf,MPI_REAL,0,idzx,MPI_COMM_WORLD,istatus)
            call mpi_send(zyr,nbuf,MPI_REAL,0,idzy,MPI_COMM_WORLD,istatus)
            call mpi_send(zzr,nbuf,MPI_REAL,0,idzz,MPI_COMM_WORLD,istatus)
        endif!iturb
        if (thermo .eq. 1) then
            call mpi_send(ttr,nbuf,MPI_REAL,0,idth,MPI_COMM_WORLD,istatus)
            call mpi_send(qvr,nbuf,MPI_REAL,0,idqv,MPI_COMM_WORLD,istatus)
        endif
        print*, 'mype', mype, 'done sending data for out.ncf'
     endif
     call mpi_barrier(MPI_COMM_WORLD,ierror)
  enddo ! end of i-loop
  endif !nrecvs

  if (mype .eq. 0 .and. nrecvs .ne.0) then
     if (ntdump .eq. nrsp) then
        istatus = nf_close(idnc)
        if (istatus .ne. nf_noerr) then
           print*,'Error closing out.ncf at ntdump = ',ntdump
           call handle_err(istatus)
        endif
     endif
  endif

 

end subroutine ncdumprsp


subroutine ncdumprst(zx,zy,zz,tt,qv,junk,ntdump,time)

! dump restart into out.ncf
  use netcdf_mod
  implicit none
  include 'param.inc'
  include 'netcdf.inc'
  include 'mpif.h'

  real :: time
  complex, dimension(iktx,ikty,iktzp) :: zx,zy,zz,tt,qv
  complex, dimension(iktx,ikty,iktzp) :: junk
  integer :: iproc,nbuf,nsends,nrecvs,ntdump,i,j
!  integer :: nbuf2
  integer :: istatus,ip,idvar


!  complex, dimension(iktx,ikty,iktz) :: junkall
  complex, allocatable, dimension(:,:,:) :: junkall 
!  complex, allocatable, dimension(:) :: zxall, junkall 
  integer, parameter :: opt = 0			! Opt = 1 doesn't work for the moment, so keep opt = 0

!if ( MPI == 1 ) then
  integer :: status(MPI_STATUS_SIZE)
!endif!mpi
  integer :: mype,ierror
  common/mpi/mype

  allocate( junkall(iktx,ikty,iktz) )
!  allocate( zxall(iktx,ikty,iktz) )

!if ( MPI == 1 ) then
  nbuf = iktx*ikty*iktzp
     nsends     = (npe-1)*(3*iturb+2*thermo)
     nrecvs     = 3*iturb+2*thermo
!endif!mpi

  ncstartk(5) = 1
  if (mype.eq.0 .and. nrecvs .ne. 0) then
     print*,'Writing to Zk.out at time ',time
     istatus    = nf_put_vara_real(idnck(ntdump),idtimesk,1,1,time)
     if (istatus .ne. nf_noerr) then
        print*, 'Error in putting time in Zk.out.ncf, ntdump = ', ntdump
        call handle_err(istatus)
     endif
  endif

  if ( opt .eq. 0 .and. nrecvs .ne. 0) then 

     if (mype .eq. 0) print*, 'in opt = 0 at start'



  if (mype.eq.0 ) then
     ncstartk(3) = 1
     ncstartk(4) = 1 ! (real part)
     if (iturb .eq. 1) then
        istatus = nf_put_vara_real(idnck(ntdump),idzxk,ncstartk,nccountk,real(zx))
        if (istatus .ne. nf_noerr) then
           print*,'Error in putting real(zx) by mype 0'
           call handle_err(istatus)
        endif
        istatus = nf_put_vara_real(idnck(ntdump),idzyk,ncstartk,nccountk,real(zy))
        if (istatus .ne. nf_noerr) then
           call handle_err(istatus)
           print*,'Error in putting real(zy) by mype 0'
        endif
        istatus = nf_put_vara_real(idnck(ntdump),idzzk,ncstartk,nccountk,real(zz))
        if (istatus .ne. nf_noerr) then
           call handle_err(istatus)
           print*,'Error in putting real(zz) by mype 0'
        endif
  endif !iturb

     if (thermo .eq. 1) then
        istatus = nf_put_vara_real(idnck(ntdump),idttk,ncstartk,nccountk,real(tt))
        if (istatus .ne. nf_noerr) then
           print*,'Error in putting real(tt) by mype 0'
           call handle_err(istatus)
        endif
        istatus = nf_put_vara_real(idnck(ntdump),idqvk,ncstartk,nccountk,real(qv))
        if (istatus .ne. nf_noerr) then
           print*,'Error in putting real(qv) by mype 0'
           call handle_err(istatus)
        endif
     endif!thermo

     ncstartk(4) = 2 ! (imag part)

  if (iturb .eq. 1) then
        istatus = nf_put_vara_real(idnck(ntdump),idzxk,ncstartk,nccountk,aimag(zx))
        if (istatus .ne. nf_noerr) then
           print*,'Error in putting imag(zx) by mype 0'
           call handle_err(istatus)
        endif

        istatus = nf_put_vara_real(idnck(ntdump),idzyk,ncstartk,nccountk,aimag(zy))
        if (istatus .ne. nf_noerr) then
           print*,'Error in putting imag(zy) by mype 0'
           call handle_err(istatus)
        endif

        istatus = nf_put_vara_real(idnck(ntdump),idzzk,ncstartk,nccountk,aimag(zz))
        if (istatus .ne. nf_noerr) then
           print*,'Error in putting imag(zz) by mype 0'
           call handle_err(istatus)
        endif
     endif!iturb
     if (thermo .eq. 1) then
        istatus = nf_put_vara_real(idnck(ntdump),idttk,ncstartk,nccountk,aimag(tt))
        if (istatus .ne. nf_noerr) then
           print*,'Error in putting imag(tt) by mype 0'
           call handle_err(istatus)
        endif
        istatus = nf_put_vara_real(idnck(ntdump),idqvk,ncstartk,nccountk,aimag(qv))
        if (istatus .ne. nf_noerr) then
           print*,'Error in putting imag(qv) by mype 0'
           call handle_err(istatus)
        endif
     endif!thermo
  endif ! mype

  do i = 1,npe-1
     if (mype .eq. 0) then
     ! receive from mype=i
        do j = 1,nrecvs
           call mpi_recv(junk,nbuf,MPI_COMPLEX,i,MPI_ANY_TAG,MPI_COMM_WORLD,status,istatus)
           ncstartk(3) = i*iktzp+1
           idvar = status(MPI_TAG)
           ncstartk(4) = 1 ! (real part)
           istatus = nf_put_vara_real(idnck(ntdump),idvar,ncstartk,nccountk,real(junk))
           if (istatus .ne. nf_noerr) then
              print*,'Error in putting variable with ncstartk(3),idvar = ',ncstartk(3), idvar
              call handle_err(istatus)
           endif
           ncstartk(4) = 2 ! (imag part)
           istatus = nf_put_vara_real(idnck(ntdump),idvar,ncstartk,nccountk,aimag(junk))
           if (istatus .ne. nf_noerr) then
              print*,'Error in putting variable with ncstartk(3),idvar = ',ncstartk(3), idvar
              call handle_err(istatus)
           endif
        enddo
        print*, 'mype = 0 has received all data for Zk.out.ncf from mype =', i
 
     elseif (mype .eq. i) then
     ! send to mype=0
     if (iturb .eq. 1) then
        call mpi_send(zx,nbuf,MPI_COMPLEX,0,idzxk,MPI_COMM_WORLD,istatus)
        call mpi_send(zy,nbuf,MPI_COMPLEX,0,idzyk,MPI_COMM_WORLD,istatus)
        call mpi_send(zz,nbuf,MPI_COMPLEX,0,idzzk,MPI_COMM_WORLD,istatus)
     endif
        if (thermo .eq. 1) then
           call mpi_send(tt,nbuf,MPI_COMPLEX,0,idttk,MPI_COMM_WORLD,istatus)
           call mpi_send(qv,nbuf,MPI_COMPLEX,0,idqvk,MPI_COMM_WORLD,istatus)
        endif
        print*, 'mype', mype, 'done sending data for Zk.out.ncf'

     endif
     call mpi_barrier(MPI_COMM_WORLD,ierror)
  enddo

  if (mype .eq. 0 ) then   
        istatus = nf_close(idnck(ntdump))
        if (istatus .ne. nf_noerr) then
           print*,'Error closing Zk.out.ncf at ntdump = ',ntdump
           call handle_err(istatus)
        endif
  endif   

  elseif (opt*nrecvs .ne. 0) then

     ncstartk(3) = 1
 
     if (mype .eq. 0) print*, 'in opt = 1 at start'
!     nbuf2 = iktx*ikty*iktzp

!	if ( MPI == 1 ) then
     		call mpi_gather(zx,nbuf,MPI_COMPLEX,junkall,nbuf,MPI_COMPLEX,0,MPI_COMM_WORLD,ierror)
!	endif!MPI
     if (mype .eq. 0) print*, 'gather of zx successfull'
     if (mype .eq. 0) then
           ncstartk(4) = 1 ! (real part)
           istatus = nf_put_vara_real(idnck(ntdump),idzxk,ncstartk,nccountk,real(junkall))
           if (istatus.ne.0) print*,'Error in putting real of variable zx with opt = 1 '
           ncstartk(4) = 2 ! (imag part)
           istatus = nf_put_vara_real(idnck(ntdump),idzxk,ncstartk,nccountk,aimag(junkall))
           if (istatus.ne.0) print*,'Error in putting aimag of variable zx with opt = 1 '
     endif
     call mpi_barrier(MPI_COMM_WORLD,ierror)
     if (mype .eq. 0) print*, 'done putting zx in Zk.out.ncf' 

!	if ( MPI == 1 ) then
     		call mpi_gather(zy,nbuf,MPI_COMPLEX,junkall,nbuf,MPI_COMPLEX,0,MPI_COMM_WORLD,ierror)
!	endif!mpi
     if (mype .eq. 0) then
           ncstartk(4) = 1 ! (real part)
           istatus = nf_put_vara_real(idnck(ntdump),idzyk,ncstartk,nccountk,real(junkall))
           if (istatus.ne.0) print*,'Error in putting real of variable zy with opt = 1 '
           ncstartk(4) = 2 ! (imag part)
           istatus = nf_put_vara_real(idnck(ntdump),idzyk,ncstartk,nccountk,aimag(junkall))
           if (istatus.ne.0) print*,'Error in putting aimag of variable zy with opt = 1 '
     endif
     call mpi_barrier(MPI_COMM_WORLD,ierror)
     if (mype .eq. 0) print*, 'done putting zy in Zk.out.ncf' 

!	if ( MPI == 1 ) then
     		call mpi_gather(zz,nbuf,MPI_COMPLEX,junkall,nbuf,MPI_COMPLEX,0,MPI_COMM_WORLD,ierror)
!	endif
     if (mype .eq. 0) then
           ncstartk(4) = 1 ! (real part)
           istatus = nf_put_vara_real(idnck(ntdump),idzzk,ncstartk,nccountk,real(junkall))
           if (istatus.ne.0) print*,'Error in putting real of variable zz with opt = 1 '
           ncstartk(4) = 2 ! (imag part)
           istatus = nf_put_vara_real(idnck(ntdump),idzzk,ncstartk,nccountk,aimag(junkall))
           if (istatus.ne.0) print*,'Error in putting aimag of variable zz with opt = 1 '
     endif
     call mpi_barrier(MPI_COMM_WORLD,ierror)
     if (mype .eq. 0) print*, 'done putting zz in Zk.out.ncf' 

     if (thermo .eq. 1) then 
!	if ( MPI == 1 ) then
     		call mpi_gather(tt,nbuf,MPI_COMPLEX,junkall,nbuf,MPI_COMPLEX,0,MPI_COMM_WORLD,ierror)
!	endif!mpi
        if (mype .eq. 0) then
              ncstartk(4) = 1 ! (real part)
              istatus = nf_put_vara_real(idnck(ntdump),idttk,ncstartk,nccountk,real(junkall))
              if (istatus.ne.0) print*,'Error in putting real of variable tt with opt = 1 '
              ncstartk(4) = 2 ! (imag part)
              istatus = nf_put_vara_real(idnck(ntdump),idttk,ncstartk,nccountk,aimag(junkall))
              if (istatus.ne.0) print*,'Error in putting aimag of variable tt with opt = 1 '
        endif
        call mpi_barrier(MPI_COMM_WORLD,ierror)
        if (mype .eq. 0) print*, 'done putting tt in Zk.out.ncf' 

!	if ( MPI == 1 ) then
     		call mpi_gather(qv,nbuf,MPI_COMPLEX,junkall,nbuf,MPI_COMPLEX,0,MPI_COMM_WORLD,ierror)
!	endif!mpi
        if (mype .eq. 0) then
              ncstartk(4) = 1 ! (real part)
              istatus = nf_put_vara_real(idnck(ntdump),idqvk,ncstartk,nccountk,real(junkall))
              if (istatus.ne.0) print*,'Error in putting real of variable qv with opt = 1 '
              ncstartk(4) = 2 ! (imag part)
              istatus = nf_put_vara_real(idnck(ntdump),idqvk,ncstartk,nccountk,aimag(junkall))
              if (istatus.ne.0) print*,'Error in putting aimag of variable qv with opt = 1 '
        endif
        call mpi_barrier(MPI_COMM_WORLD,ierror)
        if (mype .eq. 0) print*, 'done putting qv in Zk.out.ncf' 
     endif ! thermo

     if (mype .eq. 0 ) then   
           istatus = nf_close(idnck(ntdump))
           if (istatus.ne.0) print*,'Error closing Zk.out at ntdump = ',ntdump
     endif   
 
  endif ! opt

!  endif ! mype

  deallocate( junkall )

end subroutine ncdumprst


subroutine ncdumpdrop(ndropreal,ntdump,time,nddone)

! dump droplet restart into drop.out.ncf
  use netcdf_mod
  use mic_mod 
  implicit none
  include 'param.inc'
  include 'netcdf.inc'
  include 'mpif.h'

  real :: time

  ! ---- argument ----
  real, dimension(ndroppe) :: junk
  real(8), dimension(ndroppe) :: junk1
  ! ----- local ----
  real, allocatable, dimension(:) :: idpreal
  integer :: iproc,nbuf,nsends_real,nsends_real8,nrecvs_real,nrecvs_real8,ntdump,ndropreal,i,j,nddone
  integer :: istatus,ip,idvar

  real, allocatable, dimension(:) :: xall,yall,zall,udropall,vdropall,wdropall,idprealall
  real(8),allocatable,dimension(:) :: rall,dr3all,r_ccnall
  integer, parameter :: opt = 0			! Opt = 0 doesn't work for the moment, so keep opt = 1


!if ( MPI == 1 ) then
  integer :: status(MPI_STATUS_SIZE)
!endif!mpi
  integer :: mype,ierror
  common/mpi/mype

  allocate( idpreal(ndroppe) )
  allocate( xall(ndroppe*npe), yall(ndroppe*npe), zall(ndroppe*npe), rall(ndroppe*npe), dr3all(ndroppe*npe),r_ccnall(ndroppe*npe),udropall(ndroppe*npe), vdropall(ndroppe*npe), wdropall(ndroppe*npe), idprealall(ndroppe*npe) )

!if ( MPI == 1 ) then
  nbuf = ndroppe
  nsends_real = (npe-1)*7
  nsends_real8 = (npe-1)*(1+thermo*2)
  nrecvs_real=7
  nrecvs_real8 = (1+thermo*2)
!endif!mpi
  junk=udrop
  ncstartd(2) = 1

  if (mype.eq.0) then
     print*,'Writing to drop.out at time ',time
     istatus    = nf_put_vara_real(idncd(ntdump),idtimesd,1,1,time)
     if (istatus .ne. nf_noerr) then
        print*, 'Error in putting time in drop.out.ncf, ntdump =', ntdump
        call handle_err(istatus)
     endif

     if(thermo .eq. 1) then
     istatus    = nf_put_vara_double(idncd(ntdump),idPP,1,1,PP)
     if (istatus .ne. nf_noerr) then
        print*, 'Error in putting PP in drop.out.ncf, ntdump =', ntdump
        call handle_err(istatus)
     endif
     istatus    = nf_put_vara_double(idncd(ntdump),idrhoa,1,1,rhoa)
     if (istatus .ne. nf_noerr) then
        print*, 'Error in putting rhoa in drop.out.ncf, ntdump =', ntdump
        call handle_err(istatus)
     endif
     istatus    = nf_put_vara_double(idncd(ntdump),idthetapp,1,1,thetapp)
     if (istatus .ne. nf_noerr) then
        print*, 'Error in putting thetapp in drop.out.ncf, ntdump =', ntdump
        call handle_err(istatus)
     endif
     istatus    = nf_put_vara_double(idncd(ntdump),idsp,1,1,sp)
     if (istatus .ne. nf_noerr) then
        print*, 'Error in putting sp in drop.out.ncf, ntdump =', ntdump
        call handle_err(istatus)
     endif
     istatus    = nf_put_vara_double(idncd(ntdump),idqvpp,1,1,qvpp)
     if (istatus .ne. nf_noerr) then
        print*, 'Error in putting qvpp in drop.out.ncf, ntdump =', ntdump
        call handle_err(istatus)
     endif
     istatus    = nf_put_vara_real(idncd(ntdump),idrmth,1,1,rmth)
     if (istatus .ne. nf_noerr) then
        print*, 'Error in putting rmth in drop.out.ncf, ntdump =', ntdump
        call handle_err(istatus)
     endif
     istatus    = nf_put_vara_real(idncd(ntdump),idrmqv,1,1,rmqv)
     if (istatus .ne. nf_noerr) then
        print*, 'Error in putting rmqv in drop.out.ncf, ntdump =', ntdump
        call handle_err(istatus)
     endif
     endif!thermo

  endif!mype 0

! Save idp of droplets as real. Use last element to store ndropreal
! and second to last element to store nddone (used for collisional growth)
  idpreal = real(idp)
  idpreal(ndroppe) = real(ndropreal)
  idpreal(ndroppe-1) = real(nddone)

  if (opt .eq. 0) then 

     if (mype.eq.0) then
        ncstartd(1) = 1
        nccountd(1) = ndroppe
        istatus = nf_put_vara_real(idncd(ntdump),ididp,ncstartd,nccountd,idpreal)
        if (istatus .ne. nf_noerr) then
           print*,'Error in putting idpreal by mype 0'
           call handle_err(istatus)
        endif
        istatus = nf_put_vara_real(idncd(ntdump),idxd,ncstartd,nccountd,x)
        if (istatus .ne. nf_noerr) then
           print*,'Error in putting x by mype 0'
           call handle_err(istatus)
        endif
        istatus = nf_put_vara_real(idncd(ntdump),idyd,ncstartd,nccountd,y)
        if (istatus .ne. nf_noerr) then
           print*,'Error in putting y by mype 0'
           call handle_err(istatus)
        endif
        istatus = nf_put_vara_real(idncd(ntdump),idzd,ncstartd,nccountd,z)
        if (istatus .ne. nf_noerr) then
           print*,'Error in putting z by mype 0'
           call handle_err(istatus)
        endif
        istatus = nf_put_vara_real(idncd(ntdump),idudrop,ncstartd,nccountd,udrop)
        if (istatus .ne. nf_noerr) then
           print*,'Error in putting udrop by mype 0'
           call handle_err(istatus)
        endif
        istatus = nf_put_vara_real(idncd(ntdump),idvdrop,ncstartd,nccountd,vdrop)
        if (istatus .ne. nf_noerr) then
           print*,'Error in putting vdrop by mype 0'
           call handle_err(istatus)
        endif
        istatus = nf_put_vara_real(idncd(ntdump),idwdrop,ncstartd,nccountd,wdrop)
        if (istatus .ne. nf_noerr) then
           print*,'Error in putting wdrop by mype 0'
           call handle_err(istatus)
        endif 
        istatus = nf_put_vara_double(idncd(ntdump),idr,ncstartd,nccountd,r)
        if (istatus .ne. nf_noerr) then
           print*,'Error in putting r by mype 0'
           call handle_err(istatus)
        endif
        if(thermo.eq.1) then
           istatus = nf_put_vara_double(idncd(ntdump),iddr3,ncstartd,nccountd,dr3)
           if (istatus .ne. nf_noerr) then
              print*,'Error in putting dr3 by mype 0'
              call handle_err(istatus)
           endif
           istatus = nf_put_vara_double(idncd(ntdump),idr_ccn,ncstartd,nccountd,r_ccn)
           if (istatus .ne. nf_noerr) then
              print*,'Error in putting r_ccn by mype 0'
              call handle_err(istatus)
           endif
        endif
     endif ! mype .eq. 0  


  do i = 1,npe-1
     if (mype .eq. 0) then
     ! receive all vars from mype=i and collected in junk
        do j = 1,nrecvs_real
              call mpi_recv(junk,nbuf,MPI_REAL,i,MPI_ANY_TAG,MPI_COMM_WORLD,status,istatus)
           ncstartd(1) = i*ndroppe+1
           idvar = status(MPI_TAG)
           if (idvar .lt. ididp .or. idvar .gt. idwdrop) then
              print*, 'Iproc =', iproc, 'has wrong idvar =' ,idvar, 'coming from mype', status(MPI_SOURCE)
              nsends_real = nsends_real + 1
              print*, 'New nsends =', nsends_real 
           else
              istatus = nf_put_vara_real(idncd(ntdump),idvar,ncstartd,nccountd,junk)
              if (istatus .ne. nf_noerr) then
                 print*,'Error in putting variable with ncstartd(1),idvar = ',ncstartd(1), idvar
                 call handle_err(istatus)
              endif
           endif
        enddo
        print*, 'mype = 0 has received all data for drop.out.ncf from mype =', i
 
     elseif (mype .eq. i) then
     ! send to mype=0
        call mpi_send(idpreal,nbuf,MPI_REAL,0,ididp,MPI_COMM_WORLD,istatus)
        call mpi_send(x,nbuf,MPI_REAL,0,idxd,MPI_COMM_WORLD,istatus)
        call mpi_send(y,nbuf,MPI_REAL,0,idyd,MPI_COMM_WORLD,istatus)
        call mpi_send(z,nbuf,MPI_REAL,0,idzd,MPI_COMM_WORLD,istatus)
       ! call mpi_send(r,nbuf,MPI_REAL8,0,idr,MPI_COMM_WORLD,istatus)
        !if(thermo .eq. 1) then
        !        call mpi_send(dr3,nbuf,MPI_REAL8,0,iddr3,MPI_COMM_WORLD,istatus)
        !        call mpi_send(r_ccn,nbuf,MPI_REAL8,0,idr_ccn,MPI_COMM_WORLD,istatus)
        !endif
        call mpi_send(udrop,nbuf,MPI_REAL,0,idudrop,MPI_COMM_WORLD,istatus)
        call mpi_send(vdrop,nbuf,MPI_REAL,0,idvdrop,MPI_COMM_WORLD,istatus)
        call mpi_send(wdrop,nbuf,MPI_REAL,0,idwdrop,MPI_COMM_WORLD,istatus)
        print*, 'mype', mype, 'done sending data for drop.out.ncf'

     endif
     call mpi_barrier(MPI_COMM_WORLD,ierror)
  enddo
! sending real8 data
  do i = 1,npe-1
     if (mype .eq. 0) then
     ! receive all vars from mype=i and collected in junk
        do j = 1,nrecvs_real8
              call mpi_recv(junk1,nbuf,MPI_REAL8,i,MPI_ANY_TAG,MPI_COMM_WORLD,status,istatus)
           ncstartd(1) = i*ndroppe+1
           idvar = status(MPI_TAG)
           if (idvar .lt. ididp .or. idvar .gt. idwdrop) then
              print*, 'Iproc =', iproc, 'has wrong idvar =' ,idvar, 'coming from mype', status(MPI_SOURCE)
              nsends_real8 = nsends_real8 + 1
              print*, 'New nsends =', nsends_real8
           else
              istatus = nf_put_vara_double(idncd(ntdump),idvar,ncstartd,nccountd,junk1)
              if (istatus .ne. nf_noerr) then
                 print*,'Error in putting variable with ncstartd(1),idvar = ',ncstartd(1), idvar
                 call handle_err(istatus)
              endif
           endif
        enddo
        print*, 'mype = 0 has received all data for drop.out.ncf from mype =', i
     elseif (mype .eq. i) then
     ! send to mype=0
        call mpi_send(r,nbuf,MPI_REAL8,0,idr,MPI_COMM_WORLD,istatus)
        if(thermo .eq. 1) then
                call mpi_send(dr3,nbuf,MPI_REAL8,0,iddr3,MPI_COMM_WORLD,istatus)
                call mpi_send(r_ccn,nbuf,MPI_REAL8,0,idr_ccn,MPI_COMM_WORLD,istatus)
        endif
        print*, 'mype', mype, 'done sending real8 type data for drop.out.ncf'

     endif
     call mpi_barrier(MPI_COMM_WORLD,ierror)
  enddo


  if (mype .eq. 0) then
        istatus = nf_close(idncd(ntdump))
        if (istatus .ne. nf_noerr) then
           print*,'Error closing drop.out at ntdump = ',ntdump
           call handle_err(istatus)
        endif
  endif


  elseif (opt .eq. 1) then

     call mpi_gather(idpreal,ndroppe,MPI_REAL,idprealall,ndroppe,MPI_REAL,0,MPI_COMM_WORLD,ierror)
     call mpi_gather(x,ndroppe,MPI_REAL,xall,ndroppe,MPI_REAL,0,MPI_COMM_WORLD,ierror)
     call mpi_gather(y,ndroppe,MPI_REAL,yall,ndroppe,MPI_REAL,0,MPI_COMM_WORLD,ierror)
     call mpi_gather(z,ndroppe,MPI_REAL,zall,ndroppe,MPI_REAL,0,MPI_COMM_WORLD,ierror)
     call mpi_gather(r,ndroppe,MPI_REAL8,rall,ndroppe,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
     if(thermo .eq. 1) then
             call mpi_gather(dr3,ndroppe,MPI_REAL8,dr3all,ndroppe,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
             call mpi_gather(r_ccn,ndroppe,MPI_REAL8,r_ccnall,ndroppe,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
     endif
     call mpi_gather(udrop,ndroppe,MPI_REAL,udropall,ndroppe,MPI_REAL,0,MPI_COMM_WORLD,ierror)
     call mpi_gather(vdrop,ndroppe,MPI_REAL,vdropall,ndroppe,MPI_REAL,0,MPI_COMM_WORLD,ierror)
     call mpi_gather(wdrop,ndroppe,MPI_REAL,wdropall,ndroppe,MPI_REAL,0,MPI_COMM_WORLD,ierror)

     ncstartd(1) = 1

     if (mype .eq. 0) then
        istatus = nf_put_vara_real(idncd(ntdump),ididp,ncstartd,nccountd,idprealall)
        if (istatus.ne.0) print*,'Error in putting idpreal by mype 0'
        istatus = nf_put_vara_real(idncd(ntdump),idxd,ncstartd,nccountd,xall)
        if (istatus.ne.0) print*,'Error in putting x by mype 0'
        istatus = nf_put_vara_real(idncd(ntdump),idyd,ncstartd,nccountd,yall)
        if (istatus.ne.0) print*,'Error in putting y by mype 0'
        istatus = nf_put_vara_real(idncd(ntdump),idzd,ncstartd,nccountd,zall)
        if (istatus.ne.0) print*,'Error in putting z by mype 0'
        istatus = nf_put_vara_double(idncd(ntdump),idr,ncstartd,nccountd,rall)
        if (istatus.ne.0) print*,'Error in putting r by mype 0'
        if(thermo.eq.1)then
          istatus = nf_put_vara_double(idncd(ntdump),iddr3,ncstartd,nccountd,dr3all)
          if (istatus.ne.0) print*,'Error in putting dr3 by mype 0'
          istatus = nf_put_vara_double(idncd(ntdump),idr_ccn,ncstartd,nccountd,r_ccnall)
          if (istatus.ne.0) print*,'Error in putting r_ccn by mype 0'
        endif
        istatus = nf_put_vara_real(idncd(ntdump),idudrop,ncstartd,nccountd,udropall)
        if (istatus.ne.0) print*,'Error in putting udrop by mype 0'
        istatus = nf_put_vara_real(idncd(ntdump),idvdrop,ncstartd,nccountd,vdropall)
        if (istatus.ne.0) print*,'Error in putting vdrop by mype 0'
        istatus = nf_put_vara_real(idncd(ntdump),idwdrop,ncstartd,nccountd,wdropall)
        if (istatus.ne.0) print*,'Error in putting wdrop by mype 0'
!  endif ! mype

           istatus = nf_close(idncd(ntdump))
           if (istatus.ne.0) print*,'Error closing drop.out at ntdump = ',ntdump

     endif ! mype

  endif ! opt

  deallocate( idpreal )
  deallocate( xall,yall,zall,rall,dr3all,r_ccnall,udropall,vdropall,wdropall,idprealall )

end subroutine ncdumpdrop


subroutine ncreadrst(zx,zy,zz,tt,qv,time)
  use netcdf_mod
  implicit none
  include 'param.inc'
  include 'netcdf.inc'
  include 'mpif.h'

  real :: time
  complex, dimension(iktx,ikty,iktzp) :: zx,zy,zz,tt,qv
  real, dimension(iktx,ikty,iktzp) :: wr,wi
  integer iktx1,ikty1,iktz1,nwrites,iread,nbuf
  complex zi

  integer :: idncrst,nrecvs

  integer istatus,idkx,idky,idkz,idkri,idkt
  integer idtimeskrst,idzxkrst,idzykrst,idzzkrst,idttkrst,idqvkrst
  integer, dimension(5) :: ncstartrk,ncstartik,nccountkrst

!if ( MPI == 1 ) then
  integer :: status(MPI_STATUS_SIZE)
!endif!mpi
  integer :: mype
  common/mpi/mype


  zi  = cmplx(0.,1.)
  nbuf = iktx*ikty*iktzp
  nrecvs     = iturb*6+thermo*2

  if (mype.eq.0 .and. nrecvs .ne. 0) then

     print*, 'Restarting from output data Zk.in '
     istatus = nf_open('Zk.in.ncf',NF_NOWRITE,idncrst)
     if (istatus .ne. nf_noerr) then
        print*,'Zk.in bad idncrst'
        call handle_err(istatus)
     endif 

     ! get dimension IDs
     istatus = nf_inq_dimid(idncrst,'KX',idkx)
     if (istatus .ne. nf_noerr) then
        print*,'Error getting kx ID'
        call handle_err(istatus)
     endif 
     istatus = nf_inq_dimid(idncrst,'KY',idky)
     if (istatus .ne. nf_noerr) then
        print*,'Error getting ky ID'
        call handle_err(istatus)
     endif 
     istatus = nf_inq_dimid(idncrst,'KZ',idkz)
     if (istatus .ne. nf_noerr) then
        print*,'Error getting kz ID'
        call handle_err(istatus)
     endif 
     istatus = nf_inq_dimid(idncrst,'RI',idkri)
     if (istatus .ne. nf_noerr) then
        print*,'Error getting real/imag ID'
        call handle_err(istatus)
     endif 
     istatus = nf_inq_dimid(idncrst,'T', idkt)
     if (istatus .ne. nf_noerr) then
        print*,'Error getting time ID'
        call handle_err(istatus)
     endif 

     ! get dimension sizes
     istatus = nf_inq_dimlen(idncrst,idkx,iktx1)  
     if (istatus .ne. nf_noerr) then
        print*,'Error getting idkx1'
        call handle_err(istatus)
     endif 
     istatus = nf_inq_dimlen(idncrst,idky,ikty1)  
     if (istatus .ne. nf_noerr) then
        print*,'Error getting idky1'
        call handle_err(istatus)
     endif 
     istatus = nf_inq_dimlen(idncrst,idkz,iktz1)  
     if (istatus .ne. nf_noerr) then
        print*,'Error getting idkz1'
        call handle_err(istatus)
     endif 
     istatus = nf_inq_dimlen(idncrst,idkt,nwrites)
     if (istatus .ne. nf_noerr) then
        print*,'Error getting nwrites'
        call handle_err(istatus)
     endif 

     if (iktx1*ikty1*iktz1.ne.iktx*ikty*iktz) then
        print*,'Sorry, do not know how to change resolution.  Use restartres.F90'
        stop
     endif

     ! get variable IDs
     istatus = nf_inq_varid(idncrst,'TIMES',idtimeskrst)
     if (istatus .ne. nf_noerr) then
        print*,'Error getting times ID'
        call handle_err(istatus)
     endif 
     if (iturb .eq. 1) then
        istatus = nf_inq_varid(idncrst,'ZXK',idzxkrst)
        if (istatus .ne. nf_noerr) then
           print*,'Error getting zx ID'
           call handle_err(istatus)
        endif 
        istatus = nf_inq_varid(idncrst,'ZYK',idzykrst)
        if (istatus .ne. nf_noerr) then
           print*,'Error getting zy ID'
           call handle_err(istatus)
        endif 
        istatus = nf_inq_varid(idncrst,'ZZK',idzzkrst)
        if (istatus .ne. nf_noerr) then
           print*,'Error getting zz ID'
           call handle_err(istatus)
        endif 
     endif !iturb
     if (thermo .eq. 1) then
        istatus = nf_inq_varid(idncrst,'TTK',idttkrst)
        if (istatus .ne. nf_noerr) then
           print*,'Error getting tt ID'
           call handle_err(istatus)
        endif 
        istatus = nf_inq_varid(idncrst,'QVK',idqvkrst)
        if (istatus .ne. nf_noerr) then
           print*,'Error getting qv ID'
           call handle_err(istatus)
        endif 
     endif
     
     ! prep netcdf read
     ncstartrk(1) = 1
     ncstartik(1) = 1
     ncstartrk(2) = 1
     ncstartik(2) = 1
     ncstartrk(4) = 1 ! for real part
     ncstartik(4) = 2 ! for imaginary part
     ncstartrk(5) = 1
     ncstartik(5) = 1

     nccountkrst(1) = iktx
     nccountkrst(2) = ikty
     nccountkrst(3) = iktzp
     nccountkrst(4) = 1
     nccountkrst(5) = 1

     ! read in time
     istatus = nf_get_vara_real(idncrst,idtimeskrst,1,1,time)
     if (istatus .ne. nf_noerr) then
        print*,'Error reading time'
        call handle_err(istatus)
     endif 
     print*, 'in restart time =', time

     do iread=npe-1,0,-1  ! read one block at a time; read first block last

        ncstartrk(3) = iread*iktzp+1 
        ncstartik(3) = iread*iktzp+1 
        if (iturb .eq. 1) then
           istatus = nf_get_vara_real(idncrst,idzxkrst,ncstartrk,nccountkrst,wr)
           if (istatus .ne. nf_noerr) then
              print*,'Error reading real(zx)'
              call handle_err(istatus)
           endif 
           istatus = nf_get_vara_real(idncrst,idzxkrst,ncstartik,nccountkrst,wi)
           if (istatus .ne. nf_noerr) then
              print*,'Error reading imag(zx)'
              call handle_err(istatus)
           endif 
           zx = wr + zi*wi
        
           istatus = nf_get_vara_real(idncrst,idzykrst,ncstartrk,nccountkrst,wr)
           if (istatus .ne. nf_noerr) then
              print*,'Error reading real(zy)'
              call handle_err(istatus)
           endif 
           istatus = nf_get_vara_real(idncrst,idzykrst,ncstartik,nccountkrst,wi)
           if (istatus .ne. nf_noerr) then
              print*,'Error reading imag(zy)'
              call handle_err(istatus)
           endif 
           zy = wr + zi*wi
        
           istatus = nf_get_vara_real(idncrst,idzzkrst,ncstartrk,nccountkrst,wr)
           if (istatus .ne. nf_noerr) then
              print*,'Error reading real(zz)'
              call handle_err(istatus)
           endif 
           istatus = nf_get_vara_real(idncrst,idzzkrst,ncstartik,nccountkrst,wi)
           if (istatus .ne. nf_noerr) then
              print*,'Error reading imag(zz)'
              call handle_err(istatus)
           endif 
           zz = wr + zi*wi
        endif!iturb 
       
        if (thermo .eq. 1) then
  
           istatus = nf_get_vara_real(idncrst,idttkrst,ncstartrk,nccountkrst,wr)
           if (istatus .ne. nf_noerr) then
              print*,'Error reading real(tt)'
              call handle_err(istatus)
           endif 
           istatus = nf_get_vara_real(idncrst,idttkrst,ncstartik,nccountkrst,wi)
           if (istatus .ne. nf_noerr) then
              print*,'Error reading imag(tt)'
              call handle_err(istatus)
           endif 
           tt = wr + zi*wi

           istatus = nf_get_vara_real(idncrst,idqvkrst,ncstartrk,nccountkrst,wr)
           if (istatus .ne. nf_noerr) then
              print*,'Error reading real(qv)'
              call handle_err(istatus)
           endif 
           istatus = nf_get_vara_real(idncrst,idqvkrst,ncstartik,nccountkrst,wi)
           if (istatus .ne. nf_noerr) then
              print*,'Error reading imag(qv)'
              call handle_err(istatus)
           endif 
           qv = wr + zi*wi

        endif

!	if ( MPI == 1 ) then       
        if (iread.gt.0) then  ! send to processor iread
           if(iturb .eq. 1) then
              call mpi_send(zx,nbuf,MPI_COMPLEX,iread,0,MPI_COMM_WORLD,istatus)
              if (istatus.ne.0) print*,'Error sending zx to proc', iread
              call mpi_send(zy,nbuf,MPI_COMPLEX,iread,0,MPI_COMM_WORLD,istatus)
              if (istatus.ne.0) print*,'Error sending zy to proc', iread
              call mpi_send(zz,nbuf,MPI_COMPLEX,iread,0,MPI_COMM_WORLD,istatus)
              if (istatus.ne.0) print*,'Error sending zz to proc', iread
           endif
           if (thermo .eq. 1) then
              call mpi_send(tt,nbuf,MPI_COMPLEX,iread,0,MPI_COMM_WORLD,istatus)
              if (istatus.ne.0) print*,'Error sending tt to proc', iread
              call mpi_send(qv,nbuf,MPI_COMPLEX,iread,0,MPI_COMM_WORLD,istatus)
              if (istatus.ne.0) print*,'Error sending qv to proc', iread
           endif
        endif!iread
!	endif!MPI
     enddo
  endif ! mype .eq. 0 & nrecvs

if ( nrecvs .ne. 0) then
  if (mype.gt.0) then ! get block from proc 0
     if (iturb .eq. 1) then
        call mpi_recv(zx,nbuf,MPI_COMPLEX,0,MPI_ANY_TAG,MPI_COMM_WORLD,status,istatus)
        if (istatus.ne.0) print*,mype,'Error receiving zx'
        call mpi_recv(zy,nbuf,MPI_COMPLEX,0,MPI_ANY_TAG,MPI_COMM_WORLD,status,istatus)
        if (istatus.ne.0) print*,mype,'Error receiving zy'
        call mpi_recv(zz,nbuf,MPI_COMPLEX,0,MPI_ANY_TAG,MPI_COMM_WORLD,status,istatus)
        if (istatus.ne.0) print*,mype,'Error receiving zz'
     endif
     if (thermo .eq. 1) then
        call mpi_recv(tt,nbuf,MPI_COMPLEX,0,MPI_ANY_TAG,MPI_COMM_WORLD,status,istatus)
        if (istatus.ne.0) print*,mype,'Error receiving tt'
        call mpi_recv(qv,nbuf,MPI_COMPLEX,0,MPI_ANY_TAG,MPI_COMM_WORLD,status,istatus)
        if (istatus.ne.0) print*,mype,'Error receiving qv'
     endif
  endif!mype
endif!nrecvs

  if (mype.eq.0 .and. nrecvs .ne. 0) then
     istatus = nf_close(idncrst)
     if (istatus .ne. nf_noerr) then
        print*,'Error closing Zk.in'
        call handle_err(istatus)
     endif 
  endif

end subroutine ncreadrst




subroutine ncreaddrop(ndropreal,time,nddone) 
  use mic_mod
  implicit none
  include 'param.inc'
  include 'netcdf.inc'
  include 'mpif.h'

  real :: time
  real, allocatable, dimension(:) :: idpreal
  integer :: nndroppe,nwrites,nbuf,ndropreal,iread,nddone

  integer :: istatus,idncdd,idndroppe,iddt
  integer :: idtimesd,ididp,idxd,idyd,idzd,idr,iddr3,idr_ccn,idudrop,idvdrop,idwdrop
  integer :: idPP,idrhoa,idthetapp,idsp,idqvpp,idrmth,idrmqv
  integer, dimension(2) :: ncstartd,nccountd

!if ( MPI == 1 ) then
  integer :: status(MPI_STATUS_SIZE)
!endif!mpi
  integer :: mype
  common/mpi/mype

  allocate( idpreal(ndroppe) )

  nbuf = ndroppe

  if (mype.eq.0) then

     print*, 'Restarting from output data drop.in '
     istatus = nf_open('drop.in.ncf',NF_NOWRITE,idncdd)
     if (istatus .ne. nf_noerr) then
        print*,'drop.in bad idncdd'
        call handle_err(istatus)
     endif 

     ! get dimension IDs
     istatus = nf_inq_dimid(idncdd,'NDROPPE',idndroppe)
     if (istatus .ne. nf_noerr) then
        print*,'Error getting ndroppe ID'
        call handle_err(istatus)
     endif 
     istatus = nf_inq_dimid(idncdd,'TD', iddt)
     if (istatus .ne. nf_noerr) then
        print*,'Error getting time ID'
        call handle_err(istatus)
     endif 

     ! get dimension sizes
     istatus = nf_inq_dimlen(idncdd,idndroppe,nndroppe)  
     if (istatus .ne. nf_noerr) then
        print*,'Error getting nndroppe'
        call handle_err(istatus)
     endif 
     print*, 'number of droplets per cpu = ', nndroppe/npe
     if (nndroppe/npe .ne. ndroppe ) then
        print*,'Sorry, do not know how to restart with more/less droplets'
        stop
     endif
     istatus = nf_inq_dimlen(idncdd,iddt,nwrites)
     if (istatus .ne. nf_noerr) then
        print*,'Error getting nwrites'
        call handle_err(istatus)
     endif 
     print*, 'nwrites =', nwrites

     ! get variable IDs
     istatus = nf_inq_varid(idncdd,'TIMES',idtimesd)
     if (istatus .ne. nf_noerr) then
        print*,'Error getting times ID'
        call handle_err(istatus)
     endif
    if (thermo .eq. 1) then 
     istatus = nf_inq_varid(idncdd,'PP',idPP)
     if (istatus .ne. nf_noerr) then
        print*,'Error getting PP ID'
        call handle_err(istatus)
     endif
     istatus = nf_inq_varid(idncdd,'rhoa',idrhoa)
     if (istatus .ne. nf_noerr) then
        print*,'Error getting rhoa ID'
        call handle_err(istatus)
     endif
     istatus = nf_inq_varid(idncdd,'thetapp',idthetapp)
     if (istatus .ne. nf_noerr) then
        print*,'Error getting thetapp ID'
        call handle_err(istatus)
     endif
     istatus = nf_inq_varid(idncdd,'sp',idsp)
     if (istatus .ne. nf_noerr) then
        print*,'Error getting sp ID'
        call handle_err(istatus)
     endif
     istatus = nf_inq_varid(idncdd,'qvpp',idqvpp)
     if (istatus .ne. nf_noerr) then
        print*,'Error getting qvpp ID'
        call handle_err(istatus)
     endif
     istatus = nf_inq_varid(idncdd,'rmth',idrmth)
     if (istatus .ne. nf_noerr) then
        print*,'Error getting rmth ID'
        call handle_err(istatus)
     endif
     istatus = nf_inq_varid(idncdd,'rmqv',idrmqv)
     if (istatus .ne. nf_noerr) then
        print*,'Error getting rmqv ID'
        call handle_err(istatus)
     endif
    endif !thermo
     istatus = nf_inq_varid(idncdd,'IDP',ididp)
     if (istatus .ne. nf_noerr) then
        print*,'Error getting idp ID'
        call handle_err(istatus)
     endif 
     istatus = nf_inq_varid(idncdd,'XDROP',idxd)
     if (istatus .ne. nf_noerr) then
        print*,'Error getting xdrop ID'
        call handle_err(istatus)
     endif 
     istatus = nf_inq_varid(idncdd,'YDROP',idyd)
     if (istatus .ne. nf_noerr) then
        print*,'Error getting ydrop ID'
        call handle_err(istatus)
     endif 
     istatus = nf_inq_varid(idncdd,'ZDROP',idzd)
     if (istatus .ne. nf_noerr) then
        print*,'Error getting zdrop ID'
        call handle_err(istatus)
     endif 
     istatus = nf_inq_varid(idncdd,'UDROP',idudrop)
     if (istatus .ne. nf_noerr) then
        print*,'Error getting udrop ID'
        call handle_err(istatus)
     endif 
     istatus = nf_inq_varid(idncdd,'VDROP',idvdrop)
     if (istatus .ne. nf_noerr) then
        print*,'Error getting vdrop ID'
        call handle_err(istatus)
     endif 
     istatus = nf_inq_varid(idncdd,'WDROP',idwdrop)
     if (istatus .ne. nf_noerr) then
        print*,'Error getting wdrop ID'
        call handle_err(istatus)
     endif 
     istatus = nf_inq_varid(idncdd,'R',idr)
     if (istatus .ne. nf_noerr) then
        print*,'Error getting r ID'
        call handle_err(istatus)
     endif 
     if(thermo.eq.1) then
        istatus = nf_inq_varid(idncdd,'DR3',iddr3)
        if (istatus .ne. nf_noerr) then
           print*,'Error getting dr3 ID'
           call handle_err(istatus)
        endif 
        istatus = nf_inq_varid(idncdd,'R_CCN',idr_ccn)
        if (istatus .ne. nf_noerr) then
                print*,'Error getting r_ccn ID'
                call handle_err(istatus)
        endif
     endif
     ! prep netcdf read
     ncstartd(2) = 1

     nccountd(1) = ndroppe
     nccountd(2) = 1

     ! read in time
     istatus = nf_get_vara_real(idncdd,idtimesd,1,1,time)
     if (istatus .ne. nf_noerr) then
        print*,'Error reading time'
        call handle_err(istatus)
     endif 
     print*, 'in droplet restart time =', time
       if(thermo .eq. 1) then
        istatus = nf_get_vara_double(idncdd,idPP,1,1,PP)
        if (istatus .ne. nf_noerr) then
           print*,'Error reading PP'
           call handle_err(istatus)
        endif
        istatus = nf_get_vara_double(idncdd,idrhoa,1,1,rhoa)
        if (istatus .ne. nf_noerr) then
           print*,'Error reading rhoa'
           call handle_err(istatus)
        endif
        istatus = nf_get_vara_double(idncdd,idthetapp,1,1,thetapp)
        if (istatus .ne. nf_noerr) then
           print*,'Error reading thetapp'
           call handle_err(istatus)
        endif
        istatus = nf_get_vara_double(idncdd,idsp,1,1,sp)
        if (istatus .ne. nf_noerr) then
           print*,'Error reading sp'
           call handle_err(istatus)
        endif
        istatus = nf_get_vara_double(idncdd,idqvpp,1,1,qvpp)
        if (istatus .ne. nf_noerr) then
           print*,'Error reading qvpp'
           call handle_err(istatus)
        endif
        istatus = nf_get_vara_real(idncdd,idrmth,1,1,rmth)
        if (istatus .ne. nf_noerr) then
           print*,'Error reading rmth'
           call handle_err(istatus)
        endif
        istatus = nf_get_vara_real(idncdd,idrmqv,1,1,rmqv)
        if (istatus .ne. nf_noerr) then
           print*,'Error reading rmqv'
           call handle_err(istatus)
        endif
       endif !thermo
     do iread=npe-1,0,-1  ! read one block at a time; read first block last

        ncstartd(1) = iread*ndroppe+1 

        istatus = nf_get_vara_real(idncdd,ididp,ncstartd,nccountd,idpreal)
        if (istatus .ne. nf_noerr) then
           print*,'Error reading idpreal at iread, ididp', iread, ididp
           call handle_err(istatus)
        endif 
        istatus = nf_get_vara_real(idncdd,idxd,ncstartd,nccountd,x)
        if (istatus .ne. nf_noerr) then
           print*,'Error reading x at iread', iread
           call handle_err(istatus)
        endif 
        istatus = nf_get_vara_real(idncdd,idyd,ncstartd,nccountd,y)
        if (istatus .ne. nf_noerr) then
           print*,'Error reading y'
           call handle_err(istatus)
        endif 
        istatus = nf_get_vara_real(idncdd,idzd,ncstartd,nccountd,z)
        if (istatus .ne. nf_noerr) then
           print*,'Error reading z'
           call handle_err(istatus)
        endif 
        istatus = nf_get_vara_real(idncdd,idudrop,ncstartd,nccountd,udrop)
        if (istatus .ne. nf_noerr) then
           print*,'Error reading udrop'
           call handle_err(istatus)
        endif 
        istatus = nf_get_vara_real(idncdd,idvdrop,ncstartd,nccountd,vdrop)
        if (istatus .ne. nf_noerr) then
           print*,'Error reading vdrop'
           call handle_err(istatus)
        endif 
        istatus = nf_get_vara_real(idncdd,idwdrop,ncstartd,nccountd,wdrop)
        if (istatus .ne. nf_noerr) then
           print*,'Error reading wdrop'
           call handle_err(istatus)
        endif 
        istatus = nf_get_vara_double(idncdd,idr,ncstartd,nccountd,r)
        if (istatus .ne. nf_noerr) then
           print*,'Error reading r'
           call handle_err(istatus)
        endif 
        if(thermo .eq. 1)then
          istatus = nf_get_vara_double(idncdd,iddr3,ncstartd,nccountd,dr3)
          if (istatus .ne. nf_noerr) then
             print*,'Error reading dr3'
             call handle_err(istatus)
          endif 
          istatus = nf_get_vara_double(idncdd,idr_ccn,ncstartd,nccountd,r_ccn)
          if (istatus .ne. nf_noerr) then
                  print*,'Error reading r_ccn'
                  call handle_err(istatus)
          endif
        endif
        if (iread.gt.0) then  ! send to processor iread
	   if(thermo .eq. 1) then
           call mpi_send(PP,1,MPI_REAL8,iread,0,MPI_COMM_WORLD,istatus)
           if (istatus.ne.0) print*,'Error sending pp to proc', iread
           call mpi_send(rhoa,1,MPI_REAL8,iread,0,MPI_COMM_WORLD,istatus)
           if (istatus.ne.0) print*,'Error sending rhoa to proc', iread
           call mpi_send(thetapp,1,MPI_REAL8,iread,0,MPI_COMM_WORLD,istatus)
           if (istatus.ne.0) print*,'Error sending thetapp to proc', iread
           call mpi_send(sp,1,MPI_REAL8,iread,0,MPI_COMM_WORLD,istatus)
           if (istatus.ne.0) print*,'Error sending sp to proc', iread
           call mpi_send(qvpp,1,MPI_REAL8,iread,0,MPI_COMM_WORLD,istatus)
           if (istatus.ne.0) print*,'Error sending qvpp to proc', iread
           call mpi_send(rmth,1,MPI_REAL,iread,0,MPI_COMM_WORLD,istatus)
           if (istatus.ne.0) print*,'Error sending rmth to proc', iread
           call mpi_send(rmqv,1,MPI_REAL,iread,0,MPI_COMM_WORLD,istatus)
           if (istatus.ne.0) print*,'Error sending rmqv to proc', iread
	   endif!thermo

           call mpi_send(idpreal,nbuf,MPI_REAL,iread,0,MPI_COMM_WORLD,istatus)
           if (istatus.ne.0) print*,'Error sending idpreal to proc', iread
           call mpi_send(x,nbuf,MPI_REAL,iread,0,MPI_COMM_WORLD,istatus)
           if (istatus.ne.0) print*,'Error sending x to proc', iread
           call mpi_send(y,nbuf,MPI_REAL,iread,0,MPI_COMM_WORLD,istatus)
           if (istatus.ne.0) print*,'Error sending y to proc', iread
           call mpi_send(z,nbuf,MPI_REAL,iread,0,MPI_COMM_WORLD,istatus)
           if (istatus.ne.0) print*,'Error sending z to proc', iread
           call mpi_send(udrop,nbuf,MPI_REAL,iread,0,MPI_COMM_WORLD,istatus)
           if (istatus.ne.0) print*,'Error sending udrop to proc', iread
           call mpi_send(vdrop,nbuf,MPI_REAL,iread,0,MPI_COMM_WORLD,istatus)
           if (istatus.ne.0) print*,'Error sending vdrop to proc', iread
           call mpi_send(wdrop,nbuf,MPI_REAL,iread,0,MPI_COMM_WORLD,istatus)
           if (istatus.ne.0) print*,'Error sending wdrop to proc', iread
           call mpi_send(r,nbuf,MPI_REAL8,iread,0,MPI_COMM_WORLD,istatus)
           if (istatus.ne.0) print*,'Error sending r to proc', iread
           if(thermo .eq. 1)then
              call mpi_send(dr3,nbuf,MPI_REAL8,iread,0,MPI_COMM_WORLD,istatus)
              if (istatus.ne.0) print*,'Error sending dr to proc', iread
              call mpi_send(r_ccn,nbuf,MPI_REAL8,iread,0,MPI_COMM_WORLD,istatus)
              if (istatus.ne.0) print*,'Error sending r_ccn to proc', iread

           endif
        endif !iread 0
     enddo !iread
  endif ! mype .eq 0

  if (mype.gt.0) then ! get block from proc 0
     if(thermo .eq. 1) then
     call mpi_recv(pp,1,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,status,istatus)
     if (istatus.ne.0) print*,mype,'Error receiving pp'
     call mpi_recv(rhoa,1,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,status,istatus)
     if (istatus.ne.0) print*,mype,'Error receiving rhoa'
     call mpi_recv(thetapp,1,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,status,istatus)
     if (istatus.ne.0) print*,mype,'Error receiving thetapp'
     call mpi_recv(sp,1,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,status,istatus)
     if (istatus.ne.0) print*,mype,'Error receiving sp'
     call mpi_recv(qvpp,1,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,status,istatus)
     if (istatus.ne.0) print*,mype,'Error receiving qvpp'
     call mpi_recv(rmth,1,MPI_REAL,0,MPI_ANY_TAG,MPI_COMM_WORLD,status,istatus)
     if (istatus.ne.0) print*,mype,'Error receiving rmth'
     call mpi_recv(rmqv,1,MPI_REAL,0,MPI_ANY_TAG,MPI_COMM_WORLD,status,istatus)
     if (istatus.ne.0) print*,mype,'Error receiving rmqv'
     endif!thermo

     call mpi_recv(idpreal,nbuf,MPI_REAL,0,MPI_ANY_TAG,MPI_COMM_WORLD,status,istatus)
     if (istatus.ne.0) print*,mype,'Error receiving idpreal'
     call mpi_recv(x,nbuf,MPI_REAL,0,MPI_ANY_TAG,MPI_COMM_WORLD,status,istatus)
     if (istatus.ne.0) print*,mype,'Error receiving x'
     call mpi_recv(y,nbuf,MPI_REAL,0,MPI_ANY_TAG,MPI_COMM_WORLD,status,istatus)
     if (istatus.ne.0) print*,mype,'Error receiving y'
     call mpi_recv(z,nbuf,MPI_REAL,0,MPI_ANY_TAG,MPI_COMM_WORLD,status,istatus)
     if (istatus.ne.0) print*,mype,'Error receiving z'
     call mpi_recv(udrop,nbuf,MPI_REAL,0,MPI_ANY_TAG,MPI_COMM_WORLD,status,istatus)
     if (istatus.ne.0) print*,mype,'Error receiving udrop'
     call mpi_recv(vdrop,nbuf,MPI_REAL,0,MPI_ANY_TAG,MPI_COMM_WORLD,status,istatus)
     if (istatus.ne.0) print*,mype,'Error receiving vdrop'
     call mpi_recv(wdrop,nbuf,MPI_REAL,0,MPI_ANY_TAG,MPI_COMM_WORLD,status,istatus)
     if (istatus.ne.0) print*,mype,'Error receiving wdrop'
     call mpi_recv(r,nbuf,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,status,istatus)
     if (istatus.ne.0) print*,mype,'Error receiving r'
     if(thermo.eq.1)then
        call mpi_recv(dr3,nbuf,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,status,istatus)
        if (istatus.ne.0) print*,mype,'Error receiving dr3'
        call mpi_recv(r_ccn,nbuf,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,status,istatus)
        if (istatus.ne.0) print*,mype,'Error receiving r_ccn'
     endif
  endif
!endif!mpi

  idp = int(idpreal)
  ndropreal = idp(ndroppe)
  nddone = idp(ndroppe-1)

  if (mype.eq.0) then
     istatus = nf_close(idncdd)
     if (istatus .ne. nf_noerr) then
        print*,'Error closing drop.in'
        call handle_err(istatus)
     endif 
  endif

  deallocate( idpreal )

end subroutine ncreaddrop




subroutine handle_err(istatus)

  implicit none
  include 'netcdf.inc'
  include 'mpif.h'
  
  integer :: istatus

!if ( MPI == 1 ) then
  integer :: status(MPI_STATUS_SIZE)
!endif!mpi
  integer :: mype
  common/mpi/mype

  if (istatus .ne. nf_noerr) then
     print*, 'mype', mype, 'has following error:', nf_strerror(istatus), istatus
     stop
  endif

end subroutine handle_err          




