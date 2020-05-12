module mic_mod
  implicit none
    ! for mic block
      integer :: ihydro,isolu,colcal,colgr
      real(8) :: ka1,ka2,diffvnd1,diffvnd2
      real(8) :: delt,ks
      integer :: nstop,micdtout,nout,dsdout
      real(8) :: pp,sp,rhoa,thetapp,qvpp,seq
      real(8) :: thetap,qvp,deltaqp
      real    :: cellw,slabw,oneoverh,h
      real(8) :: vol
      real    :: rmth,rmqv
      real(8) :: radmin
      integer, allocatable, dimension(:) :: ixp,iyp,izp
      real, allocatable, dimension(:)    :: dtovertaup
      real, allocatable, dimension(:)    :: uudrop,vvdrop,wwdrop,flowu,flowv,floww,fv
      real, allocatable, dimension(:)    :: yrel !ndroppe
      integer,allocatable, dimension(:)  :: idpfol
      real, allocatable, dimension(:)    :: xfol,yfol,zfol
      real(8),allocatable,dimension(:)   :: rfol !ndrop/everynd
      real, allocatable, dimension(:,:)  :: infocols, infocolr!(ndroppe,16)
      integer, allocatable, dimension(:) :: head !ncell+m*m
      real, allocatable, dimension(:)    :: maxutemp !npe
      integer, allocatable, dimension(:,:) :: dcids !(5,2)
      integer, allocatable, dimension(:) :: dsd,dsdtot !nbins dsdtot only in mype0
      integer, allocatable, dimension(:) :: dsd_log2,dsdtot_log2 !nbins dsdtot only in mype0
      integer, allocatable, dimension(:) :: idp
      real, allocatable, dimension(:)    :: x,y,z,udrop,vdrop,wdrop
      real(8),allocatable,dimension(:)   :: r,r_ccn,dr3
      contains
      real function log2(rad)
      real(8) :: rad
        log2 = log(rad)/log(2.d0)
      end function log2

end module mic_mod


