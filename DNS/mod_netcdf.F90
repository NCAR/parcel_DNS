module netcdf_mod
   implicit none
	integer, parameter :: nrstd = 10
	integer, parameter :: nrst = 10
	!for netcdfrstd block
        integer, dimension(nrstd) ::idncd
        integer, dimension(2) :: ncstartd,nccountd
        integer :: idtimesd,ididp,idxd,idyd,idzd,idr,iddr3,idr_ccn,idudrop,idvdrop,idwdrop,idradmin
        integer :: idPP,idrhoa,idthetapp,idsp,idqvpp,idrmth,idrmqv
	!for netcdfrst block
        integer, dimension(nrst) :: idnck
        integer, dimension(5) :: ncstartk,nccountk
        integer :: idtimesk,idzxk,idzyk,idzzk,idttk,idqvk
	!for netcdfrsp block
        integer :: idnc
        integer,dimension(4) :: ncstart,nccount
        integer ::  idtimes,idu,idv,idw,idth,idqv,idzx,idzy,idzz
end module netcdf_mod
