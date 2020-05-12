module hydro_mod 
!!!variables for hydrodynamic effect
!     integer,parameter :: iprint = 0, idistrad = 0, cellopen = 1 
!	idistrad = 0 disregard the hydro-dynamic effect when distance > 30r; 
!	=1 include hydro-dynamic effect for any distance. cellopen=1,include droplets of neighbouring cells in the linked list, =0, only check droplets within the cell itself.
     integer,parameter :: ireloc = 1 
!	ireloc=1 allow reallocation of droplets after collision
     integer :: nbufright,nbufleft,array_siz
     integer :: udist_iter
     integer :: idinfocols
!     integer :: alloc_err0,alloc_err1,alloc_err2,alloc_err3
     real,parameter :: distRad = 6.d-4
     real :: rmsudrop,rmsvdrop,rmswdrop
     real :: dista,distbx,distby,distbz,distc,distdx,distdy,distdz,maxudist
     real :: distaji,distbxji,distbyji,distbzji,distcji,distdxji,distdyji,distdzji
     real :: distaij,distbxij,distbyij,distbzij,distcij,distdxij,distdyij,distdzij
     real :: udisttemp,vdisttemp,wdisttemp,maxuall,maxuall_old,maxuall_old2
     real :: rxji,ryji,rzji
     real, allocatable, dimension(:) :: ucharact,udist_syn,uudist_syn
!     real, allocatable, dimension(:) :: dtovertaup
     real, allocatable, dimension(:) :: dtovertaupc,flowuc,flowvc,flowwc,udist,vdist,wdist,uudist,vvdist,wwdist
     real, allocatable, dimension(:,:) :: distsendr1, distrecvr2, distsendl2, distrecvl1
!     real, allocatable, dimension(:) :: maxutemp
end module hydro_mod
