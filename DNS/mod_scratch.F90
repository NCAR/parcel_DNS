module scratch_mod

!  logical, save :: scratch_alloue = .false.
  complex, allocatable, dimension(:,:,:) :: zkt
  complex, allocatable, dimension(:,:,:) :: zk1,zk2 
  !allocate( zkt(iktx,iktz,iktyp), zk1(iktx,ikty,iktzp), zk2(iktx,ikty,iktzp) )

end module scratch_mod
