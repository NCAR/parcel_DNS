module thermo_mod 
!  real, dimension(n1d,n3d,n2dp) :: theta, qv 
  real(8) :: temp,exner,qvs,tau,source,lwc
!  real    :: lwc
  real(8)    :: sumrp
  real(8) :: temp0, esat0, factor,  sup, esat
  real(8) :: dp,  cdp, aa11, aa22
end module thermo_mod
