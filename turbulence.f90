!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! regional model module - turbulence
! usage: include all the subroutines
! and functions to calculate eddy diffusivity
! and viscousity 
! Yun Zhang 05/01/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module turbulence
  use constant_parameter
  implicit none

! not finished assume Kb is zero
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate_eddy_viscousity_diffusivity
! usage: use turbulence model to calculate
! eddy diffusivty and viscousity for momentum
! and transport equations
! Yun Zhang 05/02/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calculate_eddy_viscousity_diffusivity(nu_t,K_t)
  real(dp),dimension(NLAT,NLONG,NVERT),intent(inout):: nu_t,K_t
  nu_t=0.0_dp
  K_t=0.0_dp
end subroutine calculate_eddy_viscousity_diffusivity



end module turbulence