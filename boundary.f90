!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! regional model module - boundary
! usage: include all the functions and 
! subroutines to define boundary condition
! Yun Zhang 04/30/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module boundary
  use constant_parameter
  implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! boundary_potential_virtual_temp
! usage: provide the PVT value at domain
! boundary 
! Yun Zhang 04/30/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function boundary_potential_virtual_temp(lat,long,nlayer,t) result(output)
  real(dp),intent(in)::lat,long,t
  integer, intent(in):: nlayer
  real(dp)::output

  ! change boundary condition for differnt problem
  output=1.0_dp

end function boundary_potential_virtual_temp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! boundary_specific_humidity
! usage: provide the qv value at domain
! boundary 
! Yun Zhang 05/02/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function boundary_specific_humidity(lat,long,nlayer,t) result(output)
  real(dp),intent(in)::lat,long,t
  integer, intent(in):: nlayer
  real(dp)::output

  ! change boundary condition for differnt problem
  output=1.0_dp

end function boundary_specific_humidity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! boundary_gas
! usage: provide the gas concentration value at domain
! boundary 
! Yun Zhang 05/02/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function boundary_gas(lat,long,nlayer,t) result(output)
  real(dp),intent(in)::lat,long,t
  integer, intent(in):: nlayer
  real(dp)::output

  ! change boundary condition for differnt problem
  output=1.0_dp

end function boundary_gas

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! boundary_surf_geopotential
! usage: provide the geopotential boundary
! value at the surface layer
! Yun Zhang 04/30/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function boundary_surf_geopotential(lat,long,t) result(output)
  real(dp),intent(in)::lat,long,t
  real(dp)::output

  ! change boundary condition for differnt problem
  output=0.0_dp

end function boundary_surf_geopotential



end module boundary