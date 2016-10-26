!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! regional model module - source
! usage: include all the functions and 
! subroutines to define source and sinks
! for heat, momentum and scalar
! Yun Zhang 04/30/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module source
  use constant_parameter
  implicit none
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! heat_source
! usage: provide the value of heat sinks 
! and source to calculate potential virtual
! temperature 
! Yun Zhang 04/30/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function heat_source(lat,long,nlayer,t) result(output)
  real(dp),intent(in)::lat,long,t
  integer, intent(in):: nlayer
  real(dp)::output   
  output=0.0_dp

end function heat_source

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! gas_source
! usage: provide the value of gas sinks 
! and source to calculate gas concentration
! Yun Zhang 04/30/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function qv_source(lat,long,nlayer,t) result(output)
  real(dp),intent(in)::lat,long,t
  integer, intent(in):: nlayer
  real(dp)::output   
  output=0.0_dp

end function qv_source

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! gas_source
! usage: provide the value of gas sinks 
! and source to calculate gas concentration
! Yun Zhang 04/30/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function gas_source(lat,long,nlayer,t) result(output)
  real(dp),intent(in)::lat,long,t
  integer, intent(in):: nlayer
  real(dp)::output   
  output=0.0_dp

end function gas_source


end module source