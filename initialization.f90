!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! regional model module - initialization
! usage: initialize all the variables
! Yun Zhang 04/29/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module initialization
  use constant_parameter
  use basic_state
  implicit none

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize_lat_long
! usage: initialize all lat and long for cell center and flux face
! Yun Zhang @Stanford
! 04/25/2015
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine initialize_lat_long_coriolis(lat_c,long_c,lat_uface,long_uface,lat_vface,long_vface,f_c)
  integer::i,j,k
  real(dp),dimension(NLAT,NLONG),intent(inout):: lat_c, long_c,f_c
  real(dp),dimension(NLAT+1,NLONG),intent(inout):: lat_vface, long_vface
  real(dp),dimension(NLAT,NLONG+1),intent(inout):: lat_uface, long_uface
  ! calculate lat and long tidu for all cell
  do i=1,NLAT
    do j=1,NLONG
      lat_c(i,j)=lat_0+(i-1)*dphi
      long_c(i,j)=long_0+(j-1)*dlamda_e
    enddo
  enddo

  ! coriolis efficients
  if(coriolis==1) then
    f_c=2*Omega*sin(lat_c)
  else
    f_c=0.0_dp
  endif

  ! uface
  lat_uface(:,1:NLONG)=lat_c
  lat_uface(:,NLONG+1)=lat_c(:,NLONG)
  long_uface(:,1:NLONG)=long_c-0.5*dlamda_e
  long_uface(:,NLONG+1)=long_c(:,NLONG)+0.5*dlamda_e

  ! vface
  lat_vface(1:NLAT,:)=lat_c-0.5*dphi
  lat_vface(NLAT+1,:)=lat_c(NLAT,:)+0.5*dphi
  long_vface(1:NLAT,:)=long_c
  long_vface(NLAT+1,:)=long_c(NLAT,:)

end subroutine initialize_lat_long_coriolis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize_vertical_sigima
! usage: calculate and initialize sigma grid for test column
! Yun Zhang @Stanford
! 04/25/2015
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine initialize_vertical_sigima(zbot_test,pa_test,sigma_bot,dsigma)
  integer::i,j,k
  real(dp), dimension(0:NVERT),intent(inout):: zbot_test, pa_test, sigma_bot ! from layer 0
  real(dp), dimension(NVERT),intent(inout):: dsigma
  do i=0,NVERT
    zbot_test(i)=z_surf+(ztop_test-z_surf)*(1-DBLE(i)/NVERT)
  enddo

  ! calculate pa_test
  pa_test=standard_atmophere_interp(zbot_test,NVERT+1,'z','p')
  pa_test(0)=Pa_top
  ! calculate sigma_bot
  sigma_bot=(pa_test-Pa_top)/(pa_test(NVERT)-Pa_top)
  dsigma=sigma_bot(1:NVERT)-sigma_bot(0:(NVERT-1))
end subroutine initialize_vertical_sigima

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize_air_pressure_surf
! usage: initialize the air pressure at surface for each cell
! Yun Zhang @stanford
! 05/01/2015
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine initialize_air_pressure_surf(Pa_surf)
  real(dp),dimension(NLAT,NLONG),intent(inout)::Pa_surf
  integer::i,j
  real(dp)::tmp
  do i=1,NLAT
    do j=1,NLONG
      tmp=((Re/1000000*cos((lat_0+(i-1)*dphi+phi_center)/2)*(long_0+(j-1)*dlamda_e-lamda_center))**2)/2
      tmp=tmp+((Re/1000000*(lat_0+(i-1)*dphi-phi_center))**2)/2
      Pa_surf(i,j)=Pa_base+dpa_peak*exp(-tmp)
    enddo
  enddo 
end subroutine initialize_air_pressure_surf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize_pressure_system
! usage: initialize the column pressure, air pressure
! at the bottom and center for each layer at each cell
! Yun Zhang @stanford
! 05/01/2015
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine initialize_pressure_system(pi_c_old,pi_c_new,pi_c_tminus1,Pa_bot,P_bot,P_c,Pa_c,Pa_surf,sigma_bot)
  real(dp),dimension(NLAT,NLONG),intent(inout)::pi_c_old,pi_c_new,pi_c_tminus1
  real(dp),dimension(NLAT,NLONG,0:NVERT),intent(inout)::Pa_bot,P_bot
  real(dp),dimension(0:NVERT),intent(in)::sigma_bot
  real(dp),dimension(NLAT,NLONG,NVERT),intent(inout)::Pa_c,P_c
  real(dp),dimension(NLAT,NLONG),intent(in)::Pa_surf
  integer::i,j

  pi_c_old=Pa_surf-Pa_top
  pi_c_tminus1=pi_c_old
  ! provide approximate value for pi_timus1 at boundaries
  pi_c_tminus1(1,:)=pi_c_old(1,:)+(pi_c_old(1,:)-pi_c_old(2,:))
  pi_c_tminus1(NLAT,:)=pi_c_old(NLAT,:)+(pi_c_old(NLAT,:)-pi_c_old(NLAT-1,:))
  pi_c_tminus1(:,1)=pi_c_old(:,1)+(pi_c_old(:,1)-pi_c_old(:,2))
  pi_c_tminus1(:,NLONG)=pi_c_old(:,NLONG)+(pi_c_old(:,NLONG)-pi_c_old(:,NLONG-1))
  pi_c_new=pi_c_old

  ! calculate pressure at different vertical layer 
  ! layer bottom
  do i=1,NLAT
    do j=1,NLONG
      Pa_bot(i,j,:)=sigma_bot*pi_c_old(i,j)+Pa_top
      P_bot(i,j,:)=(Pa_bot(i,j,:)/1000)**k_therm
    enddo 
  enddo

  ! layer center
  do i=1,NLAT
    do j=1,NLONG
      P_c(i,j,:)=1/(1+k_therm)*(P_bot(i,j,1:NVERT)*Pa_bot(i,j,1:NVERT)-Pa_bot(i,j,0:(NVERT-1))*P_bot(i,j,0:(NVERT-1)))
      P_c(i,j,:)=P_c(i,j,:)/(Pa_bot(i,j,1:NVERT)-Pa_bot(i,j,0:NVERT-1))
      Pa_c(i,j,:)=1000*(P_c(i,j,:)**(1/k_therm))
    enddo
  enddo  

end subroutine initialize_pressure_system

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize_watermass_mixingratio
! usage: initialize qv for each cell each layer
! Yun Zhang 05/02/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine initialize_watermass_mixingratio(qv_c,qv_c_new,Qqv)
  real(dp),dimension(NLAT,NLONG,NVERT),intent(inout)::qv_c
  real(dp),dimension(NLAT,NLONG,NVERT),intent(inout)::qv_c_new, Qqv
  qv_c=0.0_dp
  qv_c_new=qv_c
  Qqv=0.0_dp
end subroutine initialize_watermass_mixingratio

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize_temperature_system
! usage: first initial the southwest corner of the domain
! then set all the value for all cells and layers
! Yun Zhang @stanford
! 05/01/2015
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine initialize_temperature_system(Temp_c,Pa_c,qv_c,PVT_c_old,PVT_c_new,P_c,P_bot,Qheat)
  real(dp),dimension(NLAT,NLONG,NVERT),intent(inout)::Temp_c,PVT_c_old,PVT_c_new,Qheat
  real(dp),dimension(NLAT,NLONG,NVERT),intent(in)::qv_c,Pa_c,P_c
  real(dp),dimension(NLAT,NLONG,NVERT),intent(in)::P_bot
  real(dp)::tmp
  integer::i,j,k
  ! first choose one corner use the southwest cornder (1,1,:)
  Temp_c(1,1,:)=standard_atmophere_interp(Pa_c(1,1,:),NVERT,'p','T') 

  ! set all other cells temperature and potential virtual temperature
  do i=1,NLAT
    do j=1,NLONG
      Temp_c(i,j,:)=Temp_c(1,1,:)
    enddo
  enddo

  ! set potential temperature at cell center
  PVT_c_old=Temp_c*(1+0.608*qv_c)*((1000/Pa_c)**k_therm)
  PVT_c_new=PVT_c_old
  Qheat=0.0_dp
end subroutine initialize_temperature_system

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize_gas_concentration
! usage: set the initial value for passive gas concentration
! Yun Zhang @stanford
! 05/01/2015
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine initialize_gas_concentration(gas_c_old,gas_c_new,Pa_c,P_c,P_bot,lat_c,long_c,Qgas)
  real(dp),dimension(NLAT,NLONG,NVERT),intent(inout)::gas_c_old,gas_c_new,Qgas
  real(dp),dimension(NLAT,NLONG,NVERT),intent(in)::Pa_c,P_c
  real(dp),dimension(NLAT,NLONG,NVERT),intent(in)::P_bot
  real(dp),dimension(NLAT,NLONG),intent(in)::lat_c,long_c
  real(dp)::tmp
  integer::i,j,k
  ! set all other cells temperature and potential virtual temperature
  do i=1,NLAT
    do j=1,NLONG
      if(j>18 .and. j<22 .and. i>18 .and. i<22) then
        gas_c_old(i,j,:)=1.0_dp
      else
        gas_c_old(i,j,:)=0.0_dp !Temp_c(1,1,:)
      endif
    enddo
  enddo

  ! set potential temperature at cell center
  gas_c_new=gas_c_old
  Qgas=0.0_dp
end subroutine initialize_gas_concentration

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize_velocity_field
! usage: initialize field velocity field 
! Yun Zhang @stanford
! 05/01/2015
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine initialize_velocity_field(u,u_new,v,v_new,w_sigma_old,w_sigma_new,lat_uface,lat_vface,long_uface,long_vface)
  real(dp), dimension(NLAT,NLONG+1,NVERT),intent(inout):: u,u_new
  real(dp), dimension(NLAT+1,NLONG,NVERT),intent(inout):: v,v_new
  real(dp), dimension(NLAT,NLONG,0:NVERT),intent(inout):: w_sigma_old, w_sigma_new
  real(dp), dimension(NLAT,NLONG+1),intent(in)::lat_uface,long_uface
  real(dp), dimension(NLAT+1,NLONG),intent(in)::long_vface,lat_vface
  u=0.0_dp
  v=0.0_dp
  u_new=u
  v_new=v
  w_sigma_old=0.0_dp
  w_sigma_new=0.0_dp
end subroutine initialize_velocity_field

end module initialization