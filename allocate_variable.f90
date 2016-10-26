!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! regional model module - allocate_variable
! usage: allocate all variables in the program
! Yun Zhang 04/24/2015
! @stanford
! unit clarification
! pressure: hpa
! temperature: K
! latitude,longitude: degree
! density: kg/m^3
! altitude, length: m
! gravity: m/s/s
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module allocate_variable
  use constant_parameter
  implicit none
  
  ! Initialize parameter 
  ! lat and long for each cell
  real(dp),dimension(NLAT,NLONG):: lat_c, long_c ! latitude and longitude
  real(dp),dimension(NLAT,NLONG):: f_c ! coriolis coefficient
  real(dp),dimension(NLAT+1,NLONG):: lat_vface, long_vface ! latitude at v flux face
  real(dp),dimension(NLAT,NLONG+1):: lat_uface, long_uface ! longtitude at u flux face

  ! part 1)
  ! test column to calculate dsigma
  real(dp), dimension(0:NVERT):: zbot_test ! layer bottom elevation at test column
  real(dp), dimension(0:NVERT):: pa_test ! layer bottom pressure at test column
  real(dp), dimension(0:NVERT):: sigma_bot ! sigma value at each column, constant for all time step
  real(dp), dimension(NVERT):: dsigma ! sigma(k+1)-sigma(k) 

  ! part 2) pressure setup
  ! initialize surface pressure
  real(dp), dimension(NLAT,NLONG):: Pa_surf ! surface pressure
  real(dp), dimension(NLAT,NLONG):: pi_c_new, pi_c, pi_c_tminus1, pi_c_tmp! column pressure at t+1 t t-1
  real(dp), dimension(NLAT+1,NLONG):: pi_vface ! column pressure at v face
  real(dp), dimension(NLAT,NLONG+1):: pi_uface ! column pressure at u face


  ! calculate all pressure at different layers for different cells
  real(dp), dimension(NLAT,NLONG,0:NVERT):: Pa_bot, Pa_bot_new ! layer bottom air pressure at each layer at t t+1
  real(dp), dimension(NLAT,NLONG,0:NVERT):: P_bot, P_bot_new ! layer bottom P at each layer at t t+!
  real(dp), dimension(NLAT,NLONG,NVERT):: Pa_c, Pa_c_new ! cell center air pressure at t t+1
  real(dp), dimension(NLAT,NLONG,NVERT):: P_c, P_c_new  ! cell center P at t t+1

  ! part 3) initialize temperature, humidity and gas
  real(dp), dimension(NLAT,NLONG,NVERT):: Temp_c, Temp_c_new ! cell center Temperature at t t+1
  real(dp), dimension(NLAT,NLONG,NVERT):: PVT_c, PVT_c_new, PVT_c_tmp! cell center potential virtual temperature at t t+1
  real(dp), dimension(NLAT,NLONG,0:NVERT):: PVT_bot, qv_bot, gas_bot! layer bottom potential virtual temperature and specific humidity
  real(dp), dimension(NLAT,NLONG,NVERT):: qv_c, qv_c_new !cell center specific humidity at t t+1
  real(dp), dimension(NLAT,NLONG,NVERT):: gas_c, gas_c_new !cell center gas concentration at t t+1
  real(dp), dimension(NLAT+1,NLONG,NVERT):: PVT_vface, qv_vface, gas_vface! v face potential virtual temperature and specific humidity
  real(dp), dimension(NLAT,NLONG+1,NVERT):: PVT_uface, qv_uface, gas_uface ! u face potential virtual temperature and specific humidity

  ! part 4) initialize air density
  real(dp), dimension(NLAT,NLONG,NVERT):: rhoa_c, rhoa_c_new ! cell center air density at t t+1

  ! part 5) initialize velocity field
  real(dp), dimension(NLAT,NLONG+1,NVERT):: u, u_new ! u at t t+1
  real(dp), dimension(NLAT+1,NLONG,NVERT):: v, v_new ! v at t t+1
  real(dp), dimension(NLAT,NLONG,0:NVERT):: w_sigma, w_sigma_new ! w_sigma at t t+!
  
  ! part 6) calculate column pressure at flux face

  real(dp), dimension(NLAT,NLONG,NVERT):: fluxsum ! fluxsum(i,j,k) means the sum of flux from layer 1 to k at cell i,j
  real(dp), dimension(NLAT,NLONG+1,NVERT):: flux_uface ! the sum of flux from u direction at each layer
  real(dp), dimension(NLAT+1,NLONG,NVERT):: flux_vface ! the sum of flux from v direction at each layer
  
  ! part 7) turbulence variable
  real(dp),dimension(NLAT,NLONG,NVERT):: nu_t,K_t ! cell center eddy viscousity and diffusity
 
  ! part 8) heat/gas source variable
  real(dp),dimension(NLAT,NLONG,NVERT):: Qheat ! cell center heat source
  real(dp),dimension(NLAT,NLONG,NVERT):: Qgas ! cell center gas source
  real(dp),dimension(NLAT,NLONG,NVERT):: Qqv ! cell center specific source
  
  ! part 9) geopotential
  real(dp), dimension(NLAT,NLONG,NVERT):: geopot_c,geopot_c_tminus1 ! geopotential at t t-1
  real(dp), dimension(NLAT,NLONG,0:NVERT)::geopot_bot ! geopotential at layer bottom

end module allocate_variable
