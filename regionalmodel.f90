!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! regionalmodel
! usage: generate a basic regional 
! model that solves the equations of
! atompheric dynamics.
! Yun Zhang 04/24/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! MODEL DESCRIPTION
! The purpose of this project is to 
! develop a basic regional- or global-scale 
! model that solves the equations of atmospheric 
! dynamics, except for the water-vapor continuity 
! equation. Diabatic energy sources and sinks and 
! eddy diffusion are ignored. Winds are driven 
! primarily by pressure gradients.

program regionalmodel
  use iso_fortran_env
  use allocate_variable
  use phys
  use output
  implicit none
  integer:: nt
  real(dp):: tstart,tend

  ! initialize lat and long
  call initialize_lat_long_coriolis(lat_c,long_c,lat_uface,long_uface,lat_vface,long_vface,f_c)

  ! part 1)
  ! calculate dsigma 
  ! assume mean surface altitude 0m
  ! from table B.1
  call initialize_vertical_sigima(zbot_test,pa_test,sigma_bot,dsigma)

  ! part 2)
  ! calculate initial surface pressure
  call initialize_air_pressure_surf(Pa_surf)

  ! calculate initial column pressure
  call initialize_pressure_system(pi_c,pi_c_new,pi_c_tminus1,Pa_bot,P_bot,P_c,Pa_c,Pa_surf,sigma_bot)

  ! part 3)
  ! initialize water mass mixing ratio 
  call initialize_watermass_mixingratio(qv_c,qv_c_new,Qqv)  

  ! initialize temperature for each cell
  call initialize_temperature_system(Temp_c,Pa_c,qv_c,PVT_c,PVT_c_new,P_c,P_bot,Qheat)
   
  ! initialize gas concentration for each cell
  if(gasmodel==1) call initialize_gas_concentration(gas_c,gas_c_new,Pa_c,P_c,P_bot,lat_c,long_c,Qgas)

  ! part 4) initialize air density for each cell center
  call calculate_air_density(rhoa_c,Pa_c,qv_c,Temp_c)

  ! part 5) initialize velocity field, zero everywhere 
  call initialize_velocity_field(u,u_new,v,v_new,w_sigma,w_sigma_new,lat_uface,lat_vface,&
    long_uface,long_vface)

  ! part 6) initialize turbulence term
  call calculate_eddy_viscousity_diffusivity(nu_t,K_t)

  ! part 7) initialize heat source
  if(heatsource==1) call calculate_heat_source(Pa_c,Qheat,lat_c,long_c,0) 
    
  ! initialize specific humidity source
  if(qvsource==1) call calculate_qv_source(Pa_c,Qqv,lat_c,long_c,nt)
  
  ! gas source
  if(gasmodel==1) call calculate_gas_source(Pa_c,Qgas,lat_c,long_c,0)

  ! part 8) initialize geopotential
  call calculate_geopotential(geopot_bot,geopot_c,geopot_c_tminus1,P_c,P_bot,PVT_c,lat_c,long_c,0)

  print *, "The initialization process is finished!"

  ! store cpu time 
  call cpu_time(tstart)

  ! main loop for physical function
  do nt=1,Ndt

    ! calculate column pressure at flux face
    call calculate_fluxface_column_pressure(pi_c,pi_uface,pi_vface)
    
    ! calcualte flux sum for each cell each layer
    call calculate_fluxsum(fluxsum,flux_uface,flux_vface,u,v,pi_uface,pi_vface,lat_c,lat_vface,dsigma)
    
    ! update column pressure for each cell center
    call calculate_column_pressure(pi_c,pi_c_new,lat_c,fluxsum)
    
    ! update w_sigma
    call calculate_w_sigma(w_sigma_new,pi_c,pi_c_new,lat_c,fluxsum,sigma_bot,dsigma)  
        
    if(qvmodel==1) then
      ! calculate qv at flux face
      call calculate_fluxface_qv(qv_c,qv_uface,qv_vface,qv_bot,P_c,P_bot,&
        lat_uface,long_uface,lat_vface,long_vface,nt)

      ! update specific humidity
      call calculate_scalar_field(qv_c, qv_c_new, qv_bot, &
        pi_c,pi_c_new,qv_uface, qv_vface,flux_uface,flux_vface,lat_c, rhoa_c,&
        dsigma, w_sigma_new,K_t,Pa_c,Qqv)
    endif

    if(PVTmodel==1) then
      ! calculate PVT at flux face
      call calculate_fluxface_PVT(PVT_c,PVT_uface,PVT_vface,PVT_bot,&
        P_c,P_bot,lat_uface,long_uface,lat_vface,long_vface,nt)    
      ! update potential temperature
      call calculate_scalar_field(PVT_c, PVT_c_new, PVT_bot, &
        pi_c,pi_c_new,PVT_uface, PVT_vface,flux_uface,flux_vface,lat_c, rhoa_c,&
        dsigma, w_sigma_new,K_t,Pa_c,Qheat)
    endif
 
    if(gasmodel==1) then
      ! calculate gas at flux face
      call calculate_fluxface_gas(gas_c,gas_uface,gas_vface,gas_bot,&
        P_c,P_bot,lat_uface,long_uface,lat_vface,long_vface,nt)    
      ! update potential temperature
      call calculate_scalar_field(gas_c, gas_c_new, gas_bot, &
        pi_c,pi_c_new,gas_uface, gas_vface,flux_uface,flux_vface,lat_c, rhoa_c,&
        dsigma, w_sigma_new,K_t,Pa_c,Qgas)
    endif

    ! update velocity field
    call calculate_velocity_field(pi_c,pi_c_new,pi_c,pi_c_tminus1,lat_c, lat_vface,&
      u,u_new,u,v,v_new,v,flux_uface,flux_vface,dsigma,&
      sigma_bot,w_sigma_new,f_c,geopot_c,geopot_c_tminus1,&
      PVT_c,P_bot,P_c,nu_t,rhoa_c)

    ! update air pressure
    call calculate_pressure_field(pi_c_new,Pa_c_new,Pa_bot_new,&
      P_c_new,P_bot_new,sigma_bot)

    ! update air temperature
    call calculate_air_temperature(PVT_c_new,Temp_c_new,qv_c_new,Pa_c_new)
    
    ! update air density
    call calculate_air_density(rhoa_c_new,Pa_c_new,qv_c_new,Temp_c_new)
    
    ! update eddy diffusivity and viscousity
    if(turbmodel==1) then
      call calculate_eddy_viscousity_diffusivity(nu_t,K_t)
    endif

    ! update geopotential
    call calculate_geopotential(geopot_bot,geopot_c,geopot_c_tminus1,P_c_new,P_bot_new,&
      PVT_c_new,lat_c,long_c,nt)

    ! for Matsuno method
    if(matsuno==1) then
      ! calculate column pressure at flux face
      call calculate_fluxface_column_pressure(pi_c_new,pi_uface,pi_vface) 
      ! calcualte flux sum for each cell each layer
      call calculate_fluxsum(fluxsum,flux_uface,flux_vface,u_new,v_new,pi_uface,pi_vface,lat_c,lat_vface,dsigma) 
      ! store the first calculation of pi for velocity calculation
      pi_c_tmp=pi_c_new
      ! update column pressure for each cell center
      call calculate_column_pressure(pi_c,pi_c_new,lat_c,fluxsum)  
    
      ! update w_sigma
      call calculate_w_sigma(w_sigma_new,pi_c,pi_c_new,lat_c,fluxsum,sigma_bot,dsigma)  

      if(qvmodel==1) then
       ! calculate qv at flux face
       call calculate_fluxface_qv(qv_c_new,qv_uface,qv_vface,qv_bot,P_c_new,P_bot_new,&
        lat_uface,long_uface,lat_vface,long_vface,nt)
       ! update specific humidity
       call calculate_scalar_field(qv_c, qv_c_new, qv_bot, &
        pi_c,pi_c_new,qv_uface, qv_vface,flux_uface,flux_vface,lat_c, rhoa_c_new,&
        dsigma, w_sigma_new,K_t,Pa_c_new,Qqv)
      endif 

      if(PVTmodel==1) then
        ! calculate PVT at flux face
        call calculate_fluxface_PVT(PVT_c_new,PVT_uface,PVT_vface,PVT_bot,&
          P_c_new,P_bot_new,lat_uface,long_uface,lat_vface,long_vface,nt)    
        ! store the first calculation of pVT for velocity calculation
        PVT_c_tmp=PVT_c_new
        ! update potential temperature
        call calculate_scalar_field(PVT_c, PVT_c_new, PVT_bot, &
          pi_c,pi_c_new,PVT_uface, PVT_vface,flux_uface,flux_vface,lat_c, rhoa_c_new,&
          dsigma, w_sigma_new,K_t,Pa_c_new,Qheat)
      endif  

      if(gasmodel==1) then
        ! calculate PVT at flux face
        call calculate_fluxface_gas(gas_c_new,gas_uface,gas_vface,gas_bot,&
          P_c_new,P_bot_new,lat_uface,long_uface,lat_vface,long_vface,nt)    
        ! update potential temperature
        call calculate_scalar_field(gas_c, gas_c_new, gas_bot, &
          pi_c,pi_c_new,gas_uface, gas_vface,flux_uface,flux_vface,lat_c, rhoa_c_new,&
          dsigma, w_sigma_new,K_t,Pa_c_new,Qgas)
      endif 

      ! update velocity field
      call calculate_velocity_field(pi_c,pi_c_new,pi_c_tmp,pi_c_tminus1,lat_c, lat_vface,&
        u,u_new,u_new,v,v_new,v_new,flux_uface,flux_vface,dsigma,&
        sigma_bot,w_sigma_new,f_c,geopot_c,geopot_c_tminus1,&
        PVT_c_tmp,P_bot_new,P_c_new,nu_t,rhoa_c_new)

      ! update air pressure
      call calculate_pressure_field(pi_c_new,Pa_c_new,Pa_bot_new,&
        P_c_new,P_bot_new,sigma_bot)
      
      ! update air temperature
      call calculate_air_temperature(PVT_c_new,Temp_c_new,qv_c_new,Pa_c_new)
    
      ! update eddy diffusivity and viscousity
      if(turbmodel==1) then
        call calculate_eddy_viscousity_diffusivity(nu_t,K_t)
      endif   

      ! update geopotential
      call calculate_geopotential(geopot_bot,geopot_c,geopot_c_tminus1,P_c_new,P_bot_new,&
        PVT_c_new,lat_c,long_c,nt)

    endif

    ! update heat source
    if(heatsource==1) call calculate_heat_source(Pa_c_new,Qheat,lat_c,long_c,nt)
 
    ! update specific humidity source
    if(qvsource==1) call calculate_qv_source(Pa_c,Qqv,lat_c,long_c,nt)

    ! update gas source
    if(gasmodel==1) call calculate_gas_source(Pa_c,Qgas,lat_c,long_c,nt)

    ! replace all variable with the new variables
    ! column pressure
    pi_c_tminus1=pi_c
    pi_c=pi_c_new
    ! pressure variable
    Pa_bot=Pa_bot_new
    P_bot=P_bot_new
    Pa_c=Pa_c_new
    Pa_bot=Pa_bot_new
    ! temperature
    Temp_c=Temp_c_new
    PVT_c=PVT_c_new
    gas_c=gas_c_new
    ! specific humidity
    qv_c=qv_c_new
    ! air density
    rhoa_c=rhoa_c_new
    ! velocity field
    w_sigma=w_sigma_new
    u=u_new
    v=v_new

    ! output results: pi_c,Pa_c,PTV_c,u,v,w_sigma,nu_t,K_t,geobot_c, Temp_c, rha_c, qv
    ! at first time step open all file
    if(outputswitch==1) then
      if(nt==1) call output_open_txt_files()
      if(nt==1 .or. mod(nt,Nout)==0) then
        call output_all_variables()
      endif
      if(nt==Ndt) call output_close_files() 
    endif
  enddo
  
  ! calculate runtime 
  call cpu_time(tend)
  print *, 'The program has been finished. The runtime is ',(tend-tstart),'sec'
 end program regionalmodel










