!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! regional model module - constant_parameter
! usage: define all the paramter and
! precision setup
! Yun Zhang 04/24/2015
! @stanford
! unit clarification
! pressure: hpa
! temperature: K
! latitude,longitude: rad
! density: kg/m^3
! altitude, length: m
! gravity: m/s/s
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module constant_parameter
  use iso_fortran_env, only: real32, real64, int64
  ! basic setup
  integer, parameter:: sp=real32
  integer, parameter:: dp=real64
  integer, parameter:: li=int64
  real(dp), parameter:: pi = 4.0_dp*atan(1.0_dp) ! pi
  real(dp), parameter:: k_therm=0.286_dp 
  real(dp), parameter:: R_prime=2.8704_dp !m^3 hpa kg^-1 K^-1
  real(dp), parameter:: Cp_d=1004.67 ! J kg^-1 K^-1
  real(dp), parameter:: Omega=7.2921e-5 ! rad/s rotation frequency for earth
  real(dp), parameter:: Re=6371000.0_dp ! Earth radius use unit in million of meters
  
  ! switch
  integer, parameter:: PVTmodel=1 ! whether the air is dry (0) or not(1)    
  integer, parameter:: PTVbound=0 ! whether set exact boundary condition (1) or not (0)  
  integer, parameter:: heatsource=0 ! whether there is heat source or sink in each cell
  integer, parameter:: gasmodel=1 ! whether to simulation passive gas transport
  integer, parameter:: gasbound=0 ! whether to use specified boundary condition
  integer, parameter:: gassource=0 ! whether there is a gas source within the domain
  integer, parameter:: qvmodel=1 ! whether the air is dry (0) or not(1)    
  integer, parameter:: qvbound=0 ! whether set exact boundary condition (1) or not (0)    
  integer, parameter:: qvsource=0 ! whether there is a specific humidity source within the domain
  integer, parameter:: turbmodel=0 ! whether to turn on turbulence
  integer, parameter:: outputswitch=1 ! whether output results or not
  integer, parameter:: periodicBC=0 ! whether to use periodic boundary condition 
  integer, parameter:: matsuno=1 ! whether to use matsuno scheme
  integer, parameter:: coriolis=1 ! whether to consider coriolis effects
 
  ! numerical parameter
  integer, parameter:: Ndt=500! the numbers of time steps
  integer, parameter:: Nout=50 ! how often to output results
  integer, parameter:: Nrep=1 ! how often to report progress
  real(dp), parameter:: dt=5.0_dp! time step0

  ! grid parameter
  real(dp), parameter:: lat_0=-1/180*pi ! southwest corner latitude
  real(dp), parameter:: long_0=-1/180*pi ! southwest corner longitude
  
  integer, parameter:: NLAT=40! the numbers of latitude cells
  integer, parameter:: NLONG=40 ! the numbers of longitude cells
  integer, parameter:: NVERT=15 ! the number of vertical layers

  real(dp), parameter:: dlamda_e=0.05/180*pi ! differential longitude to represent cell size 
  real(dp), parameter:: dphi=0.05/180*pi ! differential latitude to represent cell size
  real(dp), parameter:: phi_center=lat_0+((NLAT-1)/2+DBLE(mod((NLAT-1),2))/2)*dphi ! the latitude center in the domain
  real(dp), parameter:: lamda_center=long_0+((NLONG-1)/2+DBLE(mod((NLONG-1),2))/2)*dlamda_e ! the longitude center in the domain
  real(dp), parameter:: z_surf=0.0_dp ! surface elevation

  ! Phys parameter
  real(dp), parameter:: Pa_top=250_dp  ! top pressure constant
  real(dp), parameter:: Pa_below=265_dp, rhoa_below=0.414_dp, g_below=9.7764_dp, z_below=10000_dp ! values from table B.1
  real(dp), parameter:: ztop_test=z_below+(Pa_below-Pa_top)*100/rhoa_below/g_below ! top elevation for test column
  real(dp), parameter:: Pa_base=1000.0_dp, dpa_peak=10.0_dp 

  ! output results: pi_c,Pa_c,PTV_c,u,v,w_sigma,nu_t,K_t,geobot_c, Temp_c, rhoa_c
  ! output file names
  character(*), parameter:: resultfolder="/Users/zyaj/Documents/atomsphere-regional-model/results/"
  character(*), parameter:: pi_file="pi.txt" ! file for column pressure results
  character(*), parameter:: pi_format="(1xf9.4)" 
  integer, parameter:: pi_file_no=1
  character(*), parameter:: Pa_file="Pa.txt" ! file for air pressure results 
  character(*), parameter:: Pa_format="(1xf9.4)"
  integer, parameter:: Pa_file_no=2
  character(*), parameter:: PTV_file="PVT.txt" ! file for potential virtual temperature results 
  integer, parameter:: PVT_file_no=3
  character(*), parameter:: PVT_format="(1xf9.4)"
  character(*), parameter:: u_file="u.txt" ! file for u results 
  integer, parameter:: u_file_no=4
  character(*), parameter:: u_format="(1xf9.6)"
  character(*), parameter:: v_file="v.txt" ! file for v results 
  integer, parameter:: v_file_no=5
  character(*), parameter:: v_format="(1xf9.6)"
  character(*), parameter:: w_file="w.txt" ! file for w_sigma results 
  integer, parameter:: w_file_no=6
  character(*), parameter:: w_format="(1xf9.6)"
  character(*), parameter:: nu_file="nu.txt" ! file for eddy viscousity results 
  integer, parameter:: nu_file_no=7
  character(*), parameter:: nu_t_format="(1xf9.4)"
  character(*), parameter:: K_t_file="K_t.txt" ! file for eddy diffusivity results 
  integer, parameter:: K_t_file_no=8
  character(*), parameter:: K_t_format="(1xf9.4)"
  character(*), parameter:: geopot_file="geopot.txt" ! file for geopotential results
  integer, parameter:: geopot_file_no=9
  character(*), parameter:: geopot_format="(1xf12.4)"
  character(*), parameter:: Temp_file="temp.txt" ! file for temperature results
  integer, parameter:: Temp_file_no=10
  character(*), parameter:: Temp_format="(1xf9.4)"
  character(*), parameter:: rhoa_file="rhoa.txt" ! file for air density results
  integer, parameter:: rhoa_file_no=11
  character(*), parameter:: rhoa_format="(1xf6.4)"
  character(*), parameter:: qv_file="qv.txt" ! file for specific humidity results
  integer, parameter:: qv_file_no=12
  character(*), parameter:: qv_format="(1xf6.4)"
  character(*), parameter:: gas_file="gas.txt" ! file for potential virtual temperature results 
  integer, parameter:: gas_file_no=13
  character(*), parameter:: gas_format="(1xf10.8)"
end module constant_parameter

