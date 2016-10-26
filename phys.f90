!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! regional model module - phys
! usage: include all the subroutines
! and functions to calculate physical equations
! Yun Zhang 04/29/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module phys
  use constant_parameter
  use boundary
  use turbulence
  use source
  use initialization
  implicit none
  
contains 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate_air_density
! usage: calculate air density for each cell
! each layer
! Yun Zhang 04/30/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calculate_air_density(rhoa_c,Pa_c,qv_c,Temp_c)
  real(dp),dimension(NLAT,NLONG,NVERT),intent(in)::Pa_c,qv_c,Temp_c
  real(dp),dimension(NLAT,NLONG,NVERT),intent(inout)::rhoa_c
  rhoa_c=Pa_c/R_prime/(1+0.608*qv_c)/Temp_c
end subroutine calculate_air_density

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate_mass_weighed_P
! usage: calculate P_c P_bot
! each layer
! Yun Zhang 04/30/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calculate_mass_weighed_P(Pa_bot,P_bot,P_c)
  real(dp),dimension(NLAT,NLONG,0:NVERT),intent(in)::Pa_bot
  real(dp),dimension(NLAT,NLONG,0:NVERT),intent(inout)::P_bot
  real(dp),dimension(NLAT,NLONG,NVERT),intent(inout)::P_c
  integer::i,j,k
  do i=1,NLAT
    do j=1,NLONG
      P_bot(i,j,:)=(Pa_bot(i,j,:)/1000)**k_therm
      P_c(i,j,:)=1.0_dp/(1+k_therm)*(P_bot(i,j,1:NVERT)*Pa_bot(i,j,1:NVERT)-Pa_bot(i,j,0:(NVERT-1))*Pa_bot(i,j,0:(NVERT-1)))
      P_c(i,j,:)=P_c(i,j,:)/(Pa_bot(i,j,1:NVERT)-Pa_bot(i,j,0:NVERT-1))      
    enddo 
  enddo
end subroutine calculate_mass_weighed_P

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate_fluxface_column_pressure
! usage: calculate column pressure at each edge 
! of each cell. The values will be used to calculate
! cell-centered column pressure
! Yun Zhang 04/30/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calculate_fluxface_column_pressure(pi_c_old,pi_uface,pi_vface)
  integer::i,j
  real(dp), dimension(NLAT,NLONG), intent(in)::pi_c_old
  real(dp), dimension(NLAT,NLONG+1), intent(inout):: pi_uface
  real(dp), dimension(NLAT+1,NLONG), intent(inout):: pi_vface

  ! calculate uface 
  ! use central differencing
  do i=1,NLAT
    if(periodicBC==1) then
      pi_uface(i,1)=0.5*(pi_c_old(i,1)+pi_c_old(i,NLONG))
      pi_uface(i,NLONG+1)=0.5*(pi_c_old(i,1)+pi_c_old(i,NLONG))
    else
      pi_uface(i,1)=pi_c_old(i,1)
      pi_uface(i,NLONG+1)=pi_c_old(i,NLONG)
    endif
    pi_uface(i,2:NLONG)=0.5*(pi_c_old(i,1:(NLONG-1))+pi_c_old(i,2:NLONG))
  enddo

  ! calculate vface
  ! use central differencing
  do i=1,NLONG
    if(periodicBC==1) then
      pi_vface(1,i)=0.5*(pi_c_old(1,i)+pi_c_old(NLAT,i))
      pi_vface(NLAT+1,i)=0.5*(pi_c_old(1,i)+pi_c_old(NLAT,i))
    else
      pi_vface(1,i)=pi_c_old(1,i)
      pi_vface(NLAT+1,i)=pi_c_old(NLAT,i)
    endif
    pi_vface(2:NLAT,i)=0.5*(pi_c_old(1:(NLAT-1),NLONG)+pi_c_old(2:NLAT,NLONG))
  enddo
end subroutine calculate_fluxface_column_pressure

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate_fluxsum
! usage:calculate fluxsum for each cell and each layer
! fluxsum(i,j,k) means the sum of flux of 4 edges
! from layer 1 to layer k at cell i,j
! and also calculate flux_uface and flux_vface
! Yun Zhang 04/30/2015
! @stanford 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calculate_fluxsum(fluxsum,flux_uface,flux_vface,u,v,pi_uface,pi_vface,lat_c,lat_vface,dsigma)
  integer:: k
  real(dp), dimension(NVERT),intent(in):: dsigma
  real(dp), dimension(NLAT,NLONG,NVERT),intent(inout):: fluxsum
  real(dp), dimension(NLAT,NLONG+1,NVERT), intent(in):: u
  real(dp), dimension(NLAT+1,NLONG,NVERT), intent(in):: v
  real(dp), dimension(NLAT,NLONG+1), intent(in):: pi_uface
  real(dp), dimension(NLAT+1,NLONG), intent(in):: pi_vface,lat_vface  
  real(dp), dimension(NLAT,NLONG),intent(in):: lat_c
  real(dp), dimension(NLAT,NLONG+1,NVERT),intent(inout):: flux_uface
  real(dp), dimension(NLAT+1,NLONG,NVERT),intent(inout):: flux_vface    

  ! reset all flux as zero
  fluxsum=0.0_dp
  flux_vface=0.0_dp
  flux_uface=0.0_dp

  ! set for the first layer 7.15 7.16
  do k=1,NVERT
    flux_uface(:,:,k)=u(:,:,k)*pi_uface*Re*dphi
    flux_vface(:,:,k)=v(:,:,k)*pi_vface*Re*dlamda_e*cos(lat_vface)
  enddo

  k=1
  fluxsum(:,:,k)=(flux_uface(:,2:(NLONG+1),k)-flux_uface(:,1:NLONG,k))*dsigma(k)
  fluxsum(:,:,k)=fluxsum(:,:,k)+(flux_vface(2:(NLAT+1),:,k)-flux_vface(1:NLAT,:,k))*dsigma(k)

  ! for other layers  
  do k=2,NVERT
    fluxsum(:,:,k)=fluxsum(:,:,k-1)+(flux_uface(:,2:(NLONG+1),k)-flux_uface(:,1:NLONG,k))*dsigma(k)
    fluxsum(:,:,k)=fluxsum(:,:,k)+(flux_vface(2:(NLAT+1),:,k)-flux_vface(1:NLAT,:,k))*dsigma(k)
  enddo
  !print *, fluxsum(1,1,NVERT),flux_uface(1,1,1),flux_uface(1,2,1),flux_vface(1,1,1),flux_vface(2,1,1)
end subroutine calculate_fluxsum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate_column_pressure
! usage: update column pressure for each cell 
! and call by regional model every time step
! Yun Zhang 04/30/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calculate_column_pressure(pi_c_old,pi_c_new,lat_c,fluxsum)
  integer:: i,j,k
  real(dp), dimension(NLAT,NLONG,NVERT),intent(in):: fluxsum
  real(dp), dimension(NLAT,NLONG),intent(inout):: pi_c_new
  real(dp), dimension(NLAT,NLONG),intent(in):: pi_c_old
  real(dp), dimension(NLAT,NLONG),intent(in):: lat_c

  ! calculate new column pressure
  pi_c_new=pi_c_old-dt/(Re**2)/dphi/dlamda_e/cos(lat_c)*fluxsum(:,:,NVERT)
end subroutine calculate_column_pressure

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate_w_sigma
! usage: update w_sigma after each time step
! for each layer using continuity equation for 
! each layer
! Yun Zhang 04/30/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calculate_w_sigma(w_sigma, pi_c_old,pi_c_new,lat_c,fluxsum,sigma_bot,dsigma)
  integer:: i,j,k
  real(dp),dimension(NLAT,NLONG,0:NVERT),intent(inout)::w_sigma
  real(dp),dimension(NLAT,NLONG),intent(in):: pi_c_old,pi_c_new,lat_c
  real(dp),dimension(NLAT,NLONG,NVERT),intent(in):: fluxsum
  real(dp),dimension(NVERT),intent(in):: dsigma
  real(dp),dimension(0:NVERT),intent(in):: sigma_bot

  ! assume the top and bottom w_sigma=0
  w_sigma(:,:,0)=0.0_dp
  w_sigma(:,:,NVERT)=0.0_dp

  ! calculate interior layers boundaries 7.21
  do i=1,NLAT
  	do j=1,NLONG
      w_sigma(i,j,1:NVERT)=-1.0_dp/(pi_c_new(i,j)*Re*Re*cos(lat_c(i,j))*dlamda_e*dphi)*fluxsum(i,j,1:NVERT)
      w_sigma(i,j,1:NVERT)=w_sigma(i,j,1:NVERT)-sigma_bot(1:NVERT)*(pi_c_new(i,j)&
        -pi_c_old(i,j))/dt/pi_c_new(i,j)
    enddo
  enddo

end subroutine calculate_w_sigma

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate_fluxface_PVT
! usage: calculate the potential virtual temperature
! at flux face
! Yun Zhang 05/01/2015
! @STANFORD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calculate_fluxface_PVT(PVT_c_old,PVT_uface,PVT_vface,PVT_bot,&
  P_c,P_bot,lat_uface,long_uface,lat_vface,long_vface,nt)
  integer::i,j,k
  integer,intent(in):: nt
  real(dp), dimension(NLAT,NLONG,NVERT), intent(in):: PVT_c_old
  real(dp), dimension(NLAT,NLONG+1,NVERT), intent(inout):: PVT_uface
  real(dp), dimension(NLAT+1,NLONG,NVERT), intent(inout):: PVT_vface
  real(dp),dimension(NLAT,NLONG,0:NVERT),intent(inout)::PVT_bot
  real(dp),dimension(NLAT+1,NLONG),intent(in):: lat_vface, long_vface
  real(dp),dimension(NLAT,NLONG+1),intent(in):: lat_uface, long_uface
  real(dp), dimension(NLAT,NLONG,NVERT),intent(in)::P_c
  real(dp), dimension(NLAT,NLONG,0:NVERT),intent(in)::P_bot
  ! calculate uface 
  ! use central differencing
  do i=1,NLAT
    do k=1,NVERT
      if(periodicBC==1) then
        PVT_uface(i,1,k)=0.5*(PVT_c_old(i,1,k)+PVT_c_old(i,NLONG,k))
        PVT_uface(i,NLONG+1,k)=0.5*(PVT_c_old(i,1,k)+PVT_c_old(i,NLONG,k))       
      else
        PVT_uface(i,1,k)=PVT_c_old(i,1,k)
        PVT_uface(i,NLONG+1,k)=PVT_c_old(i,NLONG,k)
      endif
      PVT_uface(i,2:NLONG,k)=0.5*(PVT_c_old(i,1:(NLONG-1),k)+PVT_c_old(i,2:NLONG,k))
    enddo
  enddo

  ! calculate vface
  ! use central differencing
  do i=1,NLONG
    do k=1,NVERT
      if(periodicBC==1) then
        PVT_vface(1,i,k)=0.5*(PVT_c_old(1,i,k)+PVT_c_old(NLAT,i,k)) 
        PVT_vface(NLAT+1,i,k)=0.5*(PVT_c_old(1,i,k)+PVT_c_old(NLAT,i,k))       
      else
        PVT_vface(1,i,k)=PVT_c_old(1,i,k)
        PVT_vface(NLAT+1,i,k)=PVT_c_old(NLAT,i,k)
      endif
      PVT_vface(2:NLAT,i,k)=0.5*(PVT_c_old(1:(NLAT-1),i,k)+PVT_c_old(2:NLAT,i,k))
    enddo
  enddo

  ! set boundary PTV value
  if(PTVbound==1) then
    ! vface value
    do j=1,NLONG
      do k=1,NVERT
        PVT_vface(1,j,k)=boundary_potential_virtual_temp(lat_vface(1,j),long_vface(1,j),k,nt*dt)
        PVT_vface(NLAT+1,j,k)=boundary_potential_virtual_temp(lat_vface(NLAT+1,j),long_vface(NLAT+1,j),k,nt*dt)
      enddo
    enddo
    ! uface value
    do i=1,NLAT
      do k=1,NVERT
        PVT_uface(i,1,k)=boundary_potential_virtual_temp(lat_vface(i,1),long_uface(i,1),k,nt*dt)
        PVT_uface(i,NLONG+1,k)=boundary_potential_virtual_temp(lat_uface(i,NLONG+1),long_uface(i,NLONG+1),k,nt*dt)
      enddo
    enddo
  endif

  ! layer bottom eqn 7.11
  do i=1,NLAT
    do j=1,NLONG
      ! for zero 0, P_bot(i,j,0)=P_c(i,j,0)
      PVT_bot(i,j,0)=PVT_c_old(i,j,1)
      PVT_bot(i,j,1:(NVERT-1))=(P_bot(i,j,1:(NVERT-1))-P_c(i,j,1:(NVERT-1)))*PVT_c_old(i,j,1:(NVERT-1))
      PVT_bot(i,j,1:(NVERT-1))=PVT_bot(i,j,1:(NVERT-1))+(P_c(i,j,2:NVERT)-P_bot(i,j,1:(NVERT-1)))*PVT_c_old(i,j,2:NVERT)
      PVT_bot(i,j,1:(NVERT-1))=PVT_bot(i,j,1:(NVERT-1))/(P_c(i,j,2:NVERT)-P_c(i,j,1:(NVERT-1)))
      PVT_bot(i,j,NVERT)=PVT_c_old(i,j,NVERT)
    enddo
  enddo

end subroutine calculate_fluxface_PVT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate_fluxface_qv
! usage: calculate specific humidity
! at flux face
! Yun Zhang 05/01/2015
! @STANFORD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calculate_fluxface_qv(qv_c_old,qv_uface,qv_vface,qv_bot,P_c,P_bot,&
  lat_uface,long_uface,lat_vface,long_vface,nt)
  integer::i,j,k
  integer,intent(in):: nt
  real(dp), dimension(NLAT,NLONG,NVERT), intent(in):: qv_c_old
  real(dp), dimension(NLAT,NLONG+1,NVERT), intent(inout):: qv_uface
  real(dp), dimension(NLAT+1,NLONG,NVERT), intent(inout):: qv_vface
  real(dp), dimension(NLAT,NLONG,0:NVERT),intent(inout)::qv_bot  
  real(dp), dimension(NLAT+1,NLONG),intent(in):: lat_vface, long_vface
  real(dp), dimension(NLAT,NLONG+1),intent(in):: lat_uface, long_uface
  real(dp), dimension(NLAT,NLONG,NVERT),intent(in)::P_c
  real(dp), dimension(NLAT,NLONG,0:NVERT),intent(in)::P_bot

  ! calculate uface 
  ! use central differencing
  do i=1,NLAT
    do k=1,NVERT
      if(periodicBC==1) then
        qv_uface(i,1,k)=0.5*(qv_c_old(i,1,k)+qv_c_old(i,NLONG,k))
        qv_uface(i,NLONG+1,k)=0.5*(qv_c_old(i,1,k)+qv_c_old(i,NLONG,k))
      else
        qv_uface(i,1,k)=qv_c_old(i,1,k)
        qv_uface(i,NLONG+1,k)=qv_c_old(i,NLONG,k)
      endif
      qv_uface(i,2:NLONG,k)=0.5*(qv_c_old(i,1:(NLONG-1),k)+qv_c_old(i,2:NLONG,k))
    enddo
  enddo

  ! calculate vface
  ! use central differencing
  do i=1,NLONG
    do k=1,NVERT
      if(periodicBC==1) then
        qv_vface(1,i,k)=0.5*(qv_c_old(1,i,k)+qv_c_old(NLAT,i,k))
        qv_vface(NLAT+1,i,k)=0.5*(qv_c_old(1,i,k)+qv_c_old(NLAT,i,k))
      else
        qv_vface(1,i,k)=qv_c_old(1,i,k)
        qv_vface(NLAT+1,i,k)=qv_c_old(NLAT,i,k)
      endif
      qv_vface(2:NLAT,i,k)=0.5*(qv_c_old(1:(NLAT-1),i,k)+qv_c_old(2:NLAT,i,k))
    enddo
  enddo

  ! set boundary PTV value
  if(qvbound==1) then
    ! vface value
    do j=1,NLONG
      do k=1,NVERT
        qv_vface(1,j,k)=boundary_specific_humidity(lat_vface(1,j),long_vface(1,j),k,nt*dt)
        qv_vface(NLAT+1,j,k)=boundary_specific_humidity(lat_vface(NLAT+1,j),long_vface(NLAT+1,j),k,nt*dt)
      enddo
    enddo
    ! uface value
    do i=1,NLAT
      do k=1,NVERT
        qv_uface(i,1,k)=boundary_specific_humidity(lat_vface(i,1),long_uface(i,1),k,nt*dt)
        qv_uface(i,NLONG+1,k)=boundary_specific_humidity(lat_uface(i,NLONG+1),long_uface(i,NLONG+1),k,nt*dt)
      enddo
    enddo
  endif
  
  ! layer bottom 7.25
  do i=1,NLAT
    do j=1,NLONG
      ! for zero 0, P_bot(i,j,0)=P_c(i,j,0)
      qv_bot(i,j,0)=qv_c_old(i,j,1)
      do k=1,(NVERT-1)
        if(qv_c_old(i,j,k)==qv_c_old(i,j,k+1)) then
          qv_bot(i,j,k)=qv_c_old(i,j,k)
        else if(qv_c_old(i,j,k)==0 .or. qv_c_old(i,j,k+1)==0) then
          qv_bot(i,j,k)=0.5*(qv_c_old(i,j,k)+qv_c_old(i,j,k+1))
        else
          qv_bot(i,j,k)=(log(qv_c_old(i,j,k))-log(qv_c_old(i,j,k+1)))&
            /(1/qv_c_old(i,j,k+1)-1/qv_c_old(i,j,k))          
        endif
      enddo
      qv_bot(i,j,NVERT)=qv_c_old(i,j,NVERT)
    enddo
  enddo

end subroutine calculate_fluxface_qv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate_fluxface_gas
! usage: calculate gas concentration
! at flux face
! Yun Zhang 05/01/2015
! @STANFORD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calculate_fluxface_gas(gas_c_old,gas_uface,gas_vface,gas_bot,P_c,P_bot,&
  lat_uface,long_uface,lat_vface,long_vface,nt)
  integer::i,j,k
  integer,intent(in):: nt
  real(dp), dimension(NLAT,NLONG,NVERT), intent(in):: gas_c_old
  real(dp), dimension(NLAT,NLONG+1,NVERT), intent(inout):: gas_uface
  real(dp), dimension(NLAT+1,NLONG,NVERT), intent(inout):: gas_vface
  real(dp), dimension(NLAT,NLONG,0:NVERT),intent(inout)::gas_bot  
  real(dp), dimension(NLAT+1,NLONG),intent(in):: lat_vface, long_vface
  real(dp), dimension(NLAT,NLONG+1),intent(in):: lat_uface, long_uface
  real(dp), dimension(NLAT,NLONG,NVERT),intent(in)::P_c
  real(dp), dimension(NLAT,NLONG,0:NVERT),intent(in)::P_bot

  ! calculate uface 
  ! use central differencing
  do i=1,NLAT
    do k=1,NVERT
      if(periodicBC==1) then
        gas_uface(i,1,k)=0.5*(gas_c_old(i,1,k)+gas_c_old(i,NLONG,k))
        gas_uface(i,NLONG+1,k)=0.5*(gas_c_old(i,1,k)+gas_c_old(i,NLONG,k))
      else
        gas_uface(i,1,k)=gas_c_old(i,1,k)
        gas_uface(i,NLONG+1,k)=gas_c_old(i,NLONG,k)
      endif
      gas_uface(i,2:NLONG,k)=0.5*(gas_c_old(i,1:(NLONG-1),k)+gas_c_old(i,2:NLONG,k))
    enddo
  enddo

  ! calculate vface
  ! use central differencing
  do i=1,NLONG
    do k=1,NVERT
      if(periodicBC==1) then
        gas_vface(1,i,k)=0.5*(gas_c_old(1,i,k)+gas_c_old(NLAT,i,k))
        gas_vface(NLAT+1,i,k)=0.5*(gas_c_old(1,i,k)+gas_c_old(NLAT,i,k))
      else
        gas_vface(1,i,k)=gas_c_old(1,i,k)
        gas_vface(NLAT+1,i,k)=gas_c_old(NLAT,i,k)
      endif
      gas_vface(2:NLAT,i,k)=0.5*(gas_c_old(1:(NLAT-1),i,k)+gas_c_old(2:NLAT,i,k))
    enddo
  enddo

  ! set boundary PTV value
  if(gasbound==1) then
    ! vface value
    do j=1,NLONG
      do k=1,NVERT
        gas_vface(1,j,k)=boundary_gas(lat_vface(1,j),long_vface(1,j),k,nt*dt)
        gas_vface(NLAT+1,j,k)=boundary_gas(lat_vface(NLAT+1,j),long_vface(NLAT+1,j),k,nt*dt)
      enddo
    enddo
    ! uface value
    do i=1,NLAT
      do k=1,NVERT
        gas_uface(i,1,k)=boundary_gas(lat_vface(i,1),long_uface(i,1),k,nt*dt)
        gas_uface(i,NLONG+1,k)=boundary_gas(lat_uface(i,NLONG+1),long_uface(i,NLONG+1),k,nt*dt)
      enddo
    enddo
  endif

  ! layer bottom 7.25
  do i=1,NLAT
    do j=1,NLONG
      gas_bot(i,j,0)=gas_c_old(i,j,1)
      gas_bot(i,j,NVERT)=gas_c_old(i,j,NVERT)
      do k=1,(NVERT-1)
        gas_bot(i,j,k)=0.5*(gas_c_old(i,j,k)+gas_c_old(i,j,k+1))
      enddo
    enddo
  enddo

end subroutine calculate_fluxface_gas


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate_heat_source
! usage: calculate heat source value for each time
! step for each cell at each layer
! the heat source is treated as explicit
! and defined in the heat_source function
! in source.f90
! Yun Zhang 05/01/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calculate_heat_source(Pa_c,Q,lat_c,long_c,nt)
  real(dp),dimension(NLAT,NLONG),intent(in)::lat_c,long_c
  real(dp),dimension(NLAT,NLONG,NVERT),intent(inout):: Q
  real(dp),dimension(NLAT,NLONG,NVERT),intent(in):: Pa_c
  integer, intent(in):: nt
  integer:: i,j,k
  do i=1,NLAT
    do j=1,NLONG
      do k=1,NVERT
        Q(i,j,k)=heat_source(lat_c(i,j),long_c(i,j),k,nt*dt)
        Q(i,j,k)=((1000/Pa_c(i,j,k))**k_therm)/Cp_d*Q(i,j,k)
      enddo
    enddo
  enddo
end subroutine calculate_heat_source


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate_gas_source
! usage: calculate gas source value for each time
! step for each cell at each layer
! the heat source is treated as explicit
! and defined in the heat_source function
! in source.f90
! Yun Zhang 05/01/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calculate_gas_source(Pa_c,Q,lat_c,long_c,nt)
  real(dp),dimension(NLAT,NLONG),intent(in)::lat_c,long_c
  real(dp),dimension(NLAT,NLONG,NVERT),intent(inout):: Q
  real(dp),dimension(NLAT,NLONG,NVERT),intent(in):: Pa_c
  integer, intent(in):: nt
  integer:: i,j,k
  do i=1,NLAT
    do j=1,NLONG
      do k=1,NVERT
        Q(i,j,k)=gas_source(lat_c(i,j),long_c(i,j),k,nt*dt)
      enddo
    enddo
  enddo
end subroutine calculate_gas_source

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate_gas_source
! usage: calculate gas source value for each time
! step for each cell at each layer
! the heat source is treated as explicit
! and defined in the heat_source function
! in source.f90
! Yun Zhang 05/01/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calculate_qv_source(Pa_c,Q,lat_c,long_c,nt)
  real(dp),dimension(NLAT,NLONG),intent(in)::lat_c,long_c
  real(dp),dimension(NLAT,NLONG,NVERT),intent(inout):: Q
  real(dp),dimension(NLAT,NLONG,NVERT),intent(in):: Pa_c
  integer, intent(in):: nt
  integer:: i,j,k
  do i=1,NLAT
    do j=1,NLONG
      do k=1,NVERT
        Q(i,j,k)=qv_source(lat_c(i,j),long_c(i,j),k,nt*dt)
      enddo
    enddo
  enddo
end subroutine calculate_qv_source

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate_scalar_field
! usage: update scalar transport for each cell
! at every time step 
! eqn 7.24 7.27
! use for PTV, qv and other components
! code use PTV as example
! Yun Zhang  04/30/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calculate_scalar_field(PVT_c_old, PVT_c_new, PVT_bot, &
  pi_c_old,pi_c_new,PVT_uface, PVT_vface,flux_uface,flux_vface,lat_c, rhoa_c,&
  dsigma, w_sigma_new,K_t,Pa_c,Q)
  real(dp),dimension(NLAT,NLONG,NVERT),intent(in):: PVT_c_old, K_t, Q
  real(dp),dimension(NLAT,NLONG,NVERT),intent(inout):: PVT_c_new
  real(dp),dimension(NLAT,NLONG),intent(in)::pi_c_old, pi_c_new, lat_c
  real(dp), dimension(NLAT+1,NLONG,NVERT),intent(in):: PVT_vface, flux_vface
  real(dp), dimension(NLAT,NLONG+1,NVERT),intent(in):: PVT_uface, flux_uface 
  real(dp), dimension(NLAT,NLONG,0:NVERT),intent(in):: w_sigma_new
  real(dp), dimension(NVERT),intent(in):: dsigma
  real(dp), dimension(NLAT,NLONG,0:NVERT),intent(in):: PVT_bot
  real(dp), dimension(NLAT,NLONG,NVERT),intent(in):: Pa_c, rhoa_c
  integer:: i,j,k

  ! old value 
  do k=1,NVERT
    PVT_c_new(:,:,k)=pi_c_old*PVT_c_old(:,:,k)/pi_c_new
  enddo 
  ! horizontal flux at boundary
  do k=1,NVERT
    PVT_c_new(:,:,k)=PVT_c_new(:,:,k)+dt/(pi_c_new*Re*Re*cos(lat_c)*dlamda_e*dphi)&
    *(flux_uface(:,1:NLONG,k)*PVT_uface(:,1:NLONG,k)-flux_uface(:,2:(NLONG+1),k)*PVT_uface(:,2:(NLONG+1),k)&
      +flux_vface(1:NLAT,:,k)*PVT_vface(1:NLAT,:,k)-flux_vface(2:(NLAT+1),:,k)*PVT_vface(2:(NLAT+1),:,k))
  enddo 

  ! vertical flux 
  do k=1,NVERT
    PVT_c_new(:,:,k)=PVT_c_new(:,:,k)+dt/dsigma(k)*(w_sigma_new(:,:,k-1)*PVT_bot(:,:,k-1)&
      -w_sigma_new(:,:,k)*PVT_bot(:,:,k))
  enddo

  ! turbulence vertical flux
  ! not finished
  if(turbmodel==1) then
    do k=1,NVERT
      PVT_c_new(:,:,k)=PVT_c_new(:,:,k)
    enddo
  endif

  ! heat source
  if(heatsource==1) then
    do k=1,NVERT
      PVT_c_new(:,:,k)=PVT_c_new(:,:,k)+dt*pi_c_old/pi_c_new*Q(:,:,k)
    enddo
  endif

end subroutine calculate_scalar_field

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate_geopotential
! usage: calculate geopotential for each cell center
! and layer bottom
! need to specify the boundary geopotential at bottom
! use boundary_surf_geopotential
! Yun 05/02/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calculate_geopotential(geopot_bot,geopot_c,geopot_c_tminus1,P_c,P_bot,PVT_c_old,lat_c,long_c,nt)
  real(dp), dimension(NLAT,NLONG,NVERT),intent(in)::P_c,PVT_c_old
  real(dp), dimension(NLAT,NLONG,NVERT),intent(inout)::geopot_c,geopot_c_tminus1
  real(dp), dimension(NLAT,NLONG,0:NVERT),intent(inout):: geopot_bot
  real(dp), dimension(NLAT,NLONG,0:NVERT),intent(in):: P_bot
  real(dp), dimension(NLAT,NLONG), intent(in):: lat_c,long_c
  integer,intent(in):: nt
  integer::i,j,k

  if(nt>0) geopot_c_tminus1=geopot_c
  do i=1,NLAT
    do j=1,NLONG
      geopot_bot(i,j,NVERT)=boundary_surf_geopotential(lat_c(i,j),long_c(i,j),nt*dt)
      geopot_c(i,j,NVERT)=geopot_bot(i,j,NVERT)-Cp_d*(PVT_c_old(i,j,NVERT)*(P_c(i,j,NVERT)-P_bot(i,j,NVERT)))  
      do k=NVERT-1,0,-1
        geopot_bot(i,j,k)=geopot_c(i,j,k+1)-Cp_d*(PVT_c_old(i,j,k+1)*(P_bot(i,j,k)-P_c(i,j,k+1)))
        if(k/=0) then
          geopot_c(i,j,k)=geopot_bot(i,j,k)-Cp_d*(PVT_c_old(i,j,k)*(P_c(i,j,k)-P_bot(i,j,k)))
        endif
      enddo
    enddo
  enddo
  if(nt==0) geopot_c_tminus1=geopot_c

end subroutine calculate_geopotential

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate_velocity_field
! usage: update velocity field for each flux face
! Yun Zhang 05/02/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calculate_velocity_field(pi_c_old,pi_c_new,pi_c_tmp,pi_c_tminus1,lat_c,lat_vface,&
  u_old,u_new,u_tmp,v_old,v_new,v_tmp,flux_uface,flux_vface,dsigma,&
  sigma_bot,w_sigma_new,f_c,geopot_c,geopot_c_tminus1,&
  PVT_c_old,P_bot,P_c,nu_t,rhoa_c)
  real(dp),dimension(NLAT,NLONG),intent(in)::pi_c_old,pi_c_new,pi_c_tmp,pi_c_tminus1,lat_c,f_c
  real(dp),dimension(NLAT,NLONG+1,NVERT),intent(in)::u_old,flux_uface,u_tmp
  real(dp),dimension(NLAT,NLONG+1,NVERT),intent(inout)::u_new
  real(dp),dimension(NLAT+1,NLONG,NVERT),intent(in)::v_old,flux_vface,v_tmp
  real(dp),dimension(NLAT+1,NLONG,NVERT),intent(inout)::v_new 
  real(dp),dimension(0:NVERT),intent(in)::sigma_bot
  real(dp),dimension(NVERT),intent(in)::dsigma
  real(dp),dimension(NLAT,NLONG,0:NVERT),intent(in)::w_sigma_new,P_bot
  real(dp),dimension(NLAT,NLONG,NVERT),intent(in)::geopot_c,geopot_c_tminus1,PVT_c_old,P_c,nu_t,rhoa_c
  real(dp),dimension(NLAT+1,NLONG),intent(in)::lat_vface
  real(dp),dimension(NLAT,NLONG+1,0:NVERT)::u_bot_old  
  real(dp),dimension(NLAT+1,NLONG,0:NVERT)::v_bot_old

  ! u field  
  real(dp),dimension(NLAT,NLONG+1)::A_old,A_new
  real(dp),dimension(NLAT,NLONG+1,0:NVERT)::F
  real(dp),dimension(NLAT,0:NLONG+1,NVERT)::B
  real(dp),dimension(NLAT+1,NLONG+1,NVERT)::C
  real(dp),dimension(NLAT+1,0:NLONG+1,NVERT)::D,E
  
  ! v field
  real(dp),dimension(NLAT+1,NLONG)::P_old,P_new
  real(dp),dimension(NLAT+1,NLONG,0:NVERT)::O
  real(dp),dimension(0:NLAT+1,NLONG,NVERT)::R
  real(dp),dimension(NLAT+1,NLONG+1,NVERT)::Q
  real(dp),dimension(0:NLAT+1,NLONG+1,NVERT)::S,T 
  
  integer::i,j,k
 
  ! prepare for the boundary values
  ! use new virtual matrix to calculate velocity to simplify equation
  real(dp),dimension(0:NLAT+1,0:NLONG+1)::v_f_c,v_pi_c_new,v_pi_c_old,v_lat_c
  real(dp),dimension(0:NLAT+1,0:NLONG+1,0:NVERT)::v_F_new
  real(dp),dimension(0:NLAT+1,0:NLONG+1,0:NVERT)::v_w_sigma_new
  !real(dp),dimension(0:NLAT+1,0:NLONG+1,NVERT):: v_PVT_c_old, v_P_c, v_geopot_c
  !real(dp),dimension(0:NLAT+1,0:NLONG+1,0:NVERT):: v_P_bot
  real(dp),dimension(0:NLAT+1,0:NLONG+2,NVERT)::v_u_old, v_flux_uface
  real(dp),dimension(0:NLAT+2,0:NLONG+1,NVERT)::v_v_old, v_flux_vface

  v_f_c(1:NLAT,1:NLONG)=f_c
  v_pi_c_old(1:NLAT,1:NLONG)=pi_c_old
  v_pi_c_new(1:NLAT,1:NLONG)=pi_c_new
  v_lat_c(1:NLAT,1:NLONG)=lat_c
  v_w_sigma_new(1:NLAT,1:NLONG,0:NVERT)=w_sigma_new

  if(periodicBC==1) then
    v_pi_c_old(0,1:NLONG)=pi_c_old(NLAT,:)
    v_pi_c_old(NLAT+1,1:NLONG)=pi_c_old(1,:)
    v_pi_c_old(1:NLAT,0)=pi_c_old(:,NLONG)
    v_pi_c_old(1:NLAT,NLONG+1)=pi_c_old(:,1)
    v_pi_c_old(0,0)=0.5*(pi_c_old(1,NLONG)+pi_c_old(NLAT,1))
    v_pi_c_old(NLAT,0)=0.5*(pi_c_old(1,NLONG)+pi_c_old(1,1))
    v_pi_c_old(0,NLONG+1)=0.5*(pi_c_old(1,1)+pi_c_old(NLAT,NLONG))
    v_pi_c_old(NLAT,NLONG+1)=0.5*(pi_c_old(1,NLONG)+pi_c_old(NLAT,1))
  else
    v_pi_c_old(0,1:NLONG)=pi_c_old(1,:)
    v_pi_c_old(NLAT+1,1:NLONG)=pi_c_old(NLAT,:)
    v_pi_c_old(:,0)=v_pi_c_old(:,1)
    v_pi_c_old(:,NLONG+1)=v_pi_c_old(:,NLONG)
  endif

  if(periodicBC==1) then
    v_pi_c_new(0,1:NLONG)=pi_c_new(NLAT,:)
    v_pi_c_new(NLAT+1,1:NLONG)=pi_c_new(1,:)
    v_pi_c_new(1:NLAT,0)=pi_c_new(:,NLONG)
    v_pi_c_new(1:NLAT,NLONG+1)=pi_c_new(:,1)
    v_pi_c_new(0,0)=0.5*(pi_c_new(1,NLONG)+pi_c_new(NLAT,1))
    v_pi_c_new(NLAT,0)=0.5*(pi_c_new(1,NLONG)+pi_c_new(1,1))
    v_pi_c_new(0,NLONG+1)=0.5*(pi_c_new(1,1)+pi_c_new(NLAT,NLONG))
    v_pi_c_new(NLAT,NLONG+1)=0.5*(pi_c_new(1,NLONG)+pi_c_new(NLAT,1))
  else
    v_pi_c_new(0,1:NLONG)=pi_c_new(1,:)
    v_pi_c_new(NLAT+1,1:NLONG)=pi_c_new(NLAT,:)
    v_pi_c_new(:,0)=v_pi_c_new(:,1)
    v_pi_c_new(:,NLONG+1)=v_pi_c_new(:,NLONG)
  endif

  v_lat_c(0,:)=lat_c(1,1)-dphi
  v_lat_c(NLAT+1,:)=lat_c(NLAT,1)+dphi
  v_lat_c(1:NLAT,0)=lat_c(1:NLAT,1)
  v_lat_c(1:NLAT,NLONG+1)=lat_c(1:NLAT,NLONG)

  if(coriolis==1) then
    v_f_c=2*Omega*sin(v_lat_c)
  else
    v_f_c=0.0_dp
  endif


  if(periodicBC==1) then
    v_w_sigma_new(0,1:NLONG,:)=w_sigma_new(NLAT,:,:)
    v_w_sigma_new(NLAT+1,1:NLONG,:)=w_sigma_new(1,:,:)
    v_w_sigma_new(1:NLAT,0,:)=w_sigma_new(:,NLONG,:)
    v_w_sigma_new(1:NLAT,NLONG+1,:)=w_sigma_new(:,1,:)
    v_w_sigma_new(0,0,:)=0.5*(w_sigma_new(1,NLONG,:)+w_sigma_new(NLAT,1,:))
    v_w_sigma_new(NLAT,0,:)=0.5*(w_sigma_new(1,NLONG,:)+w_sigma_new(1,1,:))
    v_w_sigma_new(0,NLONG+1,:)=0.5*(w_sigma_new(1,1,:)+w_sigma_new(NLAT,NLONG,:))
    v_w_sigma_new(NLAT,NLONG+1,:)=0.5*(w_sigma_new(1,NLONG,:)+w_sigma_new(NLAT,1,:))
  else
    v_w_sigma_new(0,1:NLONG,:)=w_sigma_new(1,:,:)
    v_w_sigma_new(NLAT+1,1:NLONG,:)=w_sigma_new(NLAT,:,:)
    v_w_sigma_new(:,0,:)=v_w_sigma_new(:,1,:)
    v_w_sigma_new(:,NLONG+1,:)=v_w_sigma_new(:,NLONG,:)
  endif

  v_u_old(1:NLAT,1:NLONG+1,:)=u_tmp
  v_flux_uface(1:NLAT,1:NLONG+1,:)=flux_uface
  v_v_old(1:NLAT+1,1:NLONG,:)=v_tmp
  v_flux_vface(1:NLAT+1,1:NLONG,:)=flux_vface

  v_u_old(0,1:NLONG+1,:)=u_tmp(1,:,:)
  v_u_old(NLAT+1,1:NLONG+1,:)=u_tmp(NLAT,:,:)
  v_u_old(:,0,:)=v_u_old(:,1,:)
  v_u_old(:,NLONG+2,:)=v_u_old(:,NLONG+1,:)

  v_flux_uface(0,1:NLONG+1,:)=flux_uface(1,:,:)
  v_flux_uface(NLAT+1,1:NLONG+1,:)=flux_uface(NLAT,:,:)
  v_flux_uface(:,0,:)=v_flux_uface(:,1,:)
  v_flux_uface(:,NLONG+2,:)=v_flux_uface(:,NLONG+1,:)

  v_v_old(1:NLAT+1,0,:)=v_tmp(:,1,:)
  v_v_old(1:NLAT+1,NLONG+1,:)=v_tmp(:,NLONG,:)
  v_v_old(0,:,:)=v_v_old(1,:,:)
  v_v_old(NLAT+2,:,:)=v_v_old(NLAT+1,:,:)

  v_flux_vface(1:NLAT+1,0,:)=flux_vface(:,1,:)
  v_flux_vface(1:NLAT+1,NLONG+1,:)=flux_vface(:,NLONG,:)
  v_flux_vface(0,:,:)=v_flux_vface(1,:,:)
  v_flux_vface(NLAT+2,:,:)=v_flux_vface(NLAT+1,:,:)

  ! Update u field
  ! calculate column pressure multiplied by grid-cell area 7.38
  ! interior points
  A_old=1.0_dp/8.0_dp*Re*Re*dphi*dlamda_e*(v_pi_c_old(2:NLAT+1,0:NLONG)*cos(v_lat_c(2:NLAT+1,0:NLONG))&
    +v_pi_c_old(2:NLAT+1,1:NLONG+1)*cos(v_lat_c(2:NLAT+1,1:NLONG+1))+2*v_pi_c_old(1:NLAT,0:NLONG)*cos(v_lat_c(1:NLAT,0:NLONG))&
    +2*v_pi_c_old(1:NLAT,1:NLONG+1)*cos(v_lat_c(1:NLAT,1:NLONG+1))+v_pi_c_old(0:NLAT-1,0:NLONG)*cos(v_lat_c(0:NLAT-1,0:NLONG))&
    +v_pi_c_old(0:NLAT-1,1:NLONG+1)*cos(v_lat_c(0:NLAT-1,1:NLONG+1)))

  A_new=1.0_dp/8.0_dp*Re*Re*dphi*dlamda_e*(v_pi_c_new(2:NLAT+1,0:NLONG)*cos(v_lat_c(2:NLAT+1,0:NLONG))&
    +v_pi_c_new(2:NLAT+1,1:NLONG+1)*cos(v_lat_c(2:NLAT+1,1:NLONG+1))&
    +2*v_pi_c_new(1:NLAT,0:NLONG)*cos(v_lat_c(1:NLAT,0:NLONG))&
    +2*v_pi_c_new(1:NLAT,1:NLONG+1)*cos(v_lat_c(1:NLAT,1:NLONG+1))&
    +v_pi_c_new(0:NLAT-1,0:NLONG)*cos(v_lat_c(0:NLAT-1,0:NLONG))&
    +v_pi_c_new(0:NLAT-1,1:NLONG+1)*cos(v_lat_c(0:NLAT-1,1:NLONG+1)))

  ! equation 7.41
  B=1.0_dp/12.0_dp*(v_flux_uface(0:NLAT-1,0:NLONG+1,:)+v_flux_uface(0:NLAT-1,1:NLONG+2,:)&
    +2*v_flux_uface(1:NLAT,0:NLONG+1,:)+2*v_flux_uface(1:NLAT,1:NLONG+2,:)&
    +v_flux_uface(2:NLAT+1,0:NLONG+1,:)+v_flux_uface(2:NLAT+1,1:NLONG+2,:))

  ! equation 7.42
  C=1.0_dp/12.0_dp*(v_flux_vface(0:NLAT,0:NLONG,:)+v_flux_vface(0:NLAT,1:NLONG+1,:)&
    +2*v_flux_vface(1:NLAT+1,0:NLONG,:)+2*v_flux_vface(1:NLAT+1,1:NLONG+1,:)&
    +v_flux_vface(2:NLAT+2,0:NLONG,:)+v_flux_vface(2:NLAT+2,1:NLONG+1,:))

  ! equation 7.43  
  D=1.0_dp/24.0_dp*(v_flux_vface(0:NLAT,0:NLONG+1,:)+2*v_flux_vface(1:NLAT+1,0:NLONG+1,:)+v_flux_vface(2:NLAT+2,0:NLONG+1,:)&
    +v_flux_uface(0:NLAT,0:NLONG+1,:)+v_flux_uface(1:NLAT+1,0:NLONG+1,:)+v_flux_uface(0:NLAT,1:NLONG+2,:)&
    +v_flux_uface(1:NLAT+1,1:NLONG+2,:))

  ! equation 7.44
  E=1.0_dp/24.0_dp*(v_flux_vface(0:NLAT,0:NLONG+1,:)+2*v_flux_vface(1:NLAT+1,0:NLONG+1,:)+v_flux_vface(2:NLAT+2,0:NLONG+1,:)&
    -v_flux_uface(0:NLAT,0:NLONG+1,:)-v_flux_uface(1:NLAT+1,0:NLONG+1,:)-v_flux_uface(0:NLAT,1:NLONG+2,:)&
    -v_flux_uface(1:NLAT+1,1:NLONG+2,:))

  ! equation 7.45 get velocity at layer bottom
  ! the top most and bottom most u_bot are not important because w_sigma is zero
  do i=1,NLAT
    do j=1,NLONG+1
      u_bot_old(i,j,1:NVERT-1)=(u_tmp(i,j,1:NVERT-1)*dsigma(1:NVERT-1)+u_tmp(i,j,2:NVERT)*dsigma(2:NVERT))&
      /(dsigma(1:NVERT-1)+dsigma(2:NVERT))
      u_bot_old(i,j,0)=u_tmp(i,j,1)
      u_bot_old(i,j,NVERT)=u_tmp(i,j,NVERT)
    enddo
  enddo

  ! equation 7.40
  do i=0,NLAT+1
    do j=0,NLONG+1
      v_F_new(i,j,:)=v_pi_c_new(i,j)*Re*Re*dphi*dlamda_e*cos(v_lat_c(i,j))*v_w_sigma_new(i,j,:)
    enddo
  enddo

  F=1.0_dp/8.0_dp*(v_F_new(2:NLAT+1,0:NLONG,:)+v_F_new(2:NLAT+1,1:NLONG+1,:)&
    +2*v_F_new(1:NLAT,0:NLONG,:)+2*v_F_new(1:NLAT,1:NLONG+1,:)&
    +v_F_new(0:NLAT-1,0:NLONG,:)+v_F_new(0:NLAT-1,1:NLONG+1,:))

  ! equation 7.32 time difference term
  do k=1,NVERT
    u_new(:,:,k)=A_old/A_new*u_old(:,:,k)
  enddo

  ! equation 7.33 horizontal advection B 
  do k=1,NVERT
    u_new(:,:,k)=u_new(:,:,k)+dt/A_new*(B(:,0:NLONG,k)&
      *0.5*(v_u_old(1:NLAT,0:NLONG,k)+v_u_old(1:NLAT,1:NLONG+1,k))&
      -B(:,1:NLONG+1,k)*0.5*(v_u_old(1:NLAT,1:NLONG+1,k)+v_u_old(1:NLAT,2:NLONG+2,k)))
  enddo

  ! equation 7.33 horizontal advection C 
  do k=1,NVERT
    u_new(:,:,k)=u_new(:,:,k)+dt/A_new*(C(1:NLAT,1:NLONG+1,k)&
      *0.5*(v_u_old(0:NLAT-1,1:NLONG+1,k)+v_u_old(1:NLAT,1:NLONG+1,k))&
      -C(2:NLAT+1,1:NLONG+1,k)*0.5*(v_u_old(1:NLAT,1:NLONG+1,k)+v_u_old(2:NLAT+1,1:NLONG+1,k)))
  enddo

  ! equation 7.33 horizontal advection D
  do k=1,NVERT
    u_new(:,:,k)=u_new(:,:,k)+dt/A_new*(D(1:NLAT,0:NLONG,k)&
      *0.5*(v_u_old(0:NLAT-1,0:NLONG,k)+v_u_old(1:NLAT,1:NLONG+1,k))&
      -D(2:NLAT+1,1:NLONG+1,k)*0.5*(v_u_old(1:NLAT,1:NLONG+1,k)+v_u_old(2:NLAT+1,2:NLONG+2,k)))
  enddo 

  ! equation 7.33 horizontal advection E
  do k=1,NVERT
    u_new(:,:,k)=u_new(:,:,k)+dt/A_new*(E(1:NLAT,1:NLONG+1,k)&
      *0.5*(v_u_old(0:NLAT-1,2:NLONG+2,k)+v_u_old(1:NLAT,1:NLONG+1,k))&
      -E(2:NLAT+1,0:NLONG,k)*0.5*(v_u_old(1:NLAT,1:NLONG+1,k)+v_u_old(2:NLAT+1,0:NLONG,k)))
  enddo 

  ! equation 7.34 vertical transport F
  do k=1,NVERT
    u_new(:,:,k)=u_new(:,:,k)+dt/A_new/dsigma(k)*(F(:,:,k-1)*u_bot_old(:,:,k-1)&
      -F(:,:,k)*u_bot_old(:,:,k))
  enddo

  ! equation 7.35 coriolis and sperical grid conversion
  do k=1,NVERT
    u_new(:,:,k)=u_new(:,:,k)+0.5*dt/A_new*Re*dlamda_e*dphi&
      *(v_pi_c_old(1:NLAT,0:NLONG)*0.5*(v_v_old(1:NLAT,0:NLONG,k)+v_v_old(2:NLAT+1,0:NLONG,k))&
      *(v_f_c(1:NLAT,0:NLONG)*Re*cos(v_lat_c(1:NLAT,0:NLONG))&
      +0.5*(v_u_old(1:NLAT,0:NLONG,k)+v_u_old(1:NLAT,1:NLONG+1,k))*sin(v_lat_c(1:NLAT,0:NLONG)))&
      +v_pi_c_old(1:NLAT,1:NLONG+1)*0.5*(v_v_old(1:NLAT,1:NLONG+1,k)+v_v_old(2:NLAT+1,1:NLONG+1,k))&
      *(v_f_c(1:NLAT,1:NLONG+1)*Re*cos(v_lat_c(1:NLAT,1:NLONG+1))&
      +0.5*(v_u_old(1:NLAT,1:NLONG+1,k)+v_u_old(1:NLAT,2:NLONG+2,k))*sin(v_lat_c(1:NLAT,1:NLONG+1))))
  enddo

  ! equation 7.36 pressure gradient
  do k=1,NVERT
    u_new(:,2:NLONG,k)=u_new(:,2:NLONG,k)-dt/A_new(:,2:NLONG)*Re*dphi*((geopot_c(:,2:NLONG,k)&
      -geopot_c(:,1:NLONG-1,k))*0.5*(pi_c_tmp(:,1:NLONG-1)+pi_c_tmp(:,2:NLONG))&
      +0.5*(-pi_c_tmp(:,1:NLONG-1)+pi_c_tmp(:,2:NLONG))*Cp_d&
      *(PVT_c_old(:,1:NLONG-1,k)/dsigma(k)*(sigma_bot(k)&
      *(P_bot(:,1:NLONG-1,k)-P_c(:,1:NLONG-1,k))+sigma_bot(k-1)&
      *(P_c(:,1:NLONG-1,k)-P_bot(:,1:NLONG-1,k-1)))+PVT_c_old(:,2:NLONG,k)/dsigma(k)&
      *(sigma_bot(k)*(P_bot(:,2:NLONG,k)-P_c(:,2:NLONG,k))&
      +sigma_bot(k-1)*(P_c(:,2:NLONG,k)-P_bot(:,2:NLONG,k-1)))))

    ! west boundary
    u_new(:,1,k)=u_new(:,1,k)-dt/A_new(:,1)*Re*dphi*((-geopot_c_tminus1(:,1,k)+geopot_c(:,1,k))*pi_c_tmp(:,1)&
      +(-pi_c_tminus1(:,1)+pi_c_tmp(:,1))*Cp_d*(PVT_c_old(:,1,k)/dsigma(k)*(sigma_bot(k)*(P_bot(:,1,k)&
      -P_c(:,1,k))+sigma_bot(k-1)*(P_c(:,1,k)-P_bot(:,1,k-1)))))

    ! east boundary
    if(periodicBC==1) then
      u_new(:,NLONG+1,k)=u_new(:,1,k);
    else
      u_new(:,NLONG+1,k)=u_new(:,NLONG+1,k)-dt/A_new(:,NLONG+1)*Re*dphi*((geopot_c_tminus1(:,NLONG,k)&
        -geopot_c(:,NLONG,k))*pi_c_tmp(:,NLONG)&
        +(pi_c_tminus1(:,NLONG)-pi_c_tmp(:,NLONG))*Cp_d*(PVT_c_old(:,NLONG,k)/dsigma(k)*(sigma_bot(k)*(P_bot(:,NLONG,k)&
        -P_c(:,NLONG,k))+sigma_bot(k-1)*(P_c(:,NLONG,k)-P_bot(:,NLONG,k-1)))))
    endif
  enddo

  ! equation 7.37 eddy vicousity not finished!
  if(turbmodel==1) then
    do k=1,NVERT
      u_new(:,:,k)=u_new(:,:,k)
    enddo
  endif

  ! Update v field
  ! column pressure multiplied by the area eqn 7.53
  P_old=1.0_dp/8.0_dp*Re*Re*dphi*dlamda_e*(v_pi_c_old(0:NLAT,2:NLONG+1)*cos(v_lat_c(0:NLAT,2:NLONG+1))&
    +v_pi_c_old(1:NLAT+1,2:NLONG+1)*cos(v_lat_c(1:NLAT+1,2:NLONG+1))&
    +2*v_pi_c_old(0:NLAT,1:NLONG)*cos(v_lat_c(0:NLAT,1:NLONG))&
    +2*v_pi_c_old(1:NLAT+1,1:NLONG)*cos(v_lat_c(1:NLAT+1,1:NLONG))&
    +v_pi_c_old(0:NLAT,0:NLONG-1)*cos(v_lat_c(0:NLAT,0:NLONG-1))&
    +v_pi_c_old(1:NLAT+1,0:NLONG-1)*cos(v_lat_c(1:NLAT+1,0:NLONG-1)))

  P_new=1.0_dp/8.0_dp*Re*Re*dphi*dlamda_e*(v_pi_c_new(0:NLAT,2:NLONG+1)*cos(v_lat_c(0:NLAT,2:NLONG+1))&
    +v_pi_c_new(1:NLAT+1,2:NLONG+1)*cos(v_lat_c(1:NLAT+1,2:NLONG+1))&
    +2*v_pi_c_new(0:NLAT,1:NLONG)*cos(v_lat_c(0:NLAT,1:NLONG))&
    +2*v_pi_c_new(1:NLAT+1,1:NLONG)*cos(v_lat_c(1:NLAT+1,1:NLONG))&
    +v_pi_c_new(0:NLAT,0:NLONG-1)*cos(v_lat_c(0:NLAT,0:NLONG-1))&
    +v_pi_c_new(1:NLAT+1,0:NLONG-1)*cos(v_lat_c(1:NLAT+1,0:NLONG-1)))

   ! equation 7.56
  R=1.0_dp/12.0_dp*(v_flux_vface(0:NLAT+1,0:NLONG-1,:)+v_flux_vface(1:NLAT+2,0:NLONG-1,:)&
    +2*v_flux_vface(0:NLAT+1,1:NLONG,:)+2*v_flux_vface(1:NLAT+2,1:NLONG,:)&
    +v_flux_vface(0:NLAT+1,2:NLONG+1,:)+v_flux_vface(1:NLAT+2,2:NLONG+1,:))

  ! equation 7.55
  Q=1.0_dp/12.0_dp*(v_flux_uface(0:NLAT,0:NLONG,:)+v_flux_uface(1:NLAT+1,0:NLONG,:)&
    +2*v_flux_uface(0:NLAT,1:NLONG+1,:)+2*v_flux_uface(1:NLAT+1,1:NLONG+1,:)&
    +v_flux_uface(0:NLAT,2:NLONG+2,:)+v_flux_uface(1:NLAT+1,2:NLONG+2,:))

  ! equation 7.57  
  S=1.0_dp/24.0_dp*(v_flux_uface(0:NLAT+1,0:NLONG,:)+2*v_flux_uface(0:NLAT+1,1:NLONG+1,:)+v_flux_uface(0:NLAT+1,2:NLONG+2,:)&
    +v_flux_vface(0:NLAT+1,0:NLONG,:)+v_flux_vface(0:NLAT+1,1:NLONG+1,:)+v_flux_vface(1:NLAT+2,0:NLONG,:)&
    +v_flux_vface(1:NLAT+2,1:NLONG+1,:))

  ! equation 7.58
  T=1.0_dp/24.0_dp*(-v_flux_uface(0:NLAT+1,0:NLONG,:)-2*v_flux_uface(0:NLAT+1,1:NLONG+1,:)-v_flux_uface(0:NLAT+1,2:NLONG+2,:)&
    +v_flux_vface(0:NLAT+1,0:NLONG,:)+v_flux_vface(0:NLAT+1,1:NLONG+1,:)+v_flux_vface(1:NLAT+2,0:NLONG,:)&
    +v_flux_vface(1:NLAT+2,1:NLONG+1,:))

  ! equation 7.54  
  O=1.0_dp/8.0_dp*(v_F_new(0:NLAT,2:NLONG+1,:)+v_F_new(1:NLAT+1,2:NLONG+1,:)&
    +2*v_F_new(0:NLAT,1:NLONG,:)+2*v_F_new(1:NLAT+1,1:NLONG,:)&
    +v_F_new(0:NLAT,0:NLONG-1,:)+v_F_new(1:NLAT+1,0:NLONG-1,:))

  ! equation 7.45 get velocity at layer bottom
  ! the top most and bottom most u_bot are not important because w_sigma is zero
  do i=1,NLAT+1
    do j=1,NLONG
      v_bot_old(i,j,1:NVERT-1)=(v_tmp(i,j,1:NVERT-1)*dsigma(1:NVERT-1)+v_tmp(i,j,2:NVERT)*dsigma(2:NVERT))&
      /(dsigma(1:NVERT-1)+dsigma(2:NVERT))
      v_bot_old(i,j,0)=v_tmp(i,j,1)
      v_bot_old(i,j,NVERT)=v_tmp(i,j,NVERT)
    enddo
  enddo

  ! equation 7.47 time difference term
  do k=1,NVERT
    v_new(:,:,k)=P_old/P_new*v_old(:,:,k)
  enddo

  ! equation 7.48 horizontal advection R 
  do k=1,NVERT
    v_new(:,:,k)=v_new(:,:,k)+dt/P_new*(R(0:NLAT,:,k)&
      *0.5*(v_v_old(0:NLAT,1:NLONG,k)+v_v_old(1:NLAT+1,1:NLONG,k))&
      -R(1:NLAT+1,:,k)*0.5*(v_v_old(1:NLAT+1,1:NLONG,k)+v_v_old(2:NLAT+2,1:NLONG,k)))
  enddo

  ! equation 7.48 horizontal advection Q
  do k=1,NVERT
    v_new(:,:,k)=v_new(:,:,k)+dt/P_new*(Q(1:NLAT+1,1:NLONG,k)&
      *0.5*(v_v_old(1:NLAT+1,0:NLONG-1,k)+v_v_old(1:NLAT+1,1:NLONG,k))&
      -Q(1:NLAT+1,2:NLONG+1,k)*0.5*(v_v_old(1:NLAT+1,1:NLONG,k)+v_v_old(1:NLAT+1,2:NLONG+1,k)))
  enddo

  ! equation 7.48 horizontal advection S
  do k=1,NVERT
    v_new(:,:,k)=v_new(:,:,k)+dt/P_new*(S(0:NLAT,1:NLONG,k)&
      *0.5*(v_v_old(0:NLAT,0:NLONG-1,k)+v_v_old(1:NLAT+1,1:NLONG,k))&
      -S(1:NLAT+1,2:NLONG+1,k)*0.5*(v_v_old(1:NLAT+1,1:NLONG,k)+v_v_old(2:NLAT+2,2:NLONG+1,k)))
  enddo 

  ! equation 7.48 horizontal advection T
  do k=1,NVERT
    v_new(:,:,k)=v_new(:,:,k)+dt/P_new*(-T(1:NLAT+1,1:NLONG,k)&
      *0.5*(v_v_old(2:NLAT+2,0:NLONG-1,k)+v_v_old(1:NLAT+1,1:NLONG,k))&
      +T(0:NLAT,2:NLONG+1,k)*0.5*(v_v_old(1:NLAT+1,1:NLONG,k)+v_v_old(0:NLAT,2:NLONG+1,k)))
  enddo 

  ! equation 7.49 vertical transport O
  do k=1,NVERT
    v_new(:,:,k)=v_new(:,:,k)+dt/P_new/dsigma(k)*(O(:,:,k-1)*v_bot_old(:,:,k-1)&
      -O(:,:,k)*v_bot_old(:,:,k))
  enddo

  ! equation 7.50 coriolis and sperical grid conversion
  do k=1,NVERT
    v_new(:,:,k)=v_new(:,:,k)-0.5*dt/P_new*Re*dlamda_e*dphi&
      *(v_pi_c_old(0:NLAT,1:NLONG)*0.5*(v_u_old(0:NLAT,1:NLONG,k)+v_u_old(0:NLAT,2:NLONG+1,k))&
      *(v_f_c(0:NLAT,1:NLONG)*Re*cos(v_lat_c(0:NLAT,1:NLONG))&
      +0.5*(v_u_old(0:NLAT,1:NLONG,k)+v_u_old(0:NLAT,2:NLONG+1,k))*sin(v_lat_c(0:NLAT,1:NLONG)))&
      +v_pi_c_old(1:NLAT+1,1:NLONG)*0.5*(v_u_old(1:NLAT+1,1:NLONG,k)+v_u_old(1:NLAT+1,2:NLONG+1,k))&
      *(v_f_c(1:NLAT+1,1:NLONG)*Re*cos(v_lat_c(1:NLAT+1,1:NLONG))&
      +0.5*(v_u_old(1:NLAT+1,1:NLONG,k)+v_u_old(1:NLAT+1,2:NLONG+1,k))*sin(v_lat_c(1:NLAT+1,1:NLONG))))
  enddo 

  ! equation 7.51 pressure gradient
  do k=1,NVERT
    v_new(2:NLAT,:,k)=v_new(2:NLAT,:,k)-dt/P_new(2:NLAT,:)*Re*dlamda_e*cos(lat_vface(2:NLAT,:))&
      *((geopot_c(2:NLAT,:,k)-geopot_c(1:NLAT-1,:,k))*0.5*(pi_c_old(1:NLAT-1,:)+pi_c_old(2:NLAT,:))&
      +0.5*(-pi_c_old(1:NLAT-1,:)+pi_c_old(2:NLAT,:))*Cp_d&
      *(PVT_c_old(1:NLAT-1,:,k)/dsigma(k)*(sigma_bot(k)&
      *(P_bot(1:NLAT-1,:,k)-P_c(1:NLAT-1,:,k))+sigma_bot(k-1)&
      *(P_c(1:NLAT-1,:,k)-P_bot(1:NLAT-1,:,k-1)))+PVT_c_old(2:NLAT,:,k)/dsigma(k)&
      *(sigma_bot(k)*(P_bot(2:NLAT,:,k)-P_c(2:NLAT,:,k))&
      +sigma_bot(k-1)*(P_c(2:NLAT,:,k)-P_bot(2:NLAT,:,k-1)))))

    ! south boundary
    v_new(1,:,k)=v_new(1,:,k)-dt/P_new(1,:)*Re*dlamda_e*cos(lat_vface(1,:))&
      *((-geopot_c_tminus1(1,:,k)+geopot_c(1,:,k))*pi_c_old(1,:)&
      +(-pi_c_tminus1(1,:)+pi_c_old(1,:))*Cp_d*(PVT_c_old(1,:,k)/dsigma(k)*(sigma_bot(k)*(P_bot(1,:,k)&
      -P_c(1,:,k))+sigma_bot(k-1)*(P_c(1,:,k)-P_bot(1,:,k-1)))))
    ! north boundary
    if(periodicBC==1) then
      v_new(NLAT+1,:,k)=v_new(1,:,k)
    else
      v_new(NLAT+1,:,k)=v_new(NLAT+1,:,k)-dt/P_new(NLAT+1,:)*Re*dlamda_e*cos(lat_vface(NLAT+1,:))&
        *((geopot_c_tminus1(NLAT,:,k)-geopot_c(NLAT,:,k))*pi_c_old(NLAT,:)&
        +(pi_c_tminus1(NLAT,:)-pi_c_old(NLAT,:))*Cp_d*(PVT_c_old(NLAT,:,k)/dsigma(k)*(sigma_bot(k)*(P_bot(NLAT,:,k)&
        -P_c(NLAT,:,k))+sigma_bot(k-1)*(P_c(NLAT,:,k)-P_bot(NLAT,:,k-1)))))
    endif
  enddo

  ! equation 7.52 eddy vicousity not finished!
  if(turbmodel==1) then
    do k=1,NVERT
      v_new(:,:,k)=v_new(:,:,k)
    enddo
  endif

end subroutine calculate_velocity_field

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate_pressure_field
! usage: update pressure field using pi_c_new
! Yun Zhang 05/02/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calculate_pressure_field(pi_c,Pa_c,Pa_bot,P_c,P_bot,sigma_bot)
  real(dp),dimension(NLAT,NLONG),intent(in)::pi_c
  real(dp),dimension(NLAT,NLONG,0:NVERT),intent(inout)::Pa_bot,P_bot
  real(dp),dimension(0:NVERT),intent(in)::sigma_bot
  real(dp),dimension(NLAT,NLONG,NVERT),intent(inout)::Pa_c,P_c
  integer::i,j

  ! calculate pressure at different vertical layer 
  ! layer bottom
  do i=1,NLAT
    do j=1,NLONG
      Pa_bot(i,j,:)=sigma_bot*pi_c(i,j)+Pa_top
      P_bot(i,j,:)=(Pa_bot(i,j,:)/1000)**k_therm
    enddo 
  enddo

  ! layer center
  do i=1,NLAT
    do j=1,NLONG
      P_c(i,j,:)=1.0_dp/(1+k_therm)*(P_bot(i,j,1:NVERT)*Pa_bot(i,j,1:NVERT)-P_bot(i,j,0:(NVERT-1))*Pa_bot(i,j,0:(NVERT-1)))
      P_c(i,j,:)=P_c(i,j,:)/(Pa_bot(i,j,1:NVERT)-Pa_bot(i,j,0:NVERT-1))
      Pa_c(i,j,:)=1000*(P_c(i,j,:)**(1/k_therm))
    enddo
  enddo  

end subroutine calculate_pressure_field

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate_air_temperature
! usage: update air temperature with the new 
! Potential virtual temperature
! Yun Zhang 05/03/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calculate_air_temperature(PVT_c,Temp_c,qv_c,Pa_c)
  real(dp),dimension(NLAT,NLONG,NVERT),intent(in)::PVT_c,qv_c,Pa_c
  real(dp),dimension(NLAT,NLONG,NVERT),intent(inout)::Temp_c
  Temp_c=PVT_c/(1+0.608*qv_c)*((Pa_c/1000)**k_therm)
end subroutine calculate_air_temperature


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate_horizontal_bottom_value
! usage: calculate the horizontal value (u,v)
! at the layer bottom of each edege
! Yun Zhang 05/22/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calculate_horizontal_bottom_value(u,u_bot,N1,N2,dsigma)
  real(dp),dimension(N1,N2,NVERT),intent(in)::u
  real(dp),dimension(N1,N2,0:NVERT),intent(inout)::u_bot
  integer,intent(in)::N1,N2
  real(dp),dimension(NVERT),intent(in)::dsigma
  integer:: i,j
  do i=1,N1
    do j=1,N2
      u_bot(i,j,1:NVERT-1)=(u(i,j,1:NVERT-1)*dsigma(1:NVERT-1)+u(i,j,2:NVERT)*dsigma(2:NVERT))&
      /(dsigma(1:NVERT-1)+dsigma(2:NVERT))
      u_bot(i,j,0)=u(i,j,1)
      u_bot(i,j,NVERT)=u(i,j,NVERT)
    enddo
  enddo   
end subroutine calculate_horizontal_bottom_value

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate_horizontal_center_value
! usage: calculate the horizontal value (u,v)
! at the layer bottom of each edege
! Yun Zhang 05/22/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calculate_horizontal_center_value(u,u_c,N1,N2)
  integer,intent(in)::N1,N2
  real(dp),dimension(N1,N2,NVERT),intent(in)::u
  real(dp),dimension(NLAT,NLONG,NVERT),intent(inout)::u_c
  if(N2>NLONG) u_c=0.5*(u(:,1:NLONG,:)+u(:,2:NLONG+1,:))
  if(N1>NLAT) u_c=0.5*(u(1:NLAT,:,:)+u(2:NLAT+1,:,:))
end subroutine calculate_horizontal_center_value

end module phys









