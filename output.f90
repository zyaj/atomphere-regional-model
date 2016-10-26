!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! regional model module - output
! usage: output all results
! Yun Zhang 05/03/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module output
  use allocate_variable
  use constant_parameter
  implicit none

contains 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! output_open_files
! usage: open all output files for output
! Yun Zhang 05/03/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_open_bin_files()
  open (unit = pi_file_no, file = resultfolder//pi_file, status='new',form='unformatted',recl=2*10)
  open (unit = Pa_file_no, file = resultfolder//Pa_file, status='new',form='unformatted',recl=2*10)
  if(PVTmodel==1) open (unit = PVT_file_no, file = resultfolder//PTV_file, status='new',form='unformatted',recl=2*10)
  open (unit = u_file_no, file = resultfolder//u_file, status='new',form='unformatted',recl=2*10)
  open (unit = v_file_no, file = resultfolder//v_file, status='new',form='unformatted',recl=2*10)
  open (unit = w_file_no, file = resultfolder//w_file, status='new',form='unformatted',recl=2*10)
  open (unit = nu_file_no, file = resultfolder//nu_file, status='new',form='unformatted',recl=2*10)
  open (unit = K_t_file_no, file = resultfolder//K_t_file, status='new',form='unformatted',recl=2*10)
  open (unit = geopot_file_no, file = resultfolder//geopot_file, status='new',form='unformatted',recl=2*10)
  if(PVTmodel==1) open (unit = Temp_file_no, file = resultfolder//Temp_file, status='new',form='unformatted',recl=2*10)
  open (unit = rhoa_file_no, file = resultfolder//rhoa_file, status='new',form='unformatted',recl=2*10)
  if(qvmodel==1) open (unit = qv_file_no, file = resultfolder//qv_file, status='new',form='unformatted',recl=2*10)
  if(gasmodel==1) open (unit = gas_file_no, file = resultfolder//qv_file, status='new',form='unformatted',recl=2*10)

end subroutine output_open_bin_files


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! output_open_files
! usage: open all output files for output
! Yun Zhang 05/03/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_open_txt_files()
  open (unit = pi_file_no, file = resultfolder//pi_file)
  open (unit = Pa_file_no, file = resultfolder//Pa_file)
  if(PVTmodel==1) open (unit = PVT_file_no, file = resultfolder//PTV_file)
  open (unit = u_file_no, file = resultfolder//u_file)
  open (unit = v_file_no, file = resultfolder//v_file)
  open (unit = w_file_no, file = resultfolder//w_file)
  open (unit = nu_file_no, file = resultfolder//nu_file)
  open (unit = K_t_file_no, file = resultfolder//K_t_file)
  open (unit = geopot_file_no, file = resultfolder//geopot_file)
  if(PVTmodel==1) open (unit = Temp_file_no, file = resultfolder//Temp_file)
  open (unit = rhoa_file_no, file = resultfolder//rhoa_file)
  if(qvmodel==1) open (unit = qv_file_no, file = resultfolder//qv_file)
  if(gasmodel==1) open (unit = gas_file_no, file = resultfolder//gas_file)
end subroutine output_open_txt_files

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! output_all_variables
! usage: output all variables
! Yun Zhang 05/03/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_all_variables()

  ! column pressure
  call output_2Dcellcenter_variables(pi_c,pi_file_no,pi_format)

  ! air pressure
  call output_3Dcellcenter_variables(Pa_c,Pa_file_no,Pa_format)

  ! u
  call output_3Duface_variables(u,u_file_no,u_format)

  ! v
  call output_3Dvface_variables(v,v_file_no,v_format)

  ! w_sigma
  call output_3Dvertical_variables(w_sigma,w_file_no,w_format)

  ! eddy viscousity
  call output_3Dcellcenter_variables(nu_t,nu_file_no,nu_t_format)

  ! eddy diffusivity
  call output_3Dcellcenter_variables(K_t,K_t_file_no,K_t_format)

  ! geopotential
  call output_3Dcellcenter_variables(geopot_c,geopot_file_no,geopot_format)

  if(PVTmodel==1) then
    ! Temperature
    call output_3Dcellcenter_variables(Temp_c,Temp_file_no,Temp_format) 
    ! potential virtual temperature
    call output_3Dcellcenter_variables(PVT_c,PVT_file_no,PVT_format)
  endif

  ! air density
  call output_3Dcellcenter_variables(rhoa_c,rhoa_file_no,rhoa_format)

  ! specific humidity
  if(qvmodel==1) call output_3Dcellcenter_variables(qv_c,qv_file_no,qv_format)
  
  ! gas concentration
  if(gasmodel==1) call output_3Dcellcenter_variables(gas_c,gas_file_no,gas_format)

end subroutine output_all_variables 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! output_3Dcellcenter_variables
! usage: output 3D cellcenter all variables
! Yun Zhang 05/03/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_3Dcellcenter_variables(Value,file_no,format)
  real(dp),dimension(NLAT,NLONG,NVERT),intent(in):: Value
  integer, intent(in):: file_no
  integer:: i,j,k,ios
  character(*), intent(in)::format
  do i=1,NLAT
  	do j=1,NLONG
      do k=1,NVERT-1
        write(file_no, format, iostat=ios, advance="no") Value(i,j,k)
      enddo
      write(file_no, format, iostat=ios) Value(i,j,NVERT)
      if ( ios /= 0 ) then 
       print *,file_no
       stop "Write error in output file"        
      endif
    enddo
  enddo
end subroutine output_3Dcellcenter_variables 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! output_3Duface_variables
! usage: output 3D variables at u flux face
! Yun Zhang 05/03/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_3Duface_variables(Value,file_no,format)
  real(dp),dimension(NLAT,NLONG+1,NVERT),intent(in):: Value
  integer, intent(in):: file_no
  integer:: i,j,k,ios
  character(*), intent(in)::format
  do i=1,NLAT
  	do j=1,NLONG+1
      do k=1,NVERT-1
        write(file_no, format, iostat=ios,advance="no") Value(i,j,k)
      enddo
      write(file_no, format, iostat=ios) Value(i,j,NVERT)
      if ( ios /= 0 ) then 
       print *,file_no
       stop "Write error in output file"        
      endif
    enddo
  enddo
end subroutine output_3Duface_variables 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! output_3Duface_variables
! usage: output 3D variables at v flux face
! Yun Zhang 05/03/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_3Dvface_variables(Value,file_no,format)
  real(dp),dimension(NLAT+1,NLONG,NVERT),intent(in):: Value
  integer, intent(in):: file_no
  integer:: i,j,k,ios
  character(*), intent(in)::format

  do i=1,NLAT+1
  	do j=1,NLONG
      do k=1,NVERT-1
        write(file_no, format, iostat=ios,advance="no") Value(i,j,k)
      enddo
      write(file_no, format, iostat=ios) Value(i,j,NVERT)
      if ( ios /= 0 ) then 
       print *,file_no
       stop "Write error in output file"        
      endif
    enddo
  enddo
end subroutine output_3Dvface_variables 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! output_3Dvertical_variables
! usage: output 3D variables at vertical face layer top-bottom
! Yun Zhang 05/03/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_3Dvertical_variables(Value,file_no,format)
  real(dp),dimension(NLAT,NLONG,0:NVERT),intent(in):: Value
  integer, intent(in):: file_no
  integer:: i,j,k,ios
  character(*), intent(in)::format

  do i=1,NLAT
  	do j=1,NLONG
      do k=0,NVERT-1
        write(file_no, format, iostat=ios,advance="no") Value(i,j,k)
      enddo
      write(file_no, format, iostat=ios) Value(i,j,NVERT)
      if ( ios /= 0 ) then 
       print *,file_no
       stop "Write error in output file"        
      endif
      enddo
  enddo
end subroutine output_3Dvertical_variables 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! output_2Dcellcenter_variables
! usage: output 2D cellcenter all variables
! Yun Zhang 05/03/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_2Dcellcenter_variables(Value,file_no,format)
  real(dp),dimension(NLAT,NLONG),intent(in):: Value
  integer, intent(in):: file_no
  integer:: i,j,k,ios
  character(*), intent(in)::format

  do i=1,NLAT
    do j=1,NLONG-1
      write(file_no, format, iostat=ios, advance="no") Value(i,j)
    enddo
    write(file_no, format, iostat=ios) Value(i,NLONG)
    if ( ios /= 0 ) then 
     print *,file_no
     stop "Write error in output file"        
    endif
  enddo
end subroutine output_2Dcellcenter_variables 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! output_2Duface_variables
! usage: output 2D cellcenter at uface
! Yun Zhang 05/03/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_2Duface_variables(Value,file_no,format)
  real(dp),dimension(NLAT,NLONG+1),intent(in):: Value
  integer, intent(in):: file_no
  integer:: i,j,k,ios
  character(*), intent(in)::format
  do i=1,NLAT
    do j=1,NLONG
      write(file_no, format, iostat=ios,advance="no") Value(i,j)
    enddo
    write(file_no, format, iostat=ios) Value(i,NLONG+1)

    if ( ios /= 0 ) then 
      print *,file_no
      stop "Write error in output file"        
    endif
  enddo
end subroutine output_2Duface_variables 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! output_2Dvface_variables
! usage: output 2D cellcenter at v flux face
! Yun Zhang 05/03/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_2Dvface_variables(Value,file_no,format)
  real(dp),dimension(NLAT+1,NLONG),intent(in):: Value
  integer, intent(in):: file_no
  integer:: i,j,k,ios
  character(*), intent(in)::format
  do i=1,NLAT+1
    do j=1,NLONG-1
      write(file_no, format, iostat=ios,advance="no") Value(i,j)
    enddo
    write(file_no, format, iostat=ios) Value(i,NLONG)
    if ( ios /= 0 ) then 
     print *,file_no
     stop "Write error in output file"        
    endif
  enddo
end subroutine output_2Dvface_variables 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! output_close_files
! usage: close all output files for output
! Yun Zhang 05/03/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_close_files()
  close (pi_file_no)
  close (Pa_file_no)
  close (PVT_file_no)
  close (u_file_no)
  close (v_file_no)
  close (w_file_no)
  close (nu_file_no)
  close (K_t_file_no)
  close (geopot_file_no)
  close (Temp_file_no)
  close (rhoa_file_no)
  close (qv_file_no)
end subroutine output_close_files

end module output