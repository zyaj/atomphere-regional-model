!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! regional model module - basic_state
! usage: external procedure for all air 
! state equation to transfer air, temp
! pressure, density and other related 
! aspects
! Yun Zhang 04/24/2015
! @stanford
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module basic_state
  use constant_parameter
  implicit none
  
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! standard_atmophere_interp 
! usage: use table B.1 to interpolate pressure, altitude,
! gravity, temperature, air density
! input: a vertical array for horizontal cell
! output: get the relevant data based on input
! input_mode: 'z'=altitude,'p'=pressure,'g'=gravity
! 'T'=temperature,'rho'=density
! output_mode is same as input_mode
! Yun Zhang @Stanford
! 04/25/2015
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function standard_atmophere_interp(input,N,input_mode,output_mode) result(output)
  integer, intent(in)::N
  real(dp), dimension(1:N),intent(in)::input
  real(dp), dimension(1:N)::output
  character, intent(in)::input_mode,output_mode
  real(dp), dimension(5):: tmp_data
  real(dp), dimension(79):: alt, p, g, T, rho, basedata, outdata
  integer:: i,j,k,order1,order2

  ! read tableb1
  open (unit=99, file='tableb1.txt', status='old', action='read')
  do i=1,79
    read(99,*) tmp_data
    alt(i)=tmp_data(1)*1000
    p(i)=tmp_data(3)
    g(i)=tmp_data(2)
    T(i)=tmp_data(4)
    rho(i)=tmp_data(5)
  enddo
  close(99)
  order1=1
  order2=-1
  select case(input_mode)
    case ('z')
      basedata=alt
      order1=-1
      order2=-1
    case ('p')
      basedata=p
    case ('g')
      basedata=g
    case ('T')
      basedata=T
    case ('rho')
      basedata=rho
  end select
  
  select case(output_mode)
    case ('z')
      outdata=alt
    case ('p')
      outdata=p
    case ('g')
      outdata=g
    case ('T')
      outdata=T
    case ('rho')
      outdata=rho
  end select
 
  ! find section no. at basedata for each value in input data
  ! then interpolate data based on the section number
  k=79
  do i=1,N
    do j=1,k
      if((input(i)*order1)>(basedata(j)*order1)) then
        k=j
        output(i)=outdata(j+order2*1)+(input(i)-basedata(j+order2*1))&
          /(basedata(j)-basedata(j+order2*1))*(outdata(j)-outdata(j+order2*1))
        !print *, j, input_mode, input(i), basedata(j),basedata(j+order2*1), outdata(j),outdata(j+order2*1), output(i)
        
       exit
      endif
    enddo
  enddo
end function standard_atmophere_interp

end module	basic_state