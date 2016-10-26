program test
	real,dimension(0:1):: aa
    integer:: aaa,i
    real,dimension(1,2):: a,b
    open (unit=99, file='11.txt', status='new')
    a(1,1)=1
    a(1,2)=2
    write(99,*) a(1,:)
    print *,a(1,:)
end program test