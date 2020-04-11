program particle_motion_1 
	implicit none 
	
	! Declaration of variables 
	integer , parameter :: P = 110000 
	real , dimension(P) :: T,X,V 
	real :: Xin , Vin , Tfi 
	integer :: Nt , i 
	
	! Input initial conditions 
	print*, "Enter X_0 , V_0 , t_f , Nt (t_1 = 0) " 
	read(5,*) Xin , Vin , Tfi , Nt 
	
	! Error Handling 
	if (Nt .ge. P) then 
		print*, " The value of Nt must be less than the allocated array memory size : P " , Nt , P 
		stop 
	endif 
	!__________________________________________________________________________________
	! Calling Euler Function 
	call euler(Xin , Vin , Tfi , Nt , T , X , V ) 
	open(unit = 1 , file = "Euler.dat") 
	
	! Write Results of X,V vs T 
	do i = 1,Nt 
		write(1,*) T(i), X(i), V(i) 
	enddo 
	close(1) 
	!__________________________________________________________________________________ 
	! Calling Euler_Cromer Function 
	call euler_cromer(Xin , Vin , Tfi , Nt , T , X , V ) 
	open(unit = 1 , file = "Euler_Cromer.dat") 
	
	! Write Results of X,V vs T 
	do i = 1,Nt 
		write(1,*) T(i), X(i), V(i) 
	enddo 
	close(1)
	!__________________________________________________________________________________
	! Calling Euler_Verlet Function 
	call euler_verlet(Xin , Vin , Tfi , Nt , T , X , V ) 
	open(unit = 1 , file = "Euler_Verlet.dat") 
	
	! Write Results of X,V vs T 
	do i = 1,Nt 
		write(1,*) T(i), X(i), V(i) 
	enddo 
	close(1)
	!__________________________________________________________________________________
	
end program particle_motion_1

! Define Acceleration Function 

real function acc(x) 
	implicit none 
	real x 
	acc = - 10.0 * sin(x) 
end function acc


!__________________________________________________________________________________
! Euler's method 

subroutine euler(Xin,Vin,Tfi,Nt,T,X,V) 
	implicit none 
	integer :: Nt 
	real , dimension(Nt) :: T,X,V 
	real :: Xin,Vin,Tfi 
	integer :: i 
	real :: h , acc 
	
	! initialization 
	T(1) = 0 
	X(1) = Xin
	V(1) = Vin 
	
	h = (Tfi-0)/(Nt-1) 
	do i = 2 , Nt 
		T(i) = T(i-1) + h 
		X(i) = X(i-1) + V(i-1)*h 
		V(i) = V(i-1) + acc(X(i-1)) * h 
	enddo 
	
end subroutine euler 

!__________________________________________________________________________________
! Euler_Cromer method 

subroutine euler_cromer(Xin , Vin , Tfi , Nt , T , X ,V ) 

	implicit none 
	integer :: Nt 
	real,dimension(Nt) :: T,X,V 
	real :: Xin , Vin , Tfi 
	integer :: i 
	real :: acc , h 
	
	! initialization 
	T(1) = 0 
	X(1) = Xin 
	V(1) = Vin 
	
	h = (Tfi - 0 ) / (Nt -1 ) 
	do i = 2,Nt 
		T(i) = T(i-1) + h 
		V(i) = V(i-1) + acc(X(i-1)) * h 
		X(i) = X(i-1) + V(i) * h 
	enddo 
	
end subroutine euler_cromer 

!__________________________________________________________________________________

!Euler_Verlet method 

subroutine euler_verlet ( Xin ,Vin , Tfi , Nt , T ,X , V ) 
	implicit none 
	integer :: Nt 
	real , dimension(Nt) :: T,X,V 
	real  :: Xin , Vin , Tfi 
	real :: acc
	integer :: i 
	real :: h , h2 , oh2 , X0 
	
	! initialization 
	
	T(1) = 0.0 
	X(1) = Xin 
	V(1) = Vin
	h = (Tfi - 0)/(Nt-1)
	h2 = h*h
	oh2 = 0.5/h
	
	
	X0 = X(1) - V(1) * h + acc(X(1))*h2 / 2.0
	T(2) = h 
	X(2) = 2.0*X(1) - X0 + acc(X(1))*h2
	
	do i = 3 , Nt 
		T(i) = T(i-1) + h 
		X(i) = 2.0*X(i-1) - X(i-2) + acc(X(i-1))*h2
		V(i-1) = oh2 * (X(i) - X(i-2))
	enddo 
	
	V(Nt) = (X(Nt) - X(Nt-1))/h 
	
end subroutine euler_verlet 

!__________________________________________________________________________________
	

