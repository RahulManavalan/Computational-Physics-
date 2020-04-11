program rk_solver 

!_______________________________________________________________________________
! Declaration of Variables 

	implicit none 
	integer , parameter :: P = 20000.
	real , dimension(P) :: T,X1,X2 
	real :: Ti , Tf , X10 , X20 
	integer :: Nt 
	integer :: i 	
!_______________________________________________________________________________
! Input : 
	print*, "Runge Kutta method for Second Order Ordinary differential equation" 
	print*, "Enter Nt, Ti, Tf , X10 ,X20 " 
	read*, Nt , Ti , Tf , X10 , X20 
	print*, "Nt = " , Nt 
	print*, "Time : Initial Ti = ", Ti , "Final Tf = " , Tf 
	print*, " X1(Ti) = " , X10 , "X2(Ti) = " , X20 

	if (Nt.gt.P) stop "Nt>P" 
!_______________________________________________________________________________
	call RK(T,X1,X2,Ti,Tf,X10,X20,Nt) 
!_______________________________________________________________________________
! Output 
	open(unit = 11 , file = "rk.dat") 

	do i = 1,Nt 
		write(11,*) T(i), X1(i) , X2(i) 
	enddo
	close(11) 
	
end program rk_solver 
!_______________________________________________________________________________

real function f1(t,x1,x2)	
	implicit none 
	real :: t, x1, x2 
	f1 = x2 
end function f1 
!_______________________________________________________________________________

real function f2(t,x1,x2) 
	implicit none 
	real :: t, x1 , x2 
	f2 = -10.0 * x1 
end function f2 
!_______________________________________________________________________________

subroutine RK(T,X1,X2,Ti,Tf,X10,X20,Nt) 
	implicit none 
	integer :: Nt 
	integer :: i 
	real :: Ti,Tf,X10,X20 
	real,dimension(Nt) :: T, X1, X2 
	real :: dt 
	real :: TS , X1S , X2S 
	
! Initialization 

	dt = (Tf - Ti ) / (Nt -1) 
	T(1) = Ti 
	X1(1) = X10 
	X2(1) = X20 
	TS = Ti 
	X1S = X10 
	X2s = X20
	
	do i = 2 , Nt 
		call RKSTEP(TS,X1S,X2S,dt) 
		T(i) = Ts 
		X1(i) = X1S 
		X2(i) = X2S
	enddo 
end subroutine RK 
!_______________________________________________________________________________

subroutine RKSTEP(t,x1,x2,dt) 
	implicit none 
	real :: t , x1 ,x2 , dt 
	real :: f1 , f2
	real :: k11 , k12 , k13, k14 , k21 , k22 , k23 , k24 
	real :: h , h2 , h6 
	
! Initialization
	h = dt 
	h2 = 0.5D0 * h 
	h6 = h / 6.0
	
	k11 = f1(t, x1, x2) 
	k21 = f2(t, x1, x2) 
	k12 = f1(t+h2, x1+h2*k11, x2+h2*k21)
	k22 = f2(t+h2, x1+h2*k11, x2+h*k21)
	k13 = f1(t+h2, x1+h2*k12, x2+h*k22)
	k23 = f2(t+h2, x1+h2*k12, x2+h*k22) 
	k14 = f1(t+h, x1+h*k13 , x2+h*k23)
	k24 = f2(t+h, x1+h*k13 , x2+h*k23)
	
	t = t+h 
	x1 = x1 + h6 *(k11+ 2.0*(k12+k13) + k14 ) 
	x2 = x2 + h6 *(k21+ 2.0*(k22+k23) + k24 ) 

end subroutine RKSTEP
!_______________________________________________________________________________
	
	
