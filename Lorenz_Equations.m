% Lorenz Equations of Motion 

clear 
a = 10 ; 
b = 8/3 ; 
r = 25 ; 
sigma = 10 ; 
npoints = 500000 ; 
dt = 0.0001 ; 
x = zeros(npoints,1) ; 
y = zeros(npoints,1) ; 
z = zeros(npoints,1) ; 
time = zeros(npoints,1) ; 
x(1) = 1 ; 

for step = 1 : npoints - 1 
   
    x(step+1) = x(step) + sigma*(y(step) - x(step)) *dt ; 
    y(step+1) = y(step) + (-x(step)*z(step) + r*x(step) - y(step))*dt ; 
    z(step+1) = z(step) + (x(step)*y(step) - b*z(step))*dt ; 
    time(step+1) = time(step) + dt ; 
    
end

subplot(4,1,1) ;
plot(time,z,'b') ; 
xlabel("time") ; 
ylabel("z") ; 

subplot(4,1,2) ; 
plot(time,x,'r') ; 
xlabel("time") ; 
ylabel("x") ; 

subplot(4,1,3) ;
plot(time,y,'g') ; 
xlabel("time") ; 
ylabel("y") ; 

subplot(4,1,4) ; 
plot(x,z,'y') ; 
xlabel("x") ; 
ylabel("z") ;

