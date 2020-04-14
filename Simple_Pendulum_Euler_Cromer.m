% Simple Pendulum using Euler Cromer Method

clear ; 
length  = 1 ;   % length of the pendulum in metres
g = 9.81 ;      % acceleration due to gravity 
npoints = 100 ; % Discretize time into 250 intervals 
dt = 0.04 ;     % Time step in seconds 


% Vector Assignment and initialization  

omega = zeros(npoints,1) ; 
theta = zeros(npoints,1) ; 
time = zeros(npoints,1) ; 

theta(1) = 0.2 ; 

% Numerical Solution 

for step = 1 : npoints - 1 
    omega(step+1) = omega(step) - (g/length) * theta(step) * dt ; 
    theta(step+1) = theta(step) + omega(step+1) * dt ; 
    time(step+1) = time(step) +dt ; 
end

plot(time,theta,'r') ;
xlabel("Time in seconds");
ylabel("theta in radians"); 
grid() ;


