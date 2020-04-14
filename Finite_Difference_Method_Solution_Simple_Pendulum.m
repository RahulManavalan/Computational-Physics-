% Motion of a Simple Pendulum 

clear ; 
length  = 1 ;   % length of the pendulum in metres
g = 9.81 ;      % acceleration due to gravity 
npoints = 1000 ; % Discretize time into 250 intervals 
dt = 0.04 ;     % Time step in seconds 


% Vector Assignment and initialization  

omega = zeros(npoints,1) ; 
theta = zeros(npoints,1) ; 
time = zeros(npoints,1) ; 

theta(1) = 0.2 ; 

% Numerical Solution 

for step = 1 : npoints-1 
    omega(step+1) = omega(step) - (g/length) * theta(step) * dt ;
    theta(step+1) = theta(step) + omega(step) *dt ; 
    time(step+1) = time(step) + dt ; 
end

% Plot the solution 


plot(time,theta,'r') ;
grid();
xlabel("Time is seconds"); 
ylabel("theta in radians");


