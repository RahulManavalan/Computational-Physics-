% 1D radioactive Decay 

% Solve the equation : dN/dt = -N/tau 

N_uranium_initial  = 1000 ;   % initial number of uranium atoms 
npoints = 100 ;               % discretize time into 100 intervals 
dt = 1e7 ;                    % time step in years 
tau = 4.4e9 ;                 % mean lifetime of 238 U 

N_uranium = zeros(npoints,1);  % initializes N_uranium as a vector of order npoints X 1 
time = zeros(npoints,1);       % initializes time as a vector of order npoints X 1 bearing zeros


% Numerical Solution  
 
N_uranium(1) = N_uranium_initial;
time(1) = 0 ; 

for step=1:npoints-1
    N_uranium(step+1) = N_uranium(step) - (N_uranium(step) / tau) * dt ;
    time(step+1) = time(step) + dt ; 
end

% Analytical Solution 
 
t = 0 : 1e8 : 10e9 ; 
N_analytical = N_uranium_initial * exp(-t/tau) ; 

% plot both the analytical and numerical solutions
 
plot(time,N_uranium, 'r' , t, N_analytical, 'b') 
grid()
xlabel("Time in years") 
ylabel("Number of atoms")