clear ; 
length = 1.0 ; 
g = 9.81 ; 
q = 0.5 ; 

F_Drive = 1.2 ; 
Omega_D = 1/3 ; 
npoints = 15000 ; 
dt = 0.04 ;

omega = zeros(npoints,1) ; 
theta = zeros(npoints,1) ; 
time =  zeros(npoints,1) ; 

theta(1) = 0.2 ; 
omega(1) = 0 ; 

for step = 1 : npoints -1 
    omega(step+1) = omega(step) + (-(g/length) * sin (theta(step)) - q*omega(step) + ...
                    F_Drive * sin(Omega_D * time(step))) *dt ; 
    
    temporary_theta_step_plus_1 = theta(step) + omega(step + 1 ) * dt ;
    
    if(temporary_theta_step_plus_1 < -pi ) 
        temporary_theta_step_plus_1 = temporary_theta_step_plus_1 + 2*pi ; 
    elseif(temporary_theta_step_plus_1 > pi ) 
        temporary_theta_step_plus_1 = temporary_theta_step_plus_1 - 2 * pi ;
    end
    
    theta(step+1) = temporary_theta_step_plus_1 ; 
    time(step+1) = time(step) + dt ; 
end

plot(time,omega,'r') ;
xlabel("time(seconds)") ; 
ylabel("theta(radians)") ;

