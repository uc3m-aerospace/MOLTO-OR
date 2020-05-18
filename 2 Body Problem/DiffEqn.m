%% MOLTO-OR Software 

% This program is developed at Universidad Carlos III de Madrid as part 
% a Master Thesis. 

% The software and its components are developed by Sergio de Vera Muñoz,
% tutored by David Morante. 

% 2020, Madrid.

%% Two Body Problem  

function derivs = DiffEqn(state,inputs) 

global param

% derivs(1) = dp/dt; 
% derivs(2) = df/dt; 
% derivs(3) = dg/dt;
% derivs(4) = dh/dt; 
% derivs(5) = dk/dt; 
% derivs(6) = dL/dt; 
% derivs(7) = dm/dt; 


% Required parameters 

mu = inputs.mu; 
J2 = inputs.J2; 
R_e = inputs.Re; 
g0 = inputs.gravity;  % añadir al archivo input

J2_enabled = inputs.J2_enabled; 

% Equinoctial elements
p = state(1); 
f = state(2); 
g = state(3); 
h = state(4); 
k = state(5); 
L = state(6); 
m = state(7); 

% Thrusting Parameters
alpha = param.p_alpha;  % Azimuth angle
beta = param.p_beta;    % Declination angle
F = param.p_thrust;          % Thrust magnitude 
Isp = param.p_Isp;      % Specific Impulse 

% Parameters for differential equations 

w = 1 + f*cos(L) + g*sin(L); 
r = p/w; 
v = sqrt(p/mu); 
ss = 1 + h^2 + k^2; %s^2

% Thrust components 

Fr = F/m*cos(beta)*sin(alpha); 
Ft = F/m*cos(beta)*cos(alpha); 
Fn = F/m*sin(beta); 

if J2_enabled == 1 % To consider zonal harmonics perturbation
    
   C = mu*J2*R_e^2/r^4; 
   
   Fr = Fr - 1.5*C*(1 - 12*(h*sin(L) - k*cos(L))^2)/ss^2;
   Ft = Ft - 12*C*(h*sin(L) - k*cos(L))*(h*cos(L) + k*sin(L))/ss^2;
   Fn = Fn - 6*C*(1 - h^2 - k^2)*(h*sin(L) - k*cos(L))/ss^2; 
   
end




derivs(1) = 2*v*p/w*Ft; 
derivs(2) = v*sin(L)*Fr + v/w*((w + 1)*cos(L) + f)*Ft - v/w*g*(h*sin(L) - k*cos(L))*Fn;
derivs(3) = -v*cos(L)*Fr + v/w*((w + 1)*sin(L) + g)*Ft + v/w*f*(h*sin(L) - k*cos(L))*Fn;
derivs(4) = v*ss*cos(L)/(2*w)*Fn; 
derivs(5) = v*ss*sin(L)/(2*w)*Fn; 
derivs(6) = sqrt(p*mu)*w^2/p^2 + v/w*(h*sin(L) - k*cos(L))*Fn;
derivs(7) = -F/(Isp*g0); 

derivs = derivs';
end