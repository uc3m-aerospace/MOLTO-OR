%% MOLTO-OR Software 

% This program is developed at Universidad Carlos III de Madrid as part 
% a Master Thesis. 

% The software and its components are developed by Sergio de Vera Muñoz,
% tutored by David Morante. 

% 2020, Madrid.

%% Q-law 
% This script is used to compute the optimal thrusting angles for
% a Lyapunov control method used for low-thrust orbit raising 

%% Qdot extrema golden
function [f_extr, L_abscissa] =  compute_Qdot_extrema_golden(current,target,inputs,a,c,q_law_params) 

R = sqrt(5) - 0.5;  resphi = 1 - R; 

b = a + R*(c - a); 

x0 = a; x3 = c;
tol = 10^-4;  % tolerance for convergence

x1 = b; 
x2 = b + resphi*(c - b); 

f1 = compute_Qdot_min(current,target,x1,inputs,q_law_params);
f2 = compute_Qdot_min(current,target,x2,inputs,q_law_params);

iter = 0; maximum = 10; 
while ( abs(x3 - x0) > tol*(abs(x1) + abs(x2)) && iter<maximum) 
    
    if f2<f1 
        x0 = x1; 
        x1 = x2; 
        x2 = R*x2 + resphi*x3;  
        
        f1 = f2; 
        f2 =  compute_Qdot_min(current,target,x2,inputs,q_law_params);
        
    else % (f2>f1) 
        x3 = x2; 
        x2 = x1; 
        x1 = R*x1 + resphi*x0; 
        
        f2 = f1; 
        f1 =  compute_Qdot_min(current,target,x1,inputs,q_law_params);
    end
    
    iter = iter + 1; 
        
end

    if f1 < f2
        x_extr = x1; 
        f_extr = f1; 
        
    else
        x_extr = x2; 
        f_extr = f2;
    end
    
    L_abscissa = x_extr; 
    
    L_abscissa = L_abscissa - 2*pi*floor(L_abscissa/2/pi);
    
    
end
