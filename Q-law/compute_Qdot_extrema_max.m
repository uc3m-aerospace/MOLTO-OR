%% MOLTO-OR Software 

% This program is developed at Universidad Carlos III de Madrid as part 
% a Master Thesis. 

% The software and its components are developed by Sergio de Vera Muñoz,
% tutored by David Morante. 

% 2020, Madrid.

%% Q-law 
% This script is used to compute the optimal thrusting angles for
% a Lyapunov control method used for low-thrust orbit raising 

%% Qdot extrema 

function [fb, L_abscissa] = compute_Qdot_extrema_max(current,target,inputs,q_law_params)

b=0; 
fb = compute_Qdot_min(current,target,b,inputs,q_law_params); 

res = 16;
for i=1:res

    L_ = i/res*2*pi; 
    fx = compute_Qdot_min(current,target,L_,inputs,q_law_params); 
    
    if fx>fb 
        b = L_; 
        fb = fx;
    end
end
     
a = b - 2*pi/res; 
c = b + 2*pi/res; 

R = sqrt(5)/2 - 0.5; 
b = a + R*(c - a); 

init_a = a; init_c = c; 

fa = compute_Qdot_min(current,target,a,inputs,q_law_params); 
fc = compute_Qdot_min(current,target,c,inputs,q_law_params); 

iter = 0;  
maxiter = 5;

tol = 1; 
while (tol > 10^-3 && iter<maxiter) 
    
   x = b - 0.5*((fb - fc)*(b - a)^2 - (fb - fa)*(b - c)^2)/((b - a)*(fb - fc) - (b - c)*(fb - fa)); 
   fx = compute_Qdot_min(current,target,L_,inputs,q_law_params); 
   
   if fx > fb 
       
       if x < b 
          c = b;
          fc = fb;
       else
          a = b; 
          fa = fb;
       end
       tol = abs(b - x); 
       b = x; 
       fb = fx;
       
   else
       
       if x < b
           a = x; 
           fa = fx;
       else
           c = x; 
           fc = fx; 
       end
       tol = 1; 
   end
   
   iter = iter + 1; 
       
end
L_abscissa = b; 
L_abscissa = L_abscissa - 2*pi*floor(L_abscissa/2/pi);

c1 = compute_Qdot_min(current,target,b - 10^-3,inputs,q_law_params); 
c2 = compute_Qdot_min(current,target,b + 10^-3,inputs,q_law_params); 

if (c1 > fb  || c2 > fb) 
    [fb, L_abscissa] = compute_Qdot_extrema_golden(current,target,inputs,init_a,init_c,q_law_params) ;
    golden = 1; 
end

parab = 1;
end