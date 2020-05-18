function [Q] = compute_Q(current,target,inputs,q_law_inputss)

%F = 2.2e-1; % Thrust magnitude CAMBIAR!!!!!!

% Declare weights 
Wa = q_law_inputss.weights(1) ; 
Wf = q_law_inputss.weights(2);
Wg = q_law_inputss.weights(3);
Wh = q_law_inputss.weights(4);
Wk = q_law_inputss.weights(5);

ma = q_law_inputss.scalingFunc_m;
na = q_law_inputss.scalingFunc_n;
ra = q_law_inputss.scalingFunc_r;  % Tengo que ponerlos como control
% Declare current variables 
p = current.p; 
f = current.f; 
g = current.g; 
h = current.h; 
k = current.k; 
L = current.L;  
% Computing semimajor axis 
COE = EOE_getClassicalElements(current); 
a = COE.a;             %current
[tar] = EOE_getClassicalElements(target); 
a_target = tar.a;      %target

% Auxiliary variables 
% w = 1 + f*cos(L) + g*sin(L); 
% r = p/w; 
% ss = 1 + h^2 + k^2;  

% Difference current-target 
d_a = a - a_target; 
d_f = f - target.f; 
d_g = g - target.g; 
d_h = h - target.h; 
d_k = k - target.k;  

min_perigee = inputs.min_rp_enabled; 
switch min_perigee 
    case 1
        %% Consider minimum perigee restriction 
        rp_min = q_law_inputss.minRpAltitude;  
        kp = 450; % Este inputsetro que coño hago con el 
        Wp = 1; 
%         rp = p/(1+COE.e); 
%         P = exp(kp*(1 - rp/rp_min)); 

        
        t=[];
        t(1) = f^2; 
        t(2) = g^2; 
        t(3) = -t(1) -t(2) + 1; 
        t(6) = sqrt(t(1) + t(2)); 
        t(8) = 1/(1 + t(6)); 
        t(14) = exp(kp*(1 - a*t(3)*t(8)/rp_min)); 
        t(17) = d_a; 
        t(18) = abs(t(17)); 
        t(23) = (t(18)/(ma*a_target))^(na); 
        t(26) = (1 + t(23))^(1/ra); 
        t(28) = t(17)^2; 
        t(30) = a^2; 
        t(40) = d_f^2; 
        t(42) = 1/a; 
        t(43) = 1/t(3); 
        t(45) = inputs.mu*t(42)*t(43); 
        t(49) = d_g^2; 
        t(54) = d_h^2; 
        t(57) = t(43)*inputs.mu; 
        t(58) = h^2; 
        t(59) = k^2; 
        t(61) = (t(58) + t(59) + 1)^2; 
        t(62) = 1/t(61); 
        t(64) = sqrt(-t(2) + 1); 
        t(66) = (-t(64) + f)^2;  % OJO
        t(72) = d_k^2; 
        t(76) = sqrt(-t(1) + 1); 
        t(78) = (-t(76) + g)^2;  % OJO
        
        Q = (t(14)*Wp + 1)*(t(26)*Wa*t(28)/t(30)/a*inputs.mu*t(8)*(1 - t(6))/4 + Wf*t(40)*t(45)/4 + Wg*t(49)*t(45)/4 + 4*Wh*t(54)*t(42)*t(57)*t(62)*t(66) + 4*Wk*t(72)*t(42)*t(57)*t(62)*t(78));
        
        
    case 0
        %% Without perigee restriction
        
        t = []; 
        t(1) = d_a; 
        t(2) = abs(t(1)); 
        t(7) = (t(2)/(ma*a_target))^na;
        t(10) = (1 + t(7))^(1/ra); 
        t(12) = (t(1))^2;
        t(14) = a^2; 
        t(18) = f^2; 
        t(19) = g^2; 
        t(21) = sqrt(t(18) + t(19));
        t(30) = d_f^2; 
        t(32) = 1/a; 
        t(34) = 1/(1 - t(18) - t(19)); 
        t(36) = t(32)*t(34)*inputs.mu; 
        t(40) = d_g^2; 
        t(45) = d_h^2; 
        t(48) = t(34)*inputs.mu; 
        t(49) = h^2; 
        t(50) = k^2; 
        t(52) = (t(49) + t(50) + 1)^2; 
        t(53) = 1/t(52); 
        t(55) = sqrt(-t(19) + 1); 
        t(57) = (-t(55) + f)^2;  % CUIDAO CON ESTE
        t(63) = d_k^2; 
        t(67) = sqrt(-t(18) + 1); 
        t(69) = (-t(67) + g)^2;  % CUIDAO CON ESTE
        
        Q = t(10)*Wa*t(12)/t(14)/a*inputs.mu/(1 + t(21))*(1 - t(21))/4 + Wf*t(30)*t(36)/4 + Wg*t(40)*t(36)/4 + 4*Wh*t(45)*t(32)*t(48)*t(53)*t(57) + 4*Wk*t(63)*t(32)*t(48)*t(53)*t(69); 
        
end

end