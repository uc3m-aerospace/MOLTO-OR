%% Compute D1,D2,D3 

function [D1,D2,D3] = computeD1D2D3_code(current,target,inputs,q_law_params)


% Declare q_law_params
Wa = q_law_params.weights(1); 
Wf = q_law_params.weights(2);
Wg = q_law_params.weights(3);
Wh = q_law_params.weights(4);
Wk = q_law_params.weights(5);

na = q_law_params.scalingFunc_n;
ma = q_law_params.scalingFunc_m;
ra = q_law_params.scalingFunc_r;  

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
%w = 1 + f*cos(L) + g*sin(L); 
%r = p/w; 
ss = 1 + h^2 + k^2;  

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

% Compute P: low perigee restrictition 
Wp = 1; 
rp_min =  q_law_params.minRpAltitude;
kp = 450; 
%rp = p/(1+COE.e); 
%P = exp(kp*(1 - rp/rp_min)); 


% Start calculation
t = []; 

t(1) = f^2; 
t(2) = g^2; 
t(3) = 1 - t(1) - t(2); 
t(4) = t(1) + t(2); 
t(5) = 1/sqrt(t(4));
t(4) = t(4)*t(5); 
t(6) = 1 + t(4);
t(7) = 1/rp_min; 
t(6) = 1/t(6); 
t(8) = a*t(3); 
t(9) = d_a; 
t(10) = abs(t(9));
t(11) = (t(10)/(ma*a_target))^na; 
t(12) = 1 + t(11); 
t(13) = 1/ra;
t(14) = (t(12))^t(13);
t(4) = 1 - t(4); 
t(15) = d_f;
t(16) = d_g; 
t(17) = d_h; 
t(18) = ss; 
t(2) = 1 - t(2); 
t(19) = 1/sqrt(t(2));
t(2) = -t(2)*t(19) + f;
t(20) = d_k; 
t(1) = 1 - t(1); 
t(21) = 1/sqrt(t(1)); 
t(1) = -t(1)*t(21) + g; 
t(22) = 1/a; 
t(23) = 1/t(3); 
t(24) = 1/t(18); 
t(25) = (t(20))^2; 
t(26) = (t(17))^2; 
t(27) = (t(2))^2;
t(28) = (t(1))^2; 
t(29) = Wk*t(25); 
t(30) = t(28)*t(29); 
t(31) = Wh*t(26); 
t(32) = t(27)*t(31); 
t(33) = (t(24))^2;
t(34) = (t(22))^2;
t(35) = (t(9))^2; 
t(36) = (4*t(30) + 4*t(32))*t(33); 
t(37) = Wg*(t(16))^2;
t(38) = Wf*(t(15))^2; 
t(39) = t(34)*t(6); 
t(40) = t(22)*inputs.mu; 
t(41) = t(40)*(t(23)*(t(36) + t(38)/4 + t(37)/4) + t(39)*t(14)*Wa*t(35)*t(4)/4);
t(42) = Wp*exp(kp*(1 - t(8)*t(6)*t(7))); 
t(43) = 1 + t(42); 
t(12) = 1/t(12); 
t(10) = 1/t(10); 
t(42) = kp*t(42); 
t(44) = 1/inputs.mu; 
t(45) = cos(L);
t(46) = sin(L); 
t(47) = f*t(45);
t(48) = g*t(46); 
t(49) = 1 + t(47) + t(48);
t(50) = t(3)*t(6)*t(5) + 2; 
t(51) = f*t(6); 
t(52) = f*t(23); 
t(53) = t(4)*(t(6))^2;
t(5) = t(14)*Wa*t(5)*t(34)*t(35); 
t(15) = t(43)*t(40)*(t(23)*(t(29)*(8*t(33)*f*t(1)*t(21) + 8*t(28)*t(52)*t(33)) + t(52)*(t(38)/2 + t(37)/2 + 8*t(32)*t(33)) + Wf*t(15)/2 + 8*t(31)*t(33)*t(2)) + (-t(53)*f/4 - t(51)/4)*t(5)) + t(42)*t(51)*a*t(7)*t(50)*t(41);
t(8) = t(8)*t(44);
t(8) = sqrt(t(8));
t(21) = 2 + t(47) + t(48); 
t(47) = g*t(6);
t(48) = g*t(23);
t(1) = (t(43)*t(40)*(t(23)*(t(31)*(8*t(33)*g*t(2)*t(19) + 8*t(27)*t(48)*t(33)) + t(48)*((8*t(30)*t(33)) + t(38)/2 + t(37)/2) + Wg*t(16)/2 + (8*t(29)*t(33)*t(1))) + (-t(53)*g/4 - t(47)/4)*t(5)) + t(42)*t(47)*a*t(7)*t(50)*t(41));
t(2) = 1/t(49); 
t(3) = (a*(t(43)*t(34)*inputs.mu*(t(23)*(-t(36) - t(38)/4 - t(37)/4) + Wa*((-3/4*t(39) + t(13)*t(11)*na*abs(t(9))/t(9)*t(10)*t(12)*t(22)*t(6)/4)*t(4)*t(14)*t(35) + t(14)*t(9)*t(22)*t(6)*t(4)/2)) - t(42)*t(3)*t(6)*t(7)*t(41))/sqrt(t(3))*sqrt((a*t(44))));
t(2) = t(8)*t(2); 
t(4) = t(40)*t(33)*t(23); 

D1 = 2*t(3)*t(49) + t(2)*(t(15)*(t(21)*t(45) + f) + t(1)*(t(21)*t(46) + g));
D2 = t(8)*(-t(1)*t(45) + t(15)*t(46)) + 2*t(3)*(f*t(46) - g*t(45));
D3 = t(2)*(t(43)*t(18)*(t(45)*t(4)*(-(16*t(24)*t(30)*h) + Wh*(-16*t(27)*h*t(24)*t(26) + 8*t(17)*t(27))) + t(46)*t(4)*(-16*t(24)*t(32)*k + (Wk*(-16*t(28)*k*t(24)*t(25) + 8*t(20)*t(28)))))/2 +(f*t(1) - g*t(15))*(h*t(46) - k*t(45)));








    case 0 
        
%% Not considering minimum perigee restriction  

t = []; 

t(1) = d_a; 
t(2) = abs(t(1)); 
t(3) = (t(2)/(ma*a_target))^na; 
t(4) = 1 + t(3); 
t(5) = 1/ra; 
t(6) = t(4)^t(5); 
t(7) = g^2; 
t(8) = f^2; 
t(9) = t(8) + t(7); 
t(10) = 1/sqrt(t(9)); 
t(9) = t(9)*t(10); 
t(11) = 1 + t(9); 
t(9) = 1 - t(9); 
t(12) = d_f; 
t(13) = 1 - t(8) - t(7); 
t(14) = d_g;
t(15) = d_h;
t(16) = 1 + k^2 + h^2; 
t(7) = 1 - t(7); 
t(17) = 1/sqrt(t(7)); 
t(7) = -t(7)*t(17) + f; 
t(18) = d_k;
t(8) = 1 - t(8);
t(19) = 1/sqrt(t(8));
t(8) = -t(8)*t(19) + g; 
t(2) = 1/t(2); 
t(11) = 1/t(11); 
t(20) = 1/t(13); 
t(21) = 1/t(16); 
t(4) = 1/t(4); 
t(22) = 1/a; 
t(23) = t(21)^2; 
t(24) = t(15)^2; 
t(25) = t(18)^2; 
t(26) = t(22)^2; 
t(27) = t(7)^2; 
t(28) = t(8)^2; 
t(29) = t(1)^2; 
t(30) = t(25)* Wk; 
t(31) = t(28)*t(30);
t(32) = t(24)*Wh;
t(33) = t(27)*t(32);
t(34) = t(12)^2*Wf;
t(35) = t(14)^2*Wg;
t(36) = a/inputs.mu; 
t(37) = cos(L);
t(38) = sin(L); 
t(39) = f*t(37); 
t(40) = g*t(38); 
t(41) = 1 + t(39) + t(40); 
t(42) = t(26)*t(9)*t(11)^2; 
t(43) = t(26)*t(11); 
t(10) = Wa*t(6)*t(10)*t(29); 
t(44) = t(22)*inputs.mu; 
t(12) = t(44)*(t(20)*(t(30)*(8*t(23)*f*t(8)*t(19) + 8*t(28)*t(20)*f*t(23)) + (t(34)/2 + t(35)/2 + 8*t(33)*t(23))*t(20)*f + Wf*t(12)/2 + 8*t(32)*t(23)*t(7)) + (-t(42)*f/4 - t(43)*f/4)*t(10));
t(19) = t(36)*t(13); 
t(19) = sqrt(t(19)); 
t(39) = 2 + t(39) + t(40); 
t(7) = (t(44)*(t(20)*(t(32)*(8*t(23)*g*t(7)*t(17) + 8*t(27)*t(20)*g*t(23)) + ((8*t(31)*t(23)) + t(34)/2 + t(35)/2)*t(20)*g + Wg*t(14)/2 + (8*t(30)*t(23)*t(8))) + (-t(42)*g/4 - t(43)*g/4)*t(10)));
t(8) = 1/t(41); 
t(8) = t(8)*t(19); 
t(1) = (a*t(26)*inputs.mu*(t(20)*(-((4*t(31)) + 4*t(33))*t(23) - t(34)/4 - t(35)/4) + Wa*((t(5)*t(3)*na*abs(t(1))/t(1)*t(2)*t(4)*t(22)/4 - 3/4*t(26))*t(11)*t(9)*t(6)*t(29) + t(6)*t(1)*t(22)*t(11)*t(9)/2))/sqrt(t(13))*sqrt(t(36)));
t(2) = h*t(21); 
t(3) = t(44)*t(23)*t(20); 
t(4) = k*t(21);

D1 = 2*t(1)*t(41) + t(8)*(t(12)*(t(37)*t(39) + f) + t(7)*(t(38)*t(39) + g));
D2 = t(19)*(t(12)*t(38) - t(37)*t(7)) + 2*t(1)*(f*t(38) - g*t(37)); 
D3 = t(8)*(t(16)*(t(3)*(-(16*t(31)*t(2)) + Wh*(-16*t(27)*t(2)*t(24) + 8*t(15)*t(27)))*t(37)/2 + t(3)*(-16*t(33)*t(4) + Wk*(-16*t(28)*t(4)*t(25) + (8*t(18)*t(28))))*t(38)/2) + (h*t(38) - k*t(37))*(f*t(7) - g*t(12)));

end
end