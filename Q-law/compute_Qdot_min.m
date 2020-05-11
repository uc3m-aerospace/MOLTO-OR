function [Qdotmin] = compute_Qdot_min(current,target,L_override,inputs,q_law_params)


current.L = L_override; 
[D1,D2,D3] = computeD1D2D3_code(current,target,inputs,q_law_params); 

Qdotmin = -sqrt(D1^2 + D2^2 + D3^2);  
%Qdotmax = sqrt(D1*2 + D2^2 + D3^2);  

end  



