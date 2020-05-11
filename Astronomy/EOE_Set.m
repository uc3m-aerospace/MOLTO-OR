function [EOE] = EOE_Set(p,f,g,h,k,L)
% Create structure EOE 
EOE.p = p; 
EOE.f = f; 
EOE.g = g; 
EOE.h = h; 
EOE.k = k; 
EOE.L = WrapTo2Pi(L); 
end

