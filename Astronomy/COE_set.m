function [COE] = COE_set(a,e,i,o,w,nu)
% Summary of this function goes here
%   Detailed explanation goes here
COE.a = a; 
COE.e = e; 
COE.i = i; 
COE.o = wrapTo2Pi(o); 
COE.w = wrapTo2Pi(w);
COE.nu = wrapTo2Pi(nu); 

end

