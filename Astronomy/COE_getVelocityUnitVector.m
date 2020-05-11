function [vel_unit] = COE_getVelocityUnitVector(COE)
% Pre-calculate trigonometric functions 

cos_o = cos(COE.o); 
sin_o = sin(COE.o); 
cos_w_nu = cos(COE.w + COE.nu); 
sin_w_nu = sin(COE.w + COE.nu); 
cos_i = cos(COE.i); 
sin_i = sin(COE.i);
e = COE.e; 
w = COE.w; 

% Position unit vector 

vel_unit = zeros (1,3);

vel_unit(1) = -(cos_o * (sin_w_nu + e*sin (w)) + sin_o*cos_i * (cos_w_nu + e*cos (w)));
vel_unit(2) = -(sin_o * (sin_w_nu + e*sin (w)) - cos_o*cos_i * (cos_w_nu + e*cos (w)));
vel_unit(3) = (sin_i * (cos_w_nu + e*cos (w)));
end

