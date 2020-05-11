function [pos_unit] = COE_getPositionUnitVector(COE)

% Pre-calculate trigonometric functions 

cos_o = cos(COE.o); 
sin_o = sin(COE.o); 
cos_w_nu = cos(COE.w + COE.nu); 
sin_w_nu = sin(COE.w + COE.nu); 
cos_i = cos(COE.i); 
sin_i = sin(COE.i);

% Position unit vector 

pos_unit = zeros (1,3);
pos_unit(1) = cos_o*cos_w_nu - sin_o*cos_i*sin_w_nu; 
pos_unit(2) = sin_o*cos_w_nu + cos_o*cos_i*sin_w_nu;
pos_unit(3) = sin_i*sin_w_nu; 


end

