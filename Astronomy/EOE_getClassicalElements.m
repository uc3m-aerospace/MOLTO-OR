function [COE] = EOE_getClassicalElements(EOE)
EOE.L = wrapTo2Pi(EOE.L); 

COE = struct; 

COE.a = EOE.p / (1.0 - EOE.f*EOE.f - EOE.g*EOE.g);
COE.e = sqrt (EOE.f*EOE.f + EOE.g*EOE.g);
COE.i = 2.0 * atan (sqrt(EOE.h*EOE.h + EOE.k*EOE.k));
COE.w = atan2 (EOE.g, EOE.f) - atan2 (EOE.k, EOE.h); 
COE.w = wrapTo2Pi(COE.w);
COE.o = atan2 (EOE.k, EOE.h);
COE.o = wrapTo2Pi(COE.o);
COE.nu = EOE.L - (COE.o + COE.w); 
COE.nu = wrapTo2Pi(COE.nu); 

%COE.E = 2.0 * atan2 (sqrt(1.0-COE.e)*sin(COE.nu/2.0), sqrt(1.0+COE.e)*cos(COE.nu/2.0));
%COE.M = COE.E - COE.E*cos(COE.E); 
end

