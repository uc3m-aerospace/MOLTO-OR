%% MOLTO-OR Sftware 

% This program is developed at Universidad Carlos III de Madrid as part 
% a Master Thesis. 

% The software and its components are developed by Sergio de Vera Muñoz,
% tutored by David Morante. 

% 2020, Madrid.

%% Astronomy

% Transformation between coordinates

function [vel_unit] = EOE_getVelocityUnitVector(EOE)

% Calculate composite functions 

ss = 1 + EOE.h^2 + EOE.k^2;
as = EOE.h^2 - EOE.k^2;
cosL = cos(EOE.L);
sinL = sin(EOE.L); 

vel_unit = zeros(1,3); 

vel_unit(1) = -1.0 / ss * (sinL + as*sinL - 2.0*EOE.h*EOE.k*cosL + EOE.g - 2.0*EOE.f*EOE.h*EOE.k + as*EOE.g);
vel_unit(2) = -1.0 / ss * (-cosL + as*cosL + 2.0*EOE.h*EOE.k*sinL - EOE.f + 2.0*EOE.g*EOE.h*EOE.k + as*EOE.f);
vel_unit(3) = 2.0 / ss * (EOE.h*cosL + EOE.k*sinL + EOE.f*EOE.h + EOE.g*EOE.k);
end

