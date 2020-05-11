%% MOLTO-OR Software 

% This program is developed at Universidad Carlos III de Madrid as part 
% a Master Thesis. 

% The software and its components are developed by Sergio de Vera Muñoz,
% tutored by David Morante. 

% 2020, Madrid.

%% Astronomy

% Transformation between coordinates

function [pos_unit] = EOE_getPositionUnitVector(EOE)

% Compute composite functions

ss = 1 + EOE.h^2 + EOE.k^2;
as = EOE.h^2 - EOE.k^2;
cosL = cos(EOE.L);
sinL = sin(EOE.L); 

% Position Unit Vector

pos_unit = zeros(1,3); 
pos_unit(1) = 1.0 / ss * (cosL + as*cosL + 2.0 * EOE.h * EOE.k*sinL); 
pos_unit(2) = 1.0 / ss * (sinL - as*sinL + 2.0 * EOE.h * EOE.k*cosL);
pos_unit(3) = 2.0 / ss * (EOE.h*sinL - EOE.k*cosL);


end

