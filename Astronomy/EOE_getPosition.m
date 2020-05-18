%% MOLTO-OR Software 

% This program is developed at Universidad Carlos III de Madrid as part 
% a Master Thesis. 

% The software and its components are developed by Sergio de Vera Muñoz,
% tutored by David Morante. 

% 2020, Madrid.

%% Astronomy

% Transformation between coordinates

function [position] = EOE_getPosition(EOE)

% Calculate velocity

w = 1 + EOE.f*cos(EOE.L) + EOE.g*sin(EOE.L); 
rad = EOE.p/w; 

% Position vector 

position = rad*EOE_getPositionUnitVector(EOE); 

end

