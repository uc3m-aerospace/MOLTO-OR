%% MOLTO-OR Sftware 

% This program is developed at Universidad Carlos III de Madrid as part 
% a Master Thesis. 

% The software and its components are developed by Sergio de Vera Muñoz,
% tutored by David Morante. 

% 2020, Madrid.

%% Astronomy

% Transformation between coordinates

function [velocity] = EOE_getVelocity(EOE,param)

% Calculate velocity
vel = sqrt(param.mu/EOE.p); 


% Position vector 

velocity = vel*EOE_getVelocityUnitVector(EOE); 

end

