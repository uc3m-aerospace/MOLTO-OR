%% MOLTO-OR Sftware 

% This program is developed at Universidad Carlos III de Madrid as part 
% a Master Thesis. 

% The software and its components are developed by Sergio de Vera Muñoz,
% tutored by David Morante. 

% 2020, Madrid.

%% Orbit determination

% Transfer scenario planning and simulation for near-Earth orbits.
% Set of functions designed to be ready for global optimization by evolutionary
% algorithms. Optimal control based on different formulation of
% Q-law algorithm. 

% Contains some of the core features of the software. 

%% Eclipse condition 

function eclipse = getEclipseCondition(solarVector, testPosition, radius)

if dot(solarVector,testPosition) > 0 
    eclipse = -1;
end
shadowVector = -solarVector;
shadowVector = shadowVector * dot(shadowVector,testPosition); 

radialDistanceVector = testPosition - shadowVector; 


eclipse = 1 - norm(radialDistanceVector)^2/radius^2;

end