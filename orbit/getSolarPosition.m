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

%% Sun Position 

function sunPos = getSolarPosition(JD)

% JD2000
n = JD - 2451545; 

%The mean longitude of the Sun, corrected for the aberration of light
L = (280.46 + 0.9856474*n)/57.2957795131; 
L = wrapTo2Pi(L);

% The mean anomaly of the Sun
g = (357.528 + 0.9856003*n)/57.2957795131; 
g = wrapTo2Pi(g); 

% the ecliptic longitude of the Sun 
lambda = L + (1.915*sin(g) + 0.02*sin(2*g))/57.2957795131; 

%  obliquity of the ecliptic
epsilon = (23.439 - 4e-7*n)/57.2957795131; 

% Direction of sun 

sunPos = [cos(lambda),cos(epsilon)*sin(lambda),sin(epsilon)*sin(lambda)];

end

