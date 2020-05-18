%% MOLTO-OR Software 

% This program is developed at Universidad Carlos III de Madrid as part 
% a Master Thesis. 

% The software and its components are developed by Sergio de Vera Muñoz,
% tutored by David Morante. 

% 2020, Madrid.

%% Q-law 
% This script is used to compute the optimal thrusting angles for
% a Lyapunov control method used for low-thrust orbit raising 

%% Compute alpha and beta angles 

function [alpha, beta] = compute_alpha_beta (current, target, inputs,q_law_params) 

 [D1 D2 D3] = computeD1D2D3_code(current,target,inputs,q_law_params); 
 
 alpha = atan2(-D2,-D1); 
 beta = atan(-D3/sqrt(D1^2 + D2^2));
 
end