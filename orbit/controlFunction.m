%% MOLTO-OR Sftware 

% This program is developed at Universidad Carlos III de Madrid as part 
% a Master Thesis. 

% The software and its components are developed by Sergio de Vera Muñoz,
% tutored by David Morante. 

% 2020, Madrid.

%% Orbit determination

% Transfer scenario planning and siinputslation for near-Earth orbits.
% Set of functions designed to be ready for global optimization by evolutionary
% algorithms. Optimal control based on different forinputslation of
% Q-law algorithm. 

% Contains some of the core features of the software. 

%% Control Function 

function [condition,converged,direction] = controlFunction(time, eoe, inputs, q_law_params) 

global counter
global param 
% global lastEoe
% global stepsNearConvergence 
% global stepSimilar 
% Initialize variables 
t = time; 
m = eoe(7); 

current = struct; 
current.p = eoe(1); 
current.f = eoe(2);
current.g = eoe(3);
current.h = eoe(4);
current.k = eoe(5);
current.L = eoe(6); 

target = struct; 
target.p = param.p_target_p;
target.f = param.p_target_f; 
target.g = param.p_target_g; 
target.h = param.p_target_h; 
target.k = param.p_target_k; 
target.L = 0; 

maxThrust = param.p_maxThrust; 

% Auxiliary terms for further computations 
w = 1 + current.f*cos(current.L) + current.g*sin(current.L); 
r = current.p/w; 

% Compute alpha beta and Qdot
[alpha, beta] =  compute_alpha_beta (current, target, inputs,q_law_params);
QdotCurrent = compute_Qdot_min(current,target,current.L,inputs,q_law_params);  

% Check inertial thrusting direction constraints
inertialThrustDir = [cos(beta)*sin(alpha), cos(beta)*cos(alpha), sin(beta)];


% Get state vector from equinoctial elements 

position = EOE_getPosition(current);
velocity = EOE_getVelocity(current,inputs);  

% Perform rotations
transferOmegaMod = 0; % Hardcoded, es un input 
position(1) = cos(-transferOmegaMod)*position(1) - sin(-transferOmegaMod)*position(2);
position(2) = sin(-transferOmegaMod)*position(1) + cos(-transferOmegaMod)*position(2); 

velocity(1) = cos(-transferOmegaMod)*velocity(1) - sin(-transferOmegaMod)*velocity(2);
velocity(2) = sin(-transferOmegaMod)*velocity(1) + cos(-transferOmegaMod)*velocity(2); 

rstFrame = GetRSTFrame(position,velocity);   

% en realidad no se usan 
multiplicador = inertialThrustDir;
inertialThrustDir(1) = dot(rstFrame(:,1)',multiplicador);
inertialThrustDir(2) = dot(rstFrame(:,2)',multiplicador);
inertialThrustDir(3) = dot(rstFrame(:,3)',multiplicador); 

% Get unit vector pointing to Earth 
EarthRelVec = -position;  
EarthRelVec = EarthRelVec/norm(EarthRelVec);

% Get unit vector pointing Sun
sunPos = getSolarPosition(param.p_julianDate + t/86400); 

% If integrating backwards, flip alpha and beta angles 

%referenceDeltaL = inputs.ref_delta_L; 

%if referenceDeltaL < 0
 %   alpha = alpha + pi; 
  %  beta = -beta; 
%end

% Calculate current Q-value 
Qvalue = compute_Q(current,target,inputs,q_law_params); 

% BACKWARD INTEGRATION 

% lastEoe = struct; 
% lastEoe.p = 100000000; % Para inicializarla 
% if Qvalue < 10^4
%     stepsNearConvergence = stepsNearConvergence + 1; 
% else
%     if abs(current.p - lastEoe.p) < 0.01
%         stepSimilar = stepSimilar + 1; 
%     end
%     lastEoe = current; 
% end 
        
etaCutAbsolute = q_law_params.etaCutAbsolute; % From cost Function
etaCutRelative = q_law_params.etaCutRelative; % From cost Function  

if (etaCutAbsolute > 10^-6 && etaCutRelative > 10^-6) 
     [QdotMin, L_opt_min] = compute_Qdot_extrema(current,target,inputs,q_law_params); 
     [QdotMax, L_opt_max] = compute_Qdot_extrema_max(current,target,inputs,q_law_params); 
     
     % Compute effectivity parameters 
     etaAbs = min(1,max(0,QdotCurrent/QdotMin)) - etaCutAbsolute; 
     etaRel = min(1,max(0,(QdotCurrent - QdotMax)/(QdotMin - QdotMax)));
else
    % P1 problem, always thrusting 
    etaAbs = 1.0; 
    etaRel = 1.0;
end
%     step = 0; % hardcodeao CAMBIAAAAR
%     step = step + 1; % Esto hay q meterle input
    
    
% If Q is small enough, super-convergence is achieved, reduce step size to
% avoid overshoot
referenceDeltaL = 0;
if Qvalue > 10^4 
    if referenceDeltaL > 0
        deltaL = max(0.01*referenceDeltaL, (log10(Qvalue) + 1)/5*referenceDeltaL);
    else
        deltaL = min(0.01*referenceDeltaL,(log10(Qvalue) + 1)/5)*referenceDeltaL;
    end
else
    % Delfault delta L 
    deltaL = referenceDeltaL;
end

% Optimal dt from delta L 
dt = deltaL * sqrt(r^3/inputs.mu)*sqrt(w)/(1 + sqrt(current.f^2 + current.g^2));


% Eclipse calculations 

Eclipse_enabled = inputs.eclipse_enabled; % hardcoded
if Eclipse_enabled == 1 
     eclipse = getEclipseCondition(sunPos, position, inputs.Re); 
else
    eclipse = -1;
end

% Establish switch of thruster 
% Determine whether to thrust or coast based on the switch function value
switchFunc = min(min(etaAbs,etaRel),-eclipse); 

thrust = switchFunc > 0; % boolean operator 

% If Q has reached an enough low value, stop the integrator 

%converged = 0; % boolean operator 

%convergenceRadius = inputs.convergence_threshold; 
%weightNormal = 0;  % both hardcoded  

condition1 = Qvalue < (0.1); %convergenceRadius
if condition1
    thrust = 0; 
    %converged = 1; 
end

% Thrusting angles 

param.p_alpha = alpha; 
param.p_beta = beta; 

% If thrusting, apply maximum thrust , else coast 
if thrust 
    param.p_thrust = maxThrust; 
else 
    param.p_thrust = 0; 
end

eclipseThrustLevel = 0;
if (min(min(etaAbs,etaRel),1) && eclipse > 0 && eclipseThrustLevel > 10^-6)
    param.p_thrust = eclipseThrustLevel * maxThrust; 
end 

% Delta V and dt 

param.p_total_dV = param.p_total_dV + param.p_thrust/m*dt; 
param.p_dt = dt; 

% Effectivity values 
param.p_eta_absolute = etaAbs;
param.p_eta_relative = etaRel; 

param.p_eclipse = eclipse; 

param.p_Qvalue = Qvalue; 


% Conditions to stop integration 
condition = Qvalue - (0.1); 
converged = 1; 
direction = 0;


% Store values of parameters for each time step
counter = counter + 1 ; 
% global tiempo 
% tiempo(counter) = t; 
global allData

    allData(counter).p_alpha =  param.p_alpha; 
    allData(counter).p_beta = param.p_beta; 
    allData(counter).p_thrust = param.p_thrust; 
    allData(counter).p_total_dV = param.p_total_dV;
    allData(counter).p_dt = param.p_dt;
    allData(counter).p_eta_absolute = param.p_eta_absolute;
    allData(counter).p_eta_relative = param.p_eta_relative; 
    allData(counter).p_julianDate = param.p_julianDate;
    allData(counter).p_eclipse = param.p_eclipse; 
    allData(counter).p_error = param.p_error; 
    allData(counter).condition = condition1;
    
    
end