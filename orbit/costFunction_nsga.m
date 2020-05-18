%% MOLTO-OR Software 

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

%% Cost Function 

function [cost,cons] = costFunction_nsga(x,data) 

% x = globalParam to be optimized by NSGA

inputs = data.inputs; 
initialState_coe = data.initialState_coe; 
targetState_coe = data.targetState_coe;
global param 
global allData  
startDate = inputs.transferStartJulianDate;
if inputs.problem_type == 2 % Optimization problem 
    initO = initialState_coe.o; 
    targetO = targetState_coe.o; 
    
    % Selecting all the variables which are going to be introduced to the
    % optimizer
    
weights(1) = x(1);  
weights(2) = x(2); 
weights(3) = x(3); 
weights(4) = x(4); 
weights(5) = x(5); 

etaCutAbsolute = x(6);
etaCutRelative = x(7);

scalingFunc_m = x(8);
scalingFunc_n = x(9);
scalingFunc_r = x(10); 

% Free Vars 
fv = 10; 
if inputs.freeVars(1) == 1 
    initialState_coe.e = x(fv+1); 
    fv = fv + 1;
end
if inputs.freeVars(2) == 1 
    initialState_coe.i = x(fv+1); 
    fv = fv + 1;
end

if inputs.freeVars(3) == 1
    initialState_coe.o = x(fv+1); 
    fv = fv + 1;
end

if inputs.freeVars(4) == 1 
    initialState_coe.w = x(fv+1); 
    fv = fv + 1;
end
    
if inputs.freeVars(5) == 1 
    targetState_coe.e = x(fv+1); 
    fv = fv + 1;
end

if inputs.freeVars(6) == 1 
    targetState_coe.i = x(fv+1); 
    fv = fv + 1;
end

if inputs.freeVars(7) == 1 
    targetState_coe.o = x(fv+1);
    fv = fv + 1;
end

if inputs.freeVars(8) == 1 
    targetState_coe.w = x(fv+1); 
    fv = fv + 1;
end 

if inputs.freeVars(9) == 1 
    startDate = x(fv+1); 
end 
   
% Optimization variable in RAAN  % ATENTO A ESTO 

%initialState_coe.o = initO; 
%targetState_coe.o = targetO; 
else 
    
    weights(1:5) = inputs.weights; 
    scalingFunc_m = inputs.scale_func_m;
    scalingFunc_n = inputs.scale_func_n;
    scalingFunc_r = inputs.scale_func_r;
    etaCutAbsolute = inputs.eta_cut_abs;
    etaCutRelative = inputs.eta_cut_rel; 
        
end
    
% Calculate the normal for all weights in Q-law to avoid convergence
% criteria becomes tighter

weightNormal = 0; 

 
for i = 1:5 
    weightNormal = weightNormal + weights(i)^2; 
end
weightNormal = sqrt(weightNormal); 

% Set Q-law params 

%q_law_params = struct([]); 
q_law_params.weights = inputs.weights; 
q_law_params.minRpAltitude = inputs.min_rp_threshold;  
q_law_params.scalingFunc_m = scalingFunc_m;
q_law_params.scalingFunc_n = scalingFunc_n;
q_law_params.scalingFunc_r = scalingFunc_r; 
q_law_params.weightNormal = weightNormal;
q_law_params.etaCutAbsolute = etaCutAbsolute;
q_law_params.etaCutRelative = etaCutRelative;
% Transform from COE to EOE  

initialState_eoe = COE_getEquinoctialElements(initialState_coe); 
targetState_eoe = COE_getEquinoctialElements(targetState_coe); 

% Set initial condition vector 
initialCondition = zeros(7,1); 

initialCondition(1) = initialState_eoe.p; 
initialCondition(2) = initialState_eoe.f; 
initialCondition(3) = initialState_eoe.g; 
initialCondition(4) = initialState_eoe.h; 
initialCondition(5) = initialState_eoe.k; 
initialCondition(6) = initialState_eoe.L; 
initialCondition(7) = inputs.initMass; 


%parameters = struct([]); 
parameters.p_target_p = targetState_eoe.p;
parameters.p_target_f = targetState_eoe.f;
parameters.p_target_g = targetState_eoe.g;
parameters.p_target_h = targetState_eoe.h;
parameters.p_target_k = targetState_eoe.k;
parameters.p_maxThrust = inputs.maxThrust; 
parameters.p_alpha = 0; 
parameters.p_beta = 0; 
parameters.p_thrust = 0; 
parameters.p_total_dV = 0; 
parameters.p_dt = 30; 
parameters.p_Isp = inputs.Isp;
parameters.p_eta_absolute = 0; 
parameters.p_eta_relative = 0; 
parameters.p_julianDate = startDate; 
parameters.p_eclipse = -1; 
parameters.p_Qvalue = compute_Q(initialState_eoe,targetState_eoe,inputs,q_law_params); % CAMBIAR LOS INPUTS 
parameters.p_error = 0; 
parameters.p_retake_step = -1; 


% Start performance timer 
time_solver = timer('StartFcn',@(~,~)disp('timer started.'),'TimerFcn',@(~,~)disp(rand(1)));
start(time_solver)

% Start numerical integrator 

param = parameters; 
% ode45() + apaños 

% vars = size(initialCondition); 
% order = 7; 
% step = zeros(1,vars + 1)'; 
% 
% for i = 1:order 
%     k(i) = zeros(1,vars) 
% end
% step(1) = 0; 
% for n = 1:vars 
%     step(n+1) = initialCondition(n)
% end

% Define stop functions 
global stepsNearConvergence 
stepsNearConvergence = 0;
global stepSimilar 
stepSimilar = 0;
global lastEoe 
lastEoe.p = 1; 
global counter
counter  = 0;
stop_fun = @(time,eoe) controlFunction_nsga(time, eoe, inputs, q_law_params);
options = odeset('RelTol',1e-6,'AbsTol',1e-6); 
options.Events = stop_fun; 
options.Refine = 1; 
%options.MaxStep = param.p_dt;

% Define Derivatives for ode45 
derivs = @(time,eoe) DiffEqn(eoe,inputs) ; 

% Call to integrator 
[times,states] = ode45(derivs,[0,inputs.max_transfer_time],initialCondition,options); 

% lo de backward integration que de momento no lo hago 

% Transform solution into usable data. Meter todo en la estructura
% trayectory

solverDataSize = length(times); 

trajectory = struct([]); 

for j = 1:solverDataSize 
    trajectory(j).Ti = times(j); 
    trajectory(j).Eoe.p = states(j,1); 
    trajectory(j).Eoe.f = states(j,2); 
    trajectory(j).Eoe.g = states(j,3);
    trajectory(j).Eoe.h = states(j,4); 
    trajectory(j).Eoe.k = states(j,5); 
    trajectory(j).Eoe.L = states(j,6); 
    trajectory(j).Mass = states(j,7); 
    

    trajectory(j).Thrust = allData(j).p_thrust;
    trajectory(j).DeltaV = allData(j).p_total_dV;

    
    trajectory(j).Qvalue = compute_Q(trajectory(j).Eoe,targetState_eoe,inputs,q_law_params);

    trajectory(j).condition = allData(j).condition;
    
end

% Recalculate deltaV from solution. 

trajectory(1).DeltaV = 0;

for ii = 2 : length(trajectory)
    trajectory(ii).DeltaV = trajectory(ii-1).DeltaV + (trajectory(ii).Ti - trajectory(ii-1).Ti)*trajectory(ii-1).Thrust/trajectory(ii-1).Mass;
end

% Save last state 
currentState = trajectory(end); 

% If single transfer, print some useful information 

if inputs.problem_type == 1  % if single transfer 
    fprintf('Total dV = %5.4f [km/s] \n',currentState.DeltaV*10^-3); 
    fprintf('Total fuel Mass = %f [kg] \n',trajectory(1).Mass - currentState.Mass); 
    fprintf('Total Time = %f [days]\n',currentState.Ti/3600/24); 
    fprintf('Q value = %f [-] \n',currentState.Qvalue); 
%    fprintf('Total CPU time = %f',time_solver); 
end 

cost = zeros(1,2);  

% If converged, return cost function 
Qvalue = currentState.Qvalue; 
if Qvalue < 0.1 %inputs.convergence_threshold*weightNormal 
    cost(2) = abs(currentState.DeltaV/1000); % [km/s] 
    cost(1) = abs(currentState.Ti/86400); % days 
else
    % if no convergence, set a higher value
    cost(1) = 10^10; 
    cost(2) = 10^10; 
end

cons = 0; 
clear param; 
clear allData;
end


    
    
    
