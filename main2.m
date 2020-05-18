%% MOLTO-OR Software 

% This program is developed at Universidad Carlos III de Madrid as part 
% a Master Thesis. 

% The software and its components are developed by Sergio de Vera Muñoz,
% tutored by David Morante. 

% 2020, Madrid. 

%% Main script - Multi Objective Low Thrust Optimization - Orbit Raising
% This script loads user inputs and run the program



function[cost,trajectory,result] =  main2(inputs,initialState_coe,targetState_coe,nsga,app)

% File to save results 
filename = 'David_03.mat';

% Computing gravity of the central body
G = 6.6740831e-11;
inputs.gravity = G*inputs.body_mass/inputs.Re^2; 
% Convert Problem Type input
value = inputs.problem_type;
if value == "Single Transfer"
       inputs.problem_type = 1; 
elseif value == "Optimization " 
       inputs.problem_type = 2;
end

% Convert Calendar Date to Julian Date
t1 = datetime(inputs.transferStartJulianDate); 
inputs.transferStartJulianDate = juliandate(t1);   


% Creating data structure 
data.inputs = inputs; 
data.initialState_coe = initialState_coe; 
data.targetState_coe = targetState_coe;
data.app = app;

% Choose type of problem to be solved
problem_type = inputs.problem_type; 
if problem_type == 1 % single transfer
    fprintf('Single Transfer') 
    [cost,trajectory] = costFunction(0,data);  
    result = 0; 
    save(filename,'cost','trajectory')
elseif problem_type == 2 % Optimization
    
    % Set bounds 
bounds_min.weights(1:5) = [1,0,0,0,0]; 
bounds_min.eta_abs = 0; 
bounds_min.eta_rel = 0; 
bounds_min.m = 0; 
bounds_min.n = 0; 
bounds_min.r = 0; 

bounds_max.weights(1:5) = [100,100,100,100,100]; 
bounds_max.eta_abs = 0.8; 
bounds_max.eta_rel = 0.8; 
bounds_max.m = 10; 
bounds_max.n = 10; 
bounds_max.r = 10; 
    
    
    
   % Calculate number of free variables 
   numFreeVar = 0; 
   for fv = 1:length(inputs.freeVars) 
       numFreeVar = numFreeVar + inputs.freeVars(fv);
   end
    nsga.numFreeVar = numFreeVar; 
    
    % Convert from Calendar to JD if Start Date is to be optimized
    if inputs.freeVars(9) == 1
    t2 = datetime(inputs.maxStartJD); 
    inputs.maxStartJD  = juliandate(t2); 
    end
    
    % Call Genetic Algorithm 
    result = call_nsga(nsga,bounds_max,bounds_min,data); 
    
    for i = 1: nsga.max_pop
    [cost(i,:),trajectory{i}] = costFunction(result.pops(nsga.max_gen,i).var,data);
    end
    % donde el 10 iria nsga.max_gen 
    save(filename,'cost','trajectory','result')
end
end