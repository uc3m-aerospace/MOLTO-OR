%% MOLTO-OR Software 

% This program is developed at Universidad Carlos III de Madrid as part 
% a Master Thesis. 

% The software and its components are developed by Sergio de Vera Muñoz,
% tutored by David Morante. 

% 2020, Madrid. 
%% Call_NSGA 
% This script creates the necessary variables to run the Genetic Algorithm

function result = call_nsga(nsga,bounds_max,bounds_min,data)

% Creating default options structure for NSGA
options_nsga = nsgaopt(); 
options_nsga.popsize = nsga.max_pop;
options_nsga.maxGen = nsga.max_gen; 
options_nsga.mutationfraction = nsga.mutation_rate;
options_nsga.crossoverFraction = nsga.crossover_rate; 
options_nsga.numObj = 2; 
options_nsga.numVar = 10 + nsga.numFreeVar; 
options_nsga.numCons = 0;
options_nsga.vartype = linspace(1,1,options_nsga.numVar); 
options_nsga.ub = set_bounds(bounds_max,data.inputs,1);
options_nsga.lb = set_bounds(bounds_min,data.inputs,0); 
%nsga.parallel = 1;
if nsga.parallel == 1
options_nsga.useParallel = 'yes'; 
end
options_nsga.objfun = @costFunction_nsga; 
options_nsga.plotInterval = 5; 
options_nsga.nameObj={'Time of Flight [days]','$\Delta V$ [km/s]'}; 

result = nsga2(options_nsga,data);
end