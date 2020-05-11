%% MOLTO-OR Software 

% This program is developed at Universidad Carlos III de Madrid as part 
% a Master Thesis. 

% The software and its components are developed by Sergio de Vera Muñoz,
% tutored by David Morante. 

% 2020, Madrid. 

%% Set bounds 

function bounds = set_bounds(data,inputs,max)

for i = 1:5 
    bounds(i) = data.weights(i);
end

bounds(6) = data.eta_abs; 
bounds(7) = data.eta_rel; 
bounds(8) = data.m; 
bounds(9) = data.n; 
bounds(10) = data.r; 

fv = 1; 

if max == 1  % Set maximums 
    fvbound = [1,2*pi]; 
    date = inputs.maxStartJD;
elseif max == 0 % Set minimums
    fvbound = [0, 0]; 
    date = inputs.transferStartJulianDate;
end


for j = 1:8 
    if inputs.freeVars(j) == 1 
        if j == 1 || j == 5 
            bounds(10 + fv) = fvbound(1);
        else
            bounds(10 + fv) = fvbound(2);  
        end
        fv = fv + 1; 
    end
end 

if inputs.freeVars(9) == 1 
    bounds(10 + fv) = date; 
end

end