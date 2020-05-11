function [velocity] = COE_getVelocity(COE,inputs)

% Calculate velocity length 

v = sqrt(inputs.mu/(COE.a *(1-COE.e*COE.e))); % Revisar esta ecuacion

% Velocity vector 

velocity = v * COE_getVelocityUnitVector(COE); 
end

