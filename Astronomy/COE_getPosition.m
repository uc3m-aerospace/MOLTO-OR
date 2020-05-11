function [position] = COE_getPosition(COE)
% Calculate radius 

r = COE.a * (1.0 - COE.e * COE.e) / (1.0 + COE.e * cos (COE.nu)); 

position = r*COE_getPositionUnitVector(COE);
end

