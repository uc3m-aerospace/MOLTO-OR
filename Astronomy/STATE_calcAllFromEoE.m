function [COE,state] = STATE_calcAllFromEoE(EOE,param)

EOE.L = wrapTo2Pi(EOE.L); 
COE = EOE_getClassicalElements(EOE); 
[~,state] = STATE_calcAllFromCoe(COE,param);
end

