function [EOE,state] = STATE_calcAllFromCoe(COE,param)

COE.w = wrapTo2Pi(COE.w); 
COE.o = wrapTo2Pi(COE.o); 
COE.nu = wrapTo2Pi(COE.nu); 

EOE = COE_getEquinoctialElements(COE); 

pos = COE_getPosition(COE); 
vel = COE_getVelocity(COE,param); 

state.pos = pos; 
state.vel = vel;
end

