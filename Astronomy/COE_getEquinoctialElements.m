function [EOE] = COE_getEquinoctialElements(COE)

% Calculate equinoctial elements 

if sin(COE.i) >=0
    I = 1;
else
    I = -1;
end

EOE.p = COE.a*(1 - COE.e^2); 
EOE.f = COE.e*cos(I*COE.o + COE.w); 
EOE.g = COE.e*sin(I*COE.o + COE.w);
EOE.h = tan(COE.i/2)^I*cos(COE.o); 
EOE.k = tan(COE.i/2)^I*sin(COE.o);

% Limit true anomaly to [0,2Pi]
EOE.L = wrapTo2Pi(I*COE.o + COE.w + COE.nu);   

end

