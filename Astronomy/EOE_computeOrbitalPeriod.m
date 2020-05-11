function [T] = EOE_computeOrbitalPeriod(EOE)

mu = 3.986004418*10^8; %km3s-2

a =  EOE.p / (1.0 - EOE.f * EOE.f - EOE.g * EOE.g);

T = 2 * pi * sqrt (a*a*a/mu);

end

