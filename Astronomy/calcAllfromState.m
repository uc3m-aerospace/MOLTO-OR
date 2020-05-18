function [COE,EOE] = calcAllfromState(State,param)

% State.pos = Position vector (cartesian)
% State.vel = Velocity vector (cartesian) 
mu = param.mu; %km^3/s^2
pos = State.pos; %km 
vel = State.vel; %km/sv
eps=1e-10 ; 

%% 1. Norm of R and V and radial velocity
r=norm(pos);         %distance
v=norm(vel);         %velocity
vr=dot(pos,vel)/r;   %radial velocity
%% Momemtum
H=cross(pos,vel);    %momemtum vector
h=norm(H);           %momemtum 
%% Inclination
COE.i=acos(H(3)/h);  
%% Node line and right ascension of the ascendent node angle
N=cross([0 0 1],H);   %Node line vector
n=norm(N);            %Node
if n~=0
    COE.o= acos(N(1)/n);    %rad
    if N(2)<0
        COE.o=2*pi-acos(N(1)/n); %rad
    end
else
    COE.o=0;
end
%% Eccentricity and argument of periapses w
E=1/mu*((v^2-mu/r)*pos-r*vr*vel);
COE.E = E; 
e=norm(E);
COE.e = e;
if n~=0
   if e > eps
   w = acos(dot(N,E)/n/e); %rad
      if E(3) < 0
          w = 2*pi -acos(dot(N,E)/n/e) ;
      end
   else
      w = 0;
   end
else
w = 0;
end 
COE.w = w; 
%% True anomaly nu
if e>eps
    nu=acos(dot(E,pos)/(e*r));
    if vr<0
        nu=2*pi-acos(dot(E,pos)/(e*r));
    end
else 
    cp=cross(N,R);
    if cp(3)>=0
        nu=acos(dot(N,pos)/(n*r));
    else
        nu=2*pi-acos(dot(N,pos)/(n*r));
    end
end
COE.nu = nu;
%% Semimajor axis.
COE.a=h^2/(mu*(1-e^2));
EOE = struct;

EOE = COE_getEquinoctialElements(COE);
end

