function Frame = GetRSTFrame(p,v) 

p = p/norm(p); 
v = v/norm(v); 

h = cross(p,v); 
h = h/norm(h); 

h_cross_p = cross(h,p); 
h_cross_p = h_cross_p/norm(h_cross_p); 

Frame(1,:) = p; 
Frame(2,:) = h_cross_p; 
Frame(3,:) = h; 
%Frame = [p,h_cross_p,h]; 
%Frame = Frame';

end