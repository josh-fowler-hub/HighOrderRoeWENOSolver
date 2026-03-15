function [r, u, p, a, H, e, m, s] = cons_to_prim(U, gamma)
r = U(1,:);                       % Density
u = U(2,:)./U(1,:);               % Velocity
E = U(3,:);                       % Total Energy per Unit Volume
p = (gamma-1).*(E - 0.5*r.*u.^2);
p(p==0) = 1e-6;                   % Pressure
p(p<0) = 1e-6;                    % Get rid of negetive Pressures
e = (1/(gamma-1))*(p./r);         % Internal Energy
a_temp = gamma*p./r;
a_temp(a_temp==0) = 1e-6;
a_temp(a_temp<0)=1e-6;
a = sqrt(a_temp);                 % sound speed
m = u./a;                         % Mach 
s = log(p./r.^gamma);             % Entropy
H = (E + p)./r;                   % Enthalpy
end