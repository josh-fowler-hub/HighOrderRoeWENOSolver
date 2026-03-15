function [r, u, p, a, H, e, m, s] = consToPrim(U, gamma)
% Convert conservative variables to primitive variables for compressible Euler equations
% 
% INPUTS:
%   U = [rho; rho*u; rho*E] - conservative variables (3 x N array)
%   gamma - heat capacity ratio (typically 1.4 for air)
%
% OUTPUTS:
%   r - density
%   u - velocity  
%   p - pressure (potentially floored for numerical stability)
%   a - sound speed
%   H - specific enthalpy
%   e - specific internal energy
%   m - Mach number
%   s - entropy
%
% KNOWN ISSUES:
%   Pressure can become zero or negative due to:
%   - Floating-point roundoff errors when (E - 0.5*rho*u^2) is very small
%   - Non-conservative numerical schemes
%   - CFL violations in certain flow regions
%   - Vacuum or near-vacuum regions in test cases
%
%   Current handling: Floor negative/zero pressures to 1e-6
%   Better approach: Investigate root cause in numerical scheme

r = U(1,:);                       % Density
u = U(2,:)./U(1,:);               % Velocity
E = U(3,:);                       % Total Energy per Unit Volume

% Compute internal energy: e = (E - 0.5*rho*u^2)/rho
internal_energy_times_rho = E - 0.5*r.*u.^2;  % rho*e
p = (gamma-1).*internal_energy_times_rho;      % Pressure formula

% Detect and report problematic cells
pressure_violations = find(p <= 0);
if ~isempty(pressure_violations)
    % For debugging: Document and monitor pressure violations
    % Uncomment for detailed diagnostics:
    % msg = sprintf('consToPrim: %d cells with p <= 0 (max violation: %.2e)', ...
    %     length(pressure_violations), min(p(pressure_violations)));
    % warning(msg);
    
    % Floor pressure to maintain numerical stability
    p(p <= 0) = 1e-6;
else
    % Pressure is physical in all cells - good sign
    % (Pressure should be positive in Euler equations)
end

% Sound speed computation with same flooring for consistency
e = (1/(gamma-1))*(p./r);         % Internal Energy
a_temp = gamma*p./r;
a_temp(a_temp <= 0) = 1e-6;
a = sqrt(a_temp);                 % Sound speed

m = u./a;                         % Mach number
s = log(p./r.^gamma);             % Entropy
H = (E + p)./r;                   % Enthalpy
end