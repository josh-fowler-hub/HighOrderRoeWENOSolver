close all;
addpath(genpath(fullfile(pwd, 'src')));

cells = 1000;
L = 1;
gamma = 1.4;
ghost = 3;
cfl = 0.9;
flux = 1;
lam = 1;
kap = 1;
rec = 'WENO5';
CD_Term_Order = 1;
riemann_solver = 'Roe';

xs = linspace(0, 1, cells);
ys = zeros(size(xs));
zs = zeros(size(xs));
dx = xs(2);

[W, W0] = manufacturedSolutionSetUp(L, xs, ys, zs);

Uexact = primToCons(gamma, W(1,:), W(2,:), W(5,:));
U = primToCons(gamma, W0(1,:), W0(2,:), W0(5,:));
[ULBC, URBC] = manufacturedSolutionBCs(U);
U = [ULBC, U, URBC];

% Compute exact solution with BCs once before the loop
[ULBC_exact, URBC_exact] = manufacturedSolutionBCs(Uexact);
Uexact_BCs = [ULBC_exact, Uexact, URBC_exact];

steps = 0;
time = 0;
t = 0.0052;

while time < t
    steps = steps + 1;
    Uold = U;
    
    % Calculate the time step
    dt = calculateTimeStep(cfl, dx, time, t, U, gamma);
    time = time + dt;

    % Select the time integration method
    U = SSPRK65(U,flux,lam,gamma,dx,dt,kap,rec,CD_Term_Order,riemann_solver,xs,time);
    
    % Apply boundary conditions
    [ULBC, URBC] = manufacturedSolutionBCs(U);
    U = [ULBC, U, URBC];
    
    error = max(max(abs(U - Uold)));
    if error < 1e-3
        break
    end
end

plotVars(U,Uexact,xs);
