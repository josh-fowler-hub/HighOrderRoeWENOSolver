% Main driver for the HighOrderRoeWENOSolver
%
% This script sets the configuration for one or more test cases and then
% delegates the execution to `runSolver`, which is implemented in
% `src/runSolver.m`.

close all;
addpath(genpath(fullfile(pwd, 'src')));

% NOTE: Some integrators (e.g. SSPRK65) currently use a global counter.
%       This is a small legacy artifact and generally does not affect
%       correctness for standard runs.
global time_int_counter;
time_int_counter = 0;

%% Explicit Inputs - Change These
% Choose which test cases to run (see testCase.m)
test_num = [
%             1, % Modified Sod Test 1: Toro Ch. 11 Thesis: Test 1
%             2, % 123 Test Test 2: Toro Ch. 11
            3, % Left Blast Wave, Test 3: Toro Ch. 11 Thesis: Test 2
%             4, % Shock Collision Test 4: Toro Ch. 11 Thesis: Test 3
%             5, % Stationary Contact Wave Test 5: Toro Ch. 11 Thesis: Test 4
%             6, % Sod's Original Problem Thesis: Test 5
%             7, % Left Expansion and right strong shock Thesis: Test 6
%             8, % Right Expansion and left strong shock Thesis: Test 7
%             9, % Double Shock Thesis: Test 8
%             10, % Double Expansion 
%             11, % Cavitation 
%             12, % Shocktube problem of G.A. Sod, JCP 27:1, 1978 Thesis: Test 9
%             13, % Lax test case: M. Arora and P.L. Roe: JCP 132:3-11, 1997 Thesis: Test 10
%             14, % Mach 3 test case: M. Arora and P.L. Roe: JCP 132:3-11, 1997
%             15, % Shocktube problem with supersonic zone Thesis: Test 11
%             16, % Stationary shock Thesis: Test 12
%             17, % right side of 2-d riemman case 14 Thesis: Test 13
%             18, % right side of 2-d riemman case 15 Thesis: Test 14
%             19, % User Specified Test
            ];

% Select the Riemann solver(s) to run
riemann_solver = {
                  'Roe', % Roe's Riemann Solver
%                   'Roe-Pike', % Roe-Pike Riemann Solver
%                   'HLL', % HLL Riemann Solver
%                   'HLLC', % HLLC Riemann Solver
%                   'Osher', % Osher Riemann Solver
                  };

% Configure the reconstruction scheme
recon = {         % Reconstruction Method
         'WENO5', % 5th Order WENO Reconstruction
         };

% Options for MUSCL schemes (only used if recon == 'MUSCL')
flux_limiter = [ 1 ]; % Minmod
lambda = [1];        % Scaling factor
kappa = [1];         % MUSCL kappa value

% Time integration configuration
time_int_method = { 'SSPRK5-6' };

% Domain / grid settings
cells = 105;
xstart = 0;
xend = 1;
cfl = 1;

% Optional plotting control
extra_prim_var_plot = false; % plots additional primitive fields (a, M, H, s)
plot_in_time = false;        % plot at each time step (slow)

%% Run
config = struct();
config.test_num = test_num;
config.riemann_solver = riemann_solver;
config.flux_limiter = flux_limiter;
config.lambda = lambda;
config.CD_Term_Order = [1];
config.recon = recon;
config.kappa = kappa;
config.time_int_method = time_int_method;
config.cells = cells;
config.xstart = xstart;
config.xend = xend;
config.cfl = cfl;
config.extra_prim_var_plot = extra_prim_var_plot;
config.plot_in_time = plot_in_time;

results = runSolver(config);

% Display final error norms for the last configuration
if ~isempty(results)
    last = results(end);
    L2 = L2Norm(last.U(:,last.xidx), last.Uexact, config.cells);
    Linf = LinfNorm(last.U(:,last.xidx), last.Uexact);
    fprintf('Final L2 norm: %g\n', L2);
    fprintf('Final Linf norm: %g\n', Linf);
end
