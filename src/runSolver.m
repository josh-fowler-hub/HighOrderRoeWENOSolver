function results = runSolver(config)
% runSolver Runs the 1D Euler solver with a set of configuration options.
%
%   results = runSolver(config) executes the numerical solver for each test
%   case and configuration provided in the config struct. It mirrors the
%   original logic found in the legacy main.m script, but isolates the
%   execution logic in a reusable function.
%
%   CONFIG STRUCT FIELDS (recommended):
%     - test_num           : vector of test indices (see testCase.m)
%     - riemann_solver     : cell array of solver names {'Roe','HLL',...}
%     - flux_limiter       : vector of limiter IDs (see getFluxName.m)
%     - lambda             : vector of lambda values (for MUSCL)
%     - CD_Term_Order      : vector of central differencing orders
%     - recon              : cell array of reconstruction methods
%     - kappa              : vector of kappa values (MUSCL parameter)
%     - time_int_method    : cell array of time integration method names
%     - cells              : number of grid cells
%     - xstart, xend       : domain bounds
%     - cfl                : CFL number
%     - extra_prim_var_plot: bool, whether to compute extra primitive fields
%     - plot_in_time       : bool, whether to plot during time stepping
%
%   The function returns a results struct that currently contains the last
%   conserved state computed for each configuration.

% Preallocate results for traceability (optional)
results = struct();
resultIndex = 0;

% Counter to provide a consistent plotting style across solver runs
k = 0;

for testidx = 1:numel(config.test_num)
    test = config.test_num(testidx);

    % Initialize initial condition for the selected test case
    [t,rhoR,rhoL,uL,uR,PR,PL,R,gamma,middle] = testCase(test);

    % Set up the computational grid
    [xs, dx, n2, xidx] = setUpGrid(config.cells, 3, config.xstart, config.xend);

    % Set up initial conservative variables
    U = setUpU(rhoR,rhoL,uL,uR,PR,PL,gamma,xs,middle,n2);

    % Get the exact solution for plotting/comparison
    [~,~,~,~,Uexact] = Riemann(rhoL,uL,PL,rhoR,uR,PR,t,middle,R,gamma,config.cells,config.extra_prim_var_plot);

    % Loop through all solver configurations
    for solverIdx = 1:numel(config.riemann_solver)
        solver = config.riemann_solver{solverIdx};
        for Ord = 1:numel(config.CD_Term_Order)
            Order = config.CD_Term_Order(Ord);
            for count = 1:numel(config.recon)
                rec = config.recon{count};
                ghost = 3; % number of ghost cells for boundary conditions

                for int_count = 1:numel(config.time_int_method)
                    time_int = config.time_int_method{int_count};

                    flux_count = 0;
                    for flux_int = 1:numel(config.flux_limiter)
                        flux = config.flux_limiter(flux_int);
                        flux_lim = getFluxName(flux);
                        flux_count = flux_count + 1;

                        lam_count = 0;
                        for lam = config.lambda
                            lam_count = lam_count + 1;

                            % Skip original Roe method the second time if
                            % lambda is set 1 & 2
                            if flux == 18 && lam_count > 1
                                continue
                            end

                            % Skip extra lambda runs if not using MUSCL
                            if ~strcmp(rec, 'MUSCL') && (lam_count > 1 || flux_count > 1)
                                continue
                            end

                            % Increment plot line counter (used for consistent plotting)
                            k = k + 1;

                            % Reset the conservative variables for each run
                            U = setUpU(rhoR,rhoL,uL,uR,PR,PL,gamma,xs,middle,n2);

                            % Time stepping
                            time = 0;

                            steps = 0;
                            while time < t
                                steps = steps + 1;

                                % Apply boundary conditions (zero-gradient)
                                U = applyBoundaryConditions(U, ghost, xidx);

                                % Compute next time step
                                dt = calculateTimeStep(config.cfl, dx, time, t, U, gamma);
                                time = time + dt;

                                % Advance solution using selected time integration
                                switch time_int
                                    case 'FE'
                                        U = forwardEuler(U,flux,lam,gamma,dx,dt,config.kappa,rec,Order,solver,xs,time);
                                    case 'RK2'
                                        U = RK2(U,flux,lam,gamma,dx,dt,config.kappa,rec,Order,solver,xs,time);
                                    case 'RK3'
                                        U = RK3(U,flux,lam,gamma,dx,dt,config.kappa,rec,Order,solver,xs,time);
                                    case 'RK4'
                                        U = RK4(U,flux,lam,gamma,dx,dt,config.kappa,rec,Order,solver,xs,time);
                                    case 'RK4-5'
                                        U = RK45(U,flux,lam,gamma,dx,dt,config.kappa,rec,Order,solver,xs,time);
                                    case 'SSPRK4-10'
                                        U = SSPRK104(U,flux,lam,gamma,dx,dt,config.kappa,rec,Order,solver,xs,time);
                                    case 'SSPRK5-6'
                                        U = SSPRK65(U,flux,lam,gamma,dx,dt,config.kappa,rec,Order,solver,xs,time);
                                    otherwise
                                        error('Unknown time integration method: %s', time_int);
                                end

                                % Make a plot at each time step if requested
                                if config.plot_in_time
                                    clf(1,'reset');
                                    clf(2,'reset');
                                    clf(3,'reset');
                                    clf(4,'reset');
                                    if config.extra_prim_var_plot
                                        clf(5,'reset');
                                        clf(6,'reset');
                                        clf(7,'reset');
                                        clf(8,'reset');
                                    end

                                    makePlots(U,xs,gamma,flux,lam,rec,Order,xidx,k,config.extra_prim_var_plot,solver);
                                    drawnow;
                                    cla([1,2,3,4],'reset');
                                    disp(['time = ', num2str(time)]);
                                end
                            end

                            % Plot the Results at the end of the run
                            if ~config.plot_in_time
                                makePlots(U,xs,gamma,flux,lam,rec,Order,xidx,k,config.extra_prim_var_plot,solver);
                                configurePlots(rec,test_num,cells,cfl,config.extra_prim_var_plot,solver,Order,config.plot_folder);
                            end
                            if config.plot_in_time
                                Riemann(rhoL,uL,PL,rhoR,uR,PR,t,middle,R,gamma,config.cells,config.extra_prim_var_plot);
                                makePlots(U,xs,gamma,flux,lam,rec,Order,xidx,k,config.extra_prim_var_plot,solver);
                            end

                            % Store results for this configuration (optional)
                            resultIndex = resultIndex + 1;
                            results(resultIndex).test = test;
                            results(resultIndex).solver = solver;
                            results(resultIndex).recon = rec;
                            results(resultIndex).fluxLimiter = flux_lim;
                            results(resultIndex).lambda = lam;
                            results(resultIndex).U = U;
                            results(resultIndex).Uexact = Uexact;
                            results(resultIndex).xs = xs;
                            results(resultIndex).xidx = xidx;
                        end
                    end
                end
            end
        end
    end
end
end
