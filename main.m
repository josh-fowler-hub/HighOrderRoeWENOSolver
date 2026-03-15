close all;
global time_int_counter;
time_int_counter = 0;
%% Explicit Inputs - Change These
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
%                 % NOTE: Roe Solver will not work for Test 2 or 11, as it
%                         fails for low density flows
%
riemann_solver = {
                  'Roe', % Roe's Riemann Solver
%                   'Roe-Pike', % Roe-Pike Riemann Solver
%                   'HLL', % HLL Riemann Solver (not yet added)
%                   'HLLC', % HLLC Riemann Solver (not yet added)
%                   'Osher', % Osher Riemann Solver (not yet added)
                  }; 
%                         % NOTE: 'Roe' and 'Roe-Pike' options should give
%                           same results. The only difference is in how the
%                           wave strengths are calculated, but the Roe-Pike
%                           solver should be faster.
% 
flux_limiter = [  % The flux limiter to use in the MUSCL Reconstruction:
                1, % minmod
%                 2, % van Leer
%                 3, % Barth-Jesperson
%                 4, % Superbee
%                 5, % van Albada 2
%                 6, % van Albada 1
%                 7, % CHARM
%                 8, % HCUS
%                 9, % HQUICK
%                 10, % Koren
%                 11, % MC
%                 12, % Osher
%                 13, % Ospre
%                 14, % SMART
%                 15, % Sweby
%                 16, % UMIST
%                 17, % Gen Minmod
                ];
%
lambda = [       % Scaling Factor
          1,
%           2,
          ];     % see: Quantification of numerical
%                  diffusivity due to TVD schemes in the advection 
%                  equation, Shreyas Bidadi, Sarma L. Rani ∗
%                  Journal of Computational Physics 261 (2014) 65–82
%
CD_Term_Order = [   % Central Differencing Term Order
                 1, % Roe
%                  2, % 2nd Order
%                  4, % 4th Order
%                  6, % 6th Order
                 ];
%
recon = {         % Reconstruction Method
%          'ROE',   % No Reconstruction
%          'MUSCL', % MUSCL Reconstruction
%          'WENO3', % 3rd Order WENO Reconstruction
         'WENO5', % 5th Order WENO Reconstruction
         };
%
kappa = [     % Choosing the MUSCL Reconstruction Scheme
         1,   % Central Difference
%          -1,  % Upwind
%          0,   % 2nd Order Fromm
%          1/2, % 3rd Order
%          1/3 % 3rd Order
         ];
%
time_int_method = {                % The method to use for time
                   'SSPRK5-6',   % SSP Runge-Kutta 5th Order 6 step
%                    'SSPRK4-10',  % SSP Runge-Kutta 4th Order 10 step
%                    'RK4-5',      % SSP Runge-Kutta 4th Order 5 step
%                    'RK4',        % SSP Runge-Kutta 4th Order
%                    'RK3',        % SSP Runge-Kutta 3rd Order
%                    'RK2',        % SSP Runge-Kutta 2nd Order
%                    'FE',         % Forward Euler
                  };
%
cells = 105;            % Number of Cells
xstart = 0;            % First x-Position
xend = 1;              % Last x-Position
cfl = 1;               % Courant-Friedrich-Lewy Number

extra_prim_var_plot = false; % if you want speed of sound, mach number,
%                              enthalpy, and entropy plotted use true else
%                              use false
plot_in_time = false; % You can use this to plot at each time step set to
%                       true to plot at each time step

%% Implicit Inputs - Don't Change These
for testidx = 1:length(test_num)
    test = test_num(testidx);
    % Initialize the Initial Conditions based on the test number
    [t,rhoR,rhoL,uL,uR,PR,PL,R,gamma,middle] = test_case(test);

    % Set up the computational domain
    [xs, dx, n2, xidx] = set_up_grid(cells, 3, xstart, xend);

    % Set up the Conservative Variable Array
    U = set_up_U(rhoR,rhoL,uL,uR,PR,PL,gamma,xs,middle,n2);

    %% Numerical Scheme

    % Start a counter to loop through all of the solver configurations
    k = 0;

    % Get and plot the Exact solution to the Riemann Test Problem
    [~,~,~,~,Uexact] = Riemann(rhoL,uL,PL,rhoR,uR,PR,t,middle,R,gamma,cells,extra_prim_var_plot);

    % Loop through the inputs
    for solver = riemann_solver                                             
        for Ord = 1:length(CD_Term_Order)
            Order = CD_Term_Order(Ord);
            for count = 1:length(recon)
                rec = recon{count};
                ghost = 3;
                for int_count = 1:length(time_int_method)
                    flux_count = 0;
                    for flux_int = 1:length(flux_limiter)
                        flux = flux_limiter(flux_int);
                        flux_lim = get_flux_name(flux);
                        flux_count = flux_count + 1;
                        lam_count = 0;
                        for lam = lambda
                            lam_count = lam_count + 1;

                            % Skip original Roe method the second time if
                            % lambda is set 1 & 2
                            if flux == 18 && lam_count > 1
                                continue
                            end

                            % Skip multiple runs for lambda if the MUSCL method
                            % is not used
                            if ~strcmp(rec, 'MUSCL') && (lam_count > 1 ||...
                                                         flux_count > 1)
                                continue
                            end

                            % Add one to the counter to use with next solver
                            % configuration
                            k = k+1;
                            for kap_int = 1:length(kappa)
                                kap = kappa(kap_int);
                                % Set up the Conservative Variable Arrays
                                U = set_up_U(rhoR,rhoL,uL,uR,PR,PL,gamma,xs,middle,n2);

                                % Set up the Boundary Conditions
                                time = 0;
                                left_BC_idx = 1:ghost;
                                right_BC_idx = max(xidx)+1:max(xidx)+ghost;
                                BC_idx = [left_BC_idx, right_BC_idx];

                                % Start a counter for time steps
                                steps = 0;
                                while time < t
                                    steps = steps + 1;

                                    % Apply Boundary Conditions
                                    RBC = U(:,max(xidx));
                                    LBC = U(:,ghost+1);
                                    URBC = zeros([3,length(right_BC_idx)]);
                                    ULBC = zeros([3,length(left_BC_idx)]);
                                    for l = 1:length(right_BC_idx)
                                        URBC(:,l) = RBC;
                                        ULBC(:,l) = LBC;
                                    end
                                    U(:, right_BC_idx) = URBC;
                                    U(:, left_BC_idx) = ULBC;

                                    % Calculate the time step
                                    dt = calculate_time_step(cfl, dx, time, t, U, gamma);
                                    time = time + dt;

                                    % Select the Time Integration Method based on
                                    % input
                                    switch time_int_method{int_count}
                                        case 'FE'
                                            U = Forward_Euler(U,flux,lam,gamma,dx,dt,...
                                                              kap,rec,Order,solver{1},xs,time);
                                        case 'RK2'
                                            U = RK2(U,flux,lam,gamma,dx,dt,kap,rec,Order,solver{1},xs,time);
                                        case 'RK3'
                                            U = RK3(U,flux,lam,gamma,dx,dt,kap,rec,Order,solver{1},xs,time);
                                        case 'RK4'
                                            U = RK4(U,flux,lam,gamma,dx,dt,kap,rec,Order,solver{1},xs,time);
                                        case 'RK4-5'
                                            U = RK45(U,flux,lam,gamma,dx,dt,kap,rec,Order,solver{1},xs,time);
                                        case 'SSPRK4-10'
                                            U = SSPRK104(U,flux,lam,gamma,dx,dt,kap,rec,Order,solver{1},xs,time);
                                        case 'SSPRK5-6'
                                            U = SSPRK65(U,flux,lam,gamma,dx,dt,kap,rec,Order,solver{1},xs,time);
                                    end

                                    % Make a plot at each time step if plot_in_time
                                    % is set to true
                                    if plot_in_time
                                        clf(1,'reset');
                                        clf(2,'reset');
                                        clf(3,'reset');
                                        clf(4,'reset');
                                        if extra_prim_var_plot
                                            clf(5,'reset');
                                            clf(6,'reset');
                                            clf(7,'reset');
                                            clf(8,'reset');
                                        end
                                        make_plots(U,xs,gamma,flux,lam,rec,Order,xidx,k,extra_prim_var_plot,riemann_solver);
                                        drawnow
            %                             pause(0.125);
                                        cla([1,2,3,4],'reset');
                                        display("time = ", num2str(time));
                                    end
                                end
                            end

                            % Plot the Results
    %                         clf(1,'reset');
    %                         clf(2,'reset');
    %                         clf(3,'reset');
    %                         clf(4,'reset');
                            if ~plot_in_time
                                make_plots(U,xs,gamma,flux,lam,rec,Order,xidx,k,extra_prim_var_plot,riemann_solver);
                            end
                            if plot_in_time
                                Riemann(rhoL,uL,PL,rhoR,uR,PR,t,middle,R,gamma,cells,extra_prim_var_plot);
                                make_plots(U,xs,gamma,flux,lam,rec,Order,xidx,k,extra_prim_var_plot,riemann_solver);
                            end
                        end
                    end
                end
            end
        end
    end
    

    %% Plotting

    % Adjust the plots based on test number and save them
    configure_plots(recon,test,cells,cfl,extra_prim_var_plot,riemann_solver,CD_Term_Order)
    L2 = L2Norm(U(:,xidx),Uexact,cells)
    Linf = LinfNorm(U(:,xidx),Uexact)
end

%% Functions
% These are all of the functions used only in this script. The solvers and
% reconstruction methods, as well as any function that is used in multiple
% scripts have their own scripts.

function F_vec = dUdt(U, flux, lam, gamma,dx,kap,rec,CD_Term_Order,riemann_solver,xs,time)
% dUdt: Computes the Flux vector to be added to the Conservative Variable
% Vector.
%
% EXAMPLE:
%         Unew = U + dUdt(U, flux, lam, gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time)
%
% REQUIRES:
%         % Scripts
%         MUSCL_Rec.m
%         WENO3_Rec.m
%         WENO5_Rec.m
%         Roe_solver.m
%         Roe_Pike_solver.m
%         HLL_solver.m
%         HLLC_solver.m
%         Osher_solver.m
%         % Functions
%
% INPUTS:
%         U: The Conservative Variable Vector
%         flux: The flux limiter for MUSCL Reconstruction
%         lam: The Scaling Factor for MUSCL Reconstruction
%         gamma: The ratio of specific heats
%         dx: The step size in the x-direction.
%         kappa: The Reconstruction Order for MUSCL Reconstruction
%         rec: The Reconstruction Method
%         CD_Term_Order: The Central Differencing Term Order
%         riemann_solver: The type of Riemann Solver to use
%         xs: The x-vector
%         time: The current time for the time step
% OUTPUTS:
%         F_vec: The flux vector to be added to the conservative variable
%                vector
    switch rec
        case 'ROE'
            [URi,ULi,URim1,ULim1] = MUSCL_Rec(U,18,lam,kap);
            BCs = [0 0 0; 0 0 0; 0 0 0];
        case 'MUSCL'
            [URi,ULi,URim1,ULim1] = MUSCL_Rec(U,flux,lam,kap);
            BCs = [0 0 0; 0 0 0; 0 0 0];
        case 'WENO3'
            [URi,ULi] = WENO3_Rec(U(:,2:end));
            [URim1,ULim1] = WENO3_Rec(U(:,1:end-1));
            BCs = [0 0 0; 0 0 0; 0 0 0];
        case 'WENO5'
            [URi,ULi] = WENO5_Rec(U(:,2:end));
            [URim1,ULim1] = WENO5_Rec(U(:,1:end-1));
            BCs = [0 0 0; 0 0 0; 0 0 0];
    end
    switch riemann_solver
        case 'Roe'
            [Fip12, Fim12] = Roe_solver(U,URi,ULi,URim1,ULim1,gamma,CD_Term_Order);
        case 'Roe-Pike'
            [Fip12, Fim12] = Roe_Pike_solver(U,URi,ULi,URim1,ULim1,gamma,CD_Term_Order);
        case 'HLL'
            [Fip12, Fim12] = HLL_solver(U,URi,ULi,URim1,ULim1,xs,gamma,time);
        case 'HLLC'
            [Fip12, Fim12] = HLLC_solver(U,URi,ULi,URim1,ULim1,xs,gamma,time);
        case 'Osher'
            [Fip12, Fim12] = Osher_solver(U,URi,ULi,URim1,ULim1,xs,gamma,time);
    end
    F_vec = -(1/dx)*(Fip12 - Fim12);
    F_vec = [BCs, F_vec, BCs];
end

function U = Forward_Euler(U,flux,lam,gamma,dx,dt,kappa,rec,CD_Term_Order,riemann_solver,xs,time)
for i = 1:10
    dUstar = dUdt(U,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
    U = U + dUstar*dt;
end
end

function U = RK2(U,flux,lam,gamma,dx,dt,kappa,rec,CD_Term_Order,riemann_solver,xs,time)
    k1 = U + dt*dUdt(U,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
    U = 0.5*U + 0.5*k1 + 0.5*dt*dUdt(k1,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
end

function U = RK3(U,flux,lam,gamma,dx,dt,kappa,rec,CD_Term_Order,riemann_solver,xs,time)
    k1 = U + dt*dUdt(U,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
    k2 = 0.75*U + 0.25*k1 + 0.25*dt*dUdt(k1,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
    U = (1/3)*U + (2/3)*k2 + (2/3)*dt*dUdt(k2,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
end

function U = RK4(U,flux,lam,gamma,dx,dt,kappa,rec,CD_Term_Order,riemann_solver,xs,time)
    k0 = dt*dUdt(U,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
    k1 = U + 0.5*dt*dUdt(U,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
    k2 = (649/1600)*U - (10890423/25193600)*k0 + (951/1600)*k1 +...
         (5000/7873)*dt*dUdt(k1,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
    k3 = (53989/2500000)*U - (102261/5000000)*k0 + (4806213/20000000)*k1...
         - (5121/20000)*dt*dUdt(k1,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time) + ...
         (23619/32000)*k2 +...
         (7873/10000)*dt*dUdt(k2,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
    U = (1/5)*U + (1/10)*k0 + (6127/30000)*k1 + ...
        (1/6)*dt*dUdt(k1,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time) + ...
        (7873/30000)*k2 + (1/3)*k3 + ...
        (1/6)*dt*dUdt(k3,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
end

function U = RK45(U,flux,lam,gamma,dx,dt,kappa,rec,CD_Term_Order,riemann_solver,xs,time)
    U0 = U;
    K0 = dt*dUdt(U0,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
    U1 = U0 + (13736793/35065003)*K0;
    K1 = dt*dUdt(U1,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
    U2 = (18384209/41371354)*U0 + (22987145/41371354)*U1...
         + (1106722/3004045)*K1;
    K2 = dt*dUdt(U2,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
    U3 = (16461287/26546102)*U0 + (10084815/26546102)*U2...
         + (9149709/36323969)*K2;
    K3 = dt*dUdt(U3,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
    U4 = (1818641/10212497)*U0 + (8393856/10212497)*U3...
         + (11974013/21971684)*K3;
    K4 = dt*dUdt(U4,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
    U = (12153930/23498039)*U2 + (4482614/46664871)*U3...
        + (9148849/143640986)*K3...
        + (8047189/20809438)*U4...
        + (9169579/40572015)*K4;
end

function U = SSPRK104(U,flux,lam,gamma,dx,dt,kappa,rec,CD_Term_Order,riemann_solver,xs,time)
U0 = U;
U1 = U0 + dt/6 *dUdt(U0, flux, lam, gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
U2 = U1 + dt/6 * dUdt(U1,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
U3 = U2 + dt/6 * dUdt(U2,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
U4 = U3 + dt/6 * dUdt(U3,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
U5 = U4 + dt/6 * dUdt(U4,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
Ustar = 1/25 * U0 + 9/25 * U5;
U5 = 15 * Ustar - 5 * U5;
U6 = U5 + dt/6 * dUdt(U5,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
U7 = U6 + dt/6 * dUdt(U6,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
U8 = U7 + dt/6 * dUdt(U7,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
U9 = U8 + dt/6 * dUdt(U8,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
U = 6/10 * U9 + dt/10 * dUdt(U9,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time) + Ustar;
end

function U = SSPRK65(U,flux,lam,gamma,dx,dt,kap,rec,CD_Term_Order,riemann_solver,xs,time)
global time_int_counter;
time_int_counter = time_int_counter + 1;
U0 = U;
U1 = U0 + dt/2 *dUdt(U0, flux, lam, gamma,dx,kap,rec,CD_Term_Order,riemann_solver,xs,time);
U2 = U1 + dt/2 *dUdt(U1, flux, lam, gamma,dx,kap,rec,CD_Term_Order,riemann_solver,xs,time);
U3 = U2 + dt/2 *dUdt(U2, flux, lam, gamma,dx,kap,rec,CD_Term_Order,riemann_solver,xs,time);
U4 = U3 + dt/2 *dUdt(U3, flux, lam, gamma,dx,kap,rec,CD_Term_Order,riemann_solver,xs,time);
U5 = U4 + dt/2 *dUdt(U4, flux, lam, gamma,dx,kap,rec,CD_Term_Order,riemann_solver,xs,time);
U6 = (1/9)*U0 + (2/5)*U1 + (4/9)*U3 +(2/45)*U5 + dt/45 *dUdt(U5, flux, lam, gamma,dx,kap,rec,CD_Term_Order,riemann_solver,xs,time);
U = U6;
end

function file_pre = Make_File_Prefix(Schemes, test_num, n, cfl, riemann_solver,CD_Term_Order)
if strcmp(riemann_solver, 'Roe')
    solver = 'Roe';
elseif strcmp(riemann_solver, 'Roe-Pike')
    solver = 'RoePike';
elseif strcmp(riemann_solver, 'HLL')
    solver = 'HLL';
elseif strcmp(riemann_solver, 'HLLC')
    solver = 'HLLC';
elseif strcmp(riemann_solver, 'Osher')
    solver = 'Osher';
else
    solver = '';
end
Scheme_Path = [];
for i = 1:length(Schemes)
    if i ~= length(Schemes)
        curr = [char(Schemes(i)) 'v'];
    else
        curr = [char(Schemes(i))];
    end
    Scheme_Path = [Scheme_Path curr];
end
if CD_Term_Order == 1
    order = '';
else
    order = ['CD' num2str(CD_Term_Order)];
end
file_pre = ['Test' num2str(test_num) order Scheme_Path 'CFL'...
            num2str(cfl) 'Nodes' num2str(n) solver];
end

function Folder_Path = Make_Path(Schemes, test_num, n, cfl,...
                                 flux_limiter, recon,CD_Term_Order,riemann_solver)
if length(recon) ~= 1
    Scheme_Path = [];
    for i = 1:length(recon)
        if i == length(recon)
            Scheme_Path = [Scheme_Path recon{i} '/'];
        else
            Scheme_Path = [Scheme_Path recon{i} '_v_'];
        end
    end
    for i = 1:length(CD_Term_Order)
        if i == length(CD_Term_Order)
            cdname = get_CD_name(CD_Term_Order(i));
            Scheme_Path = [Scheme_Path cdname '/'];
        else
            cdname = get_CD_name(CD_Term_Order(i));
            Scheme_Path = [Scheme_Path cdname '_v_'];
        end
    end
    
    Folder_Path = ['Figures/Test#' num2str(test_num) '/'...
                    Scheme_Path num2str(n) '_Nodes' '/'...
                    'CFL=' num2str(cfl) '/']; % You can change this as
                                                  % needed
else
    if strcmp(recon{1}, 'MUSCL')
        if ~isscalar(flux_limiter)
            flux = ['/flux_lim='];
            k = 0;
            for f = flux_limiter
                k = k+1;
                if k == length(flux_limiter)
                    flux = [flux get_flux_name(f) '/'];
                else
                    flux = [flux get_flux_name(f) '_'];
                end
            end
            flux_lim = flux;
        else
            switch flux_limiter
                case 1
                    flux_lim = 'flux_lim=minmod/';
                case 2
                    flux_lim = 'flux_lim=van_Leer/';
                case 3
                    flux_lim = 'flux_lim=Barth_Jesperson/';
                case 4
                    flux_lim = 'flux_lim=Superbee/';
                case 5
                    flux_lim = 'flux_lim=van_Albada_2/';
                case 6
                    flux_lim = 'flux_lim=van_Albada_1/';
                case 7
                    flux_lim = 'flux_lim=CHARM/';
                case 8
                    flux_lim = 'flux_lim=HCUS/';
                case 9
                    flux_lim = 'flux_lim=HQUICK/';
                case 10
                    flux_lim = 'flux_lim=Koren/';
                case 11
                    flux_lim = 'flux_lim=MC/';
                case 12
                    flux_lim = 'flux_lim=Osher/';
                case 13
                    flux_lim = 'flux_lim=Ospre/';
                case 14
                    flux_lim = 'flux_lim=Smart/';
                case 15
                    flux_lim = 'flux_lim=Sweby/';
                case 16
                    flux_lim = 'flux_lim=UMIST/';
                case 17
                    flux_lim = 'flux_lim=General_Minmod/';
                case 18
                    flux_lim = 'Roe/';
                case 19
                    flux_lim = '/MultipleFluxLims/';
            end
        end
        Scheme_Path = [];
        for i = 1:length(Schemes)
            if i == length(Schemes)
                curr = [char(Schemes(i)) '/'];
            else
                curr = [char(Schemes(i)) '_'];
            end
            Scheme_Path = [Scheme_Path curr];
        end
        
        for i = 1:length(CD_Term_Order)
            if i == length(CD_Term_Order)
                cdname = get_CD_name(CD_Term_Order(i));
                Scheme_Path = [Scheme_Path cdname '/'];
            else
                cdname = get_CD_name(CD_Term_Order(i));
                Scheme_Path = [Scheme_Path cdname '_v_'];
            end
        end

        Folder_Path = ['Figures/Test#' num2str(test_num) '/'...
                        Scheme_Path flux_lim num2str(n) '_Nodes' '/'...
                        'CFL=' num2str(cfl) '/']; % You can change this as
                                                  % needed
    else
        Scheme_Path = [recon{1} '/'];
        
        for i = 1:length(CD_Term_Order)
            if i == length(CD_Term_Order)
                cdname = get_CD_name(CD_Term_Order(i));
                Scheme_Path = [Scheme_Path cdname '/'];
            else
                cdname = get_CD_name(CD_Term_Order(i));
                Scheme_Path = [Scheme_Path cdname '_v_'];
            end
        end

        Folder_Path = ['Figures/Test#' num2str(test_num) '/'...
                        Scheme_Path num2str(n) '_Nodes' '/'...
                        'CFL=' num2str(cfl) '/']; % You can change this as
                                                  % needed
    end
end
        
    mssg = ['Results being saved in:\n\t.../', Folder_Path, '\n\n'];

    if ~exist(Folder_Path, 'dir')
        mkdir(Folder_Path)
        fprintf(mssg);
    else
        fprintf(mssg);
    end
end

function [name] = get_CD_name(n)
switch n
    case 1
        name = 'Roe';
    case 2
        name = 'CD2';
    case 4
        name = 'CD4';
    case 6
        name = 'CD6';
end
end

function [smax] = find_wave_speed(U, gamma)
UL = U(:,1:end-1);
UR = U(:,2:end);
[rL, uL, ~, ~, HL, ~, ~, ~] = cons_to_prim(UL, gamma);
[rR, uR, ~, ~, HR, ~, ~, ~] = cons_to_prim(UR, gamma);
% [~, u, ~, a, ~, ~, ~, ~] = cons_to_prim(U, gamma);
ut = (sqrt(rL).*uL + sqrt(rR).*uR)./(sqrt(rL) + sqrt(rR));
Ht = (sqrt(rL).*HL + sqrt(rR).*HR)./(sqrt(rL) + sqrt(rR));
at = sqrt((gamma-1).*(Ht - 0.5.*ut.*ut));
SL = max(abs(ut - at));
SR = max(abs(ut + at));
% smax1 = max(abs(u)+a);
% smax2 = max(abs(u));
% smax3 = max(abs(u-a));
% smax = max([smax1,smax2,smax3,SL,SR]);
smax = max(SL,SR);
end

function dt = calculate_time_step(cfl, dx, time, tfinal, U, gamma)
smax = find_wave_speed(U, gamma);

dt = (cfl*dx)/smax;
if isnan(dt)
    dt = 1e-6;
end

if time + dt > tfinal
    dt = tfinal - time;
end

end

function [xs, dx, n2, xidx] = set_up_grid(n, ghost, xstart, xend)

xs = linspace(xstart, xend, n); % x-Domain of the Tests
dx = (xs(2) - xs(1)); % The spatial step
left_bound = linspace(-ghost*dx, -dx, ghost); % The left Boundary Cells
right_bound = linspace(xs(end)+dx, xs(end)+ghost*dx, ghost); % The Right Boun-
                                                          % dary Cells
xs = [left_bound, xs, right_bound]; % Concatenate the x-domain array
n2 = n+2*ghost; % The number of cells including the ghost cells
xidx = ghost+1:n+ghost; % The indexes of the test domain, x = 0:1

end

function flux_lim = get_flux_name(flux)
        switch flux
            case 1
                flux_lim = 'MM';
            case 2
                flux_lim = 'VL';
            case 3
                flux_lim = 'BJ';
            case 4
                flux_lim = 'SB';
            case 5
                flux_lim = 'VA2';
            case 6
                flux_lim = 'VA1';
            case 7
                flux_lim = 'CHARM';
            case 8
                flux_lim = 'HCUS';
            case 9
                flux_lim = 'HQUICK';
            case 10
                flux_lim = 'Koren';
            case 11
                flux_lim = 'MC';
            case 12
                flux_lim = 'Osher';
            case 13
                flux_lim = 'Ospre';
            case 14
                flux_lim = 'Smart';
            case 15
                flux_lim = 'Sweby';
            case 16
                flux_lim = 'UMIST';
            case 17
                flux_lim = 'GMM';
            case 18
                flux_lim = '';
        end
end

function name = get_legend_entry_name(flux,lam,rec,CD_Term_Order,riemann_solver)
if strcmp(riemann_solver, 'Roe')
    solver = 'Roe';
elseif strcmp(riemann_solver, 'Roe-Pike')
    solver = 'Roe-Pike';
elseif strcmp(riemann_solver, 'HLL')
    solver = 'HLL';
else
    solver = '';
end
if strcmp(rec,'ROE')
    if CD_Term_Order == 1
        name = [solver];
    else
        cdname = get_CD_name(CD_Term_Order);
        name = [cdname '--' solver];
    end
elseif strcmp(rec,'MUSCL')
    flux_lim = get_flux_name(flux);
    if flux == 18
        if CD_Term_Order == 1
            name = [solver ' ' flux_lim];
        else
            name = ['CD' num2str(CD_Term_Order) '-' solver ' ' flux_lim];
        end
    else
        if CD_Term_Order == 1
            name = [solver '--' rec];
        else
            name = ['CD' num2str(CD_Term_Order) '--' solver '--' rec];
        end
    end
else
    if CD_Term_Order == 1
        name = [solver '--' rec];
    else
        name = ['CD' num2str(CD_Term_Order) '--' solver '--' rec];
    end
end
end

function [U] = set_up_U(rhoR,rhoL,uL,uR,PR,PL,g,x,xmid,n)
% sets left and right state values (sets up riemann problem)
r = zeros(size(x));
u = zeros(size(x));
P = zeros(size(x));
a = zeros(size(x));
H = zeros(size(x));
U = zeros([3,length(x)]);
    for i = 1:n
        if x(i) <= xmid
            r(i) = rhoL;
            u(i) = uL;
            P(i) = PL;
            a(i) = sqrt(g*P(i)/r(i));
            H(i) = 0.5*u(i)^2+a(i)^2/(g-1);

            U(1,i) = rhoL;
            U(2,i) = rhoL*uL;
            U(3,i) = 0.5*rhoL*uL^2+PL/(g-1);
        else
            r(i) = rhoR;
            u(i) = uR;
            P(i) = PR;
            a(i) = sqrt(g*P(i)/r(i));
            H(i) = 0.5*u(i)^2+a(i)^2/(g-1);

            U(1,i) = rhoR;
            U(2,i) = rhoR*uR;
            U(3,i) = 0.5*rhoR*uR^2+PR/(g-1);
        end
    end
end

function Name = get_scheme_name(flux_limiter,kappa)
    if isscalar(flux_limiter)
        if flux_limiter ~= 18
            Name = ['Roe-MUSCL_kappa=' num2str(kappa)];
        else
            Name = 'Roe';
        end
    elseif ~isscalar(flux_limiter) && any(flux_limiter==18)
        Name = ['Roe_v_Roe-MUSCL_kappa=' num2str(kappa)];
    else
        Name = ['Roe-MUSCL_kappa=' num2str(kappa)];
    end
end

function [locs] = get_locations(test,CD_Term_Order)
if test == 1
    if CD_Term_Order == 1
    locs = [0.7200 0.6900 0.1682 0.1570; % Northeast
            0.1500 0.2100 0.1751 0.1735; % Southwest
            0.7200 0.6900 0.1682 0.1570; % Northeast
            0.2000 0.6900 0.1682 0.1570; % Northwest
            0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570];
    else
    locs = [0.6750 0.6900 0.1682 0.1570; % Northeast
            0.2000 0.2100 0.1751 0.1735; % Southwest
            0.6750 0.6900 0.1682 0.1570; % Northeast
            0.2000 0.6900 0.1682 0.1570; % Northwest
            0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570];
    end
elseif test == 2
    locs = [0.7200 0.6900 0.1682 0.1570; % Northeast
            0.1500 0.2100 0.1751 0.1735; % Southwest
            0.7200 0.6900 0.1682 0.1570; % Northeast
            0.1500 0.6900 0.1682 0.1570; % Northwest
            0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570];
elseif test == 3
    locs = [0.1500 0.6900 0.1682 0.1570; % Northwest
            0.4500 0.2100 0.1682 0.1570; % SouthCenter
            0.7200 0.6900 0.1682 0.1570; % Northeast
            0.1500 0.2100 0.1751 0.1735; % Southwest
            0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570];
elseif test == 4
    if CD_Term_Order == 1
    locs = [0.1500 0.6900 0.1682 0.1570; % Northwest
            0.1500 0.2100 0.1751 0.1735; % Southwest
            0.1500 0.6900 0.1682 0.1570; % Northwest
            0.1500 0.2100 0.1751 0.1735; % Southwest
            0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570];
    else
    locs = [0.2000 0.6900 0.1682 0.1570; % Northwest
            0.2000 0.2100 0.1751 0.1735; % Southwest
            0.2000 0.6900 0.1682 0.1570; % Northwest
            0.2000 0.2100 0.1751 0.1735; % Southwest
            0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570];
    end
elseif test == 5
    locs = [0.1500 0.6900 0.1682 0.1570; % Northwest
            0.4500 0.2100 0.1682 0.1570; % SouthCenter
            0.1500 0.2100 0.1751 0.1735; % Southwest
            0.1500 0.2100 0.1751 0.1735; % Southwest
            0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570];
elseif test == 6
    locs = [0.6750 0.6900 0.1682 0.1570; % Northeast
            0.6750 0.6900 0.1682 0.1570; % Northeast
            0.6750 0.6900 0.1682 0.1570; % Northeast
            0.6750 0.6900 0.1682 0.1570; % Northeast
            0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570];
elseif test == 7
    locs = [0.6750 0.6900 0.1682 0.1570; % Northeast
            0.6750 0.6900 0.1682 0.1570; % Northeast
            0.6750 0.6900 0.1682 0.1570; % Northeast
            0.6750 0.6900 0.1682 0.1570; % Northeast
            0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570];
elseif test == 8
    locs = [0.6750 0.6900 0.1682 0.1570; % Northeast
            0.4500 0.6900 0.1682 0.1570; % NorthCenter
            0.4500 0.2100 0.1682 0.1570; % SouthCenter
            0.6750 0.2100 0.1682 0.1570; % Southeast
            0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570];
elseif test == 13
    locs = [0.16 0.6900 0.1682 0.1570; % Northwest
            0.4500 0.2100 0.1682 0.1570; % SouthCenter
            0.1570 0.2100 0.1751 0.1735; % Southwest
            0.1570 0.2100 0.1751 0.1735; % Southwest
            0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570];
elseif test == 15
    locs = [0.6750 0.6900 0.1682 0.1570; % Northeast
            0.2000 0.6900 0.1682 0.1570; % Northwest
            0.6750 0.6900 0.1682 0.1570; % Northeast
            0.2000 0.6900 0.1682 0.1570; % Northwest
            0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570];
else
    locs = [0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570;
            0.7000 0.7000 0.1682 0.1570];
end
end

function ticks = find_ticks(umin,umax)
ymin = floor(min(umin));
ymax = ceil(max(umax));
scaling = 10^(max(numel(num2str(abs(ymin))), numel(num2str(abs(ymax)))) - 1);
innerticks = scaling*0.2;
ylim1 = ymin - innerticks;
ylim2 = ymax + innerticks;
ticks = ylim1:innerticks:ylim2;
end

function make_plots(U,xs,gamma,flux,lam,rec,Order,xidx,k,extra_prim_var_plot,riemann_solver)
% Inputs for Saving Plots
line_colors = (1/255)*[
                       0,0,0; % Black
                       65,105,225; % King Blue
                       237,29,36; % Comic Book Red
                       8,255,8; % Fluorescent Green
                       188,128,189; % Pastel Light Purple
                       255,127,0; % Orange
                       51,160,44; % Green
                       178,223,138; % Light Green
                       227,26,28; % Red
                       117,112,179; % Pastel Purple
                       231,41,138; % Bright Pink
                       128,177,211; % Pastel Blue
                       202,178,214; % Light Purple
                       166,206,227; % Light Blue
                       217,95,2; % Pastel Red
                       102,166,30; % Pastel Green
                       27,158,119; % Blue-Green
                       106,61,154; % Purple
                       251,154,153; % Pink
                       253,191,111; % Light Orange
                       255,255,153; % Light Yellow
                       177,89,40; % Brown
                       141,211,199; % Light Blue-Green
                       204,235,197; % Pastel Light Blue-Green
                       31,120,180; % Blue
                       251,128,114; % Salmon
                       ];
marker_styles = {'o','o','o','o','o','o','o','o','o','o','o','o','o',...
                 'o','o','o','x','x','x','x','x','x','x','x','x''x','x',...
                 'x','x','x','x','x'};
line_styles = {'-','-','-','-','-','-','-','-','-','-','-','-','-','-',...
               '-','-','--','--','--','--','--','--','--','--','--',...
               '--','--','--','--','--','--','--'};

[r, u, p, a, H, e, m, s] = cons_to_prim(U(:, xidx), gamma);

name = get_legend_entry_name(flux,lam,rec,Order,riemann_solver);

figure(1)
set(gcf,'renderer','painters');
plot(xs(xidx),r,'Color', line_colors(k,:), 'LineStyle',...
    line_styles{k},...
    'LineWidth', 8, 'MarkerSize',4,'DisplayName',name);

figure(2)
set(gcf,'renderer','painters');
plot(xs(xidx),u,'Color', line_colors(k,:), 'LineStyle',...
    line_styles{k},...
    'LineWidth', 8, 'MarkerSize',4,'DisplayName',name);

figure(3)
set(gcf,'renderer','painters');
plot(xs(xidx),p,'Color', line_colors(k,:), 'LineStyle',...
    line_styles{k},...
    'LineWidth', 8, 'MarkerSize',4,'DisplayName',name);

figure(4)
set(gcf,'renderer','painters');
plot(xs(xidx),e,'Color', line_colors(k,:), 'LineStyle',...
    line_styles{k},...
    'LineWidth', 8, 'MarkerSize',4,'DisplayName',name);

if extra_prim_var_plot == true
    figure(5)
    set(gcf,'renderer','painters');
    plot(xs(xidx),a,'Color', line_colors(k,:),...
        'LineStyle', line_styles{k},'LineWidth', 8, 'MarkerSize',...
        12,'DisplayName',name);

    figure(6)
    set(gcf,'renderer','painters');
    plot(xs(xidx),m,'Color', line_colors(k,:),...
        'LineStyle', line_styles{k},'LineWidth', 8, 'MarkerSize',...
        12,'DisplayName',name);

    figure(7)
    set(gcf,'renderer','painters');
    plot(xs(xidx),H,'Color', line_colors(k,:),...
        'LineStyle', line_styles{k},'LineWidth', 8, 'MarkerSize',...
        12,'DisplayName',name);

    figure(8)
    set(gcf,'renderer','painters');
    plot(xs(xidx),s,'Color', line_colors(k,:),...
        'LineStyle', line_styles{k},'LineWidth', 8, 'MarkerSize',...
        12,'DisplayName',name);
end

end

function configure_plots(recon,test_num,cells,cfl,extra_prim_var_plot,riemann_solver,CD_Term_Order)
% Folder = Make_Path(recon, test_num, cells, cfl, flux_limiter, recon, CD_Term_Order);
file_pre = Make_File_Prefix(recon, test_num, cells, cfl, riemann_solver,CD_Term_Order);
locs = get_locations(test_num,CD_Term_Order);
% Folder = '../thesis/CH5/EPSFDocs/';
% Folder = '';
Folder = '../thesis/FigsforCommittee/';

figure(1)
set(gcf, 'Position', get(0, 'Screensize'));
ax = gca;
ax.YAxis.FontSize = 55;
ax.XAxis.FontSize = 55;
if CD_Term_Order == 2
    ylim1 = [0,0,0,4,0,0,0,0.85,0,0,0,0,0.15,0,0,0,0,0,0];
    ylim2 = [1.1,0,7,35,7,1.1,0,1.2,0,0,0,0,1.45,0,1.1,1.1,0,0,0];
elseif CD_Term_Order == 4
    ylim1 = [0,0,0,4,0,0,0,0.875,0,0,0,0,0.2,0,0,0,0,0,0];
    ylim2 = [1.1,0,7,35,7,1.1,0,1.16,0,0,0,0,1.4,0,1.1,1.1,0,0,0];
else
    ylim1 = [0,0,0,4,0,0,0,0.875,0,0,0,0,0.2,0,0,0,0,0,0];
    ylim2 = [1.1,0,7,35,7,1.1,0,1.16,0,0,0,0,1.4,0,1.1,1.1,0,0,0];
end
ax.YAxis.Limits = [ylim1(test_num) ylim2(test_num)];
ax.LineWidth = 10;
filename = [Folder file_pre 'Density.eps'];
legend('FontSize', 60,'Interpreter','Latex')
legend('boxoff')
leg = legend();
leg.ItemTokenSize = [80,25];
leg.Position = locs(1,:);
xlabel('$x$','interpreter','latex', 'FontSize', 100, 'FontWeight','bold');
ylabel('$\rho$','interpreter','latex', 'FontSize', 100, 'FontWeight','bold');
grid off;
saveas(gcf, filename, 'epsc');

figure(2)
set(gcf, 'Position', get(0, 'Screensize'));
ax = gca;
ax.YAxis.FontSize = 55;
ax.XAxis.FontSize = 55;
if CD_Term_Order == 2
    ylim1 = [-0.2,0,-3,-10,-24,-0.05,0,-0.55,0,0,0,0,-0.1,0,-0.2,-2.1,0,0,0];
    ylim2 = [1.6,0,23,25,4,1,0,0.05,0,0,0,0,1.7,0,1.8,-0.85,0,0,0];
elseif CD_Term_Order == 4
    ylim1 = [-0.2,0,-3,-10,-24,-0.05,0,-0.475,0,0,0,0,-0.1,0,-0.2,-2.15,0,0,0];
    ylim2 = [1.6,0,23,25,4,1,0,0.05,0,0,0,0,1.6,0,1.8,-1,0,0,0];
else
    ylim1 = [-0.2,0,-3,-10,-24,-0.05,0,-0.475,0,0,0,0,-0.1,0,-0.2,-2.1,0,0,0];
    ylim2 = [1.6,0,23,25,4,1,0,0.05,0,0,0,0,1.6,0,1.8,-1,0,0,0];
end
ax.YAxis.Limits = [ylim1(test_num) ylim2(test_num)];
ax.LineWidth = 10;
filename = [Folder file_pre 'Velocity.eps'];
legend('FontSize', 60,'Interpreter','Latex')
legend('boxoff')
leg = legend();
leg.ItemTokenSize = [80,25];
leg.Position = locs(2,:);
xlabel('$x$','interpreter','latex', 'FontSize', 100, 'FontWeight','bold');
ylabel('$u$','interpreter','latex', 'FontSize', 100, 'FontWeight','bold');
grid off;
saveas(gcf, filename, 'epsc');

figure(3)
set(gcf, 'Position', get(0, 'Screensize'));
ax = gca;
ax.YAxis.FontSize = 55;
ax.XAxis.FontSize = 55;
if CD_Term_Order == 2
    ylim1 = [0,0,-100,-100,-100,0,0,6.75,0,0,0,0,0.5,0,0,0,0,0,0];
    ylim2 = [1.1,0,1100,1950,1100,1.1,0,10.25,0,0,0,0,3.75,0,1.1,1.1,0,0,0];
elseif CD_Term_Order == 4
    ylim1 = [0,0,-100,-100,-100,0,0,6.75,0,0,0,0,0.5,0,0,0,0,0,0];
    ylim2 = [1.1,0,1100,1900,1100,1.1,0,10.25,0,0,0,0,3.75,0,1.1,1.1,0,0,0];
else
    ylim1 = [0,0,-100,-100,-100,0,0,6.75,0,0,0,0,0.5,0,0,0,0,0,0];
    ylim2 = [1.1,0,1100,1900,1100,1.1,0,10.25,0,0,0,0,3.75,0,1.1,1.1,0,0,0];
end
ax.YAxis.Limits = [ylim1(test_num) ylim2(test_num)];
ax.LineWidth = 10;
filename = [Folder file_pre 'Pressure.eps'];
legend('FontSize', 60,'Interpreter','Latex')
legend('boxoff')
leg = legend();
leg.ItemTokenSize = [80,25];
leg.Position = locs(3,:);
xlabel('$x$','interpreter','latex', 'FontSize', 100, 'FontWeight','bold');
ylabel('$p$','interpreter','latex', 'FontSize', 100, 'FontWeight','bold');
grid off;
saveas(gcf, filename, 'epsc');

figure(4)
set(gcf, 'Position', get(0, 'Screensize'));
ax = gca;
ax.YAxis.FontSize = 55;
ax.XAxis.FontSize = 55;
if CD_Term_Order == 2
    ylim1 = [1.8,0,-200,-20,-100,1.6,0,17,0,0,0,0,2,0,1,1.6,0,0,0];
    ylim2 = [3.7,0,2700,400,2600,3,0,25.5,0,0,0,0,32,0,5,3.1,0,0,0];
elseif CD_Term_Order == 4
    ylim1 = [1.8,0,-200,-20,-100,1.6,0,17,0,0,0,0,2,0,1,1.6,0,0,0];
    ylim2 = [3.7,0,2700,340,2600,3,0,25.5,0,0,0,0,20,0,5,3,0,0,0];
else
    ylim1 = [1.8,0,-200,-20,-100,1.6,0,17,0,0,0,0,2,0,1,1.75,0,0,0];
    ylim2 = [3.7,0,2700,340,2600,3,0,25.5,0,0,0,0,20,0,5,3,0,0,0];
end
ax.YAxis.Limits = [ylim1(test_num) ylim2(test_num)];
ax.LineWidth = 10;
filename = [Folder file_pre 'InternalEnergy.eps'];
legend('FontSize', 60,'Interpreter','Latex')
legend('boxoff')
leg = legend();
leg.ItemTokenSize = [80,25];
leg.Position = locs(4,:);
xlabel('$x$','interpreter','latex', 'FontSize', 100, 'FontWeight','bold');
ylabel('$e$','interpreter','latex', 'FontSize', 100, 'FontWeight','bold');
grid off;
saveas(gcf, filename, 'epsc');

if extra_prim_var_plot
    figure(5)
    ax = gca;
    ax.YAxis.FontSize = 55;
    ax.XAxis.FontSize = 55;
    ymin = ax.YAxis.Limits(1);
    ymax = ax.YAxis.Limits(2);
    yticks = round(linspace(ymin,ymax,6),1);
    ax.YAxis.TickValues = yticks;
    scale = max(numel(num2str(yticks(3))), numel(num2str(yticks(4))));
    scale = 0.1*scale + 0.01;
    xloc = (ax.XAxis.Limits(1)) - scale*ax.YAxis.FontSize*0.01;
    yloc = (ax.YAxis.Limits(2) - ax.YAxis.Limits(1))/2  + ax.YAxis.Limits(1);
    ax.YAxis.Label.Position = [xloc yloc -1];
    ax.LineWidth = 10;
    filename = [Folder file_pre 'SpeedofSound.eps'];
    legend('FontSize', 60)
    legend('boxoff')
    leg = legend();
    leg.ItemTokenSize = [80,25];
    leg.Position = locs(5,:);
    xlabel('$x$','interpreter','latex', 'FontSize', 100, 'FontWeight','bold');
    ylabel('$a$','interpreter','latex', 'FontSize', 100, 'FontWeight','bold');
    grid off;
    saveas(gcf, filename, 'epsc');

    figure(6)
    ax = gca;
    ax.YAxis.FontSize = 55;
    ax.XAxis.FontSize = 55;
    ymin = ax.YAxis.Limits(1);
    ymax = ax.YAxis.Limits(2);
    yticks = round(linspace(ymin,ymax,6),1);
    ax.YAxis.TickValues = yticks;
    scale = max(numel(num2str(yticks(3))), numel(num2str(yticks(4))));
    scale = 0.1*scale + 0.01;
    xloc = (ax.XAxis.Limits(1)) - scale*ax.YAxis.FontSize*0.01;
    yloc = (ax.YAxis.Limits(2) - ax.YAxis.Limits(1))/2  + ax.YAxis.Limits(1);
    ax.YAxis.Label.Position = [xloc yloc -1];
    ax.LineWidth = 10;
    filename = [Folder file_pre 'Mach.eps'];
    legend('FontSize', 60)
    legend('boxoff')
    leg = legend();
    leg.ItemTokenSize = [80,25];
    leg.Position = locs(6,:);
    xlabel('$x$','interpreter','latex', 'FontSize', 100, 'FontWeight','bold');
    ylabel('$M$','interpreter','latex', 'FontSize', 100, 'FontWeight','bold');
    grid off;
    saveas(gcf, filename, 'epsc');

    figure(7)
    ax = gca;
    ax.YAxis.FontSize = 55;
    ax.XAxis.FontSize = 55;
    ymin = ax.YAxis.Limits(1);
    ymax = ax.YAxis.Limits(2);
    yticks = round(linspace(ymin,ymax,6),1);
    ax.YAxis.TickValues = yticks;
    scale = max(numel(num2str(yticks(3))), numel(num2str(yticks(4))));
    scale = 0.1*scale + 0.01;
    xloc = (ax.XAxis.Limits(1)) - scale*ax.YAxis.FontSize*0.01;
    yloc = (ax.YAxis.Limits(2) - ax.YAxis.Limits(1))/2  + ax.YAxis.Limits(1);
    ax.YAxis.Label.Position = [xloc yloc -1];
    ax.LineWidth = 10;
    filename = [Folder file_pre 'Enthalpy.eps'];
    legend('FontSize', 60)
    legend('boxoff')
    leg = legend();
    leg.ItemTokenSize = [80,25];
    leg.Position = locs(7,:);
    xlabel('$x$','interpreter','latex', 'FontSize', 100, 'FontWeight','bold');
    ylabel('$h$','interpreter','latex', 'FontSize', 100, 'FontWeight','bold');
    grid off;
    saveas(gcf, filename, 'epsc');

    figure(8)
    ax = gca;
    ax.YAxis.FontSize = 55;
    ax.XAxis.FontSize = 55;
    ymin = ax.YAxis.Limits(1);
    ymax = ax.YAxis.Limits(2);
    yticks = round(linspace(ymin,ymax,6),1);
    ax.YAxis.TickValues = yticks;
    scale = max(numel(num2str(yticks(3))), numel(num2str(yticks(4))));
    scale = 0.1*scale + 0.01;
    xloc = (ax.XAxis.Limits(1)) - scale*ax.YAxis.FontSize*0.01;
    yloc = (ax.YAxis.Limits(2) - ax.YAxis.Limits(1))/2  + ax.YAxis.Limits(1);
    ax.YAxis.Label.Position = [xloc yloc -1];
    ax.LineWidth = 10;
    filename = [Folder file_pre,2,4,6 'Entropy.eps'];
    legend('FontSize', 60)
    legend('boxoff')
    leg = legend();
    leg.ItemTokenSize = [80,25];
    leg.Position = locs(8,:);
    xlabel('$x$','interpreter','latex', 'FontSize', 100, 'FontWeight','bold');
    ylabel('$s$','interpreter','latex', 'FontSize', 100, 'FontWeight','bold');
    grid off;
    saveas(gcf, filename, 'epsc');
end
end