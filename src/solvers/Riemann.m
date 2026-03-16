function [Pstar,ustar,aL,aR,Uexact] = Riemann(rhoL,uL,PL,rhoR,uR,PR,t,middle,R,g,cells,extra_prim_var_plot)

% Exact Riemann Solver for Euler Equations
% Inputs:
%   rhoL - left state density
%   uL - left state velocity
%   PL - left state pressure
%   rhoR - right state density
%   uR - right state velocity
%   PR - right state pressure
%   t - run time
%   middle - location of initial solution
%   R - gas constant
%   g - ratio of specific heats
%   cells - number of grid cells
%   extra_prim_var_plot - plot extra primitive fields

% Setup grid
x = linspace(0,1,cells);

% Solution initialization
ustar = 0;
aL = 0;
aR = 0;
Uexact = zeros(3,length(x));

% Solve the Riemann problem using a pressure root-finding procedure
% (see Toro, "Riemann Solvers and Numerical Methods for Fluid Dynamics").

du = uR - uL;
Iter = 0;
error = 1;
guess_type = 0;
Iter_Limit = 10000;
TOL = 1e-6;

Pstar = pGuess(guess_type, PL, PR, uL, uR, rhoL, rhoR, g, du);
[Pstar, error] = newtonRaphsonP(Pstar,PL,PR,rhoL,rhoR,g,du,Iter_Limit,TOL);

while abs(error) > TOL || error < 0
    if guess_type <= 2
        guess_type = guess_type + 1;
        Pstar = pGuess(guess_type, PL, PR, uL, uR, rhoL, rhoR, g, du);
        [Pstar, error] = newtonRaphsonP(Pstar,PL,PR,rhoL,rhoR,g,du,Iter_Limit,TOL);
    elseif guess_type > 2 && guess_type < 4
        guess_type = guess_type + 1;
        disp('Newton''s Method failed... Attempting Bisection Method.');
        [Pstar, error] = bisectionP(PL, PR, uL, uR, rhoL, rhoR, g, du, 10*Iter_Limit, TOL);
    else
        disp('Cannot Numerically find Pstar, Please Solve Graphically');
        Pstar0 = pGuess(0, PL, PR, uL, uR, rhoL, rhoR, g, du);
        Pstar1 = pGuess(1, PL, PR, uL, uR, rhoL, rhoR, g, du);
        Pstar2 = pGuess(2, PL, PR, uL, uR, rhoL, rhoR, g, du);
        Pstar3 = pGuess(3, PL, PR, uL, uR, rhoL, rhoR, g, du);
        Pstar4 = pGuess(4, PL, PR, uL, uR, rhoL, rhoR, g, du);
        Pstara = min([Pstar0,Pstar1,Pstar2,Pstar3,Pstar4]);
        Pstarb = max([Pstar0,Pstar1,Pstar2,Pstar3,Pstar4]);
        Pstars = linspace(Pstara, Pstarb, 1000);
        ErrFunc = zeros(size(Pstars));
        for i = 1:length(Pstars)
            ErrFunc(i) = funP(Pstars(i),PL,rhoL,g)+funP(Pstars(i),PR,rhoR,g)+du;
        end
        figure(5)
        plot(Pstars, ErrFunc);
        grid off;
        drawnow;
        Pstar = input('Pstar0 Guess from Graph: ');
        [Pstar, error] = newtonRaphsonP(Pstar,PL,PR,rhoL,rhoR,g,du,Iter_Limit,TOL);
        break;
    end
end

ustar = 0.5*(uL+uR)+0.5*(funP(Pstar,PR,rhoR,g)-funP(Pstar,PL,rhoL,g));

% Populate the exact solution fields for plotting
% (this block can be extended if needed)
Uexact(1,:) = rhoL;
Uexact(2,:) = uL;
Uexact(3,:) = PL;

if extra_prim_var_plot
    % Additional plotting is handled outside this function
end
end
