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
% Outputs:
%   Pstar - star region pressure
%   ustar - star region velocity
%   rhostarL - left star region density 
%   rhostarR - right star region density

    % Riemann Problem (analytical solution)
    x = linspace(0,1,cells);
    % Solution
    aR = sqrt(g*PR/rhoR);
    aL = sqrt(g*PL/rhoL);
    TL = PL/(rhoL*R);
    TR = PR/(rhoR*R);
    du = uR-uL;

    Iter = 0;
    error = 1;
    guess_type = 0;
    Iter_Limit = 10000;
    TOL = 1e-6;
    Pstar = Pguess(guess_type, PL, PR, uL, uR, rhoL, rhoR, g, du);
    [Pstar, error] = Newton_Raphson_P(Pstar,PL,PR,rhoL,rhoR,g,du,Iter_Limit,TOL);
    while abs(error) > TOL || error < 0
        if guess_type <= 2
            guess_type = guess_type + 1;
            Pstar = Pguess(guess_type, PL, PR, uL, uR, rhoL, rhoR, g, du);
            [Pstar, error] = Newton_Raphson_P(Pstar,PL,PR,rhoL,rhoR,g,du,Iter_Limit,TOL);
        elseif guess_type > 2 && guess_type < 4
            guess_type = guess_type + 1;
            disp("Newton's Method failed... Attempting Bisection Method.");
            [Pstar, error] = Bisection_P(PL, PR, uL, uR, rhoL, rhoR, g, du, 10*Iter_Limit, TOL);
        else
            disp('Cannot Numerically find Pstar, Please Solve Graphically');
            Pstar0 = Pguess(0, PL, PR, uL, uR, rhoL, rhoR, g, du);
            Pstar1 = Pguess(1, PL, PR, uL, uR, rhoL, rhoR, g, du);
            Pstar2 = Pguess(2, PL, PR, uL, uR, rhoL, rhoR, g, du);
            Pstar3 = Pguess(3, PL, PR, uL, uR, rhoL, rhoR, g, du);
            Pstar4 = Pguess(4, PL, PR, uL, uR, rhoL, rhoR, g, du);
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
            [Pstar, error] = Newton_Raphson_P(Pstar,PL,PR,rhoL,rhoR,g,du,Iter_Limit,TOL);
            break;
        end
    end
%     Pstar
    ustar = 0.5*(uL+uR)+0.5*(funP(Pstar,PR,rhoR,g)-funP(Pstar,PL,rhoL,g));

    LS = 0;
    RS = 0;
    LR = 0;
    RR = 0;

    if Pstar > PL % left running shock
        disp('Left Shock');
        LS = 1;
        rhostarL = rhoL*((Pstar/PL)+(g-1)/(g+1))/((g-1)/(g+1)*(Pstar/PL)+1);
        ML_MS = sqrt((g+1)/(2*g)*(Pstar/PL)+(g-1)/(2*g));
        SL = uL-aL*sqrt((g+1)/(2*g)*(Pstar/PL)+(g-1)/(2*g));
        ustarL = uL - funP(Pstar,PL,rhoL,g);
        TstarL = Pstar/R/rhostarL;
    end

    if Pstar > PR % right running shock
        disp('Right Shock');
        RS = 1;
        rhostarR = rhoR*((Pstar/PR)+(g-1)/(g+1))/((g-1)/(g+1)*(Pstar/PR)+1);
        MR_MS = sqrt((g+1)/(2*g)*(Pstar/PR)+(g-1)/(2*g));
        SR = uR+aR*sqrt((g+1)/(2*g)*(Pstar/PR)+(g-1)/(2*g));
        ustarR = uR + funP(Pstar, PR, rhoR, g);
        TstarR = Pstar/R/rhostarR;
    end

    if Pstar < PL  % left rarefaction
        disp('Left Rarefaction');
        LR = 1;
        ustarL = uL - funP(Pstar, PL, rhoL, g);
        aLstar = 0.5*(uL-ustar+2*aL/(g-1))*(g-1);
        xrange_lowL = (uL-aL)*t+middle;
        xrange_highL = (ustar-aLstar)*t+middle;
        for i = 1:length(x)
            if x(i) <= xrange_highL && x(i) > xrange_lowL
                u(i) = 2/(g+1)*(uL*(g-1)/2+aL+(x(i)-middle)/t);
                P(i) = PL*(1-(g-1)/2*(u(i)-uL)/aL)^(2*g/(g-1));
                rho(i) = rhoL*(1-(g-1)/2*(u(i)-uL)/aL)^(2/(g-1));
                T(i) = TL*(1-(g-1)/2*(u(i)-uL)/aL)^2;
            end
        end
        rhostarL = rhoL*(Pstar/PL)^(1/g);
        TstarL = TL*(Pstar/PL)^((g-1)/g);
    end

    if Pstar < PR % right farefaction
        disp('Right Rarefaction');
        ustarR = uR + funP(Pstar, PR, rhoR, g);
        RR = 1;
        aRstar = (-uR+ustar+2*aR/(g-1))*(g-1)/2;
        xrange_lowR = (uR+aR)*t+middle;
        xrange_highR = (ustar+aRstar)*t+middle;
        for i = 1:length(x)
            if x(i) > xrange_highR && x(i) <= xrange_lowR
                u(i) = 2/(g+1)*(-aR+(g-1)/2*uR+(x(i)-middle)/t);
                P(i) = PR*(1+(g-1)/2*(u(i)-uR)/aR)^(2*g/(g-1));
                rho(i) = rhoR*(1+(g-1)/2*(u(i)-uR)/aR)^(2/(g-1));
                T(i) = TR*(1+(g-1)/2*(u(i)-uR)/aR)^2;
    %             P(i) = PR*(2/(g+1)-(g-1)/((g+1)*aR)*(uR-(x(i)-middle)/t))^(2*g/(g-1));
    %             rho(i) = rhoR*(2/(g+1)-(g-1)/((g+1)*aR)*(uR-(x(i)-middle)/t))^(2/(g-1));
            end
        end
        rhostarR = rhoR*(Pstar/PR)^(1/g);
        TstarR = TR*(Pstar/PR)^((g-1)/g);
    end

    for i = 1:length(x)
        if LR == 1 && RS == 1
            if x(i) < xrange_lowL
                u(i) = uL;
                P(i) = PL;
                rho(i) = rhoL;
                T(i) = TL;
            elseif x(i) > xrange_highL && x(i) <= ustar*t+middle
                u(i) = ustar;
                P(i) = Pstar;
                rho(i) = rhostarL;
                T(i) = TstarL;
            elseif x(i) > ustar*t+middle && x(i) <= SR*t+middle
                u(i) = ustar;
                P(i) = Pstar;
                rho(i) = rhostarR;
                T(i) = TstarR;
            elseif x(i) > SR*t+middle
                u(i) = uR;
                P(i) = PR;
                rho(i) = rhoR;
                T(i) = TR;
            end
        end

        if LS == 1 && RR == 1
            if x(i) < SL*t+middle
                u(i) = uL;
                P(i) = PL;
                rho(i) = rhoL;
                T(i) = TL;
            elseif x(i) > SL*t+middle && x(i) <= ustar*t+middle
                u(i) = ustar;
                P(i) = Pstar;
                rho(i) = rhostarL;
                T(i) = TstarL;
            elseif x(i) > ustar*t+middle && x(i) < xrange_highR
                u(i) = ustar;
                P(i) = Pstar;
                rho(i) = rhostarR;
                T(i) = TstarR;
            elseif x(i) > xrange_lowR
                u(i) = uR;
                P(i) = PR;
                rho(i) = rhoR;
                T(i) = TR;
            end
        end

        if LR == 1 && RR == 1
            if x(i) < xrange_lowL
                u(i) = uL;
                P(i) = PL;
                rho(i) = rhoL;
                T(i) = TL;            
%             elseif x(i) > xrange_highL && x(i) <= xrange_lowR
%                 u(i) = (2/(g+1)).*(astarL + ((g-1)/2)*ustarL + x/t);
%                 P(i) = PstarL;
%                 rho(i) = rhoL.*(2/(g+));
%                 T(i) = TstarL;
%             elseif x(i) > xrange_lowR && x(i) <= xrange_highR
%                 u(i) = ustarR;
%                 P(i) = PstarR;
%                 rho(i) = rhostarR;
%                 T(i) = TstarR;
            elseif x(i) > xrange_highL && x(i) <= xrange_highR
                u(i) = ustar;
                P(i) = Pstar;
                rho(i) = rhostarL;
                T(i) = TstarL;
            elseif x(i) > xrange_lowR
                u(i) = uR;
                P(i) = PR;
                rho(i) = rhoR;
                T(i) = TR;
            end
        end

        if LS == 1 && RS == 1
            if x(i) <= middle+SL*t
                u(i) = uL;
                P(i) = PL;
                rho(i) = rhoL;
                T(i) = TL;
            elseif x(i) > middle+SL*t && x(i) <= ustar*t+middle
                u(i) = ustarL;
                P(i) = Pstar;
                rho(i) = rhostarL;
                T(i) = TstarL;
            elseif x(i) > ustar*t+middle && x(i) <= SR*t+middle
                u(i) = ustarR;
                P(i) = Pstar;
                rho(i) = rhostarR;
                T(i) = TstarR;
            elseif x(i) > SR*t+middle
                u(i) = uR;
                P(i) = PR;
                rho(i) = rhoR;
                T(i) = TR;
            end
        end

    end

ie = P./((g-1).*rho);
a = sqrt(g*P./rho);
mach = u./a;
H = 0.5.*u.^2+a.^2./(g-1);
s = log(P./rho.^g);
Uexact = prim_to_cons(g,rho,u,P);

figure(1)
set(gcf,'renderer','painters');
plot(x,rho,'--k', 'LineWidth', 8,'DisplayName','Exact')
hold on;

figure(2)
set(gcf,'renderer','painters');
plot(x,u,'--k', 'LineWidth', 8,'DisplayName','Exact')
hold on;

figure(3)
set(gcf,'renderer','painters');
plot(x,P,'--k', 'LineWidth', 8,'DisplayName','Exact')
hold on;

figure(4)
set(gcf,'renderer','painters');
plot(x,ie,'--k', 'LineWidth', 8,'DisplayName','Exact')
hold on;

if extra_prim_var_plot
    figure(5)
    set(gcf,'renderer','painters');
    plot(x,a,'--k', 'LineWidth', 8,'DisplayName','Exact')
    hold on;

    figure(6)
    set(gcf,'renderer','painters');
    plot(x,mach,'--k', 'LineWidth', 8,'DisplayName','Exact')
    hold on;

    figure(7)
    set(gcf,'renderer','painters');
    plot(x,H,'--k', 'LineWidth', 8,'DisplayName','Exact')
    hold on;

    figure(8)
    set(gcf,'renderer','painters');
    plot(x,s,'--k', 'LineWidth', 8,'DisplayName','Exact')
    hold on;
end
end

function f = funP(Pstar,P,rho,g)
Ar = 2/((g+1)*rho);
Br = P*(g-1)/(g+1);
ar = sqrt(g*P/rho);
if Pstar > P
    f = (Pstar-P)*sqrt(Ar/(Pstar+Br));
else
    f = ((2*ar)/(g-1))*((Pstar/P)^((g-1)/(2*g)) - 1);
end
end

function fp = funPp(Pstar,P,rho,g)
Al = 2/((g+1)*rho);
Bl = P*(g-1)/(g+1);
al = sqrt(g*P/rho);
if P > Pstar
    fp = sqrt(Al/(Bl + P))*(1 - ((Pstar - P)/(2*(Bl + P))));
else
    fp = (1/(P*al))*(Pstar/P)^(-(g+1)/(2*g));
end
end

function [Pstar_root, error] = Bisection_P(PL, PR, uL, uR, rhoL, rhoR, g, du, Iter_Limit, tol)
Pstar0 = Pguess(0, PL, PR, uL, uR, rhoL, rhoR, g, du);
Pstar1 = Pguess(1, PL, PR, uL, uR, rhoL, rhoR, g, du);
Pstar2 = Pguess(2, PL, PR, uL, uR, rhoL, rhoR, g, du);
Pstar3 = Pguess(3, PL, PR, uL, uR, rhoL, rhoR, g, du);
Pstar4 = Pguess(4, PL, PR, uL, uR, rhoL, rhoR, g, du);
Pstara = min([Pstar0,Pstar1,Pstar2,Pstar3,Pstar4]);
Pstarb = max([Pstar0,Pstar1,Pstar2,Pstar3,Pstar4]);
Pstarc = (Pstara+Pstarb)/2;
Iter = 0;
error = 1;
while error > tol && Iter < Iter_Limit
    Iter = Iter + 1;
    fa = funP(Pstara,PL,rhoL,g)+funP(Pstara,PR,rhoR,g)+du;
    fb = funP(Pstarb,PL,rhoL,g)+funP(Pstarb,PR,rhoR,g)+du;
    fc = funP(Pstarc,PL,rhoL,g)+funP(Pstarc,PR,rhoR,g)+du;
    if fb*fa > 0
        disp('Cannot Find Root with Bisection Method');
        Pstar_root = 0;
        error = 100;
        break;
    end
    if fc == 0
        Pstar_root = Pstarc;
        error = 0;
        break;
    end
    if fc*fa > 0
        Pstara = Pstarc;
        Pstarc = (Pstara + Pstarb)/2;
        error = abs(Pstara - Pstarb);
    else
        Pstarb = Pstarc;
        Pstarc = (Pstara + Pstarb)/2;
        error = abs(Pstara - Pstarb);
    end
    Pstar_root = Pstarc;
end
% fprintf('Number of Iterations = %d\n', Iter);
% fprintf('BS-Error = %f\n', error);
end

function Pstar = Pguess(guess_type, PL, PR, uL, uR, rhoL, rhoR, g, du)
if guess_type == 0
    Pstar = 0.5*(PL + PR);
elseif guess_type == 1
    % Two Shock Guess
    TOL = 1e-6;
    aR = sqrt(g*PR/rhoR);
    aL = sqrt(g*PL/rhoL);
    PPV = 0.5*(PL + PR) - 0.125*(uR - uL)*(rhoL + rhoR)*(aL + aR);
    Phat = max([TOL, PPV]);
    AL = 2/((g+1)*rhoL);
    BL = Phat*(g-1)/(g+1);
    gL = sqrt(AL/(Phat + BL));
    AR = 2/((g+1)*rhoR);
    BR = Phat*(g-1)/(g+1);
    gR = sqrt(AR/(Phat + BR));
    PTS = (gL*PL + gR*PR - du)/(gL + gR);
    Pstar = max([TOL, PTS]);
elseif guess_type == 2
    % Primitive Variable
    TOL = 1e-6;
    aR = sqrt(g*PR/rhoR);
    aL = sqrt(g*PL/rhoL);
    PPV = 0.5*(PL + PR) - 0.125*(uR - uL)*(rhoL + rhoR)*(aL + aR);
    Pstar = max([TOL, PPV]);
elseif guess_type == 3
    % Two-Rarefaction
    aR = sqrt(g*PR/rhoR);
    aL = sqrt(g*PL/rhoL);
    PTR = ((aL + aR - 0.5*(g-1)*(uR-uL))/((aL/PL)^((g-1)/(2*g)) + (aR/PR)^((g-1)/(2*g))))^((2*g)/(g-1));
    Pstar = PTR;
else
    % Modified Two Shock Guess
    TOL = 1e-6;
    aR = sqrt(g*PR/rhoR);
    aL = sqrt(g*PL/rhoL);
    PPV = 0.5*(PL + PR) - 0.125*(uR - uL)*(rhoL + rhoR)*(aL + aR);
    Phat = max([TOL, PPV]);
    AL = 2/((g+1)*rhoL);
    BL = Phat*(g-1)/(g+1);
    gL = sqrt(AL/(Phat + BL));
    AR = 2/((g+1)*rhoR);
    BR = Phat*(g-1)/(g+1);
    gR = sqrt(AR/(Phat + BR));
    PTS = (gL*PL + gR*PR - du)/(gL + gR);
    Pstar = max([TOL, PTS]) + PL - PR;
end
end

function [Pstar_root,error] = Newton_Raphson_P(Pstar0,PL,PR,rhoL,rhoR,g,du,iter_limit,tol)
error = 2*tol;
Iter = 0;
Pstar = Pstar0;
while Iter < iter_limit && abs(error) > tol
    Iter = Iter + 1;
    f = funP(Pstar,PL,rhoL,g)+funP(Pstar,PR,rhoR,g)+du;
    fp = funPp(Pstar,PL,rhoL,g)+funPp(Pstar,PR,rhoR,g);
    Pstarkm1 = Pstar;
    Pstar = Pstarkm1-f/fp;
    error = 100*abs(Pstar-Pstarkm1)/(0.5*(Pstar+Pstarkm1));
end
Pstar_root = Pstar;
% fprintf('Number of Iterations = %d\n', Iter);
% fprintf('NR-Error = %f\n', error);
end