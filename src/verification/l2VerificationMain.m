close all;

cells = 1000;
L = 1;
gamma = 1.4;
ghost = 3;
cfl = 0.9;
flux = 1;
lam = 1;
kap = 1;
rec = 'WENO5';
Order = 1;
solver = 'ROE';


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

    % Select the Time Integration Method based on input
    U = SSPRK65(U,flux,lam,gamma,dx,dt,kap,rec,Order,xs);
    
    % Apply Boundary Conditions
    [ULBC, URBC] = manufacturedSolutionBCs(U);
    U = [ULBC, U, URBC];
    
    error = max(max(abs(U - Uold)));
    if error < 1e-3
        break
    end
end



plotVars(U,Uexact,xs);

function plotVars(U,Uexact,xs,step)
gamma = 1.4;
[r, u, p, a, H, e, m, s] = consToPrim(U(:,4:end-3), gamma);
% [r, u, p, a, H, e, m, s] = consToPrim(U, gamma);
[re, ue, pe, ae, He, ee, me, se] = consToPrim(Uexact, gamma);

figure(1)
plot(xs, r, 'Color', 'k', 'DisplayName', 'simulated');
hold('on');
plot(xs, re, 'Color', 'r', 'DisplayName', 'exact');
legend();
ylim([1, 1.16]);
xlim([0, 1]);
% drawnow();
% filename = ['MMS/density_', num2str(step), '.png'];
% saveas(gcf, filename);

% figure(2)
% plot(xs, ue);
% plot(xs, u);
% 
% figure(3)
% plot(xs, pe);
% plot(xs, p);
% 
% figure(4)
% plot(xs, ae);
% plot(xs, a);
% 
% figure(5)
% plot(xs, He);
% plot(xs, H);
% 
% figure(6)
% plot(xs, ee);
% plot(xs, e);
% 
% figure(7)
% plot(xs, me);
% plot(xs, m);
% 
% figure(8)
% plot(xs, se);
% plot(xs, s);
end



function F_vec = dUdt(U, flux, lam, gamma,dx,kap,rec,CD_Term_Order,xs)
    switch rec
        case 'ROE'
            [URi,ULi,URim1,ULim1] = musclRec(U,18,lam,kap);
%             [ULBC, URBC] = manufacturedSolutionBCs(U);
            BCs = [0 0 0; 0 0 0; 0 0 0];
        case 'MUSCL'
            [URi,ULi,URim1,ULim1] = musclRec(U,flux,lam,kap);
%             [ULBC, URBC] = manufacturedSolutionBCs(U);
            BCs = [0 0 0; 0 0 0; 0 0 0];
        case 'WENO3'
            [URi,ULi] = weno3Rec(U(:,2:end));
            [URim1,ULim1] = weno3Rec(U(:,1:end-1));
%             [ULBC, URBC] = manufacturedSolutionBCs(U);
            BCs = [0 0 0; 0 0 0; 0 0 0];
        case 'WENO5'
            [URi,ULi] = weno5Rec(U(:,2:end));
            [URim1,ULim1] = weno5Rec(U(:,1:end-1));
%             [ULBC, URBC] = manufacturedSolutionBCs(U);
            BCs = [0 0 0; 0 0 0; 0 0 0];
    end
    [Fip12, Fim12] = roeSolver(U,URi,ULi,URim1,ULim1,gamma,CD_Term_Order);
    fx = sourceTerms(xs);
    F_vec = fx -(1/dx)*(Fip12 - Fim12);
    F_vec = [BCs, F_vec, BCs];
end

function U = SSPRK65(U,flux,lam,gamma,dx,dt,kap,rec,CD_Term_Order,xs)
U0 = U;
U1 = U0 + dt/2 *dUdt(U0, flux, lam, gamma,dx,kap,rec,CD_Term_Order,xs);
U2 = U1 + dt/2 *dUdt(U1, flux, lam, gamma,dx,kap,rec,CD_Term_Order,xs);
U3 = U2 + dt/2 *dUdt(U2, flux, lam, gamma,dx,kap,rec,CD_Term_Order,xs);
U4 = U3 + dt/2 *dUdt(U3, flux, lam, gamma,dx,kap,rec,CD_Term_Order,xs);
U5 = U4 + dt/2 *dUdt(U4, flux, lam, gamma,dx,kap,rec,CD_Term_Order,xs);
U6 = (1/9)*U0 + (2/5)*U1 + (4/9)*U3 +(2/45)*U5 + dt/45 *dUdt(U5, flux, lam, gamma,dx,kap,rec,CD_Term_Order,xs);
U = U6(:,4:end-3);
end

function smax = findWaveSpeed(U, gamma)
UL = U(:,1:end-1);
UR = U(:,2:end);
[rL, uL, ~, ~, HL, ~, ~, ~] = consToPrim(UL, gamma);
[rR, uR, ~, ~, HR, ~, ~, ~] = consToPrim(UR, gamma);
% [~, u, ~, a, ~, ~, ~, ~] = consToPrim(U, gamma);
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

function dt = calculateTimeStep(cfl, dx, time, tfinal, U, gamma)
smax = findWaveSpeed(U, gamma);

dt = (cfl*dx)/smax;
if isnan(dt)
    dt = 1e-6;
end

if time + dt > tfinal
    dt = tfinal - time;
end

end
