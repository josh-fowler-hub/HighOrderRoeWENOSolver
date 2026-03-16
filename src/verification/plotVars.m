function plotVars(U,Uexact,xs,step)
% plotVars - Plot selected primitive variables from a simulation and exact solution
% 
% Syntax:
%   plotVars(U, Uexact, xs)
%   plotVars(U, Uexact, xs, step)
% 
% Inputs:
%   U      - simulated conservative variable matrix (3 x N)
%   Uexact - exact/converged conservative variable matrix (3 x N)
%   xs     - grid locations (1 x N)
%   step   - (optional) current time step number (used for saving plots)

if nargin < 4
    step = [];
end

gamma = 1.4;
[r, u, p, a, H, e, m, s] = consToPrim(U(:,4:end-3), gamma);
[re, ue, pe, ae, He, ee, me, se] = consToPrim(Uexact, gamma);

figure(1)
plot(xs, r, 'Color', 'k', 'DisplayName', 'simulated');
hold('on');
plot(xs, re, 'Color', 'r', 'DisplayName', 'exact');
legend();
ylim([1, 1.16]);
xlim([0, 1]);

% Uncomment the blocks below to plot additional variables
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

if ~isempty(step)
    % Example: save figure to disk
    % filename = sprintf('MMS/density_%04d.png', step);
    % saveas(gcf, filename);
end
end