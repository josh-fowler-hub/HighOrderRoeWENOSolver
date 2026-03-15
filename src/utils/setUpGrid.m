function [xs, dx, n2, xidx] = setUpGrid(n, ghost, xstart, xend)
% Set up computational grid including ghost cells
xs = linspace(xstart, xend, n); % x-domain in physical cells
dx = (xs(2) - xs(1));
left_bound = linspace(-ghost*dx, -dx, ghost);
right_bound = linspace(xs(end)+dx, xs(end)+ghost*dx, ghost);

xs = [left_bound, xs, right_bound];
n2 = n + 2*ghost; % total number of cells (including ghosts)
xidx = ghost+1:n+ghost; % indices for physical domain
end
