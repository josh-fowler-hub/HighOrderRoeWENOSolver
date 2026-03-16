function cfg = plotConfig(test_num, CD_Term_Order)
% plotConfig - Returns plotting configuration for a given test and CD order.
%
% This provides a centralized table-driven configuration for figure
% appearance, axis limits, file names, and labels.

% Define variables to plot (order corresponds to figure numbers)
vars = {
    struct('name','Density','suffix','Density','ylabel','$\\rho$'),
    struct('name','Velocity','suffix','Velocity','ylabel','$u$'),
    struct('name','Pressure','suffix','Pressure','ylabel','$p$'),
    struct('name','InternalEnergy','suffix','InternalEnergy','ylabel','$e$'),
    struct('name','SpeedOfSound','suffix','SpeedofSound','ylabel','$a$'),
    struct('name','Mach','suffix','Mach','ylabel','$M$'),
    struct('name','Enthalpy','suffix','Enthalpy','ylabel','$H$'),
    struct('name','Entropy','suffix','Entropy','ylabel','$s$')
};

% Populate y-axis limits for each variable and CD term order.
% Each entry is a 19x2 matrix: rows correspond to test_num (1..19), cols are [ymin, ymax].
%   - idx 1: CD order 2
%   - idx 2: CD order 4
%   - idx 3: other
ylims = getYlims();

% Map CD_Term_Order to index
switch CD_Term_Order
    case 2
        cdIdx = 1;
    case 4
        cdIdx = 2;
    otherwise
        cdIdx = 3;
end

% Apply per-variable limits
for i = 1:numel(vars)
    vars{i}.ylims = ylims{i}(:,:,cdIdx);
end

% Legend locations (kept as before)
locs = getLocations(test_num, CD_Term_Order);

cfg.vars = vars;
cfg.locs = locs;
end

function ylims = getYlims()
% ylims{i} is a (19 x 2 x 3) array for variable i (figure i),
% where the third dimension selects CD term order option:
%   1 => CD_Term_Order==2
%   2 => CD_Term_Order==4
%   3 => otherwise

% For each variable we store [ymin, ymax] arrays for all test cases.
% This mirrors the previous hard-coded vectors in configurePlots.

% Density
ylims{1} = cat(3, ...
    [0,1.1; 0,0; 0,7; 4,35; 0,7; 0,1.1; 0,0; 0,1.2; 0,0; 0,0; 0,0; 0,0; 0.15,1.45; 0,0; 0,1.1; 0,1.1; 0,0; 0,0; 0,0], ... % CD2
    [0,1.1; 0,0; 0,7; 4,35; 0,7; 0,1.1; 0,0; 0,1.16; 0,0; 0,0; 0,0; 0,0; 0.2,1.4; 0,0; 0,1.1; 0,1.1; 0,0; 0,0; 0,0], ... % CD4
    [0,1.1; 0,0; 0,7; 4,35; 0,7; 0,1.1; 0,0; 0,1.16; 0,0; 0,0; 0,0; 0,0; 0.2,1.4; 0,0; 0,1.1; 0,1.1; 0,0; 0,0; 0,0] ... % other
);

% Velocity
ylims{2} = cat(3, ...
    [-0.2,1.6; 0,0; -3,23; -10,25; -24,4; -0.05,1; 0,0; -0.55,0.05; 0,0; 0,0; 0,0; 0,0; -0.1,1.7; 0,0; -0.2,1.8; -2.1,-0.85; 0,0; 0,0; 0,0], ...
    [-0.2,1.6; 0,0; -3,23; -10,25; -24,4; -0.05,1; 0,0; -0.475,0.05; 0,0; 0,0; 0,0; 0,0; -0.1,1.6; 0,0; -0.2,1.8; -2.15,-1; 0,0; 0,0; 0,0], ...
    [-0.2,1.6; 0,0; -3,23; -10,25; -24,4; -0.05,1; 0,0; -0.475,0.05; 0,0; 0,0; 0,0; 0,0; -0.1,1.6; 0,0; -0.2,1.8; -2.1,-1; 0,0; 0,0; 0,0] ...
);

% Pressure
ylims{3} = cat(3, ...
    [0,1.1; 0,0; -100,1100; -100,1950; -100,1100; 0,1.1; 0,0; 6.75,10.25; 0,0; 0,0; 0,0; 0,0; 0.5,3.75; 0,0; 0,1.1; 0,1.1; 0,0; 0,0; 0,0], ...
    [0,1.1; 0,0; -100,1100; -100,1900; -100,1100; 0,1.1; 0,0; 6.75,10.25; 0,0; 0,0; 0,0; 0,0; 0.5,3.75; 0,0; 0,1.1; 0,1.1; 0,0; 0,0; 0,0], ...
    [0,1.1; 0,0; -100,1100; -100,1900; -100,1100; 0,1.1; 0,0; 6.75,10.25; 0,0; 0,0; 0,0; 0,0; 0.5,3.75; 0,0; 0,1.1; 0,1.1; 0,0; 0,0; 0,0] ...
);

% Internal energy
ylims{4} = cat(3, ...
    [1.8,3.7; 0,0; -200,2700; -20,400; -100,2600; 1.6,3; 0,0; 17,25.5; 0,0; 0,0; 0,0; 0,0; 2,32; 0,0; 1,5; 1.6,3.1; 0,0; 0,0; 0,0], ...
    [1.8,3.7; 0,0; -200,2700; -20,340; -100,2600; 1.6,3; 0,0; 17,25.5; 0,0; 0,0; 0,0; 0,0; 2,20; 0,0; 1,5; 1.6,3; 0,0; 0,0; 0,0], ...
    [1.8,3.7; 0,0; -200,2700; -20,340; -100,2600; 1.6,3; 0,0; 17,25.5; 0,0; 0,0; 0,0; 0,0; 2,20; 0,0; 1,5; 1.75,3; 0,0; 0,0; 0,0] ...
);

% The remaining plots (5-8) have no limits pre-configured and rely on autoscale
for i = 5:8
    ylims{i} = nan(19,2,3);
end
end
