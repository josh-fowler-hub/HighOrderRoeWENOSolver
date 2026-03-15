function U = applyBoundaryConditions(U, ghost, xidx)
% Apply simple constant (zero-gradient) boundary conditions by filling ghost cells
% with the nearest interior state.
% Inputs:
%   U     - full state array including ghost cells (3 x N)
%   ghost - number of ghost cells on each side
%   xidx  - indices of interior (physical) cells
% Output:
%   U     - updated with boundary values

left_BC = U(:, xidx(1));
right_BC = U(:, xidx(end));

% Fill left and right ghost zones with the nearest interior value
U(:, 1:ghost) = repmat(left_BC, 1, ghost);
U(:, end-ghost+1:end) = repmat(right_BC, 1, ghost);
end
