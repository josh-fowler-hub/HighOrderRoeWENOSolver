function [Fip12,Fim12] = osherSolver(U,URi,ULi,URim1,ULim1,xs,gamma,time)
% osherSolver - Osher numerical flux for 1D Euler equations
%
% This implementation follows an Osher-style flux construction using a
% Roe-linearization with an entropy fix (Harten–Hyman). It builds the
% dissipation matrix D = R * |Lambda| * R^{-1} and uses it to compute
% the interface flux.
%
% Inputs (all 3xN matrices):
%   U     - full conservative state array (unused, kept for signature)
%   URi   - right reconstructed state at interface
%   ULi   - left reconstructed state at interface
%   URim1 - right reconstructed state at previous interface
%   ULim1 - left reconstructed state at previous interface
%   xs, time - unused
%   gamma - ratio of specific heats
%
% Outputs:
%   Fip12, Fim12 - the numerical flux on the right and left side of the interface

% Left and right states at the interface
UL = ULi;
UR = URi;

% Convert to primitive variables
[rhoL, uL, pL, ~, ~, ~, ~, ~] = consToPrim(UL, gamma);
[rhoR, uR, pR, ~, ~, ~, ~, ~] = consToPrim(UR, gamma);

% Physical fluxes for left and right states
FL = Fvec(rhoL, uL, pL, gamma);
FR = Fvec(rhoR, uR, pR, gamma);

% Difference in conservative variables
dU = UR - UL;

% Two-point Gauss quadrature for the Osher path integral
s = [0.211324865405187, 0.788675134594813];
D = zeros(3,3,size(UL,2));
for si = s
    % State along the straight-line path between UL and UR
    Upath = UL + si .* dU;

    % Primitive variables for the path state
    [rhoP, uP, pP, ~, ~, ~, ~, ~] = consToPrim(Upath, gamma);

    % Compute |A| for each interface
    absA = absEulerJacobian(rhoP, uP, pP, gamma);

    % Accumulate integral approximation
    D = D + absA;
end

% Average the Gauss points (approximate integral)
D = D * 0.5;

% Compute dissipation term: 0.5 * D * (UR-UL)
Fdiff = zeros(size(UL));
for i = 1:size(UL,2)
    Fdiff(:,i) = 0.5 .* (D(:,:,i) * dU(:,i));
end

% Final Osher flux
Fip12 = 0.5 .* (FL + FR) - Fdiff;
Fim12 = Fip12;
end
