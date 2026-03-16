function absA = absEulerJacobian(rho, u, p, gamma)
% absEulerJacobian - compute the absolute flux Jacobian |A| for 1D Euler
%
% Returns a 3x3xN array of |A| matrices for (rho,u,p) states.

% Sound speed and enthalpy
c = sqrt(gamma .* p ./ rho);
H = 0.5 .* u.^2 + gamma ./ (gamma - 1) .* (p ./ rho);

N = size(rho,2);
absA = zeros(3,3,N);

for i = 1:N
    ui = u(i);
    ci = c(i);
    Hi = H(i);

    % Eigenvalues
    lam = [ui-ci, ui, ui+ci];
    absLam = abs(lam);

    % Right eigenvectors
    R = [1,       1,        1;
         ui-ci,   ui,       ui+ci;
         Hi-ui*ci, 0.5*ui^2, Hi+ui*ci];

    % Dissipation matrix = R * |Lambda| * R^{-1}
    absA(:,:,i) = R * diag(absLam) / R;
end
end
