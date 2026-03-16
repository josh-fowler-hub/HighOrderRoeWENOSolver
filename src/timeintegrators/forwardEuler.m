function U = forwardEuler(U,flux,lam,gamma,dx,dt,kappa,rec,CD_Term_Order,riemann_solver,xs,time)
% forwardEuler - First-order explicit Euler time integration.
% Performs a single time step: U^{n+1} = U^n + dt * dU/dt.

% Compute the spatial derivative (flux residual)
F = computeFluxDerivative(U,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);

% Advance the solution
U = U + dt * F;
end
