function U = forwardEuler(U,flux,lam,gamma,dx,dt,kappa,rec,CD_Term_Order,riemann_solver,xs,time)
% Forward Euler time integration (1st order)
% Performs one stage of Forward Euler with the supplied flux function.
for i = 1:10
    dUstar = dUdt(U,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
    U = U + dUstar*dt;
end
end
