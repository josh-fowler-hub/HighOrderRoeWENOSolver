function U = RK2(U,flux,lam,gamma,dx,dt,kappa,rec,CD_Term_Order,riemann_solver,xs,time)
% 2nd order Runge-Kutta
k1 = U + dt*computeFluxDerivative(U,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
U = 0.5*U + 0.5*k1 + 0.5*dt*computeFluxDerivative(k1,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
end
