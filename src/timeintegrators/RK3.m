function U = RK3(U,flux,lam,gamma,dx,dt,kappa,rec,CD_Term_Order,riemann_solver,xs,time)
% 3rd order Runge-Kutta
k1 = U + dt*computeFluxDerivative(U,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
k2 = 0.75*U + 0.25*k1 + 0.25*dt*computeFluxDerivative(k1,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
U = (1/3)*U + (2/3)*k2 + (2/3)*dt*computeFluxDerivative(k2,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
end
