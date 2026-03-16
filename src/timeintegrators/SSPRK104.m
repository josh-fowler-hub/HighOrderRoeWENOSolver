function U = SSPRK104(U,flux,lam,gamma,dx,dt,kappa,rec,CD_Term_Order,riemann_solver,xs,time)
% SSPRK(10,4) strong stability preserving Runge-Kutta
U0 = U;
U1 = U0 + dt/6 *computeFluxDerivative(U0, flux, lam, gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
U2 = U1 + dt/6 * computeFluxDerivative(U1,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
U3 = U2 + dt/6 * computeFluxDerivative(U2,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
U4 = U3 + dt/6 * computeFluxDerivative(U3,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
U5 = U4 + dt/6 * computeFluxDerivative(U4,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
Ustar = 1/25 * U0 + 9/25 * U5;
U5 = 15 * Ustar - 5 * U5;
U6 = U5 + dt/6 * computeFluxDerivative(U5,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
U7 = U6 + dt/6 * computeFluxDerivative(U6,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
U8 = U7 + dt/6 * computeFluxDerivative(U7,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
U9 = U8 + dt/6 * computeFluxDerivative(U8,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
U = 6/10 * U9 + dt/10 * computeFluxDerivative(U9,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time) + Ustar;
end
