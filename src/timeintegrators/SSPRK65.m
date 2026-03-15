function U = SSPRK65(U,flux,lam,gamma,dx,dt,kap,rec,CD_Term_Order,riemann_solver,xs,time)
% SSPRK(6,5) strong stability preserving Runge-Kutta
% Maintains a global counter used for debugging/analysis
global time_int_counter;
time_int_counter = time_int_counter + 1;

U0 = U;
U1 = U0 + dt/2 *dUdt(U0, flux, lam, gamma,dx,kap,rec,CD_Term_Order,riemann_solver,xs,time);
U2 = U1 + dt/2 *dUdt(U1, flux, lam, gamma,dx,kap,rec,CD_Term_Order,riemann_solver,xs,time);
U3 = U2 + dt/2 *dUdt(U2, flux, lam, gamma,dx,kap,rec,CD_Term_Order,riemann_solver,xs,time);
U4 = U3 + dt/2 *dUdt(U3, flux, lam, gamma,dx,kap,rec,CD_Term_Order,riemann_solver,xs,time);
U5 = U4 + dt/2 *dUdt(U4, flux, lam, gamma,dx,kap,rec,CD_Term_Order,riemann_solver,xs,time);
U6 = (1/9)*U0 + (2/5)*U1 + (4/9)*U3 +(2/45)*U5 + dt/45 *dUdt(U5, flux, lam, gamma,dx,kap,rec,CD_Term_Order,riemann_solver,xs,time);
U = U6;
end
