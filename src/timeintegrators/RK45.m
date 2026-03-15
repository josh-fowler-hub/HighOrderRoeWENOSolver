function U = RK45(U,flux,lam,gamma,dx,dt,kappa,rec,CD_Term_Order,riemann_solver,xs,time)
% RK4(5) method (embedded) for time stepping
U0 = U;
K0 = dt*dUdt(U0,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
U1 = U0 + (13736793/35065003)*K0;
K1 = dt*dUdt(U1,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
U2 = (18384209/41371354)*U0 + (22987145/41371354)*U1 ...
     + (1106722/3004045)*K1;
K2 = dt*dUdt(U2,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
U3 = (16461287/26546102)*U0 + (10084815/26546102)*U2 ...
     + (9149709/36323969)*K2;
K3 = dt*dUdt(U3,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
U4 = (1818641/10212497)*U0 + (8393856/10212497)*U3 ...
     + (11974013/21971684)*K3;
K4 = dt*dUdt(U4,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time);
U = (12153930/23498039)*U2 + (4482614/46664871)*U3 ...
    + (9148849/143640986)*K3 ...
    + (8047189/20809438)*U4 ...
    + (9169579/40572015)*K4;
end
