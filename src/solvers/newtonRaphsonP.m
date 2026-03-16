function [Pstar_root,error,Iter] = newtonRaphsonP(Pstar0,PL,PR,rhoL,rhoR,g,du,iter_limit,tol)
% newtonRaphsonP - Newton-Raphson root-finder for Riemann problem pressure.

error = 2*tol;
Iter = 0;
Pstar = Pstar0;
while Iter < iter_limit && abs(error) > tol
    Iter = Iter + 1;
    f = funP(Pstar,PL,rhoL,g)+funP(Pstar,PR,rhoR,g)+du;
    fp = funPp(Pstar,PL,rhoL,g)+funPp(Pstar,PR,rhoR,g);
    Pstarkm1 = Pstar;
    Pstar = Pstarkm1 - f/fp;
    error = abs((Pstar-Pstarkm1)/(0.5*(Pstar+Pstarkm1)));
end
Pstar_root = Pstar;
end