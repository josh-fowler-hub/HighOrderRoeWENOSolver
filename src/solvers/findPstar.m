function Pstar_root = findPstar(Pstar0,PL,PR,rhoL,rhoR,uL,uR,g,du,iter_limit,tol)
% findPstar - Robust root finder for the pressure in a Riemann problem.
% Uses Newton-Raphson with fallback to bisection.

[Pstar_root,error,~] = newtonRaphsonP(Pstar0,PL,PR,rhoL,rhoR,g,du,iter_limit,tol);

if error > tol
    [Pstar_root, error,~] = bisectionP(PL, PR, uL, uR, rhoL, rhoR, g, du, iter_limit, tol);
end
if error > tol
    [Pstar_root,~,~] = newtonRaphsonP(Pstar_root,PL,PR,rhoL,rhoR,g,du,iter_limit,tol);
end
end