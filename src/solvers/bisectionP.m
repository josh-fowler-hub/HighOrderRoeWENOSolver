function [Pstar_root, error,Iter] = bisectionP(PL, PR, uL, uR, rhoL, rhoR, g, du, Iter_Limit, tol)
% bisectionP - Bisection root-finding for pressure function in Riemann solver.

Pstar0 = Pguess(0, PL, PR, uL, uR, rhoL, rhoR, g, du);
Pstar1 = Pguess(1, PL, PR, uL, uR, rhoL, rhoR, g, du);
Pstar2 = Pguess(2, PL, PR, uL, uR, rhoL, rhoR, g, du);
Pstar3 = Pguess(3, PL, PR, uL, uR, rhoL, rhoR, g, du);
Pstar4 = Pguess(4, PL, PR, uL, uR, rhoL, rhoR, g, du);
Pstara = min([Pstar0,Pstar1,Pstar2,Pstar3,Pstar4]);
Pstarb = max([Pstar0,Pstar1,Pstar2,Pstar3,Pstar4]);
Pstarc = (Pstara+Pstarb)/2;
Iter = 0;
error = 1;
while error > tol && Iter < Iter_Limit
    Iter = Iter + 1;
    fa = funP(Pstara,PL,rhoL,g)+funP(Pstara,PR,rhoR,g)+du;
    fb = funP(Pstarb,PL,rhoL,g)+funP(Pstarb,PR,rhoR,g)+du;
    fc = funP(Pstarc,PL,rhoL,g)+funP(Pstarc,PR,rhoR,g)+du;
    if fb*fa > 0
        Pstar_root = Pstar0;
        error = 100;
        break;
    end
    if fc == 0
        Pstar_root = Pstarc;
        error = 0;
        break;
    end
    if fc*fa > 0
        Pstara = Pstarc;
        Pstarc = (Pstara + Pstarb)/2;
        error = abs(Pstara - Pstarb);
    else
        Pstarb = Pstarc;
        Pstarc = (Pstara + Pstarb)/2;
        error = abs(Pstara - Pstarb);
    end
    Pstar_root = Pstarc;
end
end