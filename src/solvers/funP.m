function f = funP(Pstar,P,rho,g)
% funP - Pressure function for Riemann solver root finding.

Ar = 2/((g+1)*rho);
Br = P*(g-1)/(g+1);
ar = sqrt(g*P/rho);
if Pstar > P
    f = (Pstar-P)*sqrt(Ar/(Pstar+Br));
else
    f = ((2*ar)/(g-1))*((Pstar/P)^((g-1)/(2*g)) - 1);
end
end