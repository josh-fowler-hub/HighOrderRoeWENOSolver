function fp = funPp(Pstar,P,rho,g)
% funPp - Derivative of pressure function for Riemann solver root finding.

Al = 2/((g+1)*rho);
Bl = P*(g-1)/(g+1);
al = sqrt(g*P/rho);
if P > Pstar
    fp = sqrt(Al/(Bl + P))*(1 - ((Pstar - P)/(2*(Bl + P))));
else
    fp = (1/(rho*al))*(Pstar/P)^(-(g+1)/(2*g));
end
end