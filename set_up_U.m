function [U] = set_up_U(rhoR,rhoL,uL,uR,PR,PL,g,x,xmid,n)
% sets left and right state values (sets up riemann problem)

for i = 1:n
    if x(i) <= xmid
        rho(i) = rhoL;
        u(i) = uL;
        P(i) = PL;
        a(i) = sqrt(g*P(i)/rho(i));
        H(i) = 0.5*u(i)^2+a(i)^2/(g-1);
        
        U(1,i) = rhoL;
        U(2,i) = rhoL*uL;
        U(3,i) = 0.5*rhoL*uL^2+PL/(g-1);
        
        F(1,i) = rhoL*uL;
        F(2,i) = rhoL*uL^2+PL;
        F(3,i) = uL*(rhoL*(0.5*uL^2+PL/((g-1)*rhoL))+PL);
    else
        rho(i) = rhoR;
        u(i) = uR;
        P(i) = PR;
        a(i) = sqrt(g*P(i)/rho(i));
        H(i) = 0.5*u(i)^2+a(i)^2/(g-1);
        
        U(1,i) = rhoR;
        U(2,i) = rhoR*uR;
        U(3,i) = 0.5*rhoR*uR^2+PR/(g-1);
        
        F(1,i) = rhoR*uR;
        F(2,i) = rhoR*uR^2+PR;
        F(3,i) = uR*(rhoR*(0.5*uR^2+PR/((g-1)*rhoR))+PR);
    end
end