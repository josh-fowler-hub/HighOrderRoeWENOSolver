function Pstar = Pguess(guess_type, PL, PR, uL, uR, rhoL, rhoR, g, du)
% Pguess - Provides initial guesses for pressure solver in Riemann problems.

if guess_type == 0
    Pstar = 0.5*(PL + PR);
elseif guess_type == 1
    % Two Shock Guess
    TOL = 1e-6;
    aR = sqrt(g*PR/rhoR);
    aL = sqrt(g*PL/rhoL);
    PPV = 0.5*(PL + PR) - 0.125*(uR - uL)*(rhoL + rhoR)*(aL + aR);
    Phat = max([TOL, PPV]);
    AL = 2/((g+1)*rhoL);
    BL = Phat*(g-1)/(g+1);
    gL = sqrt(AL/(Phat + BL));
    AR = 2/((g+1)*rhoR);
    BR = Phat*(g-1)/(g+1);
    gR = sqrt(AR/(Phat + BR));
    PTS = (gL*PL + gR*PR - du)/(gL + gR);
    Pstar = max([TOL, PTS]);
elseif guess_type == 2
    % Primitive Variable
    TOL = 1e-6;
    aR = sqrt(g*PR/rhoR);
    aL = sqrt(g*PL/rhoL);
    PPV = 0.5*(PL + PR) - 0.125*(uR - uL)*(rhoL + rhoR)*(aL + aR);
    Pstar = max([TOL, PPV]);
elseif guess_type == 3
    % Two-Rarefaction
    aR = sqrt(g*PR/rhoR);
    aL = sqrt(g*PL/rhoL);
    PTR = ((aL + aR - 0.5*(g-1)*(uR-uL))/((aL/PL)^((g-1)/(2*g)) + (aR/PR)^((g-1)/(2*g))))^((2*g)/(g-1));
    Pstar = PTR;
else
    % Modified Two Shock Guess
    TOL = 1e-6;
    aR = sqrt(g*PR/rhoR);
    aL = sqrt(g*PL/rhoL);
    PPV = 0.5*(PL + PR) - 0.125*(uR - uL)*(rhoL + rhoR)*(aL + aR);
    Phat = max([TOL, PPV]);
    AL = 2/((g+1)*rhoL);
    BL = Phat*(g-1)/(g+1);
    gL = sqrt(AL/(Phat + BL));
    AR = 2/((g+1)*rhoR);
    BR = Phat*(g-1)/(g+1);
    gR = sqrt(AR/(Phat + BR));
    PTS = (gL*PL + gR*PR - du)/(gL + gR);
    Pstar = max([TOL, PTS]) + PL - PR;
end
end