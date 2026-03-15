function [Fip12,Fim12] = roeSolver(U,URi,ULi,URim1,ULim1,g, CD_Term_Order)

% Pre-compute all primitive variables for the entire conservative array at once
[rho, u, P, a, H, e, M, s] = consToPrim(U, g);

% Extract derivatives of primitive variables for the stencil
Uim3 = U(:,1:end-6);
rhoim3 = rho(:,1:end-6);
uim3 = u(:,1:end-6);
Pim3 = P(:,1:end-6);
Fim3 = Fvec(rhoim3,uim3,Pim3,g);

Uim2 = U(:,2:end-5);
rhoim2 = rho(:,2:end-5);
uim2 = u(:,2:end-5);
Pim2 = P(:,2:end-5);
Fim2 = Fvec(rhoim2,uim2,Pim2,g);

Uim1 = U(:,3:end-4);
rhoim1 = rho(:,3:end-4);
uim1 = u(:,3:end-4);
Pim1 = P(:,3:end-4);
Fim1 = Fvec(rhoim1,uim1,Pim1,g);

Ui = U(:,4:end-3);
rhoi = rho(:,4:end-3);
ui = u(:,4:end-3);
Pi = P(:,4:end-3);
Fi = Fvec(rhoi,ui,Pi,g);

Uip1 = U(:,5:end-2);
rhoip1 = rho(:,5:end-2);
uip1 = u(:,5:end-2);
Pip1 = P(:,5:end-2);
Fip1 = Fvec(rhoip1,uip1,Pip1,g);

Uip2 = U(:,6:end-1);
rhoip2 = rho(:,6:end-1);
uip2 = u(:,6:end-1);
Pip2 = P(:,6:end-1);
Fip2 = Fvec(rhoip2,uip2,Pip2,g);

Uip3 = U(:,7:end);
rhoip3 = rho(:,7:end);
uip3 = u(:,7:end);
Pip3 = P(:,7:end);
Fip3 = Fvec(rhoip3,uip3,Pip3,g);

%{
      ULi    URi
|  i-1 |  i   |  i+1 |     |
ULim1 URim1  URip1  ULip1

%}

% Setting Up Flux's and Variables for reconstructed states
% Ui - use reconstructed values
[rhoRi, uRi, PRi, aRi, HRi, ~, ~, ~] = consToPrim(URi, g);
FRi = Fvec(rhoRi,uRi,PRi,g);

[rhoLi, uLi, PLi, ~, HLi, ~, ~, ~] = consToPrim(ULi, g);
FLi = Fvec(rhoLi,uLi,PLi,g);

% Uim1 - use reconstructed values
[rhoRim1, uRim1, PRim1, ~, HRim1, ~, ~, ~] = consToPrim(URim1, g);
FRim1 = Fvec(rhoRim1,uRim1,PRim1,g);

[rhoLim1, uLim1, PLim1, aLim1, HLim1, ~, ~, ~] = consToPrim(ULim1, g);
FLim1 = Fvec(rhoLim1,uLim1,PLim1,g);

% step 1 (calculating roe averaged states)


%     RT = sqrt(rR/rL);
%     r = RT*rL;
%     u = (uL+RT*uR)/(1+RT);
%     H = (HL+RT*HR)/(1+RT);
%     a = sqrt( (g-1)*(H-u*u/2) );

Rim12 = sqrt(rhoRim1./rhoLim1);
rhotim12 = Rim12.*rhoLim1;
utim12 = (uLim1 + Rim12.*uRim1)./(1 + Rim12);
Htim12 = (HLim1 + Rim12.*HRim1)./(1 + Rim12);
atim12_temp = (g-1).*(Htim12 - 0.5.*utim12.*utim12);
atim12_temp(atim12_temp<=0) = 1e-6;
atim12 = atim12_temp.^(1/2); % Get rid of the numerical error causing the complex
% value of speed of sound, inside the sqrt on Test 5 is giving a negative
% number on the order of 1E-14, Htim12 and 0.5.*utim12.^2 are very close
% and should be giving 0 for the vacuum in Test 5.


Rip12 = sqrt(rhoRi./rhoLi);
rhotip12 = Rip12.*rhoLi;
utip12 = (uLi + Rip12.*uRi)./(1 + Rip12);
Htip12 = (HLi + Rip12.*HRi)./(1 + Rip12);
atip12_temp = (g-1).*(Htip12-0.5.*utip12.^2);
atip12_temp(atip12_temp<=0) = 1e-6;
atip12 = atip12_temp.^(1/2);


% step 2 (calculating eigenvalues)
lammt1 = utim12-atim12;
lammt2 = utim12;
lammt3 = utim12+atim12;

lampt1 = utip12-atip12;
lampt2 = utip12;
lampt3 = utip12+atip12;

% step 3 (calculating eigenvectors)

eigen_ones = ones(1, length(utim12));

Kmt1 = [eigen_ones; utim12 - atim12; Htim12 - utim12.*atim12];
Kmt2 = [eigen_ones; utim12; 0.5.*utim12.*utim12];
Kmt3 = [eigen_ones; utim12 + atim12; Htim12 + utim12.*atim12];

eigen_ones = ones(1, length(utip12));

Kpt1 = [eigen_ones; utip12-atip12; Htip12-utip12.*atip12];
Kpt2 = [eigen_ones; utip12; 0.5.*utip12.*utip12];
Kpt3 = [eigen_ones; utip12 + atip12; Htip12 + utip12.*atip12];

% step 4 (calculating wave speeds)

dum1 = URim1(1,:)-ULim1(1,:);
dum2 = URim1(2,:)-ULim1(2,:);
dum3 = URim1(3,:)-ULim1(3,:);

alphamt2 = ((g-1)./atim12.^2).*(dum1.*(Htim12-utim12.^2)+utim12.*dum2-dum3);
alphamt1 = (1./(2.*atim12)).*(dum1.*(utim12+atim12)-dum2-atim12.*alphamt2);
alphamt3 = dum1-(alphamt1+alphamt2);

dup1 = URi(1,:)-ULi(1,:);
dup2 = URi(2,:)-ULi(2,:);
dup3 = URi(3,:)-ULi(3,:);

alphapt2 = ((g-1)./atip12.^2).*(dup1.*(Htip12-utip12.^2)+utip12.*dup2-dup3);
alphapt1 = 1./(2.*atip12).*(dup1.*(utip12+atip12)-dup2-atip12.*alphapt2);
alphapt3 = dup1-(alphapt1+alphapt2);

% step 4.5 entropy fix
[ustar,astarL,astarR] = adaptiveRiemannSolver(rhoLim1,rhoRi,uLim1,uRi,aLim1,aRi,PLim1,PRi,g);
[lammt1,lampt3] = entropyFix(uLim1,uRi,aLim1,aRi,ustar,astarL,astarR,lammt1,lampt3);

% Central Differencing Term
if CD_Term_Order == 1
    CDim12 = 0.5.*(FLim1+FRim1);
    CDip12 = 0.5.*(FLi+FRi);
elseif CD_Term_Order == 2
    CDim12 = 0.5.*(Fim1+Fi);
    CDip12 = 0.5.*(Fi+Fip1);
elseif CD_Term_Order == 4
    CDim12 = (1/16).*(-Fim2 + 9*Fim1 + 9*Fi - Fip1);
    CDip12 = (1/16).*(-Fim1 + 9*Fi + 9*Fip1 - Fip2);
elseif CD_Term_Order == 6
    CDim12 = (1/256).*(3*Fip2 - 25*Fip1 + 150*Fi + 150*Fim1 - 25*Fim2 + 3*Fim3);
    CDip12 = (1/256).*(3*Fip3 - 25*Fip2 + 150*Fip1 + 150*Fi - 25*Fim1 + 3*Fim2);
end

% step 5 (assembling roe flux)
DFim12 = 0.5.*(alphamt1.*abs(lammt1).*Kmt1 + alphamt2.*abs(lammt2).*Kmt2 +alphamt3.*abs(lammt3).*Kmt3);
DFip12 = 0.5.*(alphapt1.*abs(lampt1).*Kpt1 + alphapt2.*abs(lampt2).*Kpt2 + alphapt3.*abs(lampt3).*Kpt3);
Fim12 = CDim12-DFim12;
Fip12 = CDip12-DFip12;
end

function [ustar,astarL,astarR] = adaptiveRiemannSolver(rhoLim1,rhoRi,uLim1,uRi,aLim1,aRi,PLim1,PRi,g) % Toro 9.5.2
iter_limit = 100;
tol = 1e-7;
du = uRi - uLim1;
rhobar = 0.5*(rhoLim1 + rhoRi);
abar = 0.5*(aLim1 + aRi);
Pstari = 0.5*(PLim1 + PRi) + 0.5*(uLim1 - uRi).*rhobar.*abar;
Pstari(Pstari < 0) = 0;
Pmax = PRi;
Pmax(Pmax<PLim1) = PLim1(Pmax<PLim1);
Pmin = PRi;
Pmin(Pmin>PLim1) = PLim1(Pmin>PLim1);
Q = Pmax./Pmin;
Quser = 1;

PVRSidx1 = find(Q < Quser);
PVRSidx2 = find(Pstari > Pmin);
PVRSidx3 = find(Pstari < Pmax);
PVRSidx4 = PVRSidx1(ismember(PVRSidx1,PVRSidx2));
PVRSidx = PVRSidx3(ismember(PVRSidx3,PVRSidx4));
TRRSidx = find(Pstari < Pmin);
allidxssofar = [PVRSidx, TRRSidx];
allidxs = 1:length(Q);
idx = ~ismember(allidxs,allidxssofar);
TSRSidx = allidxs(idx);

Pstar = 0*Pstari;
ustar = zeros(size(uLim1));
rhostarL = zeros(size(rhoLim1));
rhostarR = zeros(size(rhoLim1));
astarL_temp = zeros(size(aLim1));
astarR_temp = zeros(size(aLim1));
astarL = zeros(size(aLim1));
astarR = zeros(size(aLim1));
ustarR = zeros(size(uLim1));
ustarL = zeros(size(uLim1));

%PVRS Toro 9.3
Pstar(PVRSidx) = 0.5.*(PLim1(PVRSidx) + PRi(PVRSidx)) + 0.5.*(uLim1(PVRSidx) - uRi(PVRSidx)).*rhobar(PVRSidx).*abar(PVRSidx);
for i = PVRSidx
    Pstar(i) = Newton_Raphson_P(Pstar(i),PLim1(i),PRi(i),rhoLim1(i),rhoRi(i),g,du(i),iter_limit,tol);
end
ustar(PVRSidx) = real(0.5.*(uLim1(PVRSidx) + uRi(PVRSidx)) + 0.5.*(PLim1(PVRSidx) - PRi(PVRSidx))./(rhobar(PVRSidx).*abar(PVRSidx)));
rhostarL(PVRSidx) = rhoLim1(PVRSidx) + (uLim1(PVRSidx) - ustar(PVRSidx)).*(rhobar(PVRSidx)./abar(PVRSidx));
rhostarR(PVRSidx) = rhoRi(PVRSidx) + (ustar(PVRSidx) - uRi(PVRSidx)).*(rhobar(PVRSidx)./abar(PVRSidx));
astarL_temp(PVRSidx) = g.*Pstar(PVRSidx)./(rhostarL(PVRSidx));
astarL_temp(astarL_temp<=0) = 1e-6;
astarR_temp(PVRSidx) = g.*Pstar(PVRSidx)./(rhostarR(PVRSidx));
astarR_temp(astarR_temp<=0) = 1e-6;
astarL(PVRSidx) = real(sqrt(astarL_temp(PVRSidx)));
astarR(PVRSidx) = real(sqrt(astarR_temp(PVRSidx)));

%TRRS Toro 9.4.1
z = (g-1)./(2*g);
Pstar(TRRSidx) = ((aLim1(TRRSidx) + aRi(TRRSidx) - z*g*(uRi(TRRSidx) - uLim1(TRRSidx)))./((aLim1(TRRSidx)./(PLim1(TRRSidx).^z)) + (aRi(TRRSidx)./(PRi(TRRSidx).^z)))).^(1/z);
for i = TRRSidx
    Pstar(i) = Newton_Raphson_P(Pstar(i),PLim1(i),PRi(i),rhoLim1(i),rhoRi(i),g,du(i),iter_limit,tol);
end
astarL_temp(TRRSidx) = Pstar(TRRSidx)./PLim1(TRRSidx);
astarL_temp(astarL_temp<=0) = 1e-6;
astarR_temp(TRRSidx) = Pstar(TRRSidx)./PRi(TRRSidx);
astarR_temp(astarR_temp<=0) = 1e-6;
astarL(TRRSidx) = real(aLim1(TRRSidx).*(astarL_temp(TRRSidx)).^z);
astarR(TRRSidx) = real(aRi(TRRSidx).*(astarR_temp(TRRSidx)).^z);
ustarL(TRRSidx) = uLim1(TRRSidx) + (1/(g*z))*(aLim1(TRRSidx) - astarL(TRRSidx));
ustarR(TRRSidx) = uRi(TRRSidx) + (1/(g*z))*(astarR(TRRSidx) - aRi(TRRSidx));
ustar(TRRSidx) = real(0.5*(ustarL(TRRSidx) + ustarR(TRRSidx)));

TOL = 1e-6;
%TSRS Toro 9.4.2
z = (g-1)/(2*g);
AR = 2./((g+1)*rhoRi(TSRSidx));
AL = 2./((g+1)*rhoLim1(TSRSidx));
BR = ((g-1)/(g+1))*PRi(TSRSidx);
BL = ((g-1)/(g+1))*PLim1(TSRSidx);
gR = @(p) sqrt(AR./(p + BR));
gL = @(p) sqrt(AL./(p + BL));
P0(TSRSidx) = 0.5.*(PLim1(TSRSidx) + PRi(TSRSidx)) + 0.5.*(uLim1(TSRSidx) - uRi(TSRSidx)).*rhobar(TSRSidx).*abar(TSRSidx);
P0(P0<=0) = TOL;
Pstar(TSRSidx) = (gL(P0(TSRSidx)).*PLim1(TSRSidx) + gR(P0(TSRSidx)).*PRi(TSRSidx) - (uRi(TSRSidx) - uLim1(TSRSidx)))./(gL(P0(TSRSidx)) + gR(P0(TSRSidx)));
for i = TSRSidx
    Pstar(i) = Newton_Raphson_P(Pstar(i),PLim1(i),PRi(i),rhoLim1(i),rhoRi(i),g,du(i),iter_limit,tol);
end
ustar(TSRSidx) = real(0.5.*(uLim1(TSRSidx) + uRi(TSRSidx)) + 0.5.*((Pstar(TSRSidx)-PRi(TSRSidx)).*gR(P0(TSRSidx)) - (Pstar(TSRSidx)-PLim1(TSRSidx)).*gL(P0(TSRSidx))));
rhostarL(TSRSidx) = rhoLim1(TSRSidx).*(((Pstar(TSRSidx)./PLim1(TSRSidx))+((g-1)./(g+1)))./(((g-1)./(g+1)).*(Pstar(TSRSidx)./PLim1(TSRSidx)) + 1));
rhostarR(TSRSidx) = rhoRi(TSRSidx).*(((Pstar(TSRSidx)./PRi(TSRSidx))+((g-1)./(g+1)))./(((g-1)./(g+1)).*(Pstar(TSRSidx)./PRi(TSRSidx)) + 1));
astarL_temp(TSRSidx) = g.*Pstar(TSRSidx)./(rhostarL(TSRSidx));
astarL_temp(astarL_temp<=0) = 1e-6;
astarR_temp(TSRSidx) = g.*Pstar(TSRSidx)./(rhostarR(TSRSidx));
astarR_temp(astarR_temp<=0) = 1e-6;
astarL(TSRSidx) = real(sqrt(astarL_temp(TSRSidx)));
astarR(TSRSidx) = real(sqrt(astarR_temp(TSRSidx)));

end


function [lammt1,lampt3] = entropyFix(uLim1,uRi,aLim1,aRi,ustar,astarL,astarR,lammt1,lampt3)

lamm1L = uLim1 - aLim1;
lamm1R = ustar - astarL;
lamp3L = ustar + astarR;
lamp3R = uRi + aRi;

% Left Transonic Rarefaction: lambda1L < 0 < lambda1R
LTRidx1 = find(lamm1L < 0);
LTRidx2 = find(lamm1R > 0);
if length(LTRidx1) < length(LTRidx2)
    LTRidx = LTRidx1(ismember(LTRidx1,LTRidx2));
else
    LTRidx = LTRidx2(ismember(LTRidx2,LTRidx1));
end
lammt1(LTRidx) = lamm1L(LTRidx).*((lamm1R(LTRidx) - lammt1(LTRidx))./(lamm1R(LTRidx) - lamm1L(LTRidx)));

% Right Transonic Rarefaction: lambda3L < 0 < lambda3R
RTRidx1 = find(lamp3L < 0);
RTRidx2 = find(lamp3R > 0);
if length(RTRidx1) < length(RTRidx2)
    RTRidx = RTRidx1(ismember(RTRidx1,RTRidx2));
else
    RTRidx = RTRidx2(ismember(RTRidx2,RTRidx1));
end
lampt3(RTRidx) = lamp3R(RTRidx).*((lampt3(RTRidx) - lamp3L(RTRidx))./(lamp3R(RTRidx) - lamp3L(RTRidx)));

end

function Pstar_root = findPstar(Pstar0,PL,PR,rhoL,rhoR,uL,uR,g,du,iter_limit,tol)
[Pstar_root,error,Iter] = newtonRaphsonP(Pstar0,PL,PR,rhoL,rhoR,g,du,iter_limit,tol);
method = 1;

if error > tol
    [Pstar_root, error,Iter] = bisectionP(PL, PR, uL, uR, rhoL, rhoR, g, du, iter_limit, tol);
    method = 2;
end
if error > tol
    [Pstar_root,error,Iter] = newtonRaphsonP(Pstar_root,PL,PR,rhoL,rhoR,g,du,iter_limit,tol);
    method = 1;
end
% if method == 1
%     disp("Newton-Raphson Method");
% else
%     disp("Bisection Method");
% end
% fprintf('Number of Iterations = %d\n', Iter);
% fprintf('NR-Error = %f\n', error);
end

function [Pstar_root, error,Iter] = bisectionP(PL, PR, uL, uR, rhoL, rhoR, g, du, Iter_Limit, tol)
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
%         disp('Cannot Find Root with Bisection Method');
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
% fprintf('Number of Iterations = %d\n', Iter);
% fprintf('BS-Error = %f\n', error);
end

function [Pstar_root,error,Iter] = newtonRaphsonP(Pstar0,PL,PR,rhoL,rhoR,g,du,iter_limit,tol)
error = 2*tol;
Iter = 0;
Pstar = Pstar0;
while Iter < iter_limit && abs(error) > tol
    Iter = Iter + 1;
    f = funP(Pstar,PL,rhoL,g)+funP(Pstar,PR,rhoR,g)+du;
    fp = funPp(Pstar,PL,rhoL,g)+funPp(Pstar,PR,rhoR,g);
    Pstarkm1 = Pstar;
    Pstar = Pstarkm1-f/fp;
    error = abs((Pstar-Pstarkm1)/(0.5*(Pstar+Pstarkm1)));
end
Pstar_root = Pstar;
% fprintf('Number of Iterations = %d\n', Iter);
% fprintf('NR-Error = %f\n', error);
end

function f = funP(Pstar,P,rho,g)
Ar = 2/((g+1)*rho);
Br = P*(g-1)/(g+1);
ar = sqrt(g*P/rho);
if Pstar > P
    f = (Pstar-P)*sqrt(Ar/(Pstar+Br));
else
    f = ((2*ar)/(g-1))*((Pstar/P)^((g-1)/(2*g)) - 1);
end
end

function fp = funPp(Pstar,P,rho,g)
Al = 2/((g+1)*rho);
Bl = P*(g-1)/(g+1);
al = sqrt(g*P/rho);
if P > Pstar
    fp = sqrt(Al/(Bl + P))*(1 - ((Pstar - P)/(2*(Bl + P))));
else
    fp = (1/(rho*al))*(Pstar/P)^(-(g+1)/(2*g));
end
end

function Pstar = Pguess(guess_type, PL, PR, uL, uR, rhoL, rhoR, g, du)
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