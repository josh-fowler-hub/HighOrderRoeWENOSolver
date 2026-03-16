function [ustar,astarL,astarR] = adaptiveRiemannSolver(rhoLim1,rhoRi,uLim1,uRi,aLim1,aRi,PLim1,PRi,g)
% adaptiveRiemannSolver - Selects a robust Riemann solver path based on wave
% strength. Uses the Toro formulation (9.3 - 9.5).

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

% PVRS (Toro 9.3)
Pstar(PVRSidx) = 0.5.*(PLim1(PVRSidx) + PRi(PVRSidx)) + 0.5.*(uLim1(PVRSidx) - uRi(PVRSidx)).*rhobar(PVRSidx).*abar(PVRSidx);
for i = PVRSidx
    Pstar(i) = newtonRaphsonP(Pstar(i),PLim1(i),PRi(i),rhoLim1(i),rhoRi(i),g,du(i),iter_limit,tol);
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

% TRRS (Toro 9.4.1)
z = (g-1)./(2*g);
Pstar(TRRSidx) = ((aLim1(TRRSidx) + aRi(TRRSidx) - z*g*(uRi(TRRSidx) - uLim1(TRRSidx)))./((aLim1(TRRSidx)./(PLim1(TRRSidx).^z)) + (aRi(TRRSidx)./(PRi(TRRSidx).^z)))).^(1/z);
for i = TRRSidx
    Pstar(i) = newtonRaphsonP(Pstar(i),PLim1(i),PRi(i),rhoLim1(i),rhoRi(i),g,du(i),iter_limit,tol);
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

% TSRS (Toro 9.4.2)
z = (g-1)/(2*g);
AR = 2./((g+1)*rhoRi(TSRSidx));
AL = 2./((g+1)*rhoLim1(TSRSidx));
BR = ((g-1)/(g+1))*PRi(TSRSidx);
BL = ((g-1)/(g+1))*PLim1(TSRSidx);
gR = @(p) sqrt(AR./(p + BR));
gL = @(p) sqrt(AL./(p + BL));
P0(TSRSidx) = 0.5.*(PLim1(TSRSidx) + PRi(TSRSidx)) + 0.5.*(uLim1(TSRSidx) - uRi(TSRSidx)).*rhobar(TSRSidx).*abar(TSRSidx);
P0(P0<=0) = 1e-6;
Pstar(TSRSidx) = (gL(P0(TSRSidx)).*PLim1(TSRSidx) + gR(P0(TSRSidx)).*PRi(TSRSidx) - (uRi(TSRSidx) - uLim1(TSRSidx)))./(gL(P0(TSRSidx)) + gR(P0(TSRSidx)));
for i = TSRSidx
    Pstar(i) = newtonRaphsonP(Pstar(i),PLim1(i),PRi(i),rhoLim1(i),rhoRi(i),g,du(i),iter_limit,tol);
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
