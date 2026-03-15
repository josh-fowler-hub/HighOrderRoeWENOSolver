function [Fip12, Fim12] = HLL_solver(U,URi,ULi,URim1,ULim1,xs,g,t)

%{
      ULi    URi
|  i-1 |  i   |  i+1 |     |
ULim1 URim1  URip1  ULip1

%}

Uim3 = U(:,1:end-6);
[rhoim3, uim3, Pim3, ~, ~, ~, ~, ~] = cons_to_prim(Uim3, g);
Fim3 = Fvec(rhoim3,uim3,Pim3,g);
Uim2 = U(:,2:end-5);
[rhoim2, uim2, Pim2, ~, ~, ~, ~, ~] = cons_to_prim(Uim2, g);
Fim2 = Fvec(rhoim2,uim2,Pim2,g);
Uim1 = U(:,3:end-4);
[rhoim1, uim1, Pim1, ~, ~, ~, ~, ~] = cons_to_prim(Uim1, g);
Fim1 = Fvec(rhoim1,uim1,Pim1,g);
Ui = U(:,4:end-3);
[rhoi, ui, Pi, ~, ~, ~, ~, ~] = cons_to_prim(Ui, g);
Fi = Fvec(rhoi,ui,Pi,g);
Uip1 = U(:,5:end-2);
[rhoip1, uip1, Pip1, ~, ~, ~, ~, ~] = cons_to_prim(Uip1, g);
Fip1 = Fvec(rhoip1,uip1,Pip1,g);
Uip2 = U(:,6:end-1);
[rhoip2, uip2, Pip2, ~, ~, ~, ~, ~] = cons_to_prim(Uip2, g);
Fip2 = Fvec(rhoip2,uip2,Pip2,g);
Uip3 = U(:,7:end);
[rhoip3, uip3, Pip3, ~, ~, ~, ~, ~] = cons_to_prim(Uip3, g);
Fip3 = Fvec(rhoip3,uip3,Pip3,g);


% find wave speeds
% Ui
[rhoRi, uRi, PRi, aRi, HRi, ~, ~, ~] = cons_to_prim(URi, g);

FRi = Fvec(rhoRi,uRi,PRi,g);

[rhoLi, uLi, PLi, ~, HLi, ~, ~, ~] = cons_to_prim(ULi, g);

FLi = Fvec(rhoLi,uLi,PLi,g);

% Uim1
[rhoRim1, uRim1, PRim1, ~, HRim1, ~, ~, ~] = cons_to_prim(URim1, g);

FRim1 = Fvec(rhoRim1,uRim1,PRim1,g);

[rhoLim1, uLim1, PLim1, aLim1, HLim1, ~, ~, ~] = cons_to_prim(ULim1, g);

FLim1 = Fvec(rhoLim1,uLim1,PLim1,g);
[ustarL,ustarR,astarL,astarR] = adaptive_riemann_solver(rhoLim1,rhoRi,uLim1,uRi,aLim1,aRi,PLim1,PRi,g);

SL = ustarL - astarL;
SR = ustarR + astarR;

% % compare wave speeds to x/t
% speeds = xs(4:end-3)./t;
% ULidx = find(speeds <= SL);
% URidx = find(speeds >= SR);
% allidxssofar = [ULidx,URidx];
% allidxs = 1:length(speeds);
% idxs = ~ismember(allidxs,allidxssofar);
% Uhllidx = allidxs(idxs);
% 
% 
% % choose appropriate U
% Uip12 = zeros(size(U));
% Uip12(ULidx) = ULi(ULidx);
% Uip12(URidx) = URi(URidx);
% Uip12(Uhllidx) = (SR(Uhllidx).*URi(Uhllidx) - SL(Uhllidx).*ULi(Uhllidx) + FLi(Uhllidx) - FRi(Uhllidx))./(SR(Uhllidx) - SL(Uhllidx));
% 
% Uim12 = zeros(size(U));
% Uim12(ULidx) = ULim1(ULidx);
% Uim12(URidx) = URim1(URidx);
% Uim12(Uhllidx) = (SR(Uhllidx).*URim1(Uhllidx) - SL(Uhllidx).*ULim1(Uhllidx) + FLim1(Uhllidx) - FRim1(Uhllidx))./(SR(Uhllidx) - SL(Uhllidx));


% compare wave speeds to 0
speeds = 0;
ULidx = find(SL >= speeds);
URidx = find(SR <= speeds);
allidxssofar = [ULidx,URidx];
allidxs = 1:length(U(1,4:end-3));
idxs = ~ismember(allidxs,allidxssofar);
Uhllidx = allidxs(idxs);

% choose appropriate F
Fip12 = zeros(size(U));
Fip12(:,ULidx) = FLi(:,ULidx);
Fip12(:,URidx) = FRi(:,URidx);
Fip12(:,Uhllidx) = (SR(Uhllidx).*FLi(:,Uhllidx) - SL(Uhllidx).*FRi(:,Uhllidx) + SL(Uhllidx).*SR(Uhllidx).*(URi(:,Uhllidx) - ULi(:,Uhllidx)))./(SR(Uhllidx) - SL(Uhllidx));

Fim12 = zeros(size(U));
Fim12(:,ULidx) = FLim1(:,ULidx);
Fim12(:,URidx) = FRim1(:,URidx);
Fim12(:,Uhllidx) = (SR(Uhllidx).*FLim1(:,Uhllidx) - SL(Uhllidx).*FRim1(:,Uhllidx) + SL(Uhllidx).*SR(Uhllidx).*(URim1(:,Uhllidx) - ULim1(:,Uhllidx)))./(SR(Uhllidx) - SL(Uhllidx));

Fip12 = Fip12(:,4:end-3);
Fim12 = Fim12(:,4:end-3);
end

function [ustarL,ustarR,astarL,astarR] = adaptive_riemann_solver(rhoLim1,rhoRi,uLim1,uRi,aLim1,aRi,PLim1,PRi,g) % Toro 9.5.2
rhobar = 0.5*(rhoLim1 + rhoRi);
abar = 0.5*(aLim1 + aRi);
Pstari = 0.5*(PLim1 + PRi) + 0.5*(uLim1 - uRi).*rhobar.*abar;
Pstari(Pstari < 0) = 0;
Pmax = PRi;
Pmax(Pmax<PLim1) = PLim1(Pmax<PLim1);
Pmin = PRi;
Pmin(Pmin>PLim1) = PLim1(Pmin>PLim1);
Q = Pmax./Pmin;
Quser = 2;

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
astarL = zeros(size(aLim1));
astarR = zeros(size(aLim1));
ustarR = zeros(size(uLim1));
ustarL = zeros(size(uLim1));

%PVRS Toro 9.3
Pstar(PVRSidx) = 0.5.*(PLim1(PVRSidx) + PRi(PVRSidx)) + 0.5.*(uLim1(PVRSidx) - uRi(PVRSidx)).*rhobar(PVRSidx).*abar(PVRSidx);
ustar(PVRSidx) = 0.5.*(uLim1(PVRSidx) + uRi(PVRSidx)) + 0.5.*(PLim1(PVRSidx) - PRi(PVRSidx))./(rhobar(PVRSidx).*abar(PVRSidx));
rhostarL(PVRSidx) = rhoLim1(PVRSidx) + (uLim1(PVRSidx) - ustar(PVRSidx)).*(rhobar(PVRSidx)./abar(PVRSidx));
rhostarR(PVRSidx) = rhoRi(PVRSidx) + (ustar(PVRSidx) - uRi(PVRSidx)).*(rhobar(PVRSidx)./abar(PVRSidx));
astarL(PVRSidx) = sqrt((g.*Pstar(PVRSidx))./(rhostarL(PVRSidx)));
astarR(PVRSidx) = sqrt((g.*Pstar(PVRSidx))./(rhostarR(PVRSidx)));

%TRRS Toro 9.4.1
z = (g-1)./(2*g);
Pstar(TRRSidx) = ((aLim1(TRRSidx) + aRi(TRRSidx) - z*g*(uRi(TRRSidx) - uLim1(TRRSidx)))./((aLim1(TRRSidx)./(PLim1(TRRSidx).^z)) + (aRi(TRRSidx)./(PRi(TRRSidx).^z)))).^(1/z);
astarL(TRRSidx) = aLim1(TRRSidx).*(Pstar(TRRSidx)./PLim1(TRRSidx)).^z;
astarR(TRRSidx) = aRi(TRRSidx).*(Pstar(TRRSidx)./PRi(TRRSidx)).^z;
ustarL(TRRSidx) = uLim1(TRRSidx) + (1/(g*z))*(aLim1(TRRSidx) - astarL(TRRSidx));
ustarR(TRRSidx) = uRi(TRRSidx) + (1/(g*z))*(astarR(TRRSidx) - aRi(TRRSidx));
ustar(TRRSidx) = 0.5*(ustarL(TRRSidx) + ustarR(TRRSidx));

%TSRS Toro 9.4.2
z = (g-1)/(2*g);
AR = 2./((g+1)*rhoRi(TSRSidx));
AL = 2./((g+1)*rhoLim1(TSRSidx));
BR = ((g-1)/(g+1))*PRi(TSRSidx);
BL = ((g-1)/(g+1))*PLim1(TSRSidx);
gR = @(p) sqrt(AR./(p + BR));
gL = @(p) sqrt(AL./(p + BL));
P0(TSRSidx) = 0.5.*(PLim1(TSRSidx) + PRi(TSRSidx)) + 0.5.*(uLim1(TSRSidx) - uRi(TSRSidx)).*rhobar(TSRSidx).*abar(TSRSidx);
P0(P0<0) = 0;
Pstar(TSRSidx) = (gL(P0(TSRSidx)).*PLim1(TSRSidx) + gR(P0(TSRSidx)).*PRi(TSRSidx) - (uRi(TSRSidx) - uLim1(TSRSidx)))./(gL(P0(TSRSidx)) + gR(P0(TSRSidx)));
ustarR(TSRSidx) = 0.5.*(uLim1(TSRSidx) + uRi(TSRSidx)) + 0.5.*((Pstar(TSRSidx)-PRi(TSRSidx)).*gR(P0(TSRSidx)) - (Pstar(TSRSidx)-PLim1(TSRSidx)).*gL(P0(TSRSidx)));
ustarL(TSRSidx) = 0.5.*(uLim1(TSRSidx) + uRi(TSRSidx)) + 0.5.*((Pstar(TSRSidx)-PRi(TSRSidx)).*gR(P0(TSRSidx)) - (Pstar(TSRSidx)-PLim1(TSRSidx)).*gL(P0(TSRSidx)));
rhostarL(TSRSidx) = rhoLim1(TSRSidx).*(((Pstar(TSRSidx)./PLim1(TSRSidx))+((g-1)./(g+1)))./(((g-1)./(g+1)).*(Pstar(TSRSidx)./PLim1(TSRSidx)) + 1));
rhostarR(TSRSidx) = rhoRi(TSRSidx).*(((Pstar(TSRSidx)./PRi(TSRSidx))+((g-1)./(g+1)))./(((g-1)./(g+1)).*(Pstar(TSRSidx)./PRi(TSRSidx)) + 1));
astarL(TSRSidx) = sqrt((g.*Pstar(TSRSidx))./(rhostarL(TSRSidx)));
astarR(TSRSidx) = sqrt((g.*Pstar(TSRSidx))./(rhostarR(TSRSidx)));

end


function [lammt1,lampt3] = entropy_fix(uLim1,uRi,aLim1,aRi,ustar,astarL,astarR,lammt1,lampt3)

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
