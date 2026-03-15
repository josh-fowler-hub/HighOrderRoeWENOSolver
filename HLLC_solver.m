function [Fip12,Fim12] = HLLC_solver(UL,UR)
% Compute HLLC flux

    % Left state
    rL = UL(1,:);
    uL = UL(2,:)./rL;
    EL = UL(3,:)./rL;
    pL = (gamma-1).*( UL(3,:) - rL.*uL.*uL/2 );
    aL = sqrt(gamma.*pL/rL);
    
    % Right state
    rR = UR(1,:);
    uR = UR(2,:)./rR;
    ER = UR(3,:)./rR;
    pR = (gamma-1).*( UR(3,:) - rR.*uR.*uR/2 );
    aR = sqrt(gamma.*pR/rR);
    
    % Left and Right fluxes
    FL=[rL.*uL; rL.*uL.^2+pL; uL.*(rL.*EL+pL)];
    FR=[rR.*uR; rR.*uR.^2+pR; uR.*(rR.*ER+pR)];

    % Compute guess pressure from PVRS Riemann solver
    PPV  = max(0 , 0.5.*(pL+pR) + 0.5.*(uL-uR) .* (0.25.*(rL+rR).*(aL+aR)));
    pmin = min(pL,pR);
    pmax = max(pL,pR);
    Qmax = pmax/pmin;
    Quser= 2.0; % <--- parameter manually set (I don't like this!)
    
     if (Qmax <= Quser) && (pmin <= PPV) && (PPV <= pmax)
     % Select PRVS Riemann solver
         pM = PPV;
         %uM = 0.5*(uL + uR) + 0.5*(pL - pR)/CUP;
      else
         if PPV < pmin
         % Select Two-Rarefaction Riemann solver
            PQ  = (pL/pR)^(gamma - 1.0)/(2.0.*gamma);
            uM  = (PQ.*uL/aL + uR/aR + 2/(gamma-1).*(PQ-1.0))/(PQ/aL+1.0/aR);
            PTL = 1 + (gamma-1)/2.0.*(uL - uM)/aL;
            PTR = 1 + (gamma-1)/2.0.*(uM - uR)/aR;
            pM  = 0.5.*(pL.*PTL^(2.*gamma/(gamma-1)) + pR.*PTR^(2.*gamma/(gamma-1)));
         else 
         % Use Two-Shock Riemann solver with PVRS as estimate
            GEL = sqrt((2/(gamma+1)/rL)/((gamma-1)/(gamma+1).*pL + PPV));
            GER = sqrt((2/(gamma+1)/rR)/((gamma-1)/(gamma+1).*pR + PPV));
            pM  = (GEL.*pL + GER.*pR - (uR - uL))/(GEL + GER);
            %uM  = 0.5*(uL + uR) + 0.5*(GER*(pM - pR) - GEL*(pM - pL));
         end
      end

    % Estimate wave speeds: SL, SR and SM (Toro, 1994)
    if pM>pL; zL=sqrt(1+(gamma+1)/(2.*gamma).*(pM/pL-1)); else, zL=1; end    
    if pM>pR; zR=sqrt(1+(gamma+1)/(2.*gamma).*(pM/pR-1)); else, zR=1; end
  
	SL = uL - aL.*zL;
    SR = uR + aR.*zR;
    SM = (pL-pR + rR.*uR.*(SR-uR) - rL.*uL.*(SL-uL))/(rR.*(SR-uR) - rL.*(SL-uL));
    
    % Compute the HLL flux.
    if 0 <= SL  % Right-going supersonic flow
        HLLC = FL;
    elseif (SL <= 0) && (0 <= SM)	% Subsonic flow to the right
        qsL = rL.*(SL-uL)/(SL-SM).*[1; SM; UL(3,:)/rL + (SM-uL).*(SM+pL/(rL.*(SL-uL)))];
        HLLC = FL + SL.*(qsL - UL);
    elseif (SM <= 0) && (0 <= SR)	% Subsonic flow to the Left
        qsR = rR.*(SR-uR)/(SR-SM).*[1; SM; UR(3,:)/rR + (SM-uR).*(SM+pR/(rR.*(SR-uR)))];
        HLLC = FR + SR.*(qsR - UR);
    elseif  0 >= SR % Left-going supersonic flow
        HLLC = FR;
    end
end

function [ustar,astarL,astarR] = adaptive_riemann_solver(rhoLim1,rhoRi,uLim1,uRi,aLim1,aRi,PLim1,PRi,g) % Toro 9.5.2
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
ustar(TSRSidx) = 0.5.*(uLim1(TSRSidx) + uRi(TSRSidx)) + 0.5.*((Pstar(TSRSidx)-PRi(TSRSidx)).*gR(P0(TSRSidx)) - (Pstar(TSRSidx)-PLim1(TSRSidx)).*gL(P0(TSRSidx)));
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
