function [Fip12,Fim12] = Roe_Pike_solver(U,URi,ULi,URim1,ULim1,g, CD_Term_Order)

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

    %{
          ULi    URi
    |  i-1 |  i   |  i+1 |     |
    ULim1 URim1  URip1  ULip1

    %}

    % Setting Up Flux's and Variables
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

   % 1. Compute the Roe average values according to Toro (11.118).
    Rim12 = sqrt(rhoRim1./rhoLim1);
    rhotim12 = Rim12.*rhoLim1;
    utim12 = (uLim1 + Rim12.*uRim1)./(1 + Rim12);
    Htim12 = (HLim1 + Rim12.*HRim1)./(1 + Rim12);
    atim12_temp = (g-1).*(Htim12 - 0.5.*utim12.*utim12);
    atim12_temp(atim12_temp<0) = 0;
    atim12 = atim12_temp.^(1/2);
    
    Rip12 = sqrt(rhoRi./rhoLi);
    rhotip12 = Rip12.*rhoLi;
    utip12 = (uLi + Rip12.*uRi)./(1 + Rip12);
    Htip12 = (HLi + Rip12.*HRi)./(1 + Rip12);
    atip12_temp = (g-1).*(Htip12-0.5.*utip12.^2);
    atip12_temp(atip12_temp<0) = 0;
    atip12 = atip12_temp.^(1/2);
   
   % 2. Compute the eigenvalues lambda_i using the analytical expressions
   % using Toro (11.107) evaluated on the averages (11.118).
    lammt1 = utim12-atim12;
    lammt2 = utim12;
    lammt3 = utim12+atim12;

    lampt1 = utip12-atip12;
    lampt2 = utip12;
    lampt3 = utip12+atip12;

   % 3. Compute the right eigenvectors using the analytical expressions
   % Toro (11.108) evaluated on the averages (11.118).

    eigen_ones = ones(1, length(utim12));

    Kmt1 = [eigen_ones; utim12 - atim12; Htim12 - utim12.*atim12];
    Kmt2 = [eigen_ones; utim12; 0.5.*utim12.*utim12];
    Kmt3 = [eigen_ones; utim12 + atim12; Htim12 + utim12.*atim12];

    eigen_ones = ones(1, length(utip12));

    Kpt1 = [eigen_ones; utip12-atip12; Htip12-utip12.*atip12];
    Kpt2 = [eigen_ones; utip12; 0.5.*utip12.*utip12];
    Kpt3 = [eigen_ones; utip12+atip12; Htip12+utip12.*atip12];
   
   % 4. Compute the wave strengths using the analytical expressions Toro
   % (11.113) evaluated on the averages (11.118).
    dPm = PRim1 - PLim1;
    drhom = rhoRim1 - rhoLim1;
    dum = uRim1 - uLim1;

    alphamt2 = drhom - dPm./(atim12.^2);
    alphamt1 = 1./(2*atim12.^2).*(dPm - rhotim12.*atim12.*dum);
    alphamt3 = 1./(2*atim12.^2).*(dPm + rhotim12.*atim12.*dum);


    dPp = PRi - PLi;
    drhop = rhoRi - rhoLi;
    dup = uRi - uLi;

    alphapt2 = drhop - dPp./(atip12.^2);
    alphapt1 = 1./(2*atip12.^2).*(dPp - rhotip12.*atip12.*dup);
    alphapt3 = 1./(2*atip12.^2).*(dPp + rhotip12.*atip12.*dup);
    
   % 4.5 Compute the Entropy fix
    [ustar,astarL,astarR] = adaptive_riemann_solver(rhoLim1,rhoRi,uLim1,uRi,aLim1,aRi,PLim1,PRi,g);
    [lammt1,lampt3] = entropy_fix(uLim1,uRi,aLim1,aRi,ustar,astarL,astarR,lammt1,lampt3);

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

   % 5. Use all of the above 1-4 to compute Fip12, according to any of the
   % formula Toro (11.27)-(11.29).DFim12 = 0.5.*(alphamt1.*abs(lammt1).*Kmt1 + alphamt2.*abs(lammt2).*Kmt2 +alphamt3.*abs(lammt3).*Kmt3);
    DFim12 = 0.5.*(alphamt1.*abs(lammt1).*Kmt1 + alphamt2.*abs(lammt2).*Kmt2 +alphamt3.*abs(lammt3).*Kmt3);
    DFip12 = 0.5.*(alphapt1.*abs(lampt1).*Kpt1 + alphapt2.*abs(lampt2).*Kpt2 + alphapt3.*abs(lampt3).*Kpt3);
    Fim12 = CDim12-DFim12;
    Fip12 = CDip12-DFip12;
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