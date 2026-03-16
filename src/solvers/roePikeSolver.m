function [Fip12,Fim12] = roePikeSolver(U,URi,ULi,URim1,ULim1,g, CD_Term_Order)

    Uim3 = U(:,1:end-6);
    [rhoim3, uim3, Pim3, ~, ~, ~, ~, ~] = consToPrim(Uim3, g);
    Fim3 = Fvec(rhoim3,uim3,Pim3,g);
    Uim2 = U(:,2:end-5);
    [rhoim2, uim2, Pim2, ~, ~, ~, ~, ~] = consToPrim(Uim2, g);
    Fim2 = Fvec(rhoim2,uim2,Pim2,g);
    Uim1 = U(:,3:end-4);
    [rhoim1, uim1, Pim1, ~, ~, ~, ~, ~] = consToPrim(Uim1, g);
    Fim1 = Fvec(rhoim1,uim1,Pim1,g);
    Ui = U(:,4:end-3);
    [rhoi, ui, Pi, ~, ~, ~, ~, ~] = consToPrim(Ui, g);
    Fi = Fvec(rhoi,ui,Pi,g);
    Uip1 = U(:,5:end-2);
    [rhoip1, uip1, Pip1, ~, ~, ~, ~, ~] = consToPrim(Uip1, g);
    Fip1 = Fvec(rhoip1,uip1,Pip1,g);
    Uip2 = U(:,6:end-1);
    [rhoip2, uip2, Pip2, ~, ~, ~, ~, ~] = consToPrim(Uip2, g);
    Fip2 = Fvec(rhoip2,uip2,Pip2,g);
    Uip3 = U(:,7:end);
    [rhoip3, uip3, Pip3, ~, ~, ~, ~, ~] = consToPrim(Uip3, g);
    Fip3 = Fvec(rhoip3,uip3,Pip3,g);

    %{
          ULi    URi
    |  i-1 |  i   |  i+1 |     |
    ULim1 URim1  URip1  ULip1

    %}

    % Setting Up Flux's and Variables
    % Ui
    [rhoRi, uRi, PRi, aRi, HRi, ~, ~, ~] = consToPrim(URi, g);

    FRi = Fvec(rhoRi,uRi,PRi,g);

    [rhoLi, uLi, PLi, ~, HLi, ~, ~, ~] = consToPrim(ULi, g);

    FLi = Fvec(rhoLi,uLi,PLi,g);

    % Uim1
    [rhoRim1, uRim1, PRim1, ~, HRim1, ~, ~, ~] = consToPrim(URim1, g);

    FRim1 = Fvec(rhoRim1,uRim1,PRim1,g);

    [rhoLim1, uLim1, PLim1, aLim1, HLim1, ~, ~, ~] = consToPrim(ULim1, g);

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

   % 5. Use all of the above 1-4 to compute Fip12, according to any of the
   % formula Toro (11.27)-(11.29).DFim12 = 0.5.*(alphamt1.*abs(lammt1).*Kmt1 + alphamt2.*abs(lammt2).*Kmt2 +alphamt3.*abs(lammt3).*Kmt3);
    DFim12 = 0.5.*(alphamt1.*abs(lammt1).*Kmt1 + alphamt2.*abs(lammt2).*Kmt2 +alphamt3.*abs(lammt3).*Kmt3);
    DFip12 = 0.5.*(alphapt1.*abs(lampt1).*Kpt1 + alphapt2.*abs(lampt2).*Kpt2 + alphapt3.*abs(lampt3).*Kpt3);
    Fim12 = CDim12-DFim12;
    Fip12 = CDip12-DFip12;
end
