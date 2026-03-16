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
Rim12 = sqrt(rhoRim1./rhoLim1);
rhotim12 = Rim12.*rhoLim1;
utim12 = (uLim1 + Rim12.*uRim1)./(1 + Rim12);
Htim12 = (HLim1 + Rim12.*HRim1)./(1 + Rim12);
atim12_temp = (g-1).*(Htim12 - 0.5.*utim12.*utim12);
atim12_temp(atim12_temp<=0) = 1e-6;
atim12 = atim12_temp.^(1/2);

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
Kpt3 = [eigen_ones; utip12+atip12; Htip12+utip12.*atip12];

% step 4 (calculating flux differences + eigenvalues + eigenvectors)

dPm = PRim1 - PLim1;
drhom = rhoRim1 - rhoLim1;
dum = uRim1 - uLim1;

alphamt2 = drhom - dPm./(atim12.^2);
alphamt1 = 1./(2*atim12.^2).*(dPm - rhotim12.*atim12.*dum);
alphamt3 = 1./(2*atim12.^2).*(dPm + rhotim12.*atim12.*dum);

alphapt2 = (rhoRi - rhoLi) - (PRi - PLi)./(atip12.^2);
alphapt1 = 1./(2*atip12.^2).*((PRi - PLi) - rhotip12.*atip12.*(uLi - uRi));
alphapt3 = 1./(2*atip12.^2).*((PRi - PLi) + rhotip12.*atip12.*(uLi - uRi));

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
