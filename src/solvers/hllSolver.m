function [Fip12, Fim12] = hllSolver(U,URi,ULi,URim1,ULim1,xs,g,t)

% HLL approximate Riemann solver (Toro)

% Pre-compute primitive variables for the full conservative array
[rho, u, P, a, H, e, M, s] = consToPrim(U, g);

% Extract stencil points for reconstruction
Uim3 = U(:,1:end-6);
Uim2 = U(:,2:end-5);
Uim1 = U(:,3:end-4);
Ui   = U(:,4:end-3);
Uip1 = U(:,5:end-2);
Uip2 = U(:,6:end-1);
Uip3 = U(:,7:end);

% Compute fluxes at stencil points
Fim3 = Fvec(rho(:,1:end-6), u(:,1:end-6), P(:,1:end-6), g);
Fim2 = Fvec(rho(:,2:end-5), u(:,2:end-5), P(:,2:end-5), g);
Fim1 = Fvec(rho(:,3:end-4), u(:,3:end-4), P(:,3:end-4), g);
Fi   = Fvec(rho(:,4:end-3), u(:,4:end-3), P(:,4:end-3), g);
Fip1 = Fvec(rho(:,5:end-2), u(:,5:end-2), P(:,5:end-2), g);
Fip2 = Fvec(rho(:,6:end-1), u(:,6:end-1), P(:,6:end-1), g);
Fip3 = Fvec(rho(:,7:end),   u(:,7:end),   P(:,7:end),   g);

% Reconstructed values at current interface
[rhoRi, uRi, PRi, aRi, HRi, ~, ~, ~] = consToPrim(URi, g);
FRi = Fvec(rhoRi,uRi,PRi,g);
[rhoLi, uLi, PLi, ~, HLi, ~, ~, ~] = consToPrim(ULi, g);
FLi = Fvec(rhoLi,uLi,PLi,g);

% Reconstructed values at upstream interface
[rhoRim1, uRim1, PRim1, ~, HRim1, ~, ~, ~] = consToPrim(URim1, g);
FRim1 = Fvec(rhoRim1,uRim1,PRim1,g);
[rhoLim1, uLim1, PLim1, aLim1, HLim1, ~, ~, ~] = consToPrim(ULim1, g);
FLim1 = Fvec(rhoLim1,uLim1,PLim1,g);

% Wave speed estimates via Roe average with entropy fix
[ustarL,ustarR,astarL,astarR] = adaptiveRiemannSolver(rhoLim1,rhoRi,uLim1,uRi,aLim1,aRi,PLim1,PRi,g);
SL = ustarL - astarL;
SR = ustarR + astarR;

% Use HLL formula to compute numerical fluxes
maskL = SL >= 0;
maskR = SR <= 0;
maskM = ~(maskL | maskR);

Fip12 = zeros(size(U));
Fip12(:,maskL) = FLi(:,maskL);
Fip12(:,maskR) = FRi(:,maskR);
Fip12(:,maskM) = (SR(maskM).*FLi(:,maskM) - SL(maskM).*FRi(:,maskM) + SL(maskM).*SR(maskM).*(URi(:,maskM) - ULi(:,maskM)))./(SR(maskM) - SL(maskM));

Fim12 = zeros(size(U));
Fim12(:,maskL) = FLim1(:,maskL);
Fim12(:,maskR) = FRim1(:,maskR);
Fim12(:,maskM) = (SR(maskM).*FLim1(:,maskM) - SL(maskM).*FRim1(:,maskM) + SL(maskM).*SR(maskM).*(URim1(:,maskM) - ULim1(:,maskM)))./(SR(maskM) - SL(maskM));

Fip12 = Fip12(:,4:end-3);
Fim12 = Fim12(:,4:end-3);
end
