function [Fip12,Fim12] = hllcSolver(U,URi,ULi,URim1,ULim1,xs,gamma,time)
% hllcSolver - HLLC approximate Riemann solver (vectorized over interfaces)
%
% Inputs:
%   U      - full conservative state array (unused except for size)
%   URi    - right reconstructed state at interface
%   ULi    - left reconstructed state at interface
%   URim1  - right reconstructed state at previous interface
%   ULim1  - left reconstructed state at previous interface
%   xs      - grid locations (unused)
%   gamma   - ratio of specific heats
%   time    - current time (unused)
%
% Outputs:
%   Fip12, Fim12 - numerical fluxes on the right and left side of each interface

% Primitive variables at each interface state
[rhoLi, uLi, PLi, aLi, HLi, ~, ~, ~] = consToPrim(ULi, gamma);
[rhoRi, uRi, PRi, aRi, HRi, ~, ~, ~] = consToPrim(URi, gamma);

% Estimate wave speeds via Roe average + entropy fix (vectorized)
[ustarL,ustarR,astarL,astarR] = adaptiveRiemannSolver(rhoLi,rhoRi,uLi,uRi,aLi,aRi,PLi,PRi,gamma);
SL = ustarL - astarL;
SR = ustarR + astarR;

% Fluxes for left and right states
FLi = Fvec(rhoLi,uLi,PLi,gamma);
FRi = Fvec(rhoRi,uRi,PRi,gamma);

% Allocate output arrays
Fip12 = zeros(size(ULi));
Fim12 = zeros(size(ULi));

% Regions
maskL = SL >= 0;
maskR = SR <= 0;
maskM = ~(maskL | maskR);

% Left region: use left flux
Fip12(:,maskL) = FLi(:,maskL);
Fim12(:,maskL) = FLi(:,maskL);

% Right region: use right flux
Fip12(:,maskR) = FRi(:,maskR);
Fim12(:,maskR) = FRi(:,maskR);

% Middle region: HLLC star-region flux (vectorized)
if any(maskM)
    SLm = SL(maskM);
    SRm = SR(maskM);
    FLm = FLi(:,maskM);
    FRm = FRi(:,maskM);
    ULm = ULi(:,maskM);
    URm = URi(:,maskM);

    Fstar = (SRm .* FLm - SLm .* FRm + SLm .* SRm .* (URm - ULm)) ./ (SRm - SLm);
    Fip12(:,maskM) = Fstar;
    Fim12(:,maskM) = Fstar;
end
end
