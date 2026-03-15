function [URip12,ULip12,URim12,ULim12] = musclRec(U,flux_limiter,lam,k)
% Computing cell state at the left and right cell interface using MUSCL reconstruction
% Uses vectorized operations for all 3 conservative variables simultaneously
ep1 = 1e-6; % small number for smoothness factor
ep2 = 1e-6; % small number for smoothness factor

%{
i-2 = 1:end-4
i-1 = 2:end-3
i = 3:end-2
i+1 = 4:end-1
i+2 = 5:end

URi = U(:, 3:end);      | Uip1, For the Roe Scheme
ULi = U(:, 2:end-1);    | Ui
URim1 = U(:, 2:end-1);  | Ui
ULim1 = U(:, 1:end-2);  | Uim1
%}

% Extract stencil points for all 3 variables at once (vectorized indexing)
Uim3 = U(:,1:end-6);
Uim2 = U(:,2:end-5);
Uim1 = U(:,3:end-4);
Ui = U(:,4:end-3);
Uip1 = U(:,5:end-2);
Uip2 = U(:,6:end-1);
Uip3 = U(:,7:end);

% Compute gradients for all 3 variables at once
dUip12 = Uip1 - Ui;         % U_{i+1} - U_i
dUim12 = Ui - Uim1;         % U_i - U_{i-1}
dUip32 = Uip2 - Uip1;       % U_{i+2} - U_{i+1}
dUim32 = Uim1 - Uim2;       % U_{i-1} - U_{i-2}

% Compute smoothness ratios for all 3 variables at once
rip1 = (dUip12 + ep1) ./ (dUip32 + ep2);
ri = (dUim12 + ep1) ./ (dUip12 + ep2);
rim1 = (dUim32 + ep1) ./ (dUim12 + ep2);

% Compute limiter functions for all 3 variables at once
phiip1 = phiFlux(flux_limiter, rip1);
phii = phiFlux(flux_limiter, ri);
phiim1 = phiFlux(flux_limiter, rim1);

% Compute coefficients once
coeff = 1 / (4 * lam);
term1 = (1 - k);
term2 = (1 + k);

% Reconstruct all 4 interfaces for all 3 variables at once
ULip12 = Ui + coeff * phii .* (term1 * dUim12 + term2 * dUip12);
URip12 = Uip1 - coeff * phiip1 .* (term2 * dUip32 + term1 * dUip12);
ULim12 = Uim1 + coeff * phiim1 .* (term1 * dUim32 + term2 * dUim12);
URim12 = Ui - coeff * phii .* (term2 * dUip12 + term1 * dUim12);

end