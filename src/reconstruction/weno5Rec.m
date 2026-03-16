function [URip12, ULip12] = weno5Rec(U)

C0 = 0.1;
C1 = 0.6;
C2 = 0.3;
ep = 1e-6;
p = 2;

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

% Extract stencil points for all 3 variables at once
Uim2 = U(:,1:end-5);
Uim1 = U(:,2:end-4);
Ui = U(:,3:end-3);
Uip1 = U(:,4:end-2);
Uip2 = U(:,5:end-1);
Uip3 = U(:,6:end);

% ===== RIGHT RECONSTRUCTION (at interface i+1/2) =====
% Compute smoothness indicators once using component 1
IS0_R = 13./12.*(Uip1(1,:)-2.*Uip2(1,:)+Uip3(1,:)).^2 + 1./4.*(3.*Uip1(1,:)-4.*Uip2(1,:)+Uip3(1,:)).^2;
IS1_R = 13./12.*(Ui(1,:)-2.*Uip1(1,:)+Uip2(1,:)).^2 + 1./4.*(Ui(1,:)-Uip2(1,:)).^2;
IS2_R = 13./12.*(Uim1(1,:)-2.*Ui(1,:)+Uip1(1,:)).^2 + 1./4.*(Uim1(1,:)-4.*Ui(1,:)+3.*Uip1(1,:)).^2;

% Compute weights once for all components
alpha0_R = C0./(ep+IS0_R).^p;
alpha1_R = C1./(ep+IS1_R).^p;
alpha2_R = C2./(ep+IS2_R).^p;
alphatot_R = alpha0_R + alpha1_R + alpha2_R;

omega0_R = alpha0_R./alphatot_R;
omega1_R = alpha1_R./alphatot_R;
omega2_R = alpha2_R./alphatot_R;

% Apply weights to all 3 variables (vectorized)
q0 = (2 .* Uip3 - 7 .* Uip2 + 11 .* Uip1) ./ 6;
q1 = (-Uip2 + 5 .* Uip1 + 2 .* Ui) ./ 6;
q2 = (2 .* Uip1 + 5 .* Ui - Uim1) ./ 6;
URip12 = omega0_R .* q0 + omega1_R .* q1 + omega2_R .* q2;

% ===== LEFT RECONSTRUCTION (at interface i+1/2) =====
% Compute smoothness indicators once using component 1
IS0_L = 13./12.*(Uim2(1,:)-2.*Uim1(1,:)+Ui(1,:)).^2 + 1./4.*(Uim2(1,:)-4.*Uim1(1,:)+3.*Ui(1,:)).^2;
IS1_L = 13./12.*(Uim1(1,:)-2.*Ui(1,:)+Uip1(1,:)).^2 + 1./4.*(Uim1(1,:)-Uip1(1,:)).^2;
IS2_L = 13./12.*(Ui(1,:)-2.*Uip1(1,:)+Uip2(1,:)).^2 + 1./4.*(3.*Ui(1,:)-4.*Uip1(1,:)+Uip2(1,:)).^2;

% Compute weights once for all components
alpha0_L = C0./(ep+IS0_L).^p;
alpha1_L = C1./(ep+IS1_L).^p;
alpha2_L = C2./(ep+IS2_L).^p;
alphatot_L = alpha0_L + alpha1_L + alpha2_L;

omega0_L = alpha0_L./alphatot_L;
omega1_L = alpha1_L./alphatot_L;
omega2_L = alpha2_L./alphatot_L;

% Apply weights to all 3 variables (vectorized)
q0 = (2 .* Uim2 - 7 .* Uim1 + 11 .* Ui) ./ 6;
q1 = (-Uim1 + 5 .* Ui + 2 .* Uip1) ./ 6;
q2 = (2 .* Ui + 5 .* Uip1 - Uip2) ./ 6;
ULip12 = omega0_L .* q0 + omega1_L .* q1 + omega2_L .* q2;
