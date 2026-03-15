function [URip12, ULip12] = weno3Rec(U)
% WENO3 interpolation scheme for linear advection equation du/dt+a.*du/dx=0
% Inputs:
%   U - advected quantity
%   i - cell face
% Outputs:
%   URip12 - reconstructed right cell face
%   ULip12 - reconstructed left cell face

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

C1 = 1/3;
C2 = 2/3;
ep = 1e-6;
p = 2;

% ===== RIGHT RECONSTRUCTION (at interface i+1/2) =====
% Compute smoothness indicators once using component 1
IS1_R = (Uip2(1,:)-Uip1(1,:)).^2;
IS2_R = (Uip1(1,:)-Ui(1,:)).^2;

% Compute weights once for all components
alpha1_R = C1./(ep+IS1_R).^p;
alpha2_R = C2./(ep+IS2_R).^p;
alphatot_R = alpha1_R + alpha2_R;

omega1_R = alpha1_R./alphatot_R;
omega2_R = alpha2_R./alphatot_R;

% Apply weights to all 3 variables
for k = 1:3
    q1 = (3/2)*Uip1(k,:) - (1/2)*Uip2(k,:);
    q2 = (1/2)*Ui(k,:) + (1/2)*Uip1(k,:);
    URip12(k,:) = omega1_R.*q1 + omega2_R.*q2;
end

% ===== LEFT RECONSTRUCTION (at interface i+1/2) =====
% Compute smoothness indicators once using component 1
IS1_L = (Ui(1,:)-Uim1(1,:)).^2;
IS2_L = (Uip1(1,:)-Ui(1,:)).^2;

% Compute weights once for all components
alpha1_L = C1./(ep+IS1_L).^p;
alpha2_L = C2./(ep+IS2_L).^p;
alphatot_L = alpha1_L + alpha2_L;

omega1_L = alpha1_L./alphatot_L;
omega2_L = alpha2_L./alphatot_L;

% Apply weights to all 3 variables
for k = 1:3
    q1 = (3/2)*Ui(k,:) - (1/2)*Uim1(k,:);
    q2 = (1/2)*Ui(k,:) + (1/2)*Uip1(k,:);
    ULip12(k,:) = omega1_L.*q1 + omega2_L.*q2;
end

end