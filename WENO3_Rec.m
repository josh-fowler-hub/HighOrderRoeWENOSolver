function [URip12, ULip12] = WENO3_Rec(U)
% WENO5 interpolation scheme for linear advection equation du/dt+a.*du/dx=0
% Inputs:
%   U - advected quanity
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

%{
     Rim3  Rim2   Rim1  Ri   Rip1  Rip2   Rip3
|     |     |     |     |     |     |     |
| i-3 | i-2 | i-1 |  i  | i+1 | i+2 | i+3 |
|     |     |     |     |     |     |     |
Lim3  Lim2  Lim1  Li   Lip1  Lip2  Lip3
%}

% U1im3 = U(1,1:end-6);
% U2im3 = U(2,1:end-6);
% U3im3 = U(3,1:end-6);
% U1im2 = U(1,2:end-5);
% U2im2 = U(2,2:end-5);
% U3im2 = U(3,2:end-5);
% U1im1 = U(1,3:end-4);
% U2im1 = U(2,3:end-4);
% U3im1 = U(3,3:end-4);
% U1i = U(1,4:end-3);
% U2i = U(2,4:end-3);
% U3i = U(3,4:end-3);
% U1ip1 = U(1,5:end-2);
% U2ip1 = U(2,5:end-2);
% U3ip1 = U(3,5:end-2);
% U1ip2 = U(1,6:end-1);
% U2ip2 = U(2,6:end-1);
% U3ip2 = U(3,6:end-1);
% U1ip3 = U(1,7:end);
% U2ip3 = U(2,7:end);
% U3ip3 = U(3,7:end);

U1im2 = U(1,1:end-5);
U2im2 = U(2,1:end-5);
U3im2 = U(3,1:end-5);
U1im1 = U(1,2:end-4);
U2im1 = U(2,2:end-4);
U3im1 = U(3,2:end-4);
U1i = U(1,3:end-3);
U2i = U(2,3:end-3);
U3i = U(3,3:end-3);
U1ip1 = U(1,4:end-2);
U2ip1 = U(2,4:end-2);
U3ip1 = U(3,4:end-2);
U1ip2 = U(1,5:end-1);
U2ip2 = U(2,5:end-1);
U3ip2 = U(3,5:end-1);
U1ip3 = U(1,6:end);
U2ip3 = U(2,6:end);
U3ip3 = U(3,6:end);

C1 = 1/3;
C2 = 2/3;
ep = 1e-6;
p = 2;

% UR1i+(1/2)
IS1 = (U1ip2-U1ip1).^2;
IS2 = (U1ip1-U1i).^2;

alpha1 = C1./(ep+IS1).^p;
alpha2 = C2./(ep+IS2).^p;
alphatot = alpha1+alpha2;

omega1 = alpha1./alphatot;
omega2 = alpha2./alphatot;

q1 = (3/2)*U1ip1-(1/2)*U1ip2;
q2 = (1/2)*U1i+(1/2)*U1ip1;

URip12(1,:) = omega1.*q1+omega2.*q2;

% UL1i+(1/2)
IS1 = (U1i-U1im1).^2;
IS2 = (U1ip1-U1i).^2;

alpha1 = C1./(ep+IS1).^p;
alpha2 = C2./(ep+IS2).^p;
alphatot = alpha1+alpha2;

omega1 = alpha1./alphatot;
omega2 = alpha2./alphatot;

q1 = (3/2)*U1i-(1/2)*U1im1;
q2 = (1/2)*U1i+(1/2)*U1ip1;
ULip12(1,:) = omega1.*q1+omega2.*q2;

% UR2i+(1/2)
IS1 = (U2ip2-U2ip1).^2;
IS2 = (U2ip1-U2i).^2;

alpha1 = C1./(ep+IS1).^p;
alpha2 = C2./(ep+IS2).^p;
alphatot = alpha1+alpha2;

omega1 = alpha1./alphatot;
omega2 = alpha2./alphatot;

q1 = (3/2)*U2ip1-(1/2)*U2ip2;
q2 = (1/2)*U2i+(1/2)*U2ip1;

URip12(2,:) = omega1.*q1+omega2.*q2;

% UL2i+(1/2)
IS1 = (U2i-U2im1).^2;
IS2 = (U2ip1-U2i).^2;

alpha1 = C1./(ep+IS1).^p;
alpha2 = C2./(ep+IS2).^p;
alphatot = alpha1+alpha2;

omega1 = alpha1./alphatot;
omega2 = alpha2./alphatot;

q1 = (3/2)*U2i-(1/2)*U2im1;
q2 = (1/2)*U2i+(1/2)*U2ip1;
ULip12(2,:) = omega1.*q1+omega2.*q2;

% UR3i+(1/2)
IS1 = (U3ip2-U3ip1).^2;
IS2 = (U3ip1-U3i).^2;

alpha1 = C1./(ep+IS1).^p;
alpha2 = C2./(ep+IS2).^p;
alphatot = alpha1+alpha2;

omega1 = alpha1./alphatot;
omega2 = alpha2./alphatot;

q1 = (3/2)*U3ip1-(1/2)*U3ip2;
q2 = (1/2)*U3i+(1/2)*U3ip1;

URip12(3,:) = omega1.*q1+omega2.*q2;

% UL3i+(1/2)
IS1 = (U3i-U3im1).^2;
IS2 = (U3ip1-U3i).^2;

alpha1 = C1./(ep+IS1).^p;
alpha2 = C2./(ep+IS2).^p;
alphatot = alpha1+alpha2;

omega1 = alpha1./alphatot;
omega2 = alpha2./alphatot;

q1 = (3/2)*U3i-(1/2)*U3im1;
q2 = (1/2)*U3i+(1/2)*U3ip1;
ULip12(3,:) = omega1.*q1+omega2.*q2;

% if any(ULip12(1,:) <= 0) || any(URip12(1,:)<= 0)
%     ULip12 = U(:,3:end-3);
%     URip12 = U(:,4:end-2);
% end