function [URip12, ULip12] = WENO5_Rec(U)

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

% UR1i+1./2
IS0 = 13./12.*(U1ip1-2.*U1ip2+U1ip3).^2+1./4.*(3.*U1ip1-4.*U1ip2+U1ip3).^2;
IS1 = 13./12.*(U1i-2.*U1ip1+U1ip2).^2+1./4.*(U1i-U1ip2).^2;
IS2 = 13./12.*(U1im1-2.*U1i+U1ip1).^2+1./4.*(U1im1-4.*U1i+3.*U1ip1).^2;

alpha0 = C0./(ep+IS0).^p;
alpha1 = C1./(ep+IS1).^p;
alpha2 = C2./(ep+IS2).^p;
alphatot = alpha0+alpha1+alpha2;

omega0 = alpha0./alphatot;
omega1 = alpha1./alphatot;
omega2 = alpha2./alphatot;

q0 = (2.*U1ip3-7.*U1ip2+11.*U1ip1)./6;
q1 = (-U1ip2+5.*U1ip1+2.*U1i)./6;
q2 = (2.*U1ip1+5.*U1i-U1im1)./6;

URip12(1,:) = omega0.*q0+omega1.*q1+omega2.*q2;

% UL1i+1./2
IS0 = 13./12.*(U1im2-2.*U1im1+U1i).^2+1./4.*(U1im2-4.*U1im1+3.*U1i).^2;
IS1 = 13./12.*(U1im1-2.*U1i+U1ip1).^2+1./4.*(U1im1-U1ip1).^2;
IS2 = 13./12.*(U1i-2.*U1ip1+U1ip2).^2+1./4.*(3.*U1i-4.*U1ip1+U1ip2).^2;

alpha0 = C0./(ep+IS0).^p;
alpha1 = C1./(ep+IS1).^p;
alpha2 = C2./(ep+IS2).^p;
alphatot = alpha0+alpha1+alpha2;

omega0 = alpha0./alphatot;
omega1 = alpha1./alphatot;
omega2 = alpha2./alphatot;

q0 = (2.*U1im2-7.*U1im1+11.*U1i)./6;
q1 = (-U1im1+5.*U1i+2.*U1ip1)./6;
q2 = (2.*U1i+5.*U1ip1-U1ip2)./6;

ULip12(1,:) = omega0.*q0+omega1.*q1+omega2.*q2;

% UR2i+1./2
IS0 = 13./12.*(U2ip1-2.*U2ip2+U2ip3).^2+1./4.*(3.*U2ip1-4.*U2ip2+U2ip3).^2;
IS1 = 13./12.*(U2i-2.*U2ip1+U2ip2).^2+1./4.*(U2i-U2ip2).^2;
IS2 = 13./12.*(U2im1-2.*U2i+U2ip1).^2+1./4.*(U2im1-4.*U2i+3.*U2ip1).^2;

alpha0 = C0./(ep+IS0).^p;
alpha1 = C1./(ep+IS1).^p;
alpha2 = C2./(ep+IS2).^p;
alphatot = alpha0+alpha1+alpha2;

omega0 = alpha0./alphatot;
omega1 = alpha1./alphatot;
omega2 = alpha2./alphatot;

q0 = (2.*U2ip3-7.*U2ip2+11.*U2ip1)./6;
q1 = (-U2ip2+5.*U2ip1+2.*U2i)./6;
q2 = (2.*U2ip1+5.*U2i-U2im1)./6;

URip12(2,:) = omega0.*q0+omega1.*q1+omega2.*q2;

% UL2i+1./2
IS0 = 13./12.*(U2im2-2.*U2im1+U2i).^2+1./4.*(U2im2-4.*U2im1+3.*U2i).^2;
IS1 = 13./12.*(U2im1-2.*U2i+U2ip1).^2+1./4.*(U2im1-U2ip1).^2;
IS2 = 13./12.*(U2i-2.*U2ip1+U2ip2).^2+1./4.*(3.*U2i-4.*U2ip1+U2ip2).^2;

alpha0 = C0./(ep+IS0).^p;
alpha1 = C1./(ep+IS1).^p;
alpha2 = C2./(ep+IS2).^p;
alphatot = alpha0+alpha1+alpha2;

omega0 = alpha0./alphatot;
omega1 = alpha1./alphatot;
omega2 = alpha2./alphatot;

q0 = (2.*U2im2-7.*U2im1+11.*U2i)./6;
q1 = (-U2im1+5.*U2i+2.*U2ip1)./6;
q2 = (2.*U2i+5.*U2ip1-U2ip2)./6;

ULip12(2,:) = omega0.*q0+omega1.*q1+omega2.*q2;

% UR3i+1./2
IS0 = 13./12.*(U3ip1-2.*U3ip2+U3ip3).^2+1./4.*(3.*U3ip1-4.*U3ip2+U3ip3).^2;
IS1 = 13./12.*(U3i-2.*U3ip1+U3ip2).^2+1./4.*(U3i-U3ip2).^2;
IS2 = 13./12.*(U3im1-2.*U3i+U3ip1).^2+1./4.*(U3im1-4.*U3i+3.*U3ip1).^2;

alpha0 = C0./(ep+IS0).^p;
alpha1 = C1./(ep+IS1).^p;
alpha2 = C2./(ep+IS2).^p;
alphatot = alpha0+alpha1+alpha2;

omega0 = alpha0./alphatot;
omega1 = alpha1./alphatot;
omega2 = alpha2./alphatot;

q0 = (2.*U3ip3-7.*U3ip2+11.*U3ip1)./6;
q1 = (-U3ip2+5.*U3ip1+2.*U3i)./6;
q2 = (2.*U3ip1+5.*U3i-U3im1)./6;

URip12(3,:) = omega0.*q0+omega1.*q1+omega2.*q2;

% UL3i+1./2
IS0 = 13./12.*(U3im2-2.*U3im1+U3i).^2+1./4.*(U3im2-4.*U3im1+3.*U3i).^2;
IS1 = 13./12.*(U3im1-2.*U3i+U3ip1).^2+1./4.*(U3im1-U3ip1).^2;
IS2 = 13./12.*(U3i-2.*U3ip1+U3ip2).^2+1./4.*(3.*U3i-4.*U3ip1+U3ip2).^2;

alpha0 = C0./(ep+IS0).^p;
alpha1 = C1./(ep+IS1).^p;
alpha2 = C2./(ep+IS2).^p;
alphatot = alpha0+alpha1+alpha2;

omega0 = alpha0./alphatot;
omega1 = alpha1./alphatot;
omega2 = alpha2./alphatot;

q0 = (2.*U3im2-7.*U3im1+11.*U3i)./6;
q1 = (-U3im1+5.*U3i+2.*U3ip1)./6;
q2 = (2.*U3i+5.*U3ip1-U3ip2)./6;

ULip12(3,:) = omega0.*q0+omega1.*q1+omega2.*q2;