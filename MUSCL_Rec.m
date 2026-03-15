function [URip12,ULip12,URim12,ULim12] = MUSCL_Rec(U,flux_limiter,lam,k)
% calculaitng cell state at the left and right cell interface
% the primitive variables are a 3x1 matrix, 3 reconstructions are preformed
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

U1im3 = U(1,1:end-6);
U2im3 = U(2,1:end-6);
U3im3 = U(3,1:end-6);
U1im2 = U(1,2:end-5);
U2im2 = U(2,2:end-5);
U3im2 = U(3,2:end-5);
U1im1 = U(1,3:end-4);
U2im1 = U(2,3:end-4);
U3im1 = U(3,3:end-4);
U1i = U(1,4:end-3);
U2i = U(2,4:end-3);
U3i = U(3,4:end-3);
U1ip1 = U(1,5:end-2);
U2ip1 = U(2,5:end-2);
U3ip1 = U(3,5:end-2);
U1ip2 = U(1,6:end-1);
U2ip2 = U(2,6:end-1);
U3ip2 = U(3,6:end-1);
U1ip3 = U(1,7:end);
U2ip3 = U(2,7:end);
U3ip3 = U(3,7:end);

% 1 Gradients
dU1ip12 = U1ip1 - U1i;
dU1im12 = U1i - U1im1;
dU1ip32 = U1ip2 - U1ip1;
dU1im32 = U1im1 - U1im2;

% 2 Smoothness factor
rip1 = (dU1ip12+ep1)./(dU1ip32+ep2);
ri = (dU1im12+ep1)./(dU1ip12+ep2);
rim1 = (dU1im32+ep1)./(dU1im12+ep2);

% 3 Limiter Functions
phiip1 = phi_flux(flux_limiter, rip1);
phii = phi_flux(flux_limiter, ri);
phiim1 = phi_flux(flux_limiter,rim1);

% first_min = phiip1;
% second_min = phii;
% third_min = phiim1;
% first_min(first_min>second_min) = second_min(second_min<first_min);
% phiip12 = first_min;
% second_min(second_min>third_min) = third_min(third_min<second_min);
% phiim12 = second_min;

% 4 1 Reconstruction
ULip12(1,:) = U1i + (phii/(4*lam)).*((1-k)*dU1im12 + (1+k)*dU1ip12);
URip12(1,:) = U1ip1 - (phiip1/(4*lam)).*((1+k)*dU1ip32 + (1-k)*dU1ip12);
ULim12(1,:) = U1im1 + (phiim1/(4*lam)).*((1-k)*dU1im32 + (1+k)*dU1im12);
URim12(1,:) = U1i - (phii/(4*lam)).*((1+k)*dU1ip12 + (1-k)*dU1im12);


% 1 Gradients
dU2ip12 = U2ip1 - U2i;
dU2im12 = U2i - U2im1;
dU2ip32 = U2ip2 - U2ip1;
dU2im32 = U2im1 - U2im2;

% 2 Smoothness factor
rip1 = (dU2ip12+ep1)./(dU2ip32+ep2);
ri = (dU2im12+ep1)./(dU2ip12+ep2);
rim1 = (dU2im32+ep1)./(dU2im12+ep2);

% 3 Limiter Functions
phiip1 = phi_flux(flux_limiter, rip1);
phii = phi_flux(flux_limiter, ri);
phiim1 = phi_flux(flux_limiter,rim1);

% first_min = phiip1;
% second_min = phii;
% third_min = phiim1;
% first_min(first_min>second_min) = second_min(second_min<first_min);
% phiip12 = first_min;
% second_min(second_min>third_min) = third_min(third_min<second_min);
% phiim12 = second_min;

% 4 2 Reconstruction
ULip12(2,:) = U2i + (phii/(4*lam)).*((1-k)*dU2im12 + (1+k)*dU2ip12);
URip12(2,:) = U2ip1 - (phiip1/(4*lam)).*((1+k)*dU2ip32 + (1-k)*dU2ip12);
ULim12(2,:) = U2im1 + (phiim1/(4*lam)).*((1-k)*dU2im32 + (1+k)*dU2im12);
URim12(2,:) = U2i - (phii/(4*lam)).*((1+k)*dU2ip12 + (1-k)*dU2im12);

% 1 Gradients
dU3ip12 = U3ip1 - U3i;
dU3im12 = U3i - U3im1;
dU3ip32 = U3ip2 - U3ip1;
dU3im32 = U3im1 - U3im2;

% 2 Smoothness factor
rip1 = (dU3ip12+ep1)./(dU3ip32+ep2);
ri = (dU3im12+ep1)./(dU3ip12+ep2);
rim1 = (dU3im32+ep1)./(dU3im12+ep2);

% 3 Limiter Functions
phiip1 = phi_flux(flux_limiter, rip1);
phii = phi_flux(flux_limiter, ri);
phiim1 = phi_flux(flux_limiter,rim1);

% first_min = phiip1;
% second_min = phii;
% third_min = phiim1;
% first_min(first_min>second_min) = second_min(second_min<first_min);
% phiip12 = first_min;
% second_min(second_min>third_min) = third_min(third_min<second_min);
% phiim12 = second_min;

% 4 3 Reconstruction
ULip12(3,:) = U3i + (phii/(4*lam)).*((1-k)*dU3im12 + (1+k)*dU3ip12);
URip12(3,:) = U3ip1 - (phiip1/(4*lam)).*((1+k)*dU3ip32 + (1-k)*dU3ip12);
ULim12(3,:) = U3im1 + (phiim1/(4*lam)).*((1-k)*dU3im32 + (1+k)*dU3im12);
URim12(3,:) = U3i - (phii/(4*lam)).*((1+k)*dU3ip12 + (1-k)*dU3im12);

% if any(ULip12(1,:)<= 0) || any(URip12(1,:)<= 0) || any(ULim12(1,:)<= 0) || any(URim12(1,:)<= 0)
%     ULip12 = U(:,4:end-3);
%     URip12 = U(:,5:end-2);
%     ULim12 = U(:,3:end-4);
%     URim12 = U(:,4:end-3);
% end