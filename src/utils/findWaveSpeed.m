function smax = findWaveSpeed(U, gamma)
% Compute max characteristic speed for CFL condition
UL = U(:,1:end-1);
UR = U(:,2:end);
[rL, uL, ~, ~, HL, ~, ~, ~] = consToPrim(UL, gamma);
[rR, uR, ~, ~, HR, ~, ~, ~] = consToPrim(UR, gamma);
ut = (sqrt(rL).*uL + sqrt(rR).*uR)./(sqrt(rL) + sqrt(rR));
Ht = (sqrt(rL).*HL + sqrt(rR).*HR)./(sqrt(rL) + sqrt(rR));
at = sqrt((gamma-1).*(Ht - 0.5.*ut.^2));
SL = max(abs(ut - at));
SR = max(abs(ut + at));
smax = max(SL,SR);
end
