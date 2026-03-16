function [lammt1,lampt3] = entropyFix(uLim1,uRi,aLim1,aRi,ustar,astarL,astarR,lammt1,lampt3)
% entropyFix - Apply entropy fix to characteristic speeds (Toro 9.5.2)

lamm1L = uLim1 - aLim1;
lamm1R = ustar - astarL;
lamp3L = ustar + astarR;
lamp3R = uRi + aRi;

% Left Transonic Rarefaction: lambda1L < 0 < lambda1R
maskLTR = (lamm1L < 0) & (lamm1R > 0);
lammt1(maskLTR) = lamm1L(maskLTR) .* ((lamm1R(maskLTR) - lammt1(maskLTR)) ./ (lamm1R(maskLTR) - lamm1L(maskLTR)));

% Right Transonic Rarefaction: lambda3L < 0 < lambda3R
maskRTR = (lamp3L < 0) & (lamp3R > 0);
lampt3(maskRTR) = lamp3R(maskRTR) .* ((lampt3(maskRTR) - lamp3L(maskRTR)) ./ (lamp3R(maskRTR) - lamp3L(maskRTR)));
end