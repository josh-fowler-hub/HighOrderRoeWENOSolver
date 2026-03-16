function lambda = entropyFixScalar(lambda, delta)
% entropyFixScalar - apply Harten-Hyman entropy fix
%
% This helper performs the Harten-Hyman entropy fix on the eigenvalues
% used in a Roe linearization. It preserves the correct amount of
% dissipation near sonic points.
%
% Inputs:
%   lambda - array of eigenvalues
%   delta  - entropy-fix threshold (same size as lambda)
%
% Output:
%   lambda - modified eigenvalues after applying the entropy fix

absLambda = abs(lambda);
idx = absLambda < delta;
lambda(idx) = (lambda(idx).^2 + delta(idx).^2) ./ (2 .* delta(idx));
end
