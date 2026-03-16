function phi = limitFlux(flux_limiter, r)
% limitFlux - Flux limiter functions for MUSCL reconstruction.
%
% Inputs:
%   flux_limiter - integer limiter selector
%   r            - smoothness ratio (vector)
%
% Output:
%   phi - limiter value (same size as r)

switch flux_limiter
    case 1 % Minmod
        phi = r;
        phi(phi>1) = 1;
        phi(phi<0) = 0;
    case 2 % van Leer
        phi = (abs(r) + r) ./ (1 + abs(r));
    case 3 % Barth-Jesperson
        first_min = 4 ./ (r + 1);
        first_min(first_min > 1) = 1;
        second_min = 4 .* r ./ (r + 1);
        second_min(second_min > 1) = 1;
        third_min = second_min;
        third_min(third_min > first_min) = first_min(first_min < third_min);
        phi = 0.5 .* (r + 1) .* third_min;
    case 4 % Superbee
        first_min = r;
        first_min(first_min > 2) = 2;
        second_min = 2 .* r;
        second_min(second_min > 1) = 1;
        third_min = first_min;
        third_min(third_min < second_min) = second_min(second_min > third_min);
        third_min(third_min < 0) = 0;
        phi = third_min;
    case 5 % van Albada 2, not 2nd order TVD
        phi = 2 .* r ./ (r.^2 + 1);
    case 6 % van Albada 1
        phi = (r + r.^2) ./ (1 + r.^2);
    case 7 % CHARM, not 2nd order TVD
        phi = r;
        phi = (phi .* (3 * phi + 1)) ./ ((phi + 1).^2);
        phi(r < 0) = 0;
    case 8 % HCUS, not 2nd order TVD
        phi = (1.5 .* (r + abs(r))) ./ (r + 2);
    case 9 % HQUICK, not 2nd order TVD
        phi = (2 .* (r + abs(r))) ./ (r + 3);
    case 10 % Koren
        first_min = (1 + 2 .* r) ./ 3;
        first_min(first_min > 2) = 2;
        second_min = 2 .* r;
        first_min(first_min > second_min) = second_min(second_min < first_min);
        first_min(first_min < 0) = 0;
        phi = first_min;
    case 11 % monotonized central (MC)
        first_min = 2 .* r;
        second_min = 0.5 .* (1 + r);
        first_min(first_min > second_min) = second_min(second_min < first_min);
        first_min(first_min > 2) = 2;
        first_min(first_min < 0) = 0;
        phi = first_min;
    case 12 % Osher
        beta = 1.5;
        phi = r;
        phi(phi > beta) = beta;
        phi(phi < 0) = 0;
    case 13 % ospre
        phi = (1.5 .* (r.^2 + r)) ./ (r.^2 + r + 1);
    case 14 % smart
        first_min = 2 .* r;
        second_min = (0.25 + 0.75 .* r);
        first_min(first_min > second_min) = second_min(second_min < first_min);
        first_min(first_min > 4) = 4;
        first_min(first_min < 0) = 0;
        phi = first_min;
    case 15 % Sweby
        beta = 1.5;
        first_min = r;
        first_min(first_min > beta) = beta;
        second_min = beta .* r;
        second_min(second_min > 1) = 1;
        first_min(first_min < second_min) = second_min(second_min > first_min);
        first_min(first_min < 0) = 0;
        phi = first_min;
    case 16 % UMIST
        first_min = 2 .* r;
        second_min = (0.25 + 0.75 .* r);
        third_min = (0.75 + 0.25 .* r);
        first_min(first_min > second_min) = second_min(second_min < first_min);
        first_min(first_min > third_min) = third_min(third_min < first_min);
        first_min(first_min > 2) = 2;
        first_min(first_min < 0) = 0;
        phi = first_min;
    case 17 % Generalized Minmod
        beta = 2; % beta is between 1 and 2
        first_min = beta .* r;
        second_min = (1 + r) ./ 2;
        first_min(first_min > second_min) = second_min(second_min < first_min);
        first_min(first_min > beta) = beta;
        first_min(first_min < 0) = 0;
        phi = first_min;
    case 18
        phi = zeros(size(r));
    otherwise
        error('Unknown flux limiter: %d', flux_limiter);
end
