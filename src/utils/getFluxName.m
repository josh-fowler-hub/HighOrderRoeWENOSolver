function fluxLim = getFluxName(flux)
% Map flux-limiter ID to a string name
switch flux
    case 1
        flux_lim = 'MM';
    case 2
        flux_lim = 'VL';
    case 3
        flux_lim = 'BJ';
    case 4
        flux_lim = 'SB';
    case 5
        flux_lim = 'VA2';
    case 6
        flux_lim = 'VA1';
    case 7
        flux_lim = 'CHARM';
    case 8
        flux_lim = 'HCUS';
    case 9
        flux_lim = 'HQUICK';
    case 10
        flux_lim = 'Koren';
    case 11
        flux_lim = 'MC';
    case 12
        flux_lim = 'Osher';
    case 13
        flux_lim = 'Ospre';
    case 14
        flux_lim = 'Smart';
    case 15
        flux_lim = 'Sweby';
    case 16
        flux_lim = 'UMIST';
    case 17
        flux_lim = 'GMM';
    case 18
        flux_lim = '';
    otherwise
        flux_lim = '';
end
end
