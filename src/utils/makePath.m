function folderPath = makePath(Schemes, test_num, n, cfl, flux_limiter, recon, CD_Term_Order, riemann_solver)
% Build folder path for saving figures based on configuration
if length(recon) ~= 1
    Scheme_Path = [];
    for i = 1:length(recon)
        if i == length(recon)
            Scheme_Path = [Scheme_Path recon{i} '/'];
        else
            Scheme_Path = [Scheme_Path recon{i} '_v_'];
        end
    end
    for i = 1:length(CD_Term_Order)
        if i == length(CD_Term_Order)
            cdname = getCDName(CD_Term_Order(i));
            Scheme_Path = [Scheme_Path cdname '/'];
        else
            cdname = getCDName(CD_Term_Order(i));
            Scheme_Path = [Scheme_Path cdname '_v_'];
        end
    end
    Folder_Path = ['Figures/Test#' num2str(test_num) '/'...
                    Scheme_Path num2str(n) '_Nodes' '/'...
                    'CFL=' num2str(cfl) '/'];
else
    if strcmp(recon{1}, 'MUSCL')
        if ~isscalar(flux_limiter)
            flux = ['/flux_lim='];
            k = 0;
            for f = flux_limiter
                k = k+1;
                if k == length(flux_limiter)
                    flux = [flux getFluxName(f) '/'];
                else
                    flux = [flux getFluxName(f) '_'];
                end
            end
            flux_lim = flux;
        else
            switch flux_limiter
                case 1
                    flux_lim = 'flux_lim=minmod/';
                case 2
                    flux_lim = 'flux_lim=van_Leer/';
                case 3
                    flux_lim = 'flux_lim=Barth_Jesperson/';
                case 4
                    flux_lim = 'flux_lim=Superbee/';
                case 5
                    flux_lim = 'flux_lim=van_Albada_2/';
                case 6
                    flux_lim = 'flux_lim=van_Albada_1/';
                case 7
                    flux_lim = 'flux_lim=CHARM/';
                case 8
                    flux_lim = 'flux_lim=HCUS/';
                case 9
                    flux_lim = 'flux_lim=HQUICK/';
                case 10
                    flux_lim = 'flux_lim=Koren/';
                case 11
                    flux_lim = 'flux_lim=MC/';
                case 12
                    flux_lim = 'flux_lim=Osher/';
                case 13
                    flux_lim = 'flux_lim=Ospre/';
                case 14
                    flux_lim = 'flux_lim=Smart/';
                case 15
                    flux_lim = 'flux_lim=Sweby/';
                case 16
                    flux_lim = 'flux_lim=UMIST/';
                case 17
                    flux_lim = 'flux_lim=General_Minmod/';
                case 18
                    flux_lim = 'Roe/';
                case 19
                    flux_lim = '/MultipleFluxLims/';
                otherwise
                    flux_lim = '';
            end
        end
        Scheme_Path = [];
        for i = 1:length(Schemes)
            if i == length(Schemes)
                curr = [char(Schemes(i)) '/'];
            else
                curr = [char(Schemes(i)) '_'];
            end
            Scheme_Path = [Scheme_Path curr];
        end
        for i = 1:length(CD_Term_Order)
            if i == length(CD_Term_Order)
                cdname = getCDName(CD_Term_Order(i));
                Scheme_Path = [Scheme_Path cdname '/'];
            else
                cdname = getCDName(CD_Term_Order(i));
                Scheme_Path = [Scheme_Path cdname '_v_'];
            end
        end
        Folder_Path = ['Figures/Test#' num2str(test_num) '/'...
                        Scheme_Path flux_lim num2str(n) '_Nodes' '/'...
                        'CFL=' num2str(cfl) '/'];
    else
        Scheme_Path = [recon{1} '/'];
        for i = 1:length(CD_Term_Order)
            if i == length(CD_Term_Order)
                cdname = getCDName(CD_Term_Order(i));
                Scheme_Path = [Scheme_Path cdname '/'];
            else
                cdname = getCDName(CD_Term_Order(i));
                Scheme_Path = [Scheme_Path cdname '_v_'];
            end
        end
        Folder_Path = ['Figures/Test#' num2str(test_num) '/'...
                        Scheme_Path num2str(n) '_Nodes' '/'...
                        'CFL=' num2str(cfl) '/'];
    end
end

mssg = ['Results being saved in:\n\t.../', Folder_Path, '\n\n'];

if ~exist(Folder_Path, 'dir')
    mkdir(Folder_Path)
    fprintf(mssg);
else
    fprintf(mssg);
end
end
