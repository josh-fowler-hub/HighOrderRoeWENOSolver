function name = getLegendEntryName(flux,lam,rec,CD_Term_Order,riemann_solver)
% Generate legend entry label for plotting based on scheme options
if strcmp(riemann_solver, 'Roe')
    solver = 'Roe';
elseif strcmp(riemann_solver, 'Roe-Pike')
    solver = 'Roe-Pike';
elseif strcmp(riemann_solver, 'HLL')
    solver = 'HLL';
else
    solver = '';
end

if strcmp(rec,'ROE')
    if CD_Term_Order == 1
        name = [solver];
    else
        cdname = getCDName(CD_Term_Order);
        name = [cdname '--' solver];
    end
elseif strcmp(rec,'MUSCL')
    flux_lim = getFluxName(flux);
    if flux == 18
        if CD_Term_Order == 1
            name = [solver ' ' flux_lim];
        else
            name = ['CD' num2str(CD_Term_Order) '-' solver ' ' flux_lim];
        end
    else
        if CD_Term_Order == 1
            name = [solver '--' rec];
        else
            name = ['CD' num2str(CD_Term_Order) '--' solver '--' rec];
        end
    end
else
    if CD_Term_Order == 1
        name = [solver '--' rec];
    else
        name = ['CD' num2str(CD_Term_Order) '--' solver '--' rec];
    end
end
end
