function makePlots(U,xs,gamma,flux,lam,rec,Order,xidx,k,extra_prim_var_plot,riemann_solver)
% Generate plots for density, velocity, pressure, internal energy, and optional vars

% Inputs for Saving Plots
line_colors = (1/255)*[
                       0,0,0; % Black
                       65,105,225; % King Blue
                       237,29,36; % Comic Book Red
                       8,255,8; % Fluorescent Green
                       188,128,189; % Pastel Light Purple
                       255,127,0; % Orange
                       51,160,44; % Green
                       178,223,138; % Light Green
                       227,26,28; % Red
                       117,112,179; % Pastel Purple
                       231,41,138; % Bright Pink
                       128,177,211; % Pastel Blue
                       202,178,214; % Light Purple
                       166,206,227; % Light Blue
                       217,95,2; % Pastel Red
                       102,166,30; % Pastel Green
                       27,158,119; % Blue-Green
                       106,61,154; % Purple
                       251,154,153; % Pink
                       253,191,111; % Light Orange
                       255,255,153; % Light Yellow
                       177,89,40; % Brown
                       141,211,199; % Light Blue-Green
                       204,235,197; % Pastel Light Blue-Green
                       31,120,180; % Blue
                       251,128,114; % Salmon
                       ];
marker_styles = {'o','o','o','o','o','o','o','o','o','o','o','o','o',...
                 'o','o','o','x','x','x','x','x','x','x','x','x''x','x',...
                 'x','x','x','x','x'};
line_styles = {'-','-','-','-','-','-','-','-','-','-','-','-','-','-',...
               '-','-','--','--','--','--','--','--','--','--','--',...
               '--','--','--','--','--','--','--'};

[r, u, p, a, H, e, m, s] = consToPrim(U(:, xidx), gamma);

name = getLegendEntryName(flux,lam,rec,Order,riemann_solver);

figure(1)
set(gcf,'renderer','painters');
plot(xs(xidx),r,'Color', line_colors(k,:), 'LineStyle',...
    line_styles{k},...
    'LineWidth', 8, 'MarkerSize',4,'DisplayName',name);

figure(2)
set(gcf,'renderer','painters');
plot(xs(xidx),u,'Color', line_colors(k,:), 'LineStyle',...
    line_styles{k},...
    'LineWidth', 8, 'MarkerSize',4,'DisplayName',name);

figure(3)
set(gcf,'renderer','painters');
plot(xs(xidx),p,'Color', line_colors(k,:), 'LineStyle',...
    line_styles{k},...
    'LineWidth', 8, 'MarkerSize',4,'DisplayName',name);

figure(4)
set(gcf,'renderer','painters');
plot(xs(xidx),e,'Color', line_colors(k,:), 'LineStyle',...
    line_styles{k},...
    'LineWidth', 8, 'MarkerSize',4,'DisplayName',name);

if extra_prim_var_plot == true
    figure(5)
    set(gcf,'renderer','painters');
    plot(xs(xidx),a,'Color', line_colors(k,:),...
        'LineStyle', line_styles{k},'LineWidth', 8, 'MarkerSize',...
        12,'DisplayName',name);

    figure(6)
    set(gcf,'renderer','painters');
    plot(xs(xidx),m,'Color', line_colors(k,:),...
        'LineStyle', line_styles{k},'LineWidth', 8, 'MarkerSize',...
        12,'DisplayName',name);

    figure(7)
    set(gcf,'renderer','painters');
    plot(xs(xidx),H,'Color', line_colors(k,:),...
        'LineStyle', line_styles{k},'LineWidth', 8, 'MarkerSize',...
        12,'DisplayName',name);

    figure(8)
    set(gcf,'renderer','painters');
    plot(xs(xidx),s,'Color', line_colors(k,:),...
        'LineStyle', line_styles{k},'LineWidth', 8, 'MarkerSize',...
        12,'DisplayName',name);
end
end
