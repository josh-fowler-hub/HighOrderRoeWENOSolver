function configurePlots(recon,test_num,cells,cfl,extra_prim_var_plot,riemann_solver,CD_Term_Order)
% Configure plot formatting, axis limits, and save EPS files.
file_pre = makeFilePrefix(recon, test_num, cells, cfl, riemann_solver, CD_Term_Order);
locs = getLocations(test_num,CD_Term_Order);
% Folder = '../thesis/CH5/EPSFDocs/';
% Folder = '';
Folder = '../thesis/FigsforCommittee/';

figure(1)
set(gcf, 'Position', get(0, 'Screensize'));
ax = gca;
ax.YAxis.FontSize = 55;
ax.XAxis.FontSize = 55;
if CD_Term_Order == 2
    ylim1 = [0,0,0,4,0,0,0,0.85,0,0,0,0,0.15,0,0,0,0,0,0];
    ylim2 = [1.1,0,7,35,7,1.1,0,1.2,0,0,0,0,1.45,0,1.1,1.1,0,0,0];
elseif CD_Term_Order == 4
    ylim1 = [0,0,0,4,0,0,0,0.875,0,0,0,0,0.2,0,0,0,0,0,0];
    ylim2 = [1.1,0,7,35,7,1.1,0,1.16,0,0,0,0,1.4,0,1.1,1.1,0,0,0];
else
    ylim1 = [0,0,0,4,0,0,0,0.875,0,0,0,0,0.2,0,0,0,0,0,0];
    ylim2 = [1.1,0,7,35,7,1.1,0,1.16,0,0,0,0,1.4,0,1.1,1.1,0,0,0];
end
ax.YAxis.Limits = [ylim1(test_num) ylim2(test_num)];
ax.LineWidth = 10;
filename = [Folder file_pre 'Density.eps'];
legend('FontSize', 60,'Interpreter','Latex')
legend('boxoff')
leg = legend();
leg.ItemTokenSize = [80,25];
leg.Position = locs(1,:);
xlabel('$x$','interpreter','latex', 'FontSize', 100, 'FontWeight','bold');
ylabel('$\rho$','interpreter','latex', 'FontSize', 100, 'FontWeight','bold');
grid off;
saveas(gcf, filename, 'epsc');

figure(2)
set(gcf, 'Position', get(0, 'Screensize'));
ax = gca;
ax.YAxis.FontSize = 55;
ax.XAxis.FontSize = 55;
if CD_Term_Order == 2
    ylim1 = [-0.2,0,-3,-10,-24,-0.05,0,-0.55,0,0,0,0,-0.1,0,-0.2,-2.1,0,0,0];
    ylim2 = [1.6,0,23,25,4,1,0,0.05,0,0,0,0,1.7,0,1.8,-0.85,0,0,0];
elseif CD_Term_Order == 4
    ylim1 = [-0.2,0,-3,-10,-24,-0.05,0,-0.475,0,0,0,0,-0.1,0,-0.2,-2.15,0,0,0];
    ylim2 = [1.6,0,23,25,4,1,0,0.05,0,0,0,0,1.6,0,1.8,-1,0,0,0];
else
    ylim1 = [-0.2,0,-3,-10,-24,-0.05,0,-0.475,0,0,0,0,-0.1,0,-0.2,-2.1,0,0,0];
    ylim2 = [1.6,0,23,25,4,1,0,0.05,0,0,0,0,1.6,0,1.8,-1,0,0,0];
end
ax.YAxis.Limits = [ylim1(test_num) ylim2(test_num)];
ax.LineWidth = 10;
filename = [Folder file_pre 'Velocity.eps'];
legend('FontSize', 60,'Interpreter','Latex')
legend('boxoff')
leg = legend();
leg.ItemTokenSize = [80,25];
leg.Position = locs(2,:);
xlabel('$x$','interpreter','latex', 'FontSize', 100, 'FontWeight','bold');
ylabel('$u$','interpreter','latex', 'FontSize', 100, 'FontWeight','bold');
grid off;
saveas(gcf, filename, 'epsc');

figure(3)
set(gcf, 'Position', get(0, 'Screensize'));
ax = gca;
ax.YAxis.FontSize = 55;
ax.XAxis.FontSize = 55;
if CD_Term_Order == 2
    ylim1 = [0,0,-100,-100,-100,0,0,6.75,0,0,0,0,0.5,0,0,0,0,0,0];
    ylim2 = [1.1,0,1100,1950,1100,1.1,0,10.25,0,0,0,0,3.75,0,1.1,1.1,0,0,0];
elseif CD_Term_Order == 4
    ylim1 = [0,0,-100,-100,-100,0,0,6.75,0,0,0,0,0.5,0,0,0,0,0,0];
    ylim2 = [1.1,0,1100,1900,1100,1.1,0,10.25,0,0,0,0,3.75,0,1.1,1.1,0,0,0];
else
    ylim1 = [0,0,-100,-100,-100,0,0,6.75,0,0,0,0,0.5,0,0,0,0,0,0];
    ylim2 = [1.1,0,1100,1900,1100,1.1,0,10.25,0,0,0,0,3.75,0,1.1,1.1,0,0,0];
end
ax.YAxis.Limits = [ylim1(test_num) ylim2(test_num)];
ax.LineWidth = 10;
filename = [Folder file_pre 'Pressure.eps'];
legend('FontSize', 60,'Interpreter','Latex')
legend('boxoff')
leg = legend();
leg.ItemTokenSize = [80,25];
leg.Position = locs(3,:);
xlabel('$x$','interpreter','latex', 'FontSize', 100, 'FontWeight','bold');
ylabel('$p$','interpreter','latex', 'FontSize', 100, 'FontWeight','bold');
grid off;
saveas(gcf, filename, 'epsc');

figure(4)
set(gcf, 'Position', get(0, 'Screensize'));
ax = gca;
ax.YAxis.FontSize = 55;
ax.XAxis.FontSize = 55;
if CD_Term_Order == 2
    ylim1 = [1.8,0,-200,-20,-100,1.6,0,17,0,0,0,0,2,0,1,1.6,0,0,0];
    ylim2 = [3.7,0,2700,400,2600,3,0,25.5,0,0,0,0,32,0,5,3.1,0,0,0];
elseif CD_Term_Order == 4
    ylim1 = [1.8,0,-200,-20,-100,1.6,0,17,0,0,0,0,2,0,1,1.6,0,0,0];
    ylim2 = [3.7,0,2700,340,2600,3,0,25.5,0,0,0,0,20,0,5,3,0,0,0];
else
    ylim1 = [1.8,0,-200,-20,-100,1.6,0,17,0,0,0,0,2,0,1,1.75,0,0,0];
    ylim2 = [3.7,0,2700,340,2600,3,0,25.5,0,0,0,0,20,0,5,3,0,0,0];
end
ax.YAxis.Limits = [ylim1(test_num) ylim2(test_num)];
ax.LineWidth = 10;
filename = [Folder file_pre 'InternalEnergy.eps'];
legend('FontSize', 60,'Interpreter','Latex')
legend('boxoff')
leg = legend();
leg.ItemTokenSize = [80,25];
leg.Position = locs(4,:);
xlabel('$x$','interpreter','latex', 'FontSize', 100, 'FontWeight','bold');
ylabel('$e$','interpreter','latex', 'FontSize', 100, 'FontWeight','bold');
grid off;
saveas(gcf, filename, 'epsc');

if extra_prim_var_plot
    figure(5)
    ax = gca;
    ax.YAxis.FontSize = 55;
    ax.XAxis.FontSize = 55;
    ymin = ax.YAxis.Limits(1);
    ymax = ax.YAxis.Limits(2);
    yticks = round(linspace(ymin,ymax,6),1);
    ax.YAxis.TickValues = yticks;
    scale = max(numel(num2str(yticks(3))), numel(num2str(yticks(4))));
    scale = 0.1*scale + 0.01;
    xloc = (ax.XAxis.Limits(1)) - scale*ax.YAxis.FontSize*0.01;
    yloc = (ax.YAxis.Limits(2) - ax.YAxis.Limits(1))/2  + ax.YAxis.Limits(1);
    ax.YAxis.Label.Position = [xloc yloc -1];
    ax.LineWidth = 10;
    filename = [Folder file_pre 'SpeedofSound.eps'];
    legend('FontSize', 60)
    legend('boxoff')
    leg = legend();
    leg.ItemTokenSize = [80,25];
    leg.Position = locs(5,:);
    xlabel('$x$','interpreter','latex', 'FontSize', 100, 'FontWeight','bold');
    ylabel('$a$','interpreter','latex', 'FontSize', 100, 'FontWeight','bold');
    grid off;
    saveas(gcf, filename, 'epsc');

    figure(6)
    ax = gca;
    ax.YAxis.FontSize = 55;
    ax.XAxis.FontSize = 55;
    ymin = ax.YAxis.Limits(1);
    ymax = ax.YAxis.Limits(2);
    yticks = round(linspace(ymin,ymax,6),1);
    ax.YAxis.TickValues = yticks;
    scale = max(numel(num2str(yticks(3))), numel(num2str(yticks(4))));
    scale = 0.1*scale + 0.01;
    xloc = (ax.XAxis.Limits(1)) - scale*ax.YAxis.FontSize*0.01;
    yloc = (ax.YAxis.Limits(2) - ax.YAxis.Limits(1))/2  + ax.YAxis.Limits(1);
    ax.YAxis.Label.Position = [xloc yloc -1];
    ax.LineWidth = 10;
    filename = [Folder file_pre 'Mach.eps'];
    legend('FontSize', 60)
    legend('boxoff')
    leg = legend();
    leg.ItemTokenSize = [80,25];
    leg.Position = locs(6,:);
    xlabel('$x$','interpreter','latex', 'FontSize', 100, 'FontWeight','bold');
    ylabel('$M$','interpreter','latex', 'FontSize', 100, 'FontWeight','bold');
    grid off;
    saveas(gcf, filename, 'epsc');

    figure(7)
    ax = gca;
    ax.YAxis.FontSize = 55;
    ax.XAxis.FontSize = 55;
    ymin = ax.YAxis.Limits(1);
    ymax = ax.YAxis.Limits(2);
    yticks = round(linspace(ymin,ymax,6),1);
    ax.YAxis.TickValues = yticks;
    scale = max(numel(num2str(yticks(3))), numel(num2str(yticks(4))));
    scale = 0.1*scale + 0.01;
    xloc = (ax.XAxis.Limits(1)) - scale*ax.YAxis.FontSize*0.01;
    yloc = (ax.YAxis.Limits(2) - ax.YAxis.Limits(1))/2  + ax.YAxis.Limits(1);
    ax.YAxis.Label.Position = [xloc yloc -1];
    ax.LineWidth = 10;
    filename = [Folder file_pre 'Enthalpy.eps'];
    legend('FontSize', 60)
    legend('boxoff')
    leg = legend();
    leg.ItemTokenSize = [80,25];
    leg.Position = locs(7,:);
    xlabel('$x$','interpreter','latex', 'FontSize', 100, 'FontWeight','bold');
    ylabel('$H$','interpreter','latex', 'FontSize', 100, 'FontWeight','bold');
    grid off;
    saveas(gcf, filename, 'epsc');

    figure(8)
    ax = gca;
    ax.YAxis.FontSize = 55;
    ax.XAxis.FontSize = 55;
    ymin = ax.YAxis.Limits(1);
    ymax = ax.YAxis.Limits(2);
    yticks = round(linspace(ymin,ymax,6),1);
    ax.YAxis.TickValues = yticks;
    scale = max(numel(num2str(yticks(3))), numel(num2str(yticks(4))));
    scale = 0.1*scale + 0.01;
    xloc = (ax.XAxis.Limits(1)) - scale*ax.YAxis.FontSize*0.01;
    yloc = (ax.YAxis.Limits(2) - ax.YAxis.Limits(1))/2  + ax.YAxis.Limits(1);
    ax.YAxis.Label.Position = [xloc yloc -1];
    ax.LineWidth = 10;
    filename = [Folder file_pre 'Entropy.eps'];
    legend('FontSize', 60)
    legend('boxoff')
    leg = legend();
    leg.ItemTokenSize = [80,25];
    leg.Position = locs(8,:);
    xlabel('$x$','interpreter','latex', 'FontSize', 100, 'FontWeight','bold');
    ylabel('$s$','interpreter','latex', 'FontSize', 100, 'FontWeight','bold');
    grid off;
    saveas(gcf, filename, 'epsc');
end
end
