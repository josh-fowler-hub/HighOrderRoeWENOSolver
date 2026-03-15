function ticks = findTicks(umin,umax)
% Generate tick values based on axis limits for consistent plotting
ymin = floor(min(umin));
ymax = ceil(max(umax));
scaling = 10^(max(numel(num2str(abs(ymin))), numel(num2str(abs(ymax)))) - 1);
innerticks = scaling*0.2;
ylim1 = ymin - innerticks;
ylim2 = ymax + innerticks;
ticks = ylim1:innerticks:ylim2;
end
