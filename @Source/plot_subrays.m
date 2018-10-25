function h = plot_subrays(source, gridR, indexRay)
% PLOT_SUBRAYS Plots the ray trajectories for the given source
%h = plot_subrays(source, grid, indexRay)
% INPUTS
% source: source to plot the rays
% grid: gridRT object that defines the domain
% nColours: number of colours to plot
%
% OUTPUTS:
% h: figure handler
% 
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke
h = figure;
hold on;
axis([0 gridR.xAxis(end) 0 gridR.yAxis(end)]);
nColours = floor(size(source.xGB, 4)/2);
colourList = cool(nColours);
plot(source.x(indexRay, :, 1), source.x(indexRay, :, 2), 'Color', colourList(1, :));
for j = 1:nColours
    plot(source.xGB(indexRay, :, 1, 2*j-1), source.xGB(indexRay, :, 2, 2*j-1), 'Color', colourList(j, :));
    plot(source.xGB(indexRay, :, 1, 2*j), source.xGB(indexRay, :, 2, 2*j), 'Color', colourList(j, :));
end
xlabel('x (m)');
ylabel('y (m)');
box on;
grid on;
title('Ray Trajectories');
