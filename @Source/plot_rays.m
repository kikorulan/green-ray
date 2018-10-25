function h = plot_rays(source, grid, nRays)
% PLOT_RAYS Plots the ray trajectories for the given sources
%function h = plot_rays(source, grid, nColours)
% INPUTS
% source: source to plot the rays
% grid: gridRT object that defines the domain
% nRays: number of rays to plot
%
% OUTPUTS:
% h: figure handler
% 
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke
totalRays = size(source.x, 1);
h = figure;
hold on;
axis([0 1e3*grid.xAxis(end) 0 1e3*grid.yAxis(end)]);
colourList = cool(min(nRays, totalRays));
factor = max(1, floor(totalRays/nRays));
for j = 1:min(nRays, totalRays)
    plot(1e3*source.x(factor*j, :, 1), 1e3*source.x(factor*j, :, 2), 'Color', colourList(j, :));
end
xlabel('x [mm]');
ylabel('y [mm]');
title('Ray Trajectories');
