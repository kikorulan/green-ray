function h = plot_rays(source, grid, nRays, scale)
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
switch scale
    case 'mm'
        factorScale = 1e3;
    case 'm'
        factorScale = 1;
    case 'km'
        factorScale = 1e-3;
end
axis([0 factorScale*grid.xAxis(end) 0 factorScale*grid.yAxis(end)]);
colourList = cool(min(nRays, totalRays));
factor = max(1, floor(totalRays/nRays));
for j = 1:min(nRays, totalRays)
    plot(factorScale*source.x(factor*j, :, 1), factorScale*source.x(factor*j, :, 2), 'Color', colourList(j, :));
end

switch scale
    case 'mm'
        xlabel('x [mm]');
        ylabel('y [mm]');
    case 'm'
        xlabel('x [m]');
        ylabel('y [m]');
    case 'km'
        xlabel('x [km]');
        ylabel('y [km]');
end     
title('Ray Trajectories');
box on;
