function h = plot_soundSpeed(grid)
% PLOT_SOUNDSPEED
% function h = plotSoundSpeed(grid)
% INPUTS
% grid: gridRT object that defines the domain
%
% OUTPUTS:
% h: figure handle
% 
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke

h = figure;
x = grid.xAxis;
y = grid.yAxis;
flipGrid = grid.c';
surf(1e3*x, 1e3*y, flipGrid, 'EdgeColor', 'none');
view(2);
box on;
xlabel('x [mm]');
ylabel('y [mm]');
%title('Sound Speed');
axis([0 1e3*x(end) 0 1e3*y(end)]);
colorbar();

