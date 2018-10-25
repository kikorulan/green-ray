function grid = copyGrid(gridOriginal)
% COPYGRID copies the grid GRIDORIGINAL into GRID
% grid = copyGrid(gridOriginal)
% INPUTS
% gridOriginal: gridRT object that contains the information
%    - Nx, dx, Ny, dy
%    - dt, tForward
%    - cFilter, delayFilter, tFilter, Filter
% 
% OUTPUTS
% grid: gridRT object that defines the domain
%
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke

% Dimensions, sound speed and pressure
grid = gridRT(gridOriginal.Nx, gridOriginal.dy, gridOriginal.Ny, gridOriginal.dy);
grid.setCMatrix(gridOriginal.c);
grid.setUMatrix(gridOriginal.u);
% Set time
grid.setTime(gridOriginal.dt, gridOriginal.tForward(end));
% Filter
grid.setCFilter(gridOriginal.cFilter);
grid.setDelayFilter(gridOriginal.delayFilter);
grid.setTFilter(gridOriginal.tFilter);
grid.setFilter(gridOriginal.Filter);

