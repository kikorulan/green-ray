function u = getU(grid, point)
% GETU returns the value of the pressure at the given point
%function u = getU(grid, point)
% INPUTS
% grid: gridRT object that defines the domain
% point: point to find the pressure (given using coordinates in the grid)
%
% OUTPUTS
% u: pressure 
% 
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke

u = grid.u(point(:, 1, 1) + (point(:, 1, 2)-1)*grid.Nx);
