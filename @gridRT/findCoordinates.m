function coord = findCoordinates(grid, point)
% FINDCOORDINATES returns the coordinates of "point" at the grid.
% Multiple points can be given along the 3rd dimension
%function coord = findCoordinates(grid, point)
% INPUTS
% grid: gridRT object that defines the domain
% point: position in the domain in [x; y; z] coordinates
%
% OUTPUTS:
% coord: coordinates of the point
% 
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke

coord(:, 1, 2) = round(point(:, 1, 2)/grid.dy) + 1;
coord(:, 1, 1) = round(point(:, 1, 1)/grid.dx) + 1;

coord(coord(:, 1, 1) > grid.Nx | coord(:, 1, 1) < 1 | coord(:, 1, 2) > grid.Ny | coord(:, 1, 2) < 1, 1, :) = inf;
