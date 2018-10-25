function source = newSource(grid, x0, angleMin, angleMax, nR, step, tMax)
% NEWSOURCE inserts a new source in the grid with the given parameters
%function source = newSource(grid, x0, angleMin, angleMax, nR, step, tMax)
% INPUTS
% grid: gridRT object that defines the domain
% x0: source point
% angleMin: minimum shooting angle
% angleMax: maximum shooting angle
% nR: number of rays shot
% step: step taken for computing the rays from this source
% tMax: maximum value for the ray parameter
%
% OUTPUTS:
% source: created source
% 
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke

nPoints = floor(tMax/step) + 1;
point = grid.findCoordinates(x0);
n = grid.getN(point);
nGn = grid.gradEta(point).*n;

source = Source(x0, n, nGn, 0, grid.deltaX, grid.deltaP, angleMin, angleMax, nR, step, nPoints);

