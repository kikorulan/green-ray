function [x, phi, n, q, p, index_ray, index_time] = trace_point(source, grid, xP)
% TRACE_POINT traces a given point to obtain the relevant information about it
% [x, q, p, index_ray, index_time] = trace_point(grid, source, xP)
%
% INPUTS
% source: source to do the conversion
% grid: object from gridRT class
% xP: point to calculate
% 
% OUTPUTS
%
% Copyright (C) Kiko Rul·lan, Marta M. Betcke

% Number of points
[nR, nPoints, dim] = size(source.x);
minDist = inf;
x = [];
phi = [];
n = [];
q = [];
p = [];
index_ray = [];
index_time = [];
% Loop over lenghtRay
for ii = 1:nPoints-1
    % Find the coordinates of the trajectories
    x0 = source.x(:, ii, :);
    xP_rep = repmat(xP, [nR, 1, 1]);
    % Compute distance to grid point
    distance = sqrt(sum((x0 - xP_rep).^2, 3));
    % Update grid
    [distanceSort, oRaysD] = sort(distance, 1, 'ascend');
    % Extract minimum
    if (distanceSort(1) < minDist)
        minDist = distanceSort(1);
        x = source.x(oRaysD(1), ii, :);
        phi = source.phi(ii);
        n = source.n(:, ii, :);
        q = source.q(oRaysD(1), ii, :);
        p = source.kIndex(oRaysD(1), ii, :);
        index_ray = oRaysD(1);
        index_time = ii;
    end
end

