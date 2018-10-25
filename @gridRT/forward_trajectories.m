function grid = forward_trajectories(grid, source)
% FORWARD_TRAJECTOREIS computes the ray trajectories of the given source nS
%function grid = forward_trajectories(grid, source)
% INPUTS
% grid: gridRT object that defines the domain
% source: source to compute the ray trajectories
%
% OUTPUTS:
% grid: gridRT object that defines the domain
% 
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke

% Number of steps of the computation
step = source.step;
[nR, nPoints, dim] = size(source.x);

% Insert time
phi = 0:step:(nPoints-1)*step;
source.insertPhiVector(phi);

% Initial points
Xn = source.x(:, 1, :);
Pn = source.p(:, 1, :);
[Xn, Pn] = grid.stepRK4Hamil(step/100, Xn, Pn);
source.insertX(Xn, 1);
source.insertP(Pn, 1);
% Insert traveled distance
coord = grid.findCoordinates(source.x(:, 1, :));
n = grid.getN(coord);
source.insertXTD(step/100./n, 1);

% Main Ray
for tIndex = 1:nPoints-1
    %[Xn, Pn] = grid.stepRK2Hamil(step, Xn, Pn);
    [Xn, Pn] = grid.stepRK4Hamil(step, Xn, Pn);
    coord = grid.findCoordinates(Xn);
    binaryCoord = coord(:, 1, 1) < inf;
    % Insert X
    source.insertX(Xn, tIndex + 1);
    % Insert P
    source.insertP(Pn, tIndex + 1);
    % Insert n
    n = coord(:, :, 1);
    n(binaryCoord, :, :) = grid.getN(coord(binaryCoord, :, :));
    source.insertN(n, tIndex + 1);
    % Insert Traveled distance
    normTrajX = step./n;
    source.insertXTD(source.xTD(:, tIndex, :) + normTrajX, tIndex+1);
    % Insert pressure
    pressure = coord(:, :, 1);
    pressure(binaryCoord, :, :) = grid.getU(coord(binaryCoord, :, :));
    source.insertPressure(pressure, tIndex + 1);
end
