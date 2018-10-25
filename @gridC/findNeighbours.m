function coord = findNeighbours(grid, point)
% FINDNEIGHBOURS returns the coordinates of the neighbours the given point
% Multiple points can be given along the 1st dimension, while the neighbors are
% given along 2nd dimension
%function coord = findNeighbours(grid, point)
% INPUTS
% grid: gridRT object that defines the domain
% point: position in the domain in [x; y; z] coordinates
%
% OUTPUTS
% coord: coordinates of the point
% 
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke

switch grid.dim
    case 2 % The first coordinate X0 contains if the point is out of the domain
        % Previous X
        X0 = double([floor(point(:, 1, 1)/grid.dx)+1]);
        X0(X0(:, 1, 1) > grid.Nx, 1, :) = inf;
        X0(X0(:, 1, 1) < 1, 1, :) = inf;
        % Following X
        X1 = double([floor(point(:, 1, 1)/grid.dx)+2]);
        X0(X1(:, 1, 1) > grid.Nx, 1, :) = inf;
        X0(X1(:, 1, 1) < 1, 1, :) = inf;
        % Previous Y
        Y0 = double([floor(point(:, 1, 2)/grid.dy)+1]);
        X0(Y0(:, 1, 1) > grid.Ny, 1, :) = inf;
        X0(Y0(:, 1, 1) < 1, 1, :) = inf;
        % Following Y
        Y1 = double([floor(point(:, 1, 2)/grid.dy)+2]);
        X0(Y1(:, 1, 1) > grid.Ny, 1, :) = inf;
        X0(Y1(:, 1, 1) < 1, 1, :) = inf;
        coord = cat(3, [X0 X1 X0 X1], [Y0 Y0 Y1 Y1]);
    case 3 %%%%%%%%% NOT IMPLEMENTED %%%%%%%%%%%%%%%%%%
        error('findNeighbors not implemented')
    otherwise
        error('Wrong dimension for grid');
end
