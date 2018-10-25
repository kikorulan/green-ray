function n = getN(grid, point)
% GETN returns the value of eta at the given point
%function n = getN(grid, point)
% INPUTS
% grid: gridRT object that defines the domain
% point: point to find eta (given using coordinates in the grid)
%
% OUTPUTS
% n: eta (inverse of sound speed)
% 
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke

switch grid.dim
    % Dimension = 2
    case 2
        n = grid.n(point(:, 1, 1) + (point(:, 1, 2)-1)*grid.Nx);
    % Dimension = 3
    case 3
        n = grid.n(point(:, 1, 1) + (point(:, 1, 2)-1)*grid.Nx ...
            + (point(:, 1, 3)-1)*grid.Nx*grid.Ny);
    otherwise
        error('Wrong dimension for grid');
end
