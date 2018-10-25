function [n, Gn] = getNGN(grid, point)
% GETNGN returns the value of eta at and its gradient the given point
%function [n, Gn] = getNGN(grid, point)
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
        coordLin = point(:, 1, 1) + (point(:, 1, 2)-1)*grid.Nx;
        % Eta
        n(:, :, 2) = grid.n(coordLin);
        n(:, :, 1) = n(:, :, 2);
        % Gradient
        Gn(:, :, 2) = grid.Gn(coordLin + grid.Nx*grid.Ny);
        Gn(:, :, 1) = grid.Gn(coordLin);
    % Dimension = 3
    case 3 %%%%%% NOT IMPLEMENTED
        error('Not implemented');        
    otherwise
        error('Wrong dimension for grid');
end
