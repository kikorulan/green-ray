function grad = gradSpeed(grid, point)
% GRADETA computes the gradient of the sound speed at the given point using central differences
%function grad = gradSpeed(grid, point)
% INPUTS
% grid: gridRT object that defines the domain
% point: point to compute the gradient (given using coordinates in the grid)
%
% OUTPUTS
% grad: gradient
% 
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke

switch grid.dim
    % Dimension = 2
    case 2
        % Y component
        Ynext = min(grid.Ny, point(:, 1, 2)+1);
        Yprev = max(1, point(:, 1, 2)-1);
        grad(:, :, 2) = (grid.c(point(:, 1, 1) + (Ynext-1)*grid.Nx) - ...
                  grid.c(point(:, 1, 1) + (Yprev-1)*grid.Nx))./double(Ynext-Yprev)./grid.dy;
        % X component
        Xnext = min(grid.Nx, point(:, 1, 1)+1);
        Xprev = max(1, point(:, 1, 1)-1);
        grad(:, :, 1) = (grid.c(Xnext + (point(:, 1, 2)-1)*grid.Nx) - ...
                  grid.c(Xprev + (point(:, 1, 2)-1)*grid.Nx))./double(Xnext-Xprev)./grid.dx;
    % Dimension = 3
    case 3
        error('Not implemented');
    otherwise
        error('Wrong dimension for grid');
end

