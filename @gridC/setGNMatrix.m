function grid = setGNMatrix(grid)
% SETGNMATRIX sets the matrix of the gradient of eta
% grid = setGNMatrix(grid)
% INPUTS
% 
% OUTPUTS
% 
% Copyright (C) 2017 Kiko Rul·lan and Marta M. Betcke
switch grid.dim
    % Dimension = 2
    case 2
        % Preallocate
        grid.Gn = zeros(grid.Nx, grid.Ny, 2);
        % Gradient of Eta - x component
        gradX = (grid.n(2:end, :) - grid.n(1:end-1, :))/grid.dx;
        gradCX = (grid.n(3:end, :) - grid.n(1:end-2, :))/(2*grid.dx);
        grid.Gn(1, :, 1) = gradX(1, :);
        grid.Gn(end, :, 1) = gradX(end, :);
        grid.Gn(2:end-1, :, 1) = gradCX;
        % Gradient of Eta - y component
        gradY = (grid.n(:, 2:end) - grid.n(:, 1:end-1))/grid.dy;
        gradCY = (grid.n(:, 3:end) - grid.n(:, 1:end-2))/(2*grid.dy);
        grid.Gn(:, 1, 2) = gradY(:, 1);
        grid.Gn(:, end, 2) = gradY(:, end);
        grid.Gn(:, 2:end-1, 2) = gradCY;
    case 3 %%%% NOT IMPLEMENTED
        %disp('Not implemented');
end

