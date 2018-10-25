function c = getC(grid, point)
% Obtain c at a given point
% grid = grid
% point = [nx; ny; nz]
    switch grid.dim
        % Dimension = 2
        case 2
            c = grid.c(point(:, 1, 1) + (point(:, 1, 2)-1)*grid.Nx);
        % Dimension = 3
        case 3
            c = grid.c(point(:, 1, 1) + (point(:, 1, 2)-1)*grid.Nx ...
                + (point(:, 1, 3)-1)*grid.Nx*grid.Ny);
        otherwise
            error('Wrong dimension for grid');
    end
end
