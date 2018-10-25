function phase = getPhase(grid, point)
% Obtain phase at a given point
% grid = grid
% point = [nx; ny; nz]
    switch grid.dim
        % Dimension = 1
        case 1
            if (length(point) == 1); phase = grid.phase(point);
            else error('Wrong dimension for "point" for getPhase - dim ~= 1');
            end
        % Dimension = 2
        case 2
            if (length(point) == 2); phase = grid.phase(point(1), point(2));
            else error('Wrong dimension for "point" for getPhase - dim ~= 2');
            end
        case 3
            if (length(point) == 3); phase = grid.phase(point(1), point(2), point(3));
            else error('Wrong dimension for "point" for getPhase - dim ~= 3');
            end
        otherwise
            error('Wrong dimension for grid');
    end
end
