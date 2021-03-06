function grid = setU(grid, point, uVal)
% Given an Alive grid point, tag as Close its given neighbour and add it to
% the Heap Loop
% grid = grid
% point = [nx; ny; nz]
    switch grid.dim
        % Dimension = 1
        case 1
            if (length(point) == 1); grid.u(point) = uVal;
            else error('Wrong dimension for "point" for tagging - dim ~= 1');
            end
        % Dimension = 2
        case 2
            if (length(point) == 2); grid.u(point(1), point(2)) = uVal;
            else error('Wrong dimension for "point" for tagging - dim ~= 2');
            end
        case 3
            if (length(point) == 3); grid.u(point(1), point(2), point(3)) = uVal;
            else error('Wrong dimension for "point" for tagging - dim ~= 3');
            end
        otherwise
            error('Wrong dimension for grid');
    end
end

