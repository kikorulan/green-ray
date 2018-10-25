function grid = setC(grid, point, cVal)
% grid = grid
% point = [nx; ny; nz]
    switch grid.dim
        % Dimension = 2
        case 2
            if (length(point) == 2); 
                grid.c(point(1), point(2)) = cVal;
                grid.n(point(1), point(2)) = 1/cVal;
            else error('Wrong dimension for "point" for tagging - dim ~= 2');
            end
        case 3
            if (length(point) == 3); grid.c(point(1), point(2), point(3)) = cVal;
            else error('Wrong dimension for "point" for tagging - dim ~= 3');
            end
        otherwise
            error('Wrong dimension for grid');
    end
end

