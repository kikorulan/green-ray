function grad = gradPhase(grid, point)
% Compute the gradient of the phase at a given point
    grad = [];
    switch grid.dim
        % Dimension = 1
        case 1
            if (length(point) == 1);
                % Compute the gradient
                if (point == 1)
                    grad_i = (grid.u(point+1) - grid.u(point))/grid.dx;
                elseif (point == grid.Nx)
                    grad_i = (grid.u(point) - grid.u(point-1))/grid.dx;
                else
                    grad_i = (grid.u(point+1) - grid.u(point-1))*0.5/grid.dx;
                end
                grad = [grad; grad_i];
            else error('Wrong dimension for "point" for getTag - dim ~= 1');
            end
        % Dimension = 2
        case 2
            if (length(point) == 2); 
                % Compute the gradient
                % X component
                if (point(1) == 1)
                    grad_i = (grid.u(point(1)+1, point(2)) - grid.u(point(1), point(2)))/grid.dx;
                elseif (point(1) == grid.Nx)
                    grad_i = (grid.u(point(1), point(2)) - grid.u(point(1)-1, point(2)))/grid.dx;
                else
                    grad_i = (grid.u(point(1)+1, point(2)) - grid.u(point(1)-1, point(2)))*0.5/grid.dx;
                end
                grad = [grad; grad_i];
                % Y component
                if (point(2) == 1)
                    grad_i = (grid.u(point(1), point(2)+1) - grid.u(point(1), point(2)))/grid.dy;
                elseif (point(2) == grid.Ny)
                    grad_i = (grid.u(point(1), point(2)) - grid.u(point(1), point(2)-1))/grid.dy;
                else
                    grad_i = (grid.u(point(1), point(2)+1) - grid.u(point(1), point(2)-1))*0.5/grid.dy;
                end
                grad = [grad; grad_i];
                
            else error('Wrong dimension for "point" for getTag - dim ~= 2');
            end
        case 3
            if (length(point) == 3); 
                % Compute the gradient
                % X component
                if (point(1) == 1)
                    grad_i = (grid.u(point(1)+1, point(2), point(3)) - grid.u(point(1), point(2), point(3)))/grid.dx;
                elseif (point(1) == grid.Nx)
                    grad_i = (grid.u(point(1), point(2), point(3)) - grid.u(point(1)-1, point(2), point(3)))/grid.dx;
                else
                    grad_i = (grid.u(point(1)+1, point(2), point(3)) - grid.u(point(1)-1, point(2), point(3)))*0.5/grid.dx;
                end
                grad = [grad; grad_i];
                % Y component
                if (point(2) == 1)
                    grad_i = (grid.u(point(1), point(2)+1, point(3)) - grid.u(point(1), point(2), point(3)))/grid.dy;
                elseif (point(2) == grid.Ny)
                    grad_i = (grid.u(point(1), point(2), point(3)) - grid.u(point(1), point(2)-1, point(3)))/grid.dy;
                else
                    grad_i = (grid.u(point(1), point(2)+1, point(3)) - grid.u(point(1), point(2)-1, point(3)))*0.5/grid.dy;
                end
                grad = [grad; grad_i];
                % Z component
                if (point(3) == 1)
                    grad_i = (grid.u(point(1), point(2), point(3)+1) - grid.u(point(1), point(2), point(3)))/grid.dz;
                elseif (point(3) == grid.Ny)
                    grad_i = (grid.u(point(1), point(2), point(3)) - grid.u(point(1), point(2), point(3)-1))/grid.dz;
                else
                    grad_i = (grid.u(point(1), point(2), point(3)+1) - grid.u(point(1), point(2), point(3)-1))*0.5/grid.dz;
                end
                grad = [grad; grad_i];
            else error('Wrong dimension for "point" for getTag - dim ~= 3');
            end
        otherwise
            error('Wrong dimension for grid');
    end
end