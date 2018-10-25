function hess = hessPhase(grid, point)
% Compute the hessian of the phase at a given point
    hess = [];
    switch grid.dim
        % Dimension = 2
        case 2
            if (length(point) == 2); 
                % Compute the hessian
                % X component
                if (point(1) == 1)
                    xMin = 1;
                    xMax = 0;
                elseif (point(1) == grid.Nx)
                    xMin = 0;
                    xMax = 1;
                else
                    xMin = 0;
                    xMax = 0;
                end
                % Y component
                if (point(2) == 1)
                    yMin = 1;
                    yMax = 0;
                elseif (point(2) == grid.Ny)
                    yMin = 0;
                    yMax = 1;
                else
                    yMin = 0;
                    yMax = 0;
                end
                pos_x = point(1) + xMin - xMax;
                pos_y = point(2) + yMin - yMax;
                hess_xx = (grid.phase(pos_x + 1, point(2)) - 2*grid.phase(pos_x, point(2)) + grid.phase(pos_x - 1, point(2)))/((grid.dx)^2);
                hess_yy = (grid.phase(point(1), pos_y + 1) - 2*grid.phase(point(1), pos_y) + grid.phase(point(1), pos_y - 1))/((grid.dy)^2);
                hess_xy = (grid.phase(point(1)+1-xMax, point(2)+1-yMax) ...
                            - grid.phase(point(1)-1-xMin, point(2)+1-yMax) ... 
                            - grid.phase(point(1)+1-xMax, point(2)-1-yMin) ...
                            + grid.phase(point(1)-1-xMin, point(2)-1-yMin))/...
                            (grid.dy*grid.dx*2^(xMin + xMax + yMin + yMax));
                hess = [hess_xx hess_xy; hess_xy hess_yy];
            else error('Wrong dimension for "point" for getTag - dim ~= 2');
            end
        case 3
            if (length(point) == 3); 
                % Compphasete the hessian
                % X component
                if (point(1) == 1)
                    xMin = 1;
                    xMax = 0;
                elseif (point(1) == grid.Nx)
                    xMin = 0;
                    xMax = 1;
                else
                    xMin = 0;
                    xMax = 0;
                end
                % Y component
                if (point(2) == 1)
                    yMin = 1;
                    yMax = 0;
                elseif (point(2) == grid.Ny)
                    yMin = 0;
                    yMax = 1;
                else
                    yMin = 0;
                    yMax = 0;
                end
                % Z component
                if (point(3) == 1)
                    zMin = 1;
                    zMax = 0;
                elseif (point(3) == grid.Nz)
                    zMin = 0;
                    zMax = 1;
                else
                    zMin = 0;
                    zMax = 0;
                end
                pos_x = point(1) + xMin - xMax;
                pos_y = point(2) + yMin - yMax;
                pos_z = point(3) + zMin - zMax;
                
                hess_xx = (grid.phase(pos_x + 1, point(2), point(3)) - 2*grid.phase(pos_x, point(2), point(3)) + grid.phase(pos_x - 1, point(2), point(3)))/((grid.dx)^2);
                hess_yy = (grid.phase(point(1), pos_y + 1, point(3)) - 2*grid.phase(point(1), pos_y, point(3)) + grid.phase(point(1), pos_y - 1, point(3)))/((grid.dy)^2);
                hess_zz = (grid.phase(point(1), point(2), pos_z + 1) - 2*grid.phase(point(1), point(2), pos_z) + grid.phase(point(1), point(2), pos_z - 1))/((grid.dy)^2);                
                
                hess_xy = (grid.phase(point(1)+1-xMax, point(2)+1-yMax, point(3)) ...
                            - grid.phase(point(1)-1-xMin, point(2)+1-yMax ,point(3)) ... 
                            - grid.phase(point(1)+1-xMax, point(2)-1-yMin, point(3)) ...
                            + grid.phase(point(1)-1-xMin, point(2)-1-yMin, point(3)))/...
                            (grid.dx*grid.dy*2^(xMin + xMax + yMin + yMax));
                hess_xz = (grid.phase(point(1)+1-xMax, point(2), point(3)+1-zMax) ...
                            - grid.phase(point(1)-1-xMin, point(2), point(3)+1-zMax) ... 
                            - grid.phase(point(1)+1-xMax, point(2), point(3)-1-zMin) ...
                            + grid.phase(point(1)-1-xMin, point(2), point(3)-1-zMin))/...
                            (grid.dx*grid.dz*2^(xMin + xMax + zMin + zMax));
                hess_yz = (grid.phase(point(1), point(2)+1-yMax, point(3)+1-zMax) ...
                            - grid.phase(point(1), point(2)+1-yMax ,point(3)-1-zMin) ... 
                            - grid.phase(point(1), point(2)-1-yMin, point(3)+1-zMax) ...
                            + grid.phase(point(1), point(2)-1-yMin, point(3)-1-zMin))/...
                            (grid.dy*grid.dz*2^(zMin + zMax + yMin + yMax));
                
                hess = [hess_xx hess_xy hess_xz; hess_xy hess_yy hess_yz; hess_xz hess_yz hess_zz];
            else error('Wrong dimension for "point" for getTag - dim ~= 3');
            end
        otherwise
            error('Wrong dimension for grid');
    end
end
