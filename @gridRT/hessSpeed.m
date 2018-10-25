function hess = hessSpeed(grid, point)
% Compute the hessian of the speed at a given point
    hess = [];
    switch grid.dim
        % Dimension = 2
        case 2
            % X component
            xMin = (point(:, :, 1) == 1);
            xMax = (point(:, :, 1) == grid.Nx);
            pos_x = point(:, :, 1) + xMin - xMax;
            % Y component
            yMin = (point(:, :, 2) == 1);
            yMax = (point(:, :, 2) == grid.Ny);
            pos_y = point(:, :, 2) + yMin - yMax;
            % Compute the Hessian
            hess_xx = (grid.getC(cat(3, pos_x + 1, point(:, :, 2))) - 2*grid.getC(cat(3, pos_x, point(:, :, 2))) + grid.getC(cat(3, pos_x - 1, point(:, :, 2))))/(grid.dx^2);
            hess_yy = (grid.getC(cat(3, point(:, :, 1), pos_y + 1)) - 2*grid.getC(cat(3, point(:, :, 1), pos_y)) + grid.getC(cat(3, point(:, :, 1), pos_y - 1)))/((grid.dy)^2);
            hess_xy = (grid.getC(cat(3, pos_x + 1, pos_y + 1)) ...
                        - grid.getC(cat(3, pos_x - 1, pos_y + 1)) ... 
                        - grid.getC(cat(3, pos_x + 1, pos_y - 1)) ...
                        + grid.getC(cat(3, pos_x - 1, pos_y - 1)))/...
                        (grid.dy*grid.dx*4);
            hess = cat(3, [hess_xx hess_xy], [hess_xy hess_yy]);
        case 3 %%%%%%%%%% NOT FULLY IMPLEMENTED!!! %%%%%%%%%%%%%%%%%%%%%%%
            if (length(point) == 3); 
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
                
                hess_xx = (grid.c(pos_x + 1, point(2), point(3)) - 2*grid.c(pos_x, point(2), point(3)) + grid.c(pos_x - 1, point(2), point(3)))/((grid.dx)^2);
                hess_yy = (grid.c(point(1), pos_y + 1, point(3)) - 2*grid.c(point(1), pos_y, point(3)) + grid.c(point(1), pos_y - 1, point(3)))/((grid.dy)^2);
                hess_zz = (grid.c(point(1), point(2), pos_z + 1) - 2*grid.c(point(1), point(2), pos_z) + grid.c(point(1), point(2), pos_z - 1))/((grid.dy)^2);                
                
                hess_xy = (grid.c(point(1)+1-xMax, point(2)+1-yMax, point(3)) ...
                            - grid.c(point(1)-1-xMin, point(2)+1-yMax ,point(3)) ... 
                            - grid.c(point(1)+1-xMax, point(2)-1-yMin, point(3)) ...
                            + grid.c(point(1)-1-xMin, point(2)-1-yMin, point(3)))/...
                            (grid.dx*grid.dy*2^(xMin + xMax + yMin + yMax));
                hess_xz = (grid.c(point(1)+1-xMax, point(2), point(3)+1-zMax) ...
                            - grid.c(point(1)-1-xMin, point(2), point(3)+1-zMax) ... 
                            - grid.c(point(1)+1-xMax, point(2), point(3)-1-zMin) ...
                            + grid.c(point(1)-1-xMin, point(2), point(3)-1-zMin))/...
                            (grid.dx*grid.dz*2^(xMin + xMax + zMin + zMax));
                hess_yz = (grid.c(point(1), point(2)+1-yMax, point(3)+1-zMax) ...
                            - grid.c(point(1), point(2)+1-yMax ,point(3)-1-zMin) ... 
                            - grid.c(point(1), point(2)-1-yMin, point(3)+1-zMax) ...
                            + grid.c(point(1), point(2)-1-yMin, point(3)-1-zMin))/...
                            (grid.dy*grid.dz*2^(zMin + zMax + yMin + yMax));
                
                hess = [hess_xx hess_xy hess_xz; hess_xy hess_yy hess_yz; hess_xz hess_yz hess_zz];
            else error('Wrong dimension for "point" for getTag - dim ~= 3');
            end
        otherwise
            error('Wrong dimension for grid');
    end
end
