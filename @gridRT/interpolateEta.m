function n = interpolateN(grid, x)
% Interpolate the value of Eta at the given point x
% Neighbours are given in the following order:
% - Dim = 1: [X0 X1]
% - Dim = 2: [X0 X1 X0 X1; Y0 Y0 Y1 Y1] 
[nR lengthX dim] = size(x); % Number of rays
n = repmat(inf, [nR 1 1]);
% Find the neighbors
neigh = grid.findNeighbours(x);
% Consider only finite coordinates
fc = (neigh(:, 1, 1) < inf);
switch dim
   case 2
        n11 = grid.getN(neigh(fc, 1, :));
        n21 = grid.getN(neigh(fc, 2, :));
        n12 = grid.getN(neigh(fc, 3, :));
        n22 = grid.getN(neigh(fc, 4, :));
        % Bilinear interpolation
        n(fc, :, :) = n11.*((neigh(fc, 4, 1)-1)*grid.dx - x(fc, :, 1)).*((neigh(fc, 4, 2) - 1)*grid.dy - x(fc, :, 2));
        n(fc, :, :) = n(fc, :, :) - n21.*((neigh(fc, 1, 1)-1)*grid.dx - x(fc, :, 1)).*((neigh(fc, 4, 2) - 1)*grid.dy - x(fc, :, 2));
        n(fc, :, :) = n(fc, :, :) - n12.*((neigh(fc, 4, 1)-1)*grid.dx - x(fc, :, 1)).*((neigh(fc, 1, 2) - 1)*grid.dy - x(fc, :, 2));
        n(fc, :, :) = n(fc, :, :) + n22.*((neigh(fc, 1, 1)-1)*grid.dx - x(fc, :, 1)).*((neigh(fc, 1, 2) - 1)*grid.dy - x(fc, :, 2));
        n(fc, :, :) = n(fc, :, :)./grid.dx/grid.dy;
    case 3
        error('Not implemented');
end
