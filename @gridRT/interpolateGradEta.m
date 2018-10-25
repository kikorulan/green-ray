function Gn = interpolateGradEta(grid, x)
% This function interpolates the gradient of eta at the given point x
% Neighbours are given in the following order:
% - Dim = 1: [X0 X1]
% - Dim = 2: [X0 X1 X0 X1; Y0 Y0 Y1 Y1] 
[nR lengthX dim] = size(x);
Gn = repmat(inf, [nR 1 dim]);
% Find the neighbors
neigh = grid.findNeighbours(x);
% Consider only the finite coordinates
fc = (neigh(:, 1, 1) < inf);
switch dim
    case 2
        n11 = grid.getN(neigh(fc, 1, :));
        n21 = grid.getN(neigh(fc, 2, :));
        n12 = grid.getN(neigh(fc, 3, :));
        n22 = grid.getN(neigh(fc, 4, :));
        % Derivative with respect to y
        GnY = - n11.*((neigh(fc, 4, 1)-1)*grid.dx - x(fc, :, 1));
        GnY = GnY + n21.*((neigh(fc, 1, 1)-1)*grid.dx - x(fc, :, 1));
        GnY = GnY + n12.*((neigh(fc, 4, 1)-1)*grid.dx - x(fc, :, 1));
        GnY = GnY - n22.*((neigh(fc, 1, 1)-1)*grid.dx - x(fc, :, 1));
        Gn(fc, :, 2) = GnY./grid.dx/grid.dy;
        % Derivative with respect to x
        GnX = - n11.*((neigh(fc, 4, 2) - 1)*grid.dy - x(fc, :, 2));
        GnX = GnX + n21.*((neigh(fc, 4, 2) - 1)*grid.dy - x(fc, :, 2));
        GnX = GnX + n12.*((neigh(fc, 1, 2) - 1)*grid.dy - x(fc, :, 2));
        GnX = GnX - n22.*((neigh(fc, 1, 2) - 1)*grid.dy - x(fc, :, 2));
        Gn(fc, :, 1) = GnX./grid.dx/grid.dy;
    case 3
        error('Not implemented');
end
    
