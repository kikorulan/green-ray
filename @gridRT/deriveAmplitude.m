function [D2x D2p] = deriveAmplitude(grid, source, index, Dx, Dp)
% DERIVEAMPLITUDE Computes the Jacobians Dx, Dp for the given index in the trajectory
% It solves one step of the ODE system for the amplitude.
%function [D2x D2p] = deriveAmplitude(grid, nS, tauIndex, Dx, Dp)
% INPUTS
% grid: gridRT object that defines the domain
% source: source to compute the derivative
% index: index in the trajectory to compute the derivative
% Dx: Jacobian of x
% Dp: Jacobian of p
%
% OUTPUTS:
% D2x: iterate for Jacobian of x
% D2p: iterate for Jacobian of p
% 
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke

%%%% Point to derive
x = source.x(:, index, :);
[nR dimY dim] = size(x);
point = grid.findCoordinates(x);
binaryPoint = point(:, 1, 1) < inf;
%%%% Update matrices
% Gradient
gradC = double(point) + 0.1;
gradC(binaryPoint, :, :) = grid.gradSpeed(point(binaryPoint, :, :));
gradCT = permute(gradC, [1 3 2]);
CMatrix = repmat(grid.getC(point(binaryPoint, :, :)), [1 dim dim]);
gradCMatrix = repmat(gradC(binaryPoint, :, :), [1 dim 1]);
gradCTMatrix = repmat(gradCT(binaryPoint, :, :), [1 1 dim]);
% Direction
p = double(point) + 0.1;
p(binaryPoint, :, :) = source.p(binaryPoint, index, :);
pT = permute(p, [1 3 2]);
pMatrix = repmat(p(binaryPoint, :, :), [1 dim 1]);
pTMatrix = repmat(pT(binaryPoint, :, :), [1 1 dim]);
% Identity
id = repmat(eye(dim), [1 1 nR]);
id_permute = permute(id, [3 2 1]);

% DxxH matrix
DxxH = repmat(inf(dim), [1 1 nR]);
DxxH_permute = grid.hessSpeed(point(binaryPoint, :, :))./CMatrix;
DxxH(:, :, binaryPoint)  = permute(DxxH_permute, [3 2 1]);

% DpxH matrix
DpxH = repmat(inf(dim), [1 1 nR]);
%DpxH_permute = CMatrix.*gradCMatrix.*pTMatrix;
DpxH_permute = CMatrix.*gradCTMatrix.*pMatrix;
DpxH(:, :, binaryPoint)  = permute(DpxH_permute, [3 2 1]);
DpxHT = permute(DpxH, [2 1 3]);

% DppH matrix
DppH = repmat(inf(dim), [1 1 nR]);
DppH_permute = CMatrix.*CMatrix.*id_permute(binaryPoint, :, :) - (CMatrix.^3).*pMatrix.*pTMatrix;
DppH(:, :, binaryPoint)  = permute(DppH_permute, [3 2 1]);

% Derive
D2x = prod3D( DpxH, Dx) + prod3D(DppH, Dp);
D2p = prod3D(-DxxH, Dx) + prod3D(-DpxHT, Dp);

end
% M = [[DpxH] [DppH]; [-DxxH] [-DpxH]];
