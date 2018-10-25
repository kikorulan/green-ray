function [Dq, Du] = deriveGB(grid, source, index, q, u)
% DERIVEGB Computes the ODE system for q
% It solves one step of the ODE system for the amplitude.
%function [D2x D2p] = deriveGauss(grid, source, index, Dx, Dp)
% INPUTS
% grid: gridRT object that defines the domain
% source: source to compute the derivative
% index: index in the trajectory to compute the derivative
% q: determinant of the jacobian
% u: pressure
%
% OUTPUTS:
% Dq: iterate for Jacobian of x
% Du: iterate for Jacobian of p
% 
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke

%%%% Point to derive
x = source.x(:, index, :);
p = source.p(:, index, :);
[nR dimY dim] = size(x);
point = grid.findCoordinates(x);
binaryPoint = point(:, 1, 1) < inf;
%%%% Update matrices
c = nan(nR, 1);
c(binaryPoint) = grid.getC(point(binaryPoint, :, :));
% Momentum
nVector = cat(3, -p(binaryPoint, :, 2), p(binaryPoint, :, 1));
NMatrix = repmat(nVector, [1 dim 1]);
NTMatrix = permute(NMatrix, [1 3 2]);
matrix = zeros(nR, dim, dim);
matrix(binaryPoint, :, :) = NTMatrix.*grid.hessSpeed(point(binaryPoint, :, :)).*NMatrix; 

% Derive
%Dq = c.*u;
%Du = -sum(sum(matrix, 2), 3).*q./c./c;
Dq = c.*c.*u;%c.*u;
Du = -sum(sum(matrix, 2), 3).*q./c;%-sum(sum(matrix, 2), 3).*q./c./c;

