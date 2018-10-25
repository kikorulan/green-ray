function [DX, DP] = deriveHamil(grid, X, P)
% DERIVEHAMIL computes the derivative of the hamiltonian for the given source
% The derivative is computed at the end of the given rays.
%function DH = deriveHamil(grid, X, P)
% INPUTS
% grid: gridRT object that defines the domain
% X: x point to compute the derivative
% P: p momentum to compute the derivative
%
% OUTPUTS:
% DH: value of the computed derivative
% 
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke

coord = grid.findCoordinates(X);
DX = coord;
DP = coord;

binaryCoord = coord(:, 1, 1) < inf;
% Stored values
[n, Gn] = grid.getNGN(coord(binaryCoord, :, :));
% Use interpolation
%n  = repmat(grid.interpolateEta(X(binaryCoord, :, :)), [1 1 2]);
%Gn = grid.interpolateGradEta(X(binaryCoord, :, :));

%Gn = grid.gradEta(coord(binaryCoord, :, :));
%n = grid.getN(coord(binaryCoord, :, :));

% Compute DX
DX(binaryCoord, :, :) = P(binaryCoord, :, :)./(n.*n);
% Compute DP
DP(binaryCoord, :, :) = Gn./n;

