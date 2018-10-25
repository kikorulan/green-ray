function A = forward_amplitudeODE(grid, source)
% AMPLITUDEODE Computes the amplitude at a given point x using the ODE method
% from source point x0.
%function A = amplitudeODE(grid, nS)
% INPUTS
% grid: gridRT object that defines the domain
% source: source to compute the amplitude
%
% OUTPUTS:
% A: amplitude
% 
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke

% Display message
[nR stepsRay dim] = size(source.x); % Number of rays
% Initial conditions
Dx = repmat(eye(grid.dim), [1 1 nR]);
Dp = grid.hessPhaseODE(source, 'real');

%==============================
% Build ray
%==============================
DxReshape = permute(reshape(Dx, dim*dim, 1, []), [3 2 1]);
DpReshape = permute(reshape(Dp, dim*dim, 1, []), [3 2 1]);
source.insertDX(DxReshape, 1);
source.insertDP(DpReshape, 1);
source.insertQ(repmat(1, [nR 1 1]), 1);
step = source.step;
nPoints = source.nPoints;
for j = 1:nPoints-1 % Last step not compatible with RK2
    % Update Dx, Dp
    [Dx Dp] = grid.stepRK4Amplitude(source, step, j, Dx, Dp);
    DxReshape = permute(reshape(Dx, dim*dim, 1, []), [3 2 1]);
    DpReshape = permute(reshape(Dp, dim*dim, 1, []), [3 2 1]);
    source.insertDX(DxReshape, j+1);
    source.insertDP(DpReshape, j+1);
    source.insertQ(permute(det3D(Dx), [3 2 1]), j+1);
    % Keller-Maslov index
    source.insertKIndex(source.kIndex(:, j) + (source.q(:, j).*source.q(:, j+1) < 0), j+1);
end
% Compute amplitude for all rays
q0 = repmat(source.q(:, 1), [1 stepsRay]);
A = 1./source.n.*sqrt(abs(q0./source.q));%.*exp(-i*source.kIndex*pi/2);
A(isnan(A)) = 0; % substitute NaN values
source.insertAVector(A);


