function A = forward_amplitudeGauss(grid, source)
% FORWARD_AMPLITUDEODE Computes the amplitude at a given point x using the ODE method
% from source point x0.
%function A = forward_amplitudeODE(grid, source)
% INPUTS
% grid: gridRT object that defines the domain
% source: source to compute the amplitude
%
% OUTPUTS:
% A: amplitude
% 
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke

source.deleteAmplitude(); % Delete amplitude in case it has been previously computed
step = source.step; % step size
[nR nPoints dim] = size(source.x); % Number of rays
% Initial conditions
Dx = repmat(eye(grid.dim), [1 1 nR]);
Dp = grid.hessPhaseODE(source, 'imag');
D2phase_permute = zeros(dim, dim, nR);
A = repmat(1, [nR 1 1]);

q = [];
%==============================
% Build ray
%==============================
Dx_reshape = permute(reshape(Dx, dim*dim, 1, []), [3 2 1]);
Dp_reshape = permute(reshape(Dp, dim*dim, 1, []), [3 2 1]);
source.insertDX(Dx_reshape, 1);
source.insertDP(Dp_reshape, 1);
source.insertD2Phase(Dp_reshape, 1);
source.insertQ(repmat(1, [nR 1 1]), 1);
source.insertA(A, 1);
for j = 1:nPoints-1 % Last step not compatible with RK2
    % Finite coordinates
    X = source.x(:, j, :);
    coord = grid.findCoordinates(X);
    binaryCoord = coord(:, 1, 1) < inf;
    c = inf(nR, 1);
    c(binaryCoord) = grid.getC(coord(binaryCoord, :, :));

    % Update Dx, Dp
    [Dx Dp] = grid.stepRK4Amplitude(source, step, j, Dx, Dp);
    %[Dx Dp] = grid.stepRK2Gauss(source, step, j, Dx, Dp);

    % Compute D2phase
    Dx_permute = permute(Dx, [2 1 3]);
    Dp_permute = permute(Dp, [2 1 3]);
    finiteCoord = ~isnan(Dx(1, 1, :));
    D2phase_permute(:) = 0;
    D2phase_permute(:, :, finiteCoord) = SliceMultiSolver(Dx_permute(:, :, finiteCoord), Dp_permute(:, :, finiteCoord));
    D2phase = permute(D2phase_permute, [2 1 3]);
    % Reshape to store
    Dx_reshape = permute(reshape(Dx, dim*dim, 1, []), [3 2 1]);
    Dp_reshape = permute(reshape(Dp, dim*dim, 1, []), [3 2 1]);
    D2phase_reshape = permute(reshape(D2phase, dim*dim, 1, []), [3 2 1]);
    source.insertDX(Dx_reshape, j+1);
    source.insertDP(Dp_reshape, j+1);
    source.insertD2Phase(D2phase_reshape, j+1);
    source.insertQ(permute(det3D(Dx), [3 2 1]), j+1);
    % Derive 
    D2phase_prev = source.D2phase(:, j, :);
    DA = -step*c.*c.*(D2phase_reshape(:, :, 1) + D2phase_reshape(:, :, 4))/2;
    DA = (DA -step*c.*c.*(D2phase_prev(:, :, 1) + D2phase_prev(:, :, 4))/2)/2;
    %DA = -step*c.*c.*(D2phase_prev(:, :, 1) + D2phase_prev(:, :, 4))/2;
    % Compute amplitude
    A = A.*exp(DA);
    %A = grid.stepRK2Amplitude(nS, j, step);
    source.insertA(A, j+1); 
end
% Compute amplitude for all rays
A = sqrt(abs(1./source.q));
A(isnan(A)) = 0; % substitute NaN values
%source.insertAVector(A);


