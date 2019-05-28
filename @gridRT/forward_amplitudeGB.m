function A = forward_amplitudeGB(grid, source)
% FORWARD_AMPLITUDEGB Computes the amplitude at a given point x using the Dynamic Ray Tracing method
% from source point x0.
%function A = forward_amplitudeGB(grid, source)
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
if (isempty(source.freq) || isempty(source.beamwidth))
    error('Please set frequency and beamwidth for GB');
end
w0 = 2*pi*source.freq;
q0 = j*w0*source.beamwidth*source.beamwidth*0.5;
q = repmat(q0, [nR 1 1]);
u = repmat(1, [nR 1 1]);

source.insertQ(q, 1);
source.insertKIndex(u, 1);
%==============================
% Gaussian Beam 
%==============================
for ii = 1:nPoints-1 % Last step not compatible with RK2
    % Update Dx, Dp    
    [q, u] = grid.stepRK2GB(source, step, ii, q, u);
    source.insertQ(q, ii+1);
    source.insertKIndex(u, ii+1);
end
% Compute amplitude for all rays
%A = sqrt(1./source.q./source.n);
qGB = source.q;
%source.insertQVector(sqrt(-2./imag(source(1).kIndex./qGB)));

%%  %==============================
%%  % Q - Proximal Ray
%%  %==============================
%%  x = source.x;
%%  % For computing the determinant
%%  difPhi1 = x(1:end-1, 1:end-1, :) - source.x(2:end, 1:end-1, :);
%%  difPhi = cat(1, difPhi1, x(end-1, 1:end-1, :) - x(end, 1:end-1, :));
%%  difTau = x(:, 2:end, :) - x(:, 1:end-1, :);
%%  flipDifTau(:, :, 2) = -difTau(:, :, 1);
%%  flipDifTau(:, :, 1) = difTau(:, :, 2); 
%%  % Q
%%  q = sum(difPhi.*flipDifTau, 3);
%%  q = [q repmat(NaN, [nR 1 1])];
%%  source.insertQVector(q);
%%  %source.insertDXVector(q);

%%  %==============================
%%  % Q - ODE method
%%  %==============================
%%  Dx = repmat(eye(grid.dim), [1 1 nR]);
%%  Dp = grid.hessPhaseODE(source, 'real');
%%  q = [];
%%  
%%  DxReshape = permute(reshape(Dx, dim*dim, 1, []), [3 2 1]);
%%  DpReshape = permute(reshape(Dp, dim*dim, 1, []), [3 2 1]);
%%  source.insertDX(DxReshape, 1);
%%  source.insertDP(DpReshape, 1);
%%  q = repmat(1, [nR 1 1]);
%%  %source.insertQ(repmat(1, [nR 1 1]), 1);
%%  for j = 1:nPoints-1 % Last step not compatible with RK2
%%      % Update Dx, Dp
%%      [Dx Dp] = grid.stepRK4Amplitude(source, step, j, Dx, Dp);
%%      DxReshape = permute(reshape(Dx, dim*dim, 1, []), [3 2 1]);
%%      DpReshape = permute(reshape(Dp, dim*dim, 1, []), [3 2 1]);
%%      source.insertDX(DxReshape, j+1);
%%      source.insertDP(DpReshape, j+1);
%%      q = [q permute(det3D(Dx), [3 2 1])];
%%      %source.insertQ(permute(det3D(Dx), [3 2 1]), j+1);
%%      % Keller-Maslov index
%%      %source.insertKIndex(source.kIndex(:, j) + (source.q(:, j).*source.q(:, j+1) < 0), j+1);
%%  end
%%  source.insertDXVector(q);

%==============================
% Amplitude
%==============================
A = sqrt(q0./qGB./source.n*source.n(1, 1));
%A = sqrt(1./qGB./source.xTD./source.n);
beamwidth = sqrt(-1./imag(source.kIndex./qGB));
%A = beamwidth.*sqrt(1./qGB./source.n);
A(isnan(A)) = 0; % substitute NaN values
source.insertAVector(A);

