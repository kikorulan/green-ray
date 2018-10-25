function grid = forward_beam(grid, source)
% FORWARD_BEAM Computes the beam for the given source
% grid = forward_beam(grid, source)
% INPUTS
% grid: gridRT object that defines the domain
% source: source to compute the beam
%
% OUTPUTS:
% grid: gridRT object that defines the domain
% 
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke

% Number of steps of the computation
step = source.step;
nPoints = source.nPoints;
nR = size(source.x, 1);

% Compute time
tBeam = 0:step:(nPoints-1)*step;

% Compute distances
rayI = source.x(2:nR, :, :);
rayF = source.x(1:nR-1, :, :);
normRay = sqrt(sum((rayI - rayF).^2, 3));
normRay(isinf(normRay) | isnan(normRay)) = 0;
dBeam0 = zeros(1, nPoints); % zero distance for first and last rays
dBeam = 0.5*(cat(1, dBeam0, normRay) + cat(1, normRay, dBeam0));
% Compute amplitude
%q = abs(source.q);
%q = source.q.*((source.q > 0) - i*(source.q < 0));
q = source.q;
%q = sqrt(source.q.*source.qGB);
q(isnan(q) | isinf(q)) = 0;
aBeam = sum(source.pressure.*q, 1); % Density correction
%aBeam = sum(source.pressure, 1);
beamwidth = sqrt(-1./imag(source.kIndex./source.q));
curvature = 1./source.n.*real(source.kIndex./source.q);
density = source.Dx;
density(isnan(density) | isinf(density)) = 0;
%aBeam = sum(source.pressure.*density, 1); % Density correction
%aBeam = abs(aBeam);
%aBeam = sum(source.pressure.*dBeam, 1); % Distance correction
%aBeam = sum(source.pressure, 1); % No distance correction

%%  aBeam = zeros(1, nPoints);
%%  for n = 1:nPoints-1
%%      condition = (source.pixelTime >= tBeam(n)) & (source.pixelTime < tBeam(n+1));
%%      aBeam(n) = sum(source.pixelPressure(condition));
%%  end
%%  aBeam(nPoints) = sum(source.pixelPressure(source.pixelTime >= tBeam(nPoints)));
% Insert values
source.setBeam(tBeam, aBeam);
