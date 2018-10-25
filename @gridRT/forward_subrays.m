function grid = forward_subrays(grid, source, nRsub, dx)
% FORWARD_SUBRAYS computes the subrays associated to a gaussian beam
%function grid = forward_subrays(grid, source)
% INPUTS
% grid: gridRT object that defines the domain
% source: source to compute the ray trajectories
% nRsub: number of subrays
% dx: separation between subrays
%
% OUTPUTS:
% grid: gridRT object that defines the domain
% 
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke

% Number of steps of the computation
step = source.step;
[nR, nPoints, dim] = size(source.x);

% Allocate rays
source.allocateGauss(nRsub);

% P orthogonal
pT = cat(3, -source.p(:, :, 2)./source.n, source.p(:, :, 1)./source.n);

% Auxiliar values for the GB
w0 = 1e6;

% Sub rays
for j = 1:nRsub
    % Trajectories
    source.insertXGB(source.x + j*pT*dx, 2*j - 1);
    source.insertXGB(source.x - j*pT*dx, 2*j);
    % Amplitudes
    attenuation = abs(exp(-i*w0*source.kIndex./(2*source.qGB)*(j*dx)^2));
    source.insertAmplitudeGB(source.amplitude.*attenuation, 2*j-1);
    source.insertAmplitudeGB(source.amplitude.*attenuation, 2*j);
end

% Pressure computation
for j = 1:2*nRsub
    for tIndex = 1:nPoints-1
        xGB = source.xGB(:, tIndex, :, j);
        coord = grid.findCoordinates(xGB);
        binaryCoord = coord(:, 1, 1) < inf;
        % Insert pressure
        pressureGB = 0*coord(:, :, 1);
        pressureGB(binaryCoord, :, :) = grid.getU(coord(binaryCoord, :, :));
        source.insertPressureGBindex(pressureGB, j, tIndex);
    end
end

% Total pressure computation
%source.amplitudeGB(isinf(source.amplitudeGB) | isnan(source.amplitudeGB)) = 0;
%source.pressureGB(isinf(source.pressureGB) | isnan(source.pressureGB)) = 0;
pressure = sum(source.pressureGB.*source.amplitudeGB, 4);
pressure(isinf(pressure) | isnan(pressure)) = 0;
pressure = pressure + source.pressure;
source.insertPressureVector(pressure);








