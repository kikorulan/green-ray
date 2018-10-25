function pixelAReverse = iterative_recon_sensor(grid, source)
% ITERATIVE_RECON_SENSOR computes an iteration of the reconstruction feeding back the data
% source = iterative_recon_sensor(grid, source)
% INPUTS
% grid: gridRT object that defines the domain
% source: source to compute the iteration
%
% OUTPUTS:
% source: vector containing the computed sources
% 
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke

step = source.step;
[nR nPoints dim] = size(source.x);
% Find pressure
for index = 1:nPoints
    coord = grid.findCoordinates(source.x(:, index, :));
    binaryCoord = coord(:, 1, 1) < inf;
    % Insert pressure
    pressure = coord(:, :, 1);
    pressure(binaryCoord, :, :) = grid.pixelAReverse(coord(binaryCoord, :, 1) + (coord(binaryCoord, :, 2)-1)*grid.Nx);
    source.insertPressure(pressure, index);
end

% Compute beam
tBeam = 0:step:(nPoints-1)*step;
q = abs(source.q);
signal = source.pressure.*q.*source.revAmplitude;
signal(isnan(signal) | isinf(signal)) = 0;
aBeam = sum(signal, 1); 
source.setBeam(tBeam, aBeam);

% Save initial data if necessary
if(isempty(source.aForward_initial))
    source.setForwardSignal_initial(source.aForward);
end

% Compute time signal
grid.forward_timeSignal(source);

% Compute difference
factorForward = max(source.aForward);
factorInitial = max(source.aForward_initial);
difference = source.aForward_initial - source.aForward/factorForward*factorInitial;
source.setForwardSignal(difference);

% Update pixelAReverse
%grid.pixelAReverse = grid.pixelAReverse - source.pixelAReverse;
grid.inverse_beam(source);
grid.pixelAReverse = grid.pixelAReverse + source.pixelAReverse;
pixelAReverse = source.pixelAReverse;
