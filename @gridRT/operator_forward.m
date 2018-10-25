function forward_data = operator_forward(grid, source, index, initial_pressure, norm_factor)
% OPERATOR_FORWARD computes the forward propagation of the given initial pressure and stores it in the source variable
% forward_data = operator_forward(grid, source, initial_pressure, norm_factor)
%
% INPUTS
%   grid               - gridRT object that defines the domain
%   source             - vector of sources
%   index              - choose between 'all' or a particular index
%   initial_pressure   - initial pressure
%   norm_factor        - normalisation factor
%
% OUTPUTS:
%   forward_data       - forward data from the given initial pressure 
%
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke

% Choose between all sources or a particular index
if (strcmp(index, 'all'))
    selectSource = source;
else 
    selectSource = source(index);
end

% Compute the forward data over all sources
for n = 1:length(selectSource)
    step = selectSource(n).step;
    [nR nPoints dim] = size(selectSource(n).x);
    % Find pressure
    for index = 1:nPoints
        coord = grid.findCoordinates(selectSource(n).x(:, index, :));
        binaryCoord = coord(:, 1, 1) < inf;
        % Insert pressure
        pressure = coord(:, :, 1);
        pressure(binaryCoord, :, :) = norm_factor*initial_pressure(coord(binaryCoord, :, 1) + (coord(binaryCoord, :, 2)-1)*grid.Nx);
        selectSource(n).insertPressure(pressure, index);
    end
    
    % Compute beam
    tBeam = 0:step:(nPoints-1)*step;
    q = abs(selectSource(n).q);
    signal = selectSource(n).pressure.*q.*selectSource(n).revAmplitude;
    signal(isnan(signal) | isinf(signal)) = 0;
    aBeam = sum(signal, 1); 
    selectSource(n).setBeam(tBeam, aBeam);
    
    % Save initial data if necessary
    if(isempty(selectSource(n).aForward_initial))
        selectSource(n).setForwardSignal_initial(selectSource(n).aForward);
    end
    
    % Compute time signal
    grid.forward_timeSignal(selectSource(n));
    % Store forward data
    forward_data(n, :) = selectSource(n).aForward;
end

