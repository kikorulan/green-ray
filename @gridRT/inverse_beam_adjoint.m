function aReverse = inverse_beam_adjoint(grid, source)
% INVERSE_BEAM_ADJOINT computes the reverse signal considering the effect of other sensors
%aReverse = inverse_beam_adjoint(grid, source)
%
% INPUTS
% grid: object from gridRT class
% source: vector of sources 
%
% OUTPUTS
% aReverse: computed signal
%
% Copyright (C) Kiko Rul·lan, Marta M. Betcke


% Filter time array
dt = grid.dt; % time signal increment

% Create variables
lengthForward = length(grid.tForward);
nSources = length(source);
aForward_flip = nan(nSources, lengthForward);
aForward_int = zeros(nSources, lengthForward);
% Find values for all sources
for n = 1:nSources
    % Find coordinates
    coord = grid.findCoordinates(source(n).x0);
    % Find attenuation
    attenuation(n) = source(n).pixelAttenuation(coord(:, :, 1), coord(:, :, 2));
    delay(n) = floor(source(n).pixelTime(coord(:, :, 1), coord(:, :, 2))/dt); 
    aForward_flip(n, :) = fliplr(source(n).aForward);
end
% Integrate signals
aForward_int(:, 1) = aForward_flip(:, 1);
for n = 2:lengthForward
    aForward_int(:, n) = aForward_flip(:, n) + aForward_int(:, n-1);
end


% Loop over all sources
aReverse = zeros(1, lengthForward);
for n = 1:nSources
    disp(n);
    aReverse(:) = 0;
    for j = 1:nSources
       % Delay signal
        aForward_delay = attenuation(j)*padarray(aForward_int(j, 1:end-delay(j)), [0 delay(j)], 'pre');
        %aForward_delay = attenuation(j)*padarray(aForward_flip(j, 1:end-delay(j)), [0 delay(j)], 'pre');
        % Sum signal
        if (j == n)
            aReverse = aReverse + aForward_delay;
        else
            aReverse = aReverse - aForward_delay;
        end
    end
    source(n).aReverse = aReverse;
end

for n = 1:nSources
    source(n).aReverse = aForward_int(n, :);
    %source(n).aReverse = conv(aForward_flip, [10 1 1 1 1 1 1 1 1], 'same');
end

% Assign output
aReverse = zeros(nSources, lengthForward);
for n = 1:nSources
    aReverse(n, :) = source(n).aReverse;
end
