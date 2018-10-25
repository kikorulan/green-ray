function source = inverse_timeSource(grid, source)
% inverse_timeSource computes the time propagation from the given source
%source = inverse_timeSource(grid, source)
%
% INPUTS
% grid: object from gridRT class
% source: source
% 
% OUTPUTS
% source:
%    - aPropagation: signal to propagate
%
% Copyright (C) Kiko Rul·lan, Marta M. Betcke

if (isempty(grid.source(nS).pixelAttenuation))
    error('Rays not converted to grid.');
end

if (isempty(grid.Filter) || isempty(grid.cFilter) || isempty(grid.tFilter))
   error('The filter has not been set.');
end

% Display initialisation message
msg = strcat({'Computing time propagation for source '}, int2str(nS), {'...'});
disp(msg{1});

% Filter time array
cMax = max(grid.c(:)); % maximum sound speed
cMin = min(grid.c(:)); % minumum sound speed
dt = grid.dt; % time signal increment
nFilters = size(grid.FilterReverse, 1);
deltaC = (cMax - cMin)/(nFilters-1);
tMax = max(source.pixelTime(source.pixelTime < inf)); % maximum time of arrival to pixel
delayPropMax = floor(tMax/dt); % maximum propagation delay

% Time reverse
%aForwardFlip = fliplr(grid.source(nS).aForward);
aForward = source.aForward;
lengthAForward = length(aForward);
% Convolve with filters
source.aReverse = [];
for n = 1:nFilters
    source.aReverse = [source.aReverse; ...
                          zeros(1, delayPropMax) conv(aForward, grid.FilterReverse(n, :))];
end

lengthAReverse = size(source.aReverse, 2);
% Initialisation
lengthPropagation = lengthAForward + delayPropMax;
pixelDelay = min(floor(source.pixelTime/dt), delayPropMax); % Find delay
pixelAttenuation = repmat(source.pixelAttenuation, [1 1 lengthPropagation]);
%pixelAttenuation(1:120, :, :) = 0;
source.pixelAPropagation = nan(grid.Nx, grid.Ny, lengthPropagation);
switch grid.dim
    case 2
        for i = grid.Nx:-1:1
            disp(i);
            for j = grid.Ny:-1:1
                c = grid.c(i, j);
                if (c == cMin) indexFilter = 1;
                else indexFilter = 1 + floor((c-cMin)/deltaC);
                end
                source.pixelAPropagation(i, j, :) = source.aReverse(indexFilter, delayPropMax-pixelDelay(i, j)+1:lengthPropagation+delayPropMax-pixelDelay(i, j));
            end
        end
        source.pixelAPropagation = source.pixelAPropagation.*pixelAttenuation;
    case 3
        error('Not implemented');
end
  
% Correct NaN
source.pixelAPropagation(isnan(source.pixelAPropagation(:))) = 0;

