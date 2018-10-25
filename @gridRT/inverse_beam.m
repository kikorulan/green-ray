function grid = inverse_beam(grid, source)
% INVERSE_BEAM computes the time reverse signal for the given source
%function grid = inverse_beam(grid, nS)
%
% INPUTS
% grid: object from gridRT class
% source: source
% 
% OUTPUTS
% grid: object from gridRT class
%   - pixelABeamFFT: FFT transform of the aBeam signal
%
% Copyright (C) Kiko Rul·lan, Marta M. Betcke

if (isempty(grid.Filter) || isempty(grid.cFilter) || isempty(grid.tFilter))
   error('The filter has not been set.');
end

% Display initialisation message
%msg = strcat({'Computing time propagation for the given source...'});
%disp(msg{1});

% Filter time array
cMax = max(grid.c(:)); % maximum sound speed
cMin = min(grid.c(:)); % minumum sound speed
dt = grid.dt; % time signal increment
nFilters = size(grid.FilterReverse, 1);
deltaC = (cMax - cMin)/(nFilters-1);
tMax = max(source.pixelTime(source.pixelTime < inf)); % maximum time of arrival to pixel
delayPropMax = floor(tMax/dt); % maximum propagation delay
%lengthFilter = length(grid.Filter);
lengthFilter = size(grid.FilterReverse, 2);

% Time reverse
aForward = source.aForward;
aForwardFlip = fliplr(aForward);
% Convolve with filters
%aForwardFlip = grid.source(nS).aReverse;
source.aReverse = [];
for n = 1:nFilters
    %convSignal = conv(aForward, grid.FilterReverse(n, :));
    %convSignal = fliplr(convSignal);
    %filter = 1:lengthFilter;
    %grid.FilterReverse(n, :) = 1./(filter.^2);
    %source.aReverse = [source.aReverse; padarray(aForwardFlip, [0 lengthFilter], 'post')];
    convSignal = conv(aForwardFlip, grid.FilterReverse(n, :));
    source.aReverse = [source.aReverse; convSignal];
end

% Initialisation
source.pixelAReverse = zeros(grid.Nx, grid.Ny);
%source.pixelTReverse = nan(grid.Nx, grid.Ny);

% Find Filter
if (cMax == cMin) 
    pixelIndexFilter = ones(grid.Nx, grid.Ny);
else 
    pixelIndexFilter = 1 + floor((grid.c-cMin)/deltaC);
end
% Find delay
pixelDelay = min(floor(source.pixelTime/dt), delayPropMax);
% Loop over filters
for n = 1:nFilters
    pixelLogicalIndexFilter = (pixelIndexFilter == n);
    source.pixelAReverse = source.pixelAReverse + ...
                                    pixelLogicalIndexFilter.*(reshape(source.aReverse(n, max(1, -pixelDelay' + grid.delayFilterReverse(n) + length(grid.tForward))), grid.Ny, grid.Nx))';
end
source.pixelAReverse = source.pixelAReverse.*source.pixelAttenuation;
%source.pixelAReverse = source.pixelAReverse.*grid.c;
source.aReverse = [];
