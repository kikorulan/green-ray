function source = forward_timeSignal(grid, source)
% FORWARDTIMESIGNAL computes the time signal corresponding to the given sensor
%grid = computeTimeSignal(grid, nS, sensor)
%
% INPUTS
% grid: object from gridRT class
% source: source in which to compute the signal
%
% OUTPUTS
% grid: object from gridRT class
%   - source(nS).tSignal: time signal
%   - source(nS).aSignal: amplitude signal
%
% Copyright (C) Kiko Rul·lan, Marta M. Betcke

% Find sound speed at the sensor
cSource = 1/grid.getN(grid.findCoordinates(source.x0));

% Generate Spline
lForward = length(grid.tForward);
% Obtain spline filter
tFilterSpline = 0:grid.dt:grid.tFilter(end)*grid.cFilter/cSource; % time array for the filter
filterSpline = spline(grid.tFilter*grid.cFilter/cSource, grid.Filter, tFilterSpline);%*cSource;  % Q: multiply by cSource?
% Convolve Signal and Filter
signalConv = conv(source.aBeam, filterSpline);
% Adjust delay
delayFilter = find(filterSpline == max(filterSpline)) - 1;
aForwardDelay = signalConv(1+delayFilter:end);
% Adjust total length
if (length(aForwardDelay) < lForward)
    aForward = padarray(aForwardDelay, [0 lForward-length(aForwardDelay)], 0, 'post');
else 
    aForward = aForwardDelay(1:lForward);
end
% Save Signal
source.setForwardSignal(aForward);

