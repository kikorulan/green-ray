function h = plot_pressure(grid, timeStamp)
% PLOTPRESSURE plots the pressure at the given timestamp
%function h = plotPressure(grid, timeStamp)
% INPUTS
% grid: gridRT object that defines the domain
% timeStamp: time when to plot the pressure
%
% OUTPUTS:
% h: figure handler
% 
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke

tMax = grid.timeSignal(end);

% Verify validity of timeStamp
if (timeStamp > tMax)
    error('Time stamp is greater than the maximum time');
end

% Find index
index = find(grid.timeSignal > timeStamp, 1);
% Plot pressure
pressure = real(grid.pixelAmplitudeSignal(:, :, index));
h = figure;
surf(grid.xAxis, grid.yAxis, pressure, 'EdgeColor', 'none');
title('Amplitude');

% Plot complex part
pressure = imag(grid.pixelAmplitudeSignal(:, :, index));
binaryPressure = (pressure ~= 0);
if (sum(binaryPressure(:)) > 0)
    disp('Warning: Non-zero complex part');
    figure;
    surf(x, y, pressure, 'EdgeColor', 'none');
    title('Complex part');
end
            
