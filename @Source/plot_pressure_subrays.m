function [h1, h2] = plot_pressure_subrays(source, indexRay)
% PLOT_AMPLITUDE_SUBRAYS Plots the Amplitude values for the given sources
% h = plot_amplitude_subrays(source, indexRay)
% INPUTS
% source: source to plot
% indexRay: choose ray index to plot subrays
%
% OUTPUTS:
% h: figure handler
% 
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke

[nR nPoints dim nRsub] = size(source.xGB);
colourList = cool(nRsub + 1);

if (nR > 0)
%==============================
% Rays
%==============================
h1 = figure;
pressure = source.pressure(indexRay, :, :);
plot(source.phi(indexRay, :), pressure, 'Color', colourList(1, :));
hold on;
for j = 1:nRsub
    pressure = source.pressureGB(indexRay, :, :, j);
    plot(source.phi(indexRay, :), pressure, 'Color', colourList(j+1, :));
end
title('Pressure');
xlabel('t (s)');
ylabel('Pressure');
grid on;

%==============================
% Imagesc
%==============================
pressure = [];
% Reorder pressure
for j = 1:nRsub
    if (mod(j, 2) == 1)
        pressure = [source.pressureGB(indexRay, :, :, j); pressure];
    else
        pressure = [pressure; source.pressureGB(indexRay, :, :, j)];
    end
end

h2 = figure;
imagesc(pressure);
title(strcat('Pressure for ray ', int2str(indexRay)));
xlabel('t (s)');
ylabel('Subray number');
colorbar();
end
