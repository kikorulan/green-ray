function h1 = plot_pressure(source)
% PLOT_AMPLITUDE Plots the Amplitude values for the given sources
% h = plot_amplitude(source, nColours, plotImag)
% INPUTS
% source: source to plot
%
% OUTPUTS:
% h: figure handler
% 
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke

nRays = size(source.x, 1);
colourList = cool(nRays);

%==============================
% Rays
%==============================
h1 = figure;
for j = 1:nRays
    pressure = source.pressure(j, :);
    plot(source.phi(j, :), pressure, 'Color', colourList(j, :));
    hold on;
end
title('Pressure');
xlabel('t (s)');
ylabel('Pressure');
grid on;

%%  %==============================
%%  % Imagesc
%%  %==============================
%%  h2 = figure;
%%  imagesc(source.pressure);
%%  title('Pressure');
%%  xlabel('t (s)');
%%  ylabel('Pressure');
%%  grid on;
