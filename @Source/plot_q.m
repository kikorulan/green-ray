function h1 = plot_q(source, nColours, plotImag)
% PLOT_Q Plots the Amplitude values for the given sources
% h = plot_q(source, nColours, plotImag)
% INPUTS
% source: source to plot
% nColours: number of colours to plot
% type: plot real and imaginary parts
%
% OUTPUTS:
% h: figure handler
% 
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke

colourList = cool(nColours);
nRays = size(source.x, 1);

%==============================
% Real part
%==============================
h1 = figure;
for j = 1:nRays
    colourNum = floor(nColours*(j-1)/nRays) + 1;
    q = real(source.q(j, :));
    plot(source.phi(j, :), q, 'Color', colourList(colourNum, :));
    hold on;
end
title('q - real');
xlabel('t (s)');
ylabel('q');
grid on;

%==============================
% Imaginary part
%==============================
if(plotImag == true)
    h2 = figure;
    for j = 1:nRays
        colourNum = floor(nColours*(j-1)/nRays) + 1;
        q = imag(source.q(j, :));
        plot(source.phi(j, :), q, 'Color', colourList(colourNum, :));
        hold on;
    end
    title('q - imag');
    xlabel('t (s)');
    ylabel('q');
    grid on;
end

