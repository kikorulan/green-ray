function h1 = plot_amplitude(source, nR, nColours, plotImag)
% PLOT_AMPLITUDE Plots the Amplitude values for the given sources
% h = plot_amplitude(source, nColours, plotImag)
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
for j = 1:nR
    index = 1 + (j-1)*floor((nRays-1)/(nR-1));
    colourNum = floor(nColours*(j-1)/nR) + 1;
    A = real(source.amplitude(index, :));
    semilogy(source.phi(:), A, 'Color', colourList(colourNum, :));
    hold on;
end
title('Amplitude - real');
xlabel('t (s)');
ylabel('Amplitude');
grid on;

%==============================
% Imaginary part
%==============================
if(plotImag == true)
    h2 = figure;
    for j = 1:nRays
        colourNum = floor(nColours*(j-1)/nRays) + 1;
        A = imag(source.amplitude(j, :));
        semilogy(source.phi(:), A, 'Color', colourList(colourNum, :));
        hold on;
    end
    title('Amplitude - imag');
    xlabel('t (s)');
    ylabel('Amplitude');
    grid on;
end
