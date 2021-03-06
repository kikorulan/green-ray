function [h1, h2] = plot_amplitude_subrays(source, indexRay, plotImag)
% PLOT_AMPLITUDE_SUBRAYS Plots the Amplitude values for the given sources
% h = plot_amplitude_subrays(source, nColours, plotImag)
% INPUTS
% source: source to plot
% indexRay: ray index to plot
% type: plot real and imaginary parts
%
% OUTPUTS:
% h: figure handler
% 
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke

[nR nPoints dim nRsub] = size(source.xGB);
colourList = cool(nRsub + 1);

if (nR > 0)
%==============================
% Real part
%==============================
h1 = figure;
A = real(source.amplitude(indexRay, :));
semilogy(source.phi(indexRay, :), A, 'Color', colourList(1, :));
hold on;
for j = 1:nRsub
    A = real(source.amplitudeGB(indexRay, :, :, j));
    semilogy(source.phi(indexRay, :), A, 'Color', colourList(j+1, :));
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
    A = imag(source.amplitude(indexRay, :));
    semilogy(source.phi(indexRay, :), A, 'Color', colourList(1, :));
    hold on;
    for j = 1:nRsub
        A = imag(source.amplitudeGB(indexRay, :, :, j));
        semilogy(source.phi(indexRay, :), A, 'Color', colourList(j+1, :));
    end
    title('Amplitude - imag');
    xlabel('t (s)');
    ylabel('Amplitude');
    grid on;
end

end
