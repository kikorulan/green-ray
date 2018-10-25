function h = plot_YNM(source, nR)
% PLOT_YNM plots the matrices for the Riccati equation
%function h = plotAmplitude(grid, nSVector, nColours)
% INPUTS
% source: source
% nR: ray number
% 
%
% OUTPUTS:
% h: figure handler
% 
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke

colourList = cool(9);
%===========
% Y Matrix
%===========
Y = source.Y(nR, :, :);

h = figure;
title('Y - real');
hold on;
plot(real(Y(:, :, 1)), 'Color', colourList(1, :));
plot(real(Y(:, :, 2)), 'Color', colourList(2, :));
plot(real(Y(:, :, 3)), 'Color', colourList(3, :));
plot(real(Y(:, :, 4)), 'Color', colourList(4, :));
plot(real(Y(:, :, 5)), 'Color', colourList(5, :));
plot(real(Y(:, :, 6)), 'Color', colourList(6, :));
plot(real(Y(:, :, 7)), 'Color', colourList(7, :));
plot(real(Y(:, :, 8)), 'Color', colourList(8, :));
plot(real(Y(:, :, 9)), 'Color', colourList(9, :));
grid on;

h = figure;
title('Y - imag');
hold on;
plot(imag(Y(:, :, 1)), 'Color', colourList(1, :));
plot(imag(Y(:, :, 2)), 'Color', colourList(2, :));
plot(imag(Y(:, :, 3)), 'Color', colourList(3, :));
plot(imag(Y(:, :, 4)), 'Color', colourList(4, :));
plot(imag(Y(:, :, 5)), 'Color', colourList(5, :));
plot(imag(Y(:, :, 6)), 'Color', colourList(6, :));
plot(imag(Y(:, :, 7)), 'Color', colourList(7, :));
plot(imag(Y(:, :, 8)), 'Color', colourList(8, :));
plot(imag(Y(:, :, 9)), 'Color', colourList(9, :));
grid on;

%===========
% N Matrix
%===========
N = source.N(nR, :, :);

h = figure;
title('N - real');
hold on;
plot(real(N(:, :, 1)), 'Color', colourList(1, :));
plot(real(N(:, :, 2)), 'Color', colourList(2, :));
plot(real(N(:, :, 3)), 'Color', colourList(3, :));
plot(real(N(:, :, 4)), 'Color', colourList(4, :));
plot(real(N(:, :, 5)), 'Color', colourList(5, :));
plot(real(N(:, :, 6)), 'Color', colourList(6, :));
plot(real(N(:, :, 7)), 'Color', colourList(7, :));
plot(real(N(:, :, 8)), 'Color', colourList(8, :));
plot(real(N(:, :, 9)), 'Color', colourList(9, :));
grid on;

h = figure;
title('N - imag');
hold on;
plot(imag(N(:, :, 1)), 'Color', colourList(1, :));
plot(imag(N(:, :, 2)), 'Color', colourList(2, :));
plot(imag(N(:, :, 3)), 'Color', colourList(3, :));
plot(imag(N(:, :, 4)), 'Color', colourList(4, :));
plot(imag(N(:, :, 5)), 'Color', colourList(5, :));
plot(imag(N(:, :, 6)), 'Color', colourList(6, :));
plot(imag(N(:, :, 7)), 'Color', colourList(7, :));
plot(imag(N(:, :, 8)), 'Color', colourList(8, :));
plot(imag(N(:, :, 9)), 'Color', colourList(9, :));
grid on;

%===========
% M Matrix
%===========
M = source.M(nR, :, :);

h = figure;
title('M - real');
hold on;
plot(real(M(:, :, 1)), 'Color', colourList(1, :));
plot(real(M(:, :, 2)), 'Color', colourList(2, :));
plot(real(M(:, :, 3)), 'Color', colourList(3, :));
plot(real(M(:, :, 4)), 'Color', colourList(4, :));
plot(real(M(:, :, 5)), 'Color', colourList(5, :));
plot(real(M(:, :, 6)), 'Color', colourList(6, :));
plot(real(M(:, :, 7)), 'Color', colourList(7, :));
plot(real(M(:, :, 8)), 'Color', colourList(8, :));
plot(real(M(:, :, 9)), 'Color', colourList(9, :));
grid on;

h = figure;
title('M - imag');
hold on;
plot(imag(M(:, :, 1)), 'Color', colourList(1, :));
plot(imag(M(:, :, 2)), 'Color', colourList(2, :));
plot(imag(M(:, :, 3)), 'Color', colourList(3, :));
plot(imag(M(:, :, 4)), 'Color', colourList(4, :));
plot(imag(M(:, :, 5)), 'Color', colourList(5, :));
plot(imag(M(:, :, 6)), 'Color', colourList(6, :));
plot(imag(M(:, :, 7)), 'Color', colourList(7, :));
plot(imag(M(:, :, 8)), 'Color', colourList(8, :));
plot(imag(M(:, :, 9)), 'Color', colourList(9, :));
grid on;
