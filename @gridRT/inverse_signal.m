function pixelAReverseSensors = inverse_signal(grid, source)
% INVERSE_SIGNAL computes the time reverse signal for all sources
%grid = inverse_signal(grid, source)
%
% INPUTS
% grid: object from gridRT class
% source: vector of sources
% 
% OUTPUTS
% pixelAReverseSensors: matrix with the time reverse
% Copyright (C) Kiko Rul·lan, Marta M. Betcke

% Number of sources
nSources = length(source);

% Allocate matrices
grid.pixelAReverse = zeros(grid.Nx, grid.Ny);
pixelAReverseSensors = zeros(grid.Nx, grid.Ny, length(source));

% Sum all time reversal matrices
for n = 1:nSources
    pixelAReverse = source(n).pixelAReverse;
    pixelAReverse(isnan(pixelAReverse)) = 0;
    % Option 1
    grid.pixelAReverse = grid.pixelAReverse + pixelAReverse;
    % Option 2
    %factor = sqrt(sum(source(n).pixelAReverse(:).^2));
    %grid.pixelAReverse = grid.pixelAReverse + source(n).pixelAReverse/(factor^2);
    % Option 3
    %grid.pixelAReverse = grid.pixelAReverse + (1-0.8*(max(20, abs(128-n)))/128)*source(n).pixelAReverse;
    
    pixelAReverseSensors(:, :, n) = source(n).pixelAReverse;
end

% Positivity condition
%grid.pixelAReverse(grid.pixelAReverse < 0) = 0;
% Normalisation
%maxA = max(grid.pixelAReverse(:));
%grid.pixelAReverse = real(grid.pixelAReverse/maxA);
