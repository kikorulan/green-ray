function pixelAReverse = operator_inverse(grid, source, index, forward_data)
% OPERATOR_INVERSE computes the inverse operator over the forward data
% pixelAReverse = operator_inverse(grid, source, forward_data)
%
% INPUTS
%   grid          - gridRT object that defines the domain
%   source        - vector of sources to compute the iteration
%   index         - choose the source between 'all' or a particular index
%   forward_data  - forward data to use in the operator
%
% OUTPUTS:
%   source        - vector containing the computed sources
% 
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke

% Initialise the pressure
pixelAReverse = zeros(grid.Nx, grid.Ny);

% Choose between all sources and a particular one
if(strcmp(index, 'all'))
    selectSource = source;
else
    selectSource = source(index);
end

% Loop over all sources
for n = 1:length(selectSource)
    selectSource(n).setForwardSignal(forward_data(n, :));
    grid.inverse_beam(selectSource(n));
    pixelAReverse = pixelAReverse + selectSource(n).pixelAReverse;
end

% Force to zero pressure at sensors
for n = 1:length(source)
    coord = grid.findCoordinates(source(n).x0);
    grid.pixelAReverse(coord(1), coord(2)) = 0;
    pixelAReverse(coord(1), coord(2)) = 0;
end

