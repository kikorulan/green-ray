function pixelAReverse = iterative_recon(grid, source, step, nIter)
% ITERATIVE_RECON computes an iteration of the reconstruction feeding back the data
% source = iterative_recon(grid, source)
% INPUTS
% grid: gridRT object that defines the domain
% source: source to compute the iteration
% step: step size for the iterations
% nIter: number of iterations
%
% OUTPUTS:
% source: vector containing the computed sources
% 
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke

% Initialise forward data
for n = 1:length(source)
    if(isempty(source(n).aForward_initial))
        source(n).setForwardSignal_initial(source(n).aForward);
    end
    b(n, :) = source(n).aForward_initial;
end

% Initialise pressure
pixelAReverse = grid.pixelAReverse;
% Compute norm
normU = max(grid.u(:));
% Compute the iterations
for n = 1:nIter
    % Compute norm
    normPixelAReverse = max(pixelAReverse(:));
    norm_factor = normU/normPixelAReverse;
    % Positivity constraint
    pixelAReverse = max(0, pixelAReverse);
    % Compute gradient    
    forward_data = grid.operator_forward(source, pixelAReverse, norm_factor);
    difference = forward_data - b;
    initial_pressure = 2*step*grid.operator_inverse(source, difference);
    % Update
    pixelAReverse = pixelAReverse - initial_pressure;
end
