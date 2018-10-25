function grid = forward_pressure(grid, source)
% FORWARD_PRESSURE computes the pressure contribution of each point
% in a Ray given an initial pressure
% grid = forward_pressure(grid, source)
% INPUTS
% grid: gridRT object that defines the domain
% source: source to compute the pressure
%
% OUTPUTS
% grid: gridRT object
% 
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke

%%  % Display message
%%  msg = strcat({'Computing pressure for source '}, int2str(nS), {'...'});
%%  disp(msg{1});

% Compute Pressure
%pressure = source.pressure.*source.revAmplitude;
pressure = source.pressure.*source.amplitude;
pressure(isinf(pressure) | isnan(pressure)) = 0;
source.insertPressureVector(pressure);
