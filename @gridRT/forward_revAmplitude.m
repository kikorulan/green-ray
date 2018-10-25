function revA = forward_revAmplitude(grid, source)
% REVAMPLITUDE Computes the attenuation at the ray trajectories by inverting the interface losses
% from the source to the ray points
%function A = revAmplitude(grid, source)
% INPUTS
% grid: gridRT object that defines the domain
% source: source to compute the reverse amplitude
%
% OUTPUTS:
% A: amplitude
% 
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke

%%  % Display message
%%  msg = strcat({'Computing the reverse amplitude for source '}, int2str(nS), {'...'});
%%  disp(msg{1});

% Compute forward amplitude
A = source.amplitude;
% Number of rays
[nR nPoints dim] = size(source.x);
% Amplitude due to spreading losses
deltaX = repmat(grid.deltaX, [nR nPoints]);
Aspread = sqrt(deltaX./(deltaX + source.xTD));
% Compute the reverse amplitude
revA = Aspread.*Aspread./A;
%revA = A;
% Insert computed reverse amplitude
source.insertRevAVector(revA);
%source.insertRevAVector(A);

