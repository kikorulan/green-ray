function grid = computeHamil(grid, source, mode, del)
% COMPUTEHAMIL computes the hamiltonian trajectories of the given source nS
%function grid = computeHamil(grid, source)
% INPUTS
% grid: gridRT object that defines the domain
% source: source to compute the ray trajectories
% mode: choose from
%    - 'p': proximal ray mode
%    - 'g': dynamical ray tracing mode
% del: true to delete info from source
%
% OUTPUTS:
% grid: gridRT object that defines the domain
% 
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke

% Number of steps of the computation
step = source.step;
[nR, ~, dim] = size(source.x);
nPoints = source.nPoints;
% Allocate rays
source.allocateRays(nR, nPoints, dim);
% Trajectory
grid.forward_trajectories(source);

% Amplitude
switch mode
    case 'p'
        disp('PR mode')
        grid.forward_amplitude(source);
    case 'o'
        disp('ODE mode')
        grid.forward_amplitudeODE(source);
    case 's'
        disp('Gauss mode')
        grid.forward_amplitudeGauss(source);
    case 'g'
        disp('GB mode')
        grid.forward_amplitudeGB(source);
end

% Reverse Amplitude
grid.forward_revAmplitude(source);

% Pressure
grid.forward_pressure(source);

% Rays to grid
grid.forward_raysToGrid(source);

% Beam
grid.forward_beam(source);

% Forward Signal
grid.forward_timeSignal(source);

% Delete auxiliary info
if nargin > 3
    if del
        source.deleteAuxInfo();
    end
end
% Display message
msg = strcat({'Computed given source with '}, int2str(nPoints), {' steps.'});
disp(msg{1});

