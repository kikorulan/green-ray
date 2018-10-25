function grid = impulse_dirichlet(grid)
% IMPULSE_DIRICHLET computes the impulse response for a Dirichlet source for the given domain
% using k-Wave
%function grid = impulse_dirichlet(grid)
% INPUTS
% grid: gridRT object that defines the domain
%
% OUTPUTS:
% grid: gridRT object that defines the domain
% 
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke

% Define the properties of the propagation medium
c0 = min(grid.c(:));  % [m/s]
% Create the computational grid
Nx = 128;
Ny = 128;
dx = grid.dx;
dy = grid.dy;
% Create the time array
xMax = 2*Nx*dx;
tMax = 2*xMax/c0;
Filter = 1;
t_array = 1;
dt = grid.dt;

%==============================
% k-Wave simulation
%==============================
medium.sound_speed = c0;
kgrid = makeGrid(Nx, dx, Ny, dy);
if (grid.dt == 0)
    error('Time step for grid not set: grid.dt = 0.');
end
kgrid.t_array = 0:dt:tMax;
% Define a single sensor point
sensor.mask = zeros(Nx, Ny);
sensor.mask(Nx/4, Ny/2) = 1;
% Define the acoustic parameters to record
sensor.record = {'p'};

% Define a single source point
source.p_mask = zeros(Nx, Ny);
source.p_mask(end - Nx/4, Ny/2) = 1;
% Define a time varying gaussian
lengthSignal = length(kgrid.t_array);
x = 1:lengthSignal;
stdX2 = lengthSignal/10;
source.p = exp((-(x-stdX2).^2)/stdX2);
% Force dirichlet mode
source.p_mode = 'dirichlet';
% Filter the source to remove high frequencies not supported by the grid
source.p = filterTimeSeries(kgrid, medium, source.p);

% Run the simulation
sensorData = kspaceFirstOrder2D(kgrid, medium, source, sensor);
% Rename data
t_array = kgrid.t_array;
Filter = sensorData.p;     

% Find filter
maxF = max(Filter);
delay = find(Filter == maxF) - 1;
flipFilter = fliplr(Filter(1:delay));
indexIni = find(flipFilter < 1e-4*maxF, 1);
flipFilter = fliplr(Filter(delay+1:end));

lengthFilter = length(t_array) - delay + indexIni;
tFilter = t_array(1:lengthFilter);
Filter = Filter(delay-indexIni+1:end);
delay = find(Filter == maxF) - 1;

% Set values to grid
grid.delayFilterReverse = [delay; delay; delay];
grid.tFilterReverse = tFilter;
grid.FilterReverse = [Filter; Filter; Filter];

