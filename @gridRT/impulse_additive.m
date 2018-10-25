function grid = impulse_additive(grid, impulseType)
% IMPULSE_ADDITIVE computes the impulse response for the given domain
% using kWave
%function grid = impulse_additive(grid, nS)
% INPUTS
% grid: gridRT object that defines the domain
% impulseType: choose between the possible impulse responses
%   - 'IV': initial value problem - solved via k-Wave
%   - 'TS': time varying source - solved via k-Wave
%   - 'GF': use an approximation of the Green's function derivative
% saveResults: flag to indicate to save results or not
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

if (strcmp(impulseType, 'IV') | strcmp(impulseType, 'TS'))
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
    
    % Switch for source type
    switch impulseType
        case 'TS'
        % Define a single source point
        source.p_mask = zeros(Nx, Ny);
        source.p_mask(end - Nx/4, Ny/2) = 1;
        % Define a time varying gaussian
        source.p = zeros(1, length(kgrid.t_array));
        source.p(floor(length(kgrid.t_array)/5)) = 1;
        % Filter the source to remove high frequencies not supported by the grid
        source.p = filterTimeSeries(kgrid, medium, source.p);
    
        case 'IV'
        % define a single source point
        source.p0 = zeros(Nx, Ny);
        source.p0(end - Nx/4, Ny/2) = 1;
    end
    
    % Run the simulation
    sensorData = kspaceFirstOrder2D(kgrid, medium, source, sensor);
    % Rename data
    t_array = kgrid.t_array;
    Filter = sensorData.p;     
elseif (strcmp(impulseType, 'GF'))
    %==============================
    % Green's function approximation
    %==============================
    deltaGaussP = @(x, a) -2*x.*a^3/sqrt(pi).*exp(-x.*x.*a^2);
    H = @(x) 0.5*sign(x) + 0.5;
    % Time array
    t_array = 2*dt:dt:tMax;
    % Green's function
    GF = 1/2/pi*H(t_array)./(c0*c0*t_array);
    GF(GF == inf) = 0;
    % Delta derivative approximation
    epsilon = 1/(dt*80);
    t_array_delta = -5e2*dt:dt:3e2*dt;
    deltaDerivative = deltaGaussP(t_array_delta, epsilon);
    figure;
    plot(t_array_delta, deltaDerivative);
    title('Delta');
    % Compute filter
    Filter = conv(GF, deltaDerivative);
    Filter = Filter(1:length(GF));
    figure, plot(Filter);
else
    error('Not known impulse response');
end


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
grid.setCFilter(c0);
grid.setDelayFilter(delay);
grid.setTFilter(tFilter);
grid.setFilter(Filter);



