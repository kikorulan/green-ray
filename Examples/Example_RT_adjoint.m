%================================================================================
% Example for gridRT class
% 
% Note: Please set the working directory to the location of this file
%===============================================================================

close all; clear all;

% LOAD k-Wave FORWARD DATA
load sensor_data.mat;
run colourMap;

%========================================
% Rgrid definition
%========================================
Nx = 128;           % number of Rgrid points in the x (row) direction
Ny = 256;           % number of Rgrid points in the y (column) direction
dx = 2e-4;          % Rgrid point spacing in the x direction [m]
dy = 2e-4;          % Rgrid point spacing in the y direction [m]
Rgrid = gridRT(Nx, dx, Ny, dy);

% Build domain
Rgrid.setCMatrix(medium.sound_speed);

% Non-smoothed initial pressure
inputIm = imread('../Phantoms/Veins_modified.jpg');
u0 = double(255-inputIm)/255;
Rgrid.setUMatrix(u0);

%========================================
% Impulse Response
%========================================
% Set time
dt = 5e-8;
tMax = 4e-5;
Rgrid.setTime(dt, tMax);

% Compute impulse response
Rgrid.impulse_additive('IV');
%========================================
% Ray Shooting Parameters
%========================================
cMax = max(Rgrid.c(:));
dt = 5e-8;

% Measure computational time
tic;
start_time = clock;

% Number of rays & sources
nRays = 2000;% 800
nSources = 764;%256

% Parametrisation
tMax = 4e-5;
tStep = dt;

%================================================================================
% FORWARD PROBLEM
%================================================================================
% Sources locations
clear x;
for n = 1:Rgrid.Nx
    x{n} = cat(3, (n-1)*Rgrid.dx, 0);
end
for n = 1:Rgrid.Ny-2
    x{2*n-1+Rgrid.Nx}   = cat(3,                   0, n*Rgrid.dy);
    x{2*n  +Rgrid.Nx}   = cat(3, (Rgrid.Nx-1)*Rgrid.dx, n*Rgrid.dy);
end
for n = 1:Rgrid.Nx
    x{n+Rgrid.Nx+2*(Rgrid.Ny-2)} = cat(3, (n-1)*Rgrid.dx, (Rgrid.Ny-1)*Rgrid.dy);
end
source = Rgrid.computeForwardParallel(x, 0, 2*pi-0.01, nRays, tStep, tMax, 'p', true);
signalRT_nonsmooth = zeros(nSources, length(source(1).aForward));
for n = 1:nSources
    signalRT_nonsmooth(n, :) = source(n).aForward;
end

save signalRT_nonsmooth signalRT_nonsmooth;

%==============================
% Measure time
%==============================
end_time = clock;
% Measure computational time
disp(['  FORWARD PROBLEM: total computation time  ' num2str(etime(end_time, start_time))]);

%================================================================================
% ADJOINT PROBLEM
%================================================================================
load recon_data_adjoint.mat;
% Measure computational time
tic;
start_time = clock;
 
%========================================
% Compute filters
%========================================
nFilters = 100;
Rgrid.inverse_filter(nFilters);

%========================================
% Compute reverse signal
%========================================
nSources = 764;
for n = 1:nSources
    disp(n)
    Rgrid.inverse_beam(source(n));
end
%Rgrid.computeAdjointParallel(source);

pixelAReverseSensors = Rgrid.inverse_signal(source);

adjointRT_nonsmooth = Rgrid.pixelAReverse;
save adjointRT_nonsmooth adjointRT_nonsmooth;

% Measure computational time
end_time = clock;
disp(['  ADJOINT PROBLEM: total computation time ' num2str(etime(end_time, start_time))]);


%================================================================================
% MIX 
%================================================================================
load sensor_data;
load signalRT_nonsmooth;
%=========================================
% kWave reconstruction using HG data
%=========================================
sensor_data_RT = signalRT_nonsmooth;
% run the time reversal reconstruction
input_args = {'PMLInside', false, 'PlotPML', false, 'Smooth', false};
% Build sensor
sensor_adjoint.mask = ones(kgrid.Nx, kgrid.Ny);
sensor_adjoint.record = {'p_final'};
source_adjoint.p_mask = zeros(kgrid.Nx, kgrid.Ny);
source_adjoint.p_mask(1, :) = 1;
source_adjoint.p_mask(end, :) = 1;
source_adjoint.p_mask(:, 1) = 1;
source_adjoint.p_mask(:, end) = 1;
source_adjoint.p = fliplr(sensor_data);
p0_recon_adjoint = kspaceFirstOrder2D(kgrid, medium, source_adjoint, sensor_adjoint, input_args{:});

%=========================================
% HG reconstruction using kWave data
%=========================================
% Compute filters
nFilters = 100;
Rgrid.inverse_filter(nFilters);
%========================================
% Compute reverse signal
%========================================
nSources = 764;
for n = 1:nSources
    source(n).setForwardSignal(sensor_data(n, :));
end
for n = 1:nSources
    disp(n)
    Rgrid.inverse_beam(source(n));
end
pixelAReverseSensors = Rgrid.inverse_signal(source);

adjointKWaveForward_RT = Rgrid.pixelAReverse;
adjointRTForward_kWave = p0_recon_adjoint;
save adjointKWaveForward_RT adjointKWaveForward_RT;
save adjointRTForward_kWave adjointRTForward_kWave;

