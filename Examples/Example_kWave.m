%================================================================================
% Example for gridRT class - kWave simulation for comparison purposes
% 
% Note: Please set the working directory to the location of this file
%===============================================================================
clear all;
close all;

run colourMap;
%=========================================================================
% DOMAIN DEFINITION
%=========================================================================
% create the computational grid
Nx = 128;           % number of grid points in the x (row) direction
Ny = 256;           % number of grid points in the y (column) direction
dx = 2e-4;          % grid point spacing in the x direction [m]
dy = 2e-4;          % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy);

%==============================
% define the properties of the propagation medium
%==============================
% Build domain
c0 = 1500;
factor = 0.1;
% Build Peaks
frame = 2.5;
x = -frame:(2*frame)/(Ny-1):frame;
y = -frame:(2*frame)/(Nx-1):frame;
[X, Y] = meshgrid(x, y);
p = peaks(X, Y);
medium.sound_speed = c0*(ones(Nx, Ny) + factor*p/max(p(:)));
medium.density = 1;
    
% compute time
dt = 5e-8;
tMax = 4e-5;
kgrid.t_array = 0:dt:tMax;
%[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);

% Build initial pressure
inputIm = imread('../Phantoms/Veins_modified.jpg');
sourceKW.p0 = double(255-inputIm)/255;
% smooth the initial pressure distribution and restore the magnitude
sourceKW.p0 = smooth(kgrid, sourceKW.p0, true);


%=========================================================================
% SIMULATION
%=========================================================================
% Define the sensors
sensor.mask = zeros(Nx, Ny);
sensor.mask(1, :) = 1;
sensor.mask(end, :) = 1;
sensor.mask(:, 1) = 1;
sensor.mask(:, end) = 1;

% set the input arguements: force the PML to be outside the computational
% grid; switch off p0 smoothing within kspaceFirstOrder2D
input_args = {'PMLInside', false, 'PlotPML', false, 'Smooth', false};
% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, sourceKW, sensor, input_args{:});

save sensor_data.mat kgrid sensor sourceKW medium sensor_data input_args;
% reset the initial pressure
sensor.time_reversal_boundary_data = sensor_data;

%==============================
% Adjoint
%==============================
sensor_adjoint.mask = ones(kgrid.Nx, kgrid.Ny);
sensor_adjoint.record = {'p_final'};
source_adjoint.p_mask = zeros(kgrid.Nx, kgrid.Ny);
source_adjoint.p_mask(1, :) = 1;
source_adjoint.p_mask(end, :) = 1;
source_adjoint.p_mask(:, 1) = 1;
source_adjoint.p_mask(:, end) = 1;
source_adjoint.p = fliplr(sensor_data);
p0_recon_adjoint = kspaceFirstOrder2D(kgrid, medium, source_adjoint, sensor_adjoint, input_args{:});

% Save data
save recon_data_adjoint.mat kgrid medium source_adjoint sensor_adjoint p0_recon_adjoint input_args;
