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

% Smoothed initial pressure
Rgrid.setUMatrix(sourceKW.p0);

%========================================
% Impulse Response
%========================================
% Set time
dt = 5e-8;
tMax = 4e-5;
Rgrid.setTime(dt, tMax);

% Compute impulse response
Rgrid.impulse_additive('IV');
save gridRT_impulse.mat Rgrid;
%========================================
% Ray Shooting Parameters
%========================================
load gridRT_impulse;
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
signalRT_smooth = zeros(nSources, length(source(1).aForward));
for n = 1:nSources
    signalRT_smooth(n, :) = source(n).aForward;
end

save signalRT_smooth signalRT_smooth;
%========================================
% Sensor Selection
%========================================
n1 = round(Rgrid.Nx + 2*(round(Rgrid.Ny/3) - 2)   + 2);      % 1st sensor: Nx + 2*(Ny/3-2)   + 2
n2 = round(Rgrid.Nx + 2*(round(2*Rgrid.Ny/3) - 2) + 1);      % 2nd sensor: Nx + 2*(2*Ny/3-2) + 1
n3 = round(Rgrid.Nx + 2*(Rgrid.Ny - 2) + round(Rgrid.Nx/2)); % 3rd sensor: Nx + 2*(Ny-2)     + Nx/2
sensor1 = x{n1};
sensor2 = x{n2};
sensor3 = x{n3};
sourceSel(1) = Rgrid.newSource(sensor1, pi/2, 3*pi/2, nRays, tStep, tMax);
sourceSel(2) = Rgrid.newSource(sensor2, -pi/2, pi/2, nRays, tStep, tMax);
sourceSel(3) = Rgrid.newSource(sensor3, pi, 2*pi, nRays, tStep, tMax);
Rgrid.computeHamil(sourceSel(1), 'p');
Rgrid.computeHamil(sourceSel(2), 'p');
Rgrid.computeHamil(sourceSel(3), 'p');


%==============================
% Measure time
%==============================
end_time = clock;
% Measure computational time
disp(['  FORWARD PROBLEM: total computation time  ' num2str(etime(end_time, start_time))]);

%==================================================================================
% Plot results
%==================================================================================
position       = [700 700 320 630];
positionNoY    = [700 700 300 600];
positionNoYBar = [700 700 363 600];
positionYBar   = [700 700 390 630];

%==============================
% Sound Speed
%==============================
figure;
surf(1e3*Rgrid.xAxis, 1e3*Rgrid.yAxis, medium.sound_speed', 'EdgeColor', 'none');
axis image;
xlabel('x [mm]');
ylabel('y [mm]');
colorbar();
view(2);

%==============================
% Initial Pressure
%==============================
figure;
surf(1e3*Rgrid.xAxis, 1e3*Rgrid.yAxis, Rgrid.u', 'EdgeColor', 'none');
axis image;
xlabel('x [mm]');
ylabel('y [mm]');
colorbar();
view(2);

%==============================
% Sources and their respective Rays
%==============================
axisRgrid = [0 1e3*Rgrid.xAxis(end) 0 1e3*Rgrid.yAxis(end)];
nRaysPlot = 100;
colorList = cool(nRaysPlot);
% Sensor 1
figure;
hold on;
axis(axisRgrid);
for n = 1:nRaysPlot
    index = n*floor(nRays/nRaysPlot);
    plot(1e3*sourceSel(1).x(index, :, 1), 1e3*sourceSel(1).x(index, :, 2), 'Color', colorList(n, :), 'LineWidth', 1.5);
end
plot(1e3*sensor1(1), 1e3*sensor1(2), 'ok', 'MarkerSize', 8, 'LineWidth', 1, 'MarkerFaceColor','r');
box on;
set(gcf, 'pos', position);
xlabel('x [mm]');
ylabel('y [mm]');
set(gca,'FontSize',13)

% Sensor 2
figure;
hold on;
axis(axisRgrid);
for n = 1:nRaysPlot
    index = n*floor(nRays/nRaysPlot);
    plot(1e3*sourceSel(2).x(index, :, 1), 1e3*sourceSel(2).x(index, :, 2), 'Color', colorList(n, :), 'LineWidth', 1.5);
end
plot(1e3*sensor2(1), 1e3*sensor2(2), 'ok', 'MarkerSize', 8, 'LineWidth', 1, 'MarkerFaceColor','g');
box on;
set(gcf, 'pos', position);
xlabel('x [mm]');
set(gca,'FontSize',13)


% Sensor 3
figure;
hold on;
axis(axisRgrid);
for n = 1:nRaysPlot
    index = n*floor(nRays/nRaysPlot);
    plot(1e3*sourceSel(3).x(index, :, 1), 1e3*sourceSel(3).x(index, :, 2), 'Color', colorList(n, :), 'LineWidth', 1.5);
end
plot(1e3*sensor3(1), 1e3*sensor3(2), 'ok', 'MarkerSize', 8, 'LineWidth', 1, 'MarkerFaceColor','b');
box on;
set(gcf, 'pos', position);
xlabel('x [mm]');
set(gca,'FontSize',13)


%==============================
% Time Signals - ONLY RT
%==============================
posForward = [700 700 600 400];
set(0,'DefaultFigurePaperPositionMode','auto');


normRT = max(signalRT_smooth(:));
normKWave = max(sensor_data(:));

% Signals
signalRT1 = sourceSel(1).aForward/normRT;
signalRT2 = sourceSel(2).aForward/normRT;
signalRT3 = sourceSel(3).aForward/normRT;
signalKW1 = sensor_data(n1, :)/normKWave;
signalKW2 = sensor_data(n2, :)/normKWave;
signalKW3 = sensor_data(n3, :)/normKWave;

% Sensor 1
figure; hold on;
axis([0 40 -.4 .8]);
plot(1e6*Rgrid.tForward, signalRT1, 'Color', colourMapV(1), 'LineWidth', 3);
plot(1e6*kgrid.t_array,  signalKW1, 'Color',           'k', 'LineWidth', 2);
box on; grid on;
legend('Sensor 1 - HG', 'Sensor 1 - kWave');
xlabel('t [\mus]');
ylabel('Amplitude');
set(gca,'FontSize',13);
set(gcf, 'pos', posForward);
% Sensor 2;
figure; hold on;
axis([0 40 -.4 .8]);
plot(1e6*Rgrid.tForward, signalRT2, 'Color', colourMapV(2), 'LineWidth', 3);
plot(1e6*kgrid.t_array,  signalKW2, 'Color',           'k', 'LineWidth', 2);
box on; grid on;
legend('Sensor 2 - HG', 'Sensor 2 - kWave');
xlabel('t [\mus]');
ylabel('Amplitude');
set(gca,'FontSize',13)
set(gcf, 'pos', posForward);
% Sensor 3
figure; hold on;
axis([0 40 -.4 .8]);
plot(1e6*Rgrid.tForward, signalRT3, 'Color', colourMapV(3), 'LineWidth', 3);
plot(1e6*kgrid.t_array,  signalKW3, 'Color',           'k', 'LineWidth', 2);
box on; grid on;
legend('Sensor 3 - HG', 'Sensor 3 - kWave');
xlabel('t [\mus]');
ylabel('Amplitude');
set(gca,'FontSize',13)
set(gcf, 'pos', posForward);
% Errors
error1 = signalRT1 - signalKW1;
error2 = signalRT2 - signalKW2;
error3 = signalRT3 - signalKW3;

figure; 
set(gcf, 'pos', posForward);
axis([0 40 -.2 .3]);
hold on;
plot(1e6*Rgrid.tForward, error1, 'Color', colourMapV(1), 'LineWidth', 2);
plot(1e6*Rgrid.tForward, error2, 'Color', colourMapV(2), 'LineWidth', 2);
plot(1e6*Rgrid.tForward, error3, 'Color', colourMapV(3), 'LineWidth', 2);
legend('Sensor 1', 'Sensor 2', 'Sensor 3');
box on; grid on;
xlabel('t [\mus]');
ylabel('Error HG');
set(gca,'FontSize',13);



