function grid = forward_raysToGrid(grid, source)
% FORWARD_RAYSTOGRID converts the information from the rays from the given source
% into the grid format
%grid = forward_raysToGrid(grid, source)
%
% INPUTS
% grid: object from gridRT class
% source: source to do the conversion
% 
% OUTPUTS
% grid: object from gridRT class
%   - pixelDistance: distance matrix from the source to the pixel
%   - pixelTime: propagation time from the source to the pixel
%   - pixelAttenuation: attenuation from the source to the pixel
%
% Copyright (C) Kiko Rul·lan, Marta M. Betcke

% Number of points
nPoints = source.nPoints;

%%  % Display initialisation message
%%  msg = strcat({'Convert info from source '}, int2str(nS), {' rays to the grid...'});
%%  disp(msg{1});

% Matrix initialisation
source.pixelDistance = inf(grid.Nx, grid.Ny);
source.pixelTime = inf(grid.Nx, grid.Ny);
source.pixelAttenuation = zeros(grid.Nx, grid.Ny, nPoints);
source.pixelPressure = zeros(grid.Nx, grid.Ny, nPoints);

% Loop over lenghtRay
for j = 1:nPoints-1
    % Find the coordinates of the trajectories
    x = source.x(:, j, :);
    xCoord = grid.findCoordinates(x);
    binaryCoord = xCoord(:, :, 1) < inf;
    xFin = x(binaryCoord, :, :);
    xCoordFin = xCoord(binaryCoord, :, :);
    xCoordLin = xCoordFin(:, 1, 1) + (xCoordFin(:, 1, 2) - 1)*grid.Nx;
    xCoordLinZ = xCoordLin + (j-1)*grid.Nx*grid.Ny;
    xGrid = cat(3, (xCoordFin(:, 1, 1) - 1)*grid.dx, (xCoordFin(:, 1, 2) - 1)*grid.dy);
    % Compute distance to grid point
    distance = sqrt(sum((xFin - xGrid).^2, 3));
    % Compute time vector
    time = source.phi(binaryCoord, j, :);
    %%%%% Update grid
    [distanceSort, oRaysD] = sort(distance, 1, 'descend');
    [timeSort, oRaysT] = sort(time, 1, 'descend');
    % Distance
    source.pixelDistance(xCoordLin(oRaysD)) = min(distanceSort, source.pixelDistance(xCoordLin(oRaysD)));
    modifiedRays = (source.pixelDistance(xCoordLin) - distance == 0);
    %source.pixelDistance(xCoordLin(modifiedRays)) = distance(modifiedRays);
    % Time
    %source.pixelTime(xCoordLin(oRaysT)) = min(timeSort, source.pixelTime(xCoordLin(oRaysT)));
    %modifiedRays = (source.pixelTime(xCoordLin) - time == 0);
    source.pixelTime(xCoordLin(modifiedRays)) = time(modifiedRays);;
    % Attenuation                
    %attenuation = source.amplitude(binaryCoord, j, :).*sqrt(abs(source.q(binaryCoord, j, :)));
    %attenuation = source.amplitude(binaryCoord, j, :).*sqrt(abs(source.q(binaryCoord, j, :))).*correctionAngleBin;
    %attenuation = source.amplitude(binaryCoord, j, :).*abs(source.q(binaryCoord, j, :)).*correctionAngleBin.*correctionAngleBin;
    %attenuation = 0*source.q(binaryCoord, j, :) + 1;
    %attenuation = source.amplitude(binaryCoord, j, :).*source.xTD(binaryCoord, j, :);
    attenuation = source.amplitude(binaryCoord, j, :);
    %attenuation = source.amplitude(binaryCoord, j, :).*source.q(binaryCoord, j, :);
    source.pixelAttenuation(xCoordLinZ(modifiedRays)) = attenuation(modifiedRays);
    % Pressure 
    pressure = source.pressure(binaryCoord, j, :);
    source.pixelPressure(xCoordLinZ(modifiedRays)) = pressure(modifiedRays);
    %%%%%%%%%%%%%% DEBUG %%%%%%%%%%%%%%%
    %%  binaryAttenuation = double(isnan(source.pixelAttenuation));
    %%  if(sum(binaryAttenuation(:)) > 0)
    %%      save x.mat xCoord xFin xCoordFin xCoordLin xGrid modifiedRays j;
    %%      error('NaN attenuation introduced');
    %%  end
end

source.pixelPressure(isnan(source.pixelPressure)) = 0;

 
