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
source.pixelDistance = inf(grid.Nx, grid.Ny, grid.Nz);
source.pixelTime = inf(grid.Nx, grid.Ny, grid.Nz);
source.pixelAttenuation = zeros(grid.Nx, grid.Ny, grid.Nz);
source.pixelPressure = zeros(grid.Nx, grid.Ny, grid.Nz);
%source.pixelAngleCorrection = zeros(grid.Nx, grid.Ny, grid.Nz);

% Loop over lenghtRay
for j = 1:nPoints-1
    % Find the coordinates of the trajectories
    x = source.x(:, j, :);
    xCoord = grid.findCoordinates(x);
    binaryCoord = xCoord(:, :, 1) < inf;
    xFin = x(binaryCoord, :, :);
    xCoordFin = xCoord(binaryCoord, :, :);
    xCoordLin = xCoordFin(:, 1, 1) + (xCoordFin(:, 1, 2) - 1)*grid.Nx;
    xGrid = cat(3, (xCoordFin(:, 1, 1) - 1)*grid.dx, (xCoordFin(:, 1, 2) - 1)*grid.dy);
    % Find the angles 
    direction = source.x(:, j+1, :) - x;
    angle = atan(abs(direction(:, :, 2)./direction(:, :, 1)));
    angle = min(angle, pi/2 - angle);
    correctionAngle = 1./cos(angle);
    correctionAngleBin = correctionAngle(binaryCoord);
    % Compute distance to grid point
    distance = sqrt(sum((xFin - xGrid).^2, 3));
    % Compute time vector
    time = source.phi(1, j);
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
    source.pixelTime(xCoordLin(modifiedRays)) = repmat(time, sum(modifiedRays), 1);
    % Attenuation                
    %attenuation = source.amplitude(binaryCoord, j, :).*sqrt(abs(source.q(binaryCoord, j, :)));
    %attenuation = source.amplitude(binaryCoord, j, :).*sqrt(abs(source.q(binaryCoord, j, :))).*correctionAngleBin;
    %attenuation = source.amplitude(binaryCoord, j, :).*abs(source.q(binaryCoord, j, :)).*correctionAngleBin.*correctionAngleBin;
    %attenuation = 0*source.q(binaryCoord, j, :) + 1;
    %attenuation = source.amplitude(binaryCoord, j, :).*source.xTD(binaryCoord, j, :);
    %attenuation = source.amplitude(binaryCoord, j, :);
    %attenuation = source.amplitude(binaryCoord, j, :).*source.q(binaryCoord, j, :);
    attenuation = source.revAmplitude(binaryCoord, j, :);
    attenuation(isinf(attenuation) | isnan(attenuation)) = 0;
    source.pixelAttenuation(xCoordLin(modifiedRays)) = attenuation(modifiedRays);
    % Pressure 
    pressure = source.pressure(binaryCoord, j, :);
    source.pixelPressure(xCoordLin(modifiedRays)) = pressure(modifiedRays);
    % Angle correction
    source.pixelAngleCorrection(xCoordLin(modifiedRays)) = correctionAngleBin(modifiedRays);
    %%%%%%%%%%%%%% DEBUG %%%%%%%%%%%%%%%
    %%  binaryAttenuation = double(isnan(source.pixelAttenuation));
    %%  if(sum(binaryAttenuation(:)) > 0)
    %%      save x.mat xCoord xFin xCoordFin xCoordLin xGrid modifiedRays j;
    %%      error('NaN attenuation introduced');
    %%  end
end

source.pixelPressure(isnan(source.pixelPressure)) = 0;

 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MULTI RAY IMPLEMENTATION 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  function pixelNRay = raysToGrid(grid, nS)
%%  % Compute the time propagation from source nS given the time signal
%%      [dim lengthRay nR] = size(source.x);
%%      
%%      if (isempty(grid.Filter) || isempty(grid.cFilter) || isempty(grid.tFilter))
%%         error('The filter has not been set.');
%%      end
%%      % Display initialisation message
%%      msg = strcat({'Computing time propagation for source '}, int2str(nS), {'...'});
%%      disp(msg{1});
%%      
%%      % Minimum tau increment to guarantee different ray sources
%%      deltaTau = 10*sqrt(grid.dx^2 + grid.dy^2 + grid.dz^2)*min(grid.c(:));
%%      switch dim
%%          case 2
%%              nRays = 10; % number of superposed rays
%%              % Matrix initialisation
%%              source.pixelDistance = inf(grid.Nx, grid.Ny, nRays);
%%              source.pixelTime = inf(grid.Nx, grid.Ny, nRays);
%%              source.pixelAttenuation = zeros(grid.Nx, grid.Ny, nRays);
%%              pixelTau = -inf(grid.Nx, grid.Ny);
%%              pixelNRay = zeros(grid.Nx, grid.Ny);
%%              % Loop over lenghtRay
%%              for tau = 1:lengthRay
%%                  % Find the coordinates of the trajectories
%%                  x = source.x(:, tau, :);
%%                  xCoord = grid.findCoordinates(x);
%%                  xFin = x(:, :, xCoord(1, :, :) < inf);
%%                  xCoordFin = xCoord(:, :, xCoord(1, :, :) < inf); % finite coordinates
%%                  xCoordLin = xCoordFin(1, 1, :) + (xCoordFin(2, 1, :) - 1)*grid.Nx; % conversion to linear coordinates
%%                  xGrid = [(xCoordFin(1, 1, :) - 1)*grid.dx; (xCoordFin(2, 1, :) - 1)*grid.dy];
%%                  % Compute tau condition
%%                  tauCondition = double(tau - pixelTau(xCoordLin) > deltaTau);
%%                  xCoordLin = xCoordLin(tauCondition == 1);
%%                  pixelTau(xCoordLin) = tau;
%%                  pixelNRay(xCoordLin) = pixelNRay(xCoordLin) + 1;
%%                  xCoordLinZ = xCoordLin + min((pixelNRay(xCoordLin)-1), nRays-1)*grid.Nx*grid.Ny; % consider number of rays
%%                  % Compute time vector
%%                  time = source.phi(:, tau, xCoord(1, :, :) < inf);
%%                  timeTauC = time(tauCondition == 1);
%%                  % Compute distance to grid point
%%                  distance = sqrt(sum((xFin - xGrid).^2, 1));
%%                  distanceTauC = distance(tauCondition == 1);
%%                  %%%%% Update grid
%%                  [timeSort, oRaysT] = sort(timeTauC, 3, 'descend');
%%                  [distanceSort, oRaysD] = sort(distanceTauC, 3, 'descend');
%%                  xCoordLinSort = xCoordLinZ(oRaysT); % reorder
%%                  % Time
%%                  source.pixelTime(xCoordLinSort) = min(timeSort, source.pixelTime(xCoordLinSort));
%%                  modifiedRays = (source.pixelTime(xCoordLinZ) - timeTauC == 0);
%%                  % Distance
%%                  source.pixelDistance(xCoordLinZ(modifiedRays)) = distance(modifiedRays);;
%%                  % Attenuation                
%%                  attenuation = source.amplitude(:, tau, xCoord(1, :, :) < inf);
%%                  attenuationTauC = attenuation(tauCondition == 1);
%%                  source.pixelAttenuation(xCoordLinZ(modifiedRays)) = attenuationTauC(modifiedRays);
%%                  %%%%%%%%%%%%%% DEBUG %%%%%%%%%%%%%%%
%%                  binaryAttenuation = double(isnan(source.pixelAttenuation));
%%                  if(sum(binaryAttenuation(:)) > 0)
%%                      save x.mat xCoord xFin xCoordFin xCoordLin xGrid modifiedRays tau;
%%                      error('NaN attenuation introduced');
%%                  end
%%              end
%%          case 3
%%              error('Not implemented');
%%      end
%%  end   
