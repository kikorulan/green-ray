function grid = inverse_filter_fft(grid, nFilters)
% INVERSE_FILTER_FFT computes FFT of the filters at each pixel
% grid = inverse_filter_fft(grid, nFilters)
%
% INPUTS
% grid: object from gridRT class
% nFilters: number of filters to compute
% tMaxReverse: maximum propagation reverse time
% 
% OUTPUTS
% grid: object from gridRT class
%   - pixelFilterFFT: FFT of the filters at each pixel
%
% Copyright (C) Kiko Rul·lan, Marta M. Betcke

if (isempty(grid.Filter) || isempty(grid.cFilter) || isempty(grid.tFilter))
   error('The filter has not been set.');
end
% Display initialisation message
msg = strcat({'Computing pixel FFT filters...'});
disp(msg{1});

% Filter time array
cMax = max(grid.c(:)); % maximum sound speed
cMin = min(grid.c(:)); % minumum sound speed
dt = grid.dt; % time signal increment
delayPropMax = length(grid.tForward); % maximum propagation delay
delayFilterMax = floor(grid.delayFilter*grid.cFilter/cMin); % maximum delay for the filter
tFilter = grid.tFilter; % time array for the filter
lengthFilter = length(tFilter);
lengthForward = length(grid.tForward);

grid.tReverse = 0:dt:dt*(lengthForward + lengthFilter + delayPropMax - 1);
lengthReverse = length(grid.tReverse);
switch grid.dim
    case 2
        if (cMin == cMax)
            FilterPad = padarray(grid.Filter, [0 lengthForward + delayPropMax + delayFilterMax], 0, 'post');
            grid.pixelFilterFFT = repmat(permute(fft(FilterPad), [1 3 2]), [grid.Nx grid.Ny 1]);
            grid.pixelDelayFilter = repmat(grid.delayFilter, [grid.Nx grid.Ny 1]);
        else
            % Generate Filter time array
            grid.pixelFilterFFT = zeros(grid.Nx, grid.Ny, lengthReverse + delayFilterMax);
            grid.pixelDelayFilter = zeros(grid.Nx, grid.Ny);
            % Vector of filters
            c = [cMin:(cMax - cMin)/nFilters:cMax + (cMax - cMin)/nFilters];
            for j = 1:length(c)-1
                msg = strcat({'Computing pixel FFT filter '}, int2str(j));
                disp(msg{1});
                % Compute splines for filter
                splineFilter = c(j)*spline(grid.tFilter*grid.cFilter/c(j), grid.Filter, tFilter);
                FilterPad = padarray(splineFilter, [0 lengthForward + delayPropMax + delayFilterMax], 0, 'post');
                FilterFFT = permute(fft(FilterPad), [1 3 2]);
                pixelFilterFFT = repmat(FilterFFT, [grid.Nx grid.Ny 1]);
                binaryC = (c(j+1) > grid.c) & (grid.c >= c(j));
                binaryCMatrix = repmat(binaryC, [1 1 lengthReverse+delayFilterMax]);
                grid.pixelFilterFFT = grid.pixelFilterFFT + pixelFilterFFT.*binaryCMatrix;
                % Compute delay for filter
                delayFilter = find(splineFilter == max(splineFilter)) - 1;
                grid.pixelDelayFilter = grid.pixelDelayFilter + delayFilter.*binaryC;
            end
        end
    case 3
        error('Not implemented');
end
  
