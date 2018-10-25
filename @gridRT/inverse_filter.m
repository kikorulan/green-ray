function grid = inverse_filter(grid, nFilters)
% INVERSE_FILTER computes the filters corresponding to the given domain
%grid = inverse_filter(grid, nFilters)
%
% INPUTS
% grid: object from gridRT class
% nFilters: number of filters to compute
% 
% OUTPUTS
% grid: object from gridRT class
%   - pixelFilterFFT: filters at each pixel of the domain
%   - pixelDelayFilter: delay of the corresponding filter at each pixel of the domain
%
% Copyright (C) Kiko Rul·lan, Marta M. Betcke

if (isempty(grid.Filter) || isempty(grid.cFilter) || isempty(grid.tFilter))
   error('The filter has not been set.');
end
% Display initialisation message
msg = strcat({'Computing filters...'});
disp(msg{1});

% Filter time array
cMax = max(grid.c(:)); % maximum sound speed
cMin = min(grid.c(:)); % minumum sound speed
dt = grid.dt;

lengthFilter = length(grid.Filter);
grid.FilterReverse = nan(nFilters, lengthFilter);
grid.delayFilterReverse = nan(nFilters, 1);
% Convert information from rays to grid
switch grid.dim
    case 2
        % Generate Filter time array
        if (cMin == cMax)
            grid.FilterReverse = grid.Filter;
            grid.delayFilterReverse = grid.delayFilter;
        else
            c = [cMin:(cMax - cMin)/nFilters:cMax];
            for j = 1:length(c)
                % Compute splines for filter
                splineFilter = spline(grid.tFilter*grid.cFilter/c(j), grid.Filter, grid.tFilter); %c(j)*
                grid.FilterReverse(j, :) = splineFilter; 
                % Compute delay for filter
                delayFilter = find(splineFilter == max(splineFilter)) - 1;
                grid.delayFilterReverse(j) = delayFilter;
            end
        end
        lengthForward = length(grid.tForward);
        grid.tReverse = 0:dt:(lengthFilter + lengthForward - 2)*dt;
    case 3
        error('Not implemented');
end
  
