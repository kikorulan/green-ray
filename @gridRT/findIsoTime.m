function x = findIsoTime(grid, source, nCurves, tMin, tMax)
% FINDISOTIME computes the isotime surface with the given numbers of curves nCurves
%function x = findIsoTime(grid, nS, nCurves)
% INPUTS
% grid: gridRT object that defines the domain
% source: source to compute the pressure
% nCurves: number of iso-time curves to return
% tMin: minimum value for the iso-time curves
% tMax: maximum value for the iso-time curves
%
% OUTPUTS
% x: iso-time curves
% 
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke


% Find the dimensions
[nR lengthRay dim] = size(source.x);
% Loop over the number of curves
tVector = tMin:(tMax-tMin)/nCurves:tMax;
x = nan(nCurves, nR, dim);
for i = 1:length(tVector)-1
    for j = 1:nR
        % Find the position where time is greater than the given isotime curve
        position = find(source.phi(j, :, :) > tVector(i) & source.phi(j, :, :) < tVector(i+1), 1);
        if(isempty(position)); x(i, j, :) = repmat(inf, [1 1 dim]);
        else x(i, j, :) = source.x(j, position, :); end;
    end
end          

