function [ray_index, prop_time] = ssestimation_rays(grid, source_vec, nS)
% SSESTIMATION_RAYS computes and selects those rays from sensor nS that are closest to 
% the vector of sensors source_vec
%function grid = ssestimation_rays(grid, nS, source_vec)
% INPUTS
% grid: gridRT object that defines the domain
% nS: selected source
% source_vec: source to compute the ray trajectories
%
% OUTPUTS:
% grid: gridRT object that defines the domain
% 
% Copyright (C) 2018 Kiko Rul·lan, Marta M. Betcke

% Number of steps of the computation
step = source_vec(nS).step;
[nR, ~, dim] = size(source_vec(nS).x);
nPoints = source_vec(nS).nPoints;

% Allocate rays
source_vec(nS).allocateRays(nR, nPoints, dim);
% Trajectory
grid.forward_trajectories(source_vec(nS));

% Number of sources
lSources = length(source_vec);
% Select origin
source_origin = zeros(lSources, 1, dim);
for i = 1:lSources
    source_origin(i, :, :) = source_vec(i).x0;
end
% Minimum distance
dif_dist = bsxfun(@minus, source_origin, source_vec(nS).x0);
min_dist_vec = sqrt(sum(dif_dist.*dif_dist, 3));
% Ray index
ray_index = ones(lSources, 1);
t_index = ones(lSources, 1);
% Loop over sensors
for i = 1:lSources
    dif_dist = bsxfun(@minus, source_vec(nS).x(:, :, :), source_origin(i, :, :));
    dist = sqrt(sum(dif_dist.*dif_dist, 3));
    min_dist = min(dist(:));
    if (min_dist < min_dist_vec(i)) 
        index = find(dist == min_dist, 1);
        ray_index(i) = mod(index-1, nR) + 1;
        t_index(i) = floor((index-1)/nR);
        min_dist_vec(i) = min_dist;
    end
end
% Select rays
source_vec(nS).selectRays(ray_index, t_index);
% Propagation time
prop_time = source_vec(nS).phi(t_index)';
% Delete aux info 
source_vec(nS).deleteAuxInfo();

% Display message
msg = strcat({'Computed given source with '}, int2str(nPoints), {' steps.'});
disp(msg{1});
