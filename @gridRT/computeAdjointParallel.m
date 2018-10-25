function source = computeAdjointParallel(grid, source)
% COMPUTEADJOINTPARALLEL computes the reverse propagation of all the given sources parallely
% grid = computeTimeReverseParallel(grid, source)
% INPUTS
% grid: gridRT object that defines the domain
% source: vector of sources
%
% OUTPUTS:
% grid: gridRT object that defines the domain
% 
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke

myCluster = parcluster('local');
nW = myCluster.NumWorkers;
try
    myPool = parpool(nW);
catch exception
end
nS = length(source);
parfor n = 1:nS
    grid.inverse_beam(source(n));
end
