function source = computeForwardParallel(grid, x0, angleMin, angleMax, nRays, tStep, tMax, mode, deleteRays)
% COMPUTEFORWARDPARALLEL computes the forward propagation of all the given sources parallely
% grid = computeForwardParallel(grid, numWorkers)
% INPUTS
% grid: gridRT object that defines the domain
% x0: shooting points for a vector of sources
% angleMin: minimum angle for shooting the rays
% angleMax: maximum angle for shooting the rays
% nRays: number of rays shot
% tStep: step taken to compute the Hamiltonian
% tMax: maximum value for the parameter
% mode: choose between amplitude computation modes
% deleteRays: boolean that indicates if we wish to delete the rays after computing the trajectories
%
% OUTPUTS:
% source: vector containing the computed sources
% 
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke

% Create workers
myCluster = parcluster('local');
nW = myCluster.NumWorkers;
try
    myPool = parpool(nW);
catch exception
end

% Parallelise multiple
nS = length(x0);
parforSize = 3*nW;
numParfor = ceil(nS/parforSize);
for j = 1:numParfor
    vectorParfor = (1+parforSize*(j-1)):min(parforSize*j, nS);
    parfor n = vectorParfor
        source(n) = grid.newSource(x0{n}, angleMin, angleMax, nRays, tStep, tMax);
        grid.computeHamil(source(n), mode);
        source(n).deleteAuxInfo();
        if deleteRays
            source(n).deleteRays();
        end
    end
end

% Parallelise execution
%%  nS = length(x0);
%%  parfor n = 1:nS
%%      msg = strcat({'Source '}, int2str(n), {'...'});
%%      disp(msg{1});
%%      source(n) = grid.newSource(x0{n}, angleMin, angleMax, nRays, tStep, tMax);
%%      grid.computeHamil(source(n));
%%      source(n).deleteAuxInfo();
%%      if deleteRays
%%          source(n).deleteRays();
%%      end
%%  end

