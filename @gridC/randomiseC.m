function  grid = randomiseC(grid, spatialSDV, inc)
% RANDOMISEC randomises the speed of sound of the given grid.
%function  grid = randomiseC(grid, spatialSDV, inc)
% STEPS
% - Take the mean value of the current sound speed
% - Generate a gaussian random spatial vector
% - Sum the Gaussian random vector with the mean sound speed
% 
% INPUTS
% grid: gridRT object that defines the domain
% spatialSDV: the spatial standard deviation of the gaussian kernel as a percentage of the size grid.Nx
% inc: maximum increase of the speed of sound 
% 
% OUTPUTS
% grid: gridRT object that defines the domain
%
% Copyright (C) 2017 Kiko Rul·lan, Felix Lucka, Marta M. Betcke

% Get the current state of the random generator (will be set to that state
% after constructing the noIC setting)
currRandState = rng;
% Obtain mean sound speed
meanSoundSpeed = mean(grid.c(:));

% Set the random generator to the state specified in the noise model of the
% input setting
rng(floor(100*inc)+2);

soundSpeed = meanSoundSpeed*ones(grid.Nx, grid.Ny, grid.Nz);
maxSpeedVar = inc*meanSoundSpeed;
% Generate vectors for construction
stdGaussian = spatialSDV*grid.Nx;
xVector = repmat(permute(-(grid.Nx-1)/2:1:(grid.Nx-1)/2, [2 1]), [1 grid.Ny grid.Nz]);
yVector = repmat(-(grid.Ny-1)/2:1:(grid.Ny-1)/2, [grid.Nx 1 grid.Nz]);
zVector = repmat(permute(-(grid.Nz-1)/2:1:(grid.Nz-1)/2, [1 3 2]), [grid.Nx grid.Ny 1]);
% Generate 3D gaussian
Gaussian3D = exp(- 1/stdGaussian^2 * (xVector.^2 + yVector.^2 + zVector.^2));
padSize = 20;
Gaussian3DFFT = fftn(padarray(Gaussian3D,[padSize,padSize,padSize],0));
mrfSpeedVariationFFT = Gaussian3DFFT .* fftn(padarray(randn(grid.Nx, grid.Ny, grid.Nz),[padSize,padSize,padSize],0));
mrfSpeedVariation = real(ifftshift(ifftn(mrfSpeedVariationFFT)));
clear mrfSpeedVariationFFT;
% Correct size, maximum and mean
mrfSpeedVariation = mrfSpeedVariation(padSize+1:end-padSize,padSize+1:end-padSize,padSize+1:end-padSize);
mrfSpeedVariation = mrfSpeedVariation/max(mrfSpeedVariation(:));
mrfSpeedVariation = mrfSpeedVariation - mean(mrfSpeedVariation(:));

soundSpeed = soundSpeed + maxSpeedVar * (mrfSpeedVariation);
plot_cube(mrfSpeedVariation);

clear mrfSpeedVariation;
grid.setCMatrix(soundSpeed);
% set the internal state of the random generator to what it was before
rng(currRandState);
