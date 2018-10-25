function [Dxn, Dpn] = stepRK2Gauss(grid, source, step, index, Dx, Dp)
% STEPRK2GAUSS Computes the Jacobians Dx, Dp for the given index in the trajectory
% It solves one step of the ODE system for the amplitude.
%function [Dxn, Dpn] = stepRK2Amplitude(grid, nS, tauIndex, Dx, Dp)
% INPUTS
% grid: gridRT object that defines the domain
% source: source to compute the derivative
% index: index to compute the derivative
% h: step size
%
% OUTPUTS:
% Dxn: iterate for Jacobian of x
% Dpn: iterate for Jacobian of p
% 
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke

% Compute RK
%[k1x k1p] = grid.deriveAmplitude(source, index, Dx, Dp);
%[k2x k2p] = grid.deriveAmplitude(source, index, Dx + k1x*h, Dp + k1p*h);
[k1x k1p] = grid.deriveGauss(source, index, Dx, Dp);
[k2x k2p] = grid.deriveGauss(source, index, Dx + k1x*step, Dp + k1p*step);
Dxn = Dx + step*0.5*(k1x + k2x);
Dpn = Dp + step*0.5*(k1p + k2p);
