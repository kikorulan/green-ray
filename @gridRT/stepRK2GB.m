function [qn, un] = stepRK2GB(grid, source, step, index, q, u)
% STEPRK2GB Computes the Jacobians Dx, Dp for the given index in the trajectory
% It solves one step of the ODE system for the amplitude.
%function [Dxn, Dpn] = stepRK2Amplitude(grid, nS, tauIndex, Dx, Dp)
% INPUTS
% grid: gridRT object that defines the domain
% source: source to compute the derivative
% index: index to compute the derivative
% h: step size
%
% OUTPUTS:
% qn: iterate for q
% un: iterate for pressure
% 
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke

% Compute RK
[k1q k1u] = grid.deriveGB(source, index, q, u);
[k2q k2u] = grid.deriveGB(source, index, q + k1q*step, u + k1u*step);
qn = q + step*0.5*(k1q + k2q);
un = u + step*0.5*(k1u + k2u);
