function [Xn, Pn] = stepRK2Hamil(grid, h, X, P)
% STEPRK2HAMIL Computes the Runge-Kutta step for the given ray point using step h
%function Xn = stepRK2Hamil(grid, h, X, P)
% INPUTS
% grid: gridRT object that defines the domain
% h: step size
% X: current x point in the hamiltonian trajectory
% P: current p point in the hamiltonian trajectory
% source: source to compute the ray step
%
% OUTPUTS:
% Xn: next x point in the ray trajectory
% Pn: next p point in the ray trajectory
% 
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke

[k1X, k1P] = grid.deriveHamil(X,         P        );
[k2X, k2P] = grid.deriveHamil(X + h*k1X, P + h*k1P);

% Compute RK4
Xn = X + 1/2*h.*double(k1X + k2X);
Pn = P + 1/2*h.*double(k1P + k2P);
