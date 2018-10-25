function [Xn, Pn] = stepRK4Hamil(grid, h, X, P)
% STEPRK4HAMIL Computes the Runge-Kutta step for the given ray point using step h
%function Xn = stepRK4Hamil(grid, tau, X, h, nS)
% INPUTS
% grid: gridRT object that defines the domain
% h: step size
% X: current x point in the hamiltonian trajectory
% P: current p point in the hamiltonian trajectory
%
% OUTPUTS:
% Xn: next x point in the ray trajectory
% Pn: next p point in the ray trajectory
% 
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke

[k1X, k1P] = grid.deriveHamil(X,             P            );
[k2X, k2P] = grid.deriveHamil(X + 0.5*h*k1X, P + 0.5*h*k1P);
[k3X, k3P] = grid.deriveHamil(X + 0.5*h*k2X, P + 0.5*h*k2P);
[k4X, k4P] = grid.deriveHamil(X + h*k3X,     P + h*k3P    );

% Compute RK4
Xn = X + 1/6*h.*double(k1X + 2*k2X + 2*k3X + k4X);
Pn = P + 1/6*h.*double(k1P + 2*k2P + 2*k3P + k4P);
