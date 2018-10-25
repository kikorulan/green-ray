function hess = hessPhaseODE(grid, source, typeSource)
% HESSPHASEODE computes the hessian of the phase either at the start or the end
% of the trajectories at the given source
%function hess = hessPhaseODE(grid, nS, point)
% INPUTS
% grid: gridRT object that defines the domain
% source: source to compute the hessian
% typeSource: choose between real and imag
%
% OUTPUTS
% hess: hessian
% 
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke

sourceCoord = grid.findCoordinates(source.x0);
c = grid.getC(sourceCoord);
n = source.n(1, 1, 1);
delta = grid.deltaX;
% Compute the coordinates of the source with respect to the shooting point
x1 = permute(delta*source.p(:, 1, 1)/n/n, [3 2 1]);
x2 = permute(delta*source.p(:, 1, 2)/n/n, [3 2 1]);

% Compute the hessian
denom = (x1.^2 + x2.^2).^(1.5);
hess_xx = (x2.^2)./denom./c; %./c;
hess_yy = (x1.^2)./denom./c;
hess_xy = (-x1.*x2)./denom./c;

%%  hess_xx = (x2.^2)./denom./c./c; %./c;
%%  hess_yy = (x1.^2)./denom./c./c;
%%  hess_xy = (-x1.*x2)./denom./c./c;

hess = [hess_xx hess_xy; hess_xy hess_yy];



% Compute the imaginary part
hess_xx_imag = (x2.^2); %./c;
hess_yy_imag = (x1.^2);
hess_xy_imag = (-x1.*x2);
hess_imag = [hess_xx_imag hess_xy_imag; hess_xy_imag hess_yy_imag];

% Rescale hessian
nR = size(source.x, 1);
scaleFactor = permute(10.^(0:log10(nR)/(nR-1):log10(nR)), [1 3 2]);
%scaleFactor = permute(1:1:nR, [1 3 2]);
scaleMatrix = repmat(scaleFactor, [2 2 1]);

% Compute imaginary part
if (typeSource == 'imag')
    %hess = (hess + i*hess/100);
    hess = i*hess;
    %hess = (hess/1e2 + i*hess/200);
end


%%  % Alternative initialisation
%%  hess_xx(:) = 1;
%%  hess_xy(:) = 0;
%%  hess_yy(:) = 1;
%%  hess = 1e-5*[hess_xx hess_xy; hess_xy hess_yy];

