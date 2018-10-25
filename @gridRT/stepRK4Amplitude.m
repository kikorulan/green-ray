function [Dxn, Dpn] = stepRK2Amplitude(grid, source, h, index, Dx, Dp)
% Function for computing the RK2 method for deriving the x and p over initial conditions
    [k1x, k1p] = grid.deriveAmplitude(source, index, Dx,             Dp);
    [k2x, k2p] = grid.deriveAmplitude(source, index, Dx + 0.5*h*k1x, Dp + 0.5*h*k1p);
    [k3x, k3p] = grid.deriveAmplitude(source, index, Dx + 0.5*h*k2x, Dp + 0.5*h*k2p);
    [k4x, k4p] = grid.deriveAmplitude(source, index, Dx + h*k3x,     Dp + h*k3p);

    Dxn = Dx + 1/6*h.*double(k1x + 2*k2x + 2*k3x + k4x);
    Dpn = Dp + 1/6*h.*double(k1p + 2*k2p + 2*k3p + k4p);
end
