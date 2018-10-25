function [Dxn, Dpn] = stepRK2Amplitude(grid, source, step, index, Dx, Dp)
% Function for computing the RK2 method for deriving the x and p over initial conditions
    [k1x k1p] = grid.deriveAmplitude(source, index, Dx, Dp);
    [k2x k2p] = grid.deriveAmplitude(source, index, Dx + k1x*step, Dp + k1p*step);
    Dxn = Dx + step*0.5*(k1x + k2x);
    Dpn = Dp + step*0.5*(k1p + k2p);
end
