%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEPRECATED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s0 = flipTrajectories(s0, deltaX, deltaP)
% This function flips the trajectories x of the given source s0, stored in xR
    [dim, lengthX, nR] = size(s0.x);
    flipX = fliplr(s0.x);
    flipP = fliplr(s0.p);
    % Initialization
    s0.xR = [];
    s0.pR = [];
    % Phi - time
    %maxPhi = max(max(s0.phi(s0.phi<inf)));
    %flipPhi = -fliplr(s0.phi) + maxPhi;
    %flipPhi(flipPhi == -inf) = inf;
    %s0.phi = [];
    for n = 1:nR
        s0.xR(:, :, n) = [flipX(:, flipX(1, :, n) < inf, n) flipX(:, flipX(1, :, n) == inf, n)];
        s0.pR(:, :, n) = [-flipP(:, flipX(1, :, n) < inf, n) flipP(:, flipX(1, :, n) == inf, n)];
        %s0.phi(:, :, n) = [flipPhi(:, flipPhi(1, :, n) < inf, n) flipPhi(:, flipPhi(1, :, n) == inf, n)];
    end
    xR = s0.xR(:, 1, :);
    pR = s0.pR(:, 1, :);
    % Phi and Theta rays
    switch dim
        case 1
        case 2
            % Phi Ray
            npR = sqrt(sum(pR.*pR, 1));
            s0.pPhi = [pR(1, 1, :)*cos(deltaP) - pR(2, 1, :)*sin(deltaP); ...
                        pR(2, 1, :)*cos(deltaP) + pR(1, 1, :)*sin(deltaP)];
            s0.xPhi = xR - deltaX*pR./repmat(npR, [dim 1 1]) + deltaX*s0.pPhi./repmat(npR, [dim 1 1]);
        case 3 %%%%%%%%%%%%%% NOT IMPLEMENTED %%%%%%%%%%%%%%%%%%%%%%%%%
            error('Not implemented');
   end
end 
        
