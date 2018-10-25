function A = forward_amplitude(grid, source)
% FORWARD_AMPLITUDE Computes the amplitude at a given point x using the proximal ray method
% from source point x0.
% A = forward_amplitude(grid, source)
% INPUTS
% grid: gridRT object that defines the domain
% source: source to compute the amplitude
%
% OUTPUTS:
% A: amplitude
% 
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke

%%  % Display message
%%  msg = strcat({'Computing amplitude for source '}, int2str(nS), {'...'});
%%  disp(msg{1});

% Number of rays
[nR nPoints dim] = size(source.x);

switch grid.dim % Switch over the dimension
    case 2
        x = source.x;
        % For computing the determinant
        difPhi1 = x(1:end-1, 1:end-1, :) - source.x(2:end, 1:end-1, :);
        difPhi = cat(1, difPhi1, x(end-1, 1:end-1, :) - x(end, 1:end-1, :));
        difTau = x(:, 2:end, :) - x(:, 1:end-1, :);
        flipDifTau(:, :, 2) = -difTau(:, :, 1);
        flipDifTau(:, :, 1) = difTau(:, :, 2); 
        %%%% Amplitude & Q
        DxVector(:, :, 3:4) = difTau;
        DxVector(:, :, 1:2) = difPhi;
        source.insertDXVector(DxVector);
        q = sum(difPhi.*flipDifTau, 3);
        q = [q repmat(NaN, [nR 1 1])];
        %A = sqrt(abs(repmat(q(:, 1, :), [1 nPoints 1])./q)); % All q(1) are the same for the nR rays
        %A = repmat(source.n(:, 1, :), [1 nPoints 1])./source.n.* ...
        %    sqrt(abs(repmat(q(:, 1, :), [1 nPoints 1])./q))./source.n./source.n; % All q(1) are the same for the nR rays
        A = repmat(source.n(:, 1, :), [1 nPoints 1])./source.n.* ...
            sqrt(abs(repmat(q(:, 1, :), [1 nPoints 1])./q)); % All q(1) are the same for the nR rays
        %q(isnan(q) | isinf(q)) = 0; % substitute NaN values
        A(isnan(A)) = 0; % substitute NaN values
        A(A > 0.1) = 0.1;
        source.insertQVector(q);
        source.insertAVector(A);
    case 3 %%%%%%%%%%%%%% NOT IMPLEMENTED %%%%%%%%%%%%%
        error('Compute amplitude not implemented');
end



