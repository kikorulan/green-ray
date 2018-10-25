function ssestimation_length(grid, source_vec, nS)
% SSESTIMATION_LENGTH computes the hamiltonian trajectories of the given source nS
%function grid = ssestimation_length(grid, source)
% INPUTS
% grid: gridRT object that defines the domain
% source: source to compute the ray trajectories
%
% OUTPUTS:
% grid: gridRT object that defines the domain
% 
% Copyright (C) 2018 Kiko Rul·lan, Marta M. Betcke

% Number of steps of the computation
[nR, nPoints, dim] = size(source_vec(nS).x);
nSources = length(source_vec);

if nS == 1
    grid.LMatrix(:) = 0;
end

% Run index
for index = 1:nPoints-1
    % Coordinates
    P1 = source_vec(nS).x(:, index, :);
    P2 = source_vec(nS).x(:, index+1, :);
    c1 = grid.findCoordinates(P1)-1;
    c2 = grid.findCoordinates(P2)-1;
    coordFin = find(c2(:, 1, 1) < inf);
    % Subset of rays - with finite coordinates
    P1 = P1(coordFin, :, :);
    P2 = P2(coordFin, :, :);
    c1 = c1(coordFin, :, :);
    c2 = c2(coordFin, :, :);
    % Triangle parameters
    B = abs(P1(:, :, 1) - P2(:, :, 1));
    H = abs(P1(:, :, 2) - P2(:, :, 2));
    L = sqrt(B.*B + H.*H);
    b1 = abs(grid.dx/2*(c1(:, :, 1)+c2(:, :, 1)) - P1(:, :, 1));
    b2 = abs(grid.dx/2*(c1(:, :, 1)+c2(:, :, 1)) - P2(:, :, 1));
    h1 = abs(grid.dx/2*(c1(:, :, 2)+c2(:, :, 2)) - P1(:, :, 2));
    h2 = abs(grid.dx/2*(c1(:, :, 2)+c2(:, :, 2)) - P2(:, :, 2));
    cross_bh_1 = b1.*h2;
    cross_bh_2 = b2.*h1;
    % Linear coordinates inside matrix
    coordL1 = (c1(:, :, 1) + c1(:, :, 2)*grid.Nx)*nR*nSources + (nS-1)*nR + coordFin;
    coordL2 = (c2(:, :, 1) + c2(:, :, 2)*grid.Nx)*nR*nSources + (nS-1)*nR + coordFin;
    coordL3 = ((c2(:, :, 1) + c1(:, :, 2)*grid.Nx)*nR*nSources + (nS-1)*nR + coordFin).*(cross_bh_1 <  cross_bh_2) + ...
              ((c1(:, :, 1) + c2(:, :, 2)*grid.Nx)*nR*nSources + (nS-1)*nR + coordFin).*(cross_bh_1 >= cross_bh_2);
    % Distiguish between 5 cases
    case1 = ((c1(:, :, 1) == c2(:, :, 1)) & (c1(:, :, 2) == c2(:, :, 2)));                               % x1 == x2 & y1 == y2
    case2 = ((c1(:, :, 1) ~= c2(:, :, 1)) & (c1(:, :, 2) == c2(:, :, 2)));                               % x1 ~= x2 & y1 == y2
    case3 = ((c1(:, :, 1) == c2(:, :, 1)) & (c1(:, :, 2) ~= c2(:, :, 2)));                               % x1 == x2 & y1 ~= y2
    case4 = ((c1(:, :, 1) ~= c2(:, :, 1)) & (c1(:, :, 2) ~= c2(:, :, 2)) & (cross_bh_1 <  cross_bh_2));  % x1 ~= x2 & y1 ~= y2 & b1*h2 <  b2*h1
    case5 = ((c1(:, :, 1) ~= c2(:, :, 1)) & (c1(:, :, 2) ~= c2(:, :, 2)) & (cross_bh_1 >= cross_bh_2));  % x1 ~= x2 & y1 ~= y2 & b1*h2 >= b2*h1
    % Compute lengths
    l1 = zeros(length(coordFin), 1);
    l2 = zeros(length(coordFin), 1);
    l3 = zeros(length(coordFin), 1);
    % Case 1
    l1(case1) = L(case1);
    % Case 2
    l1(case2) = b1(case2).*L(case2)./B(case2);
    l2(case2) = b2(case2).*L(case2)./B(case2);
    % Case 3
    l1(case3) = h1(case3).*L(case3)./H(case3);
    l2(case3) = h2(case3).*L(case3)./H(case3);
    % Case 4
    l1(case4) = b1(case4).*L(case4)./B(case4);
    l2(case4) = h2(case4).*L(case4)./H(case4);
    l3(case4) = L(case4) - l1(case4) - l2(case4);
    % Case 5
    l1(case5) = h1(case5).*L(case5)./H(case5);
    l2(case5) = b2(case5).*L(case5)./B(case5);
    l3(case5) = L(case5) - l1(case5) - l2(case5);
    % Modify matrix
    grid.LMatrix(coordL1) =  grid.LMatrix(coordL1) + l1; 
    grid.LMatrix(coordL2) =  grid.LMatrix(coordL2) + l2; 
    grid.LMatrix(coordL3) =  grid.LMatrix(coordL3) + l3; 


    % OUTPUT
    %%  if max(L) > 1.1*grid.c(1,1)*source_vec(nS).step
    %%      maxL = max(L)
    %%  end
    %%  
    %%  if max(l1) > 1.1*grid.c(1,1)*source_vec(nS).step
    %%      maxL1 = max(l1)
    %%  end
    %%  if max(l2) > 1.1*grid.c(1,1)*source_vec(nS).step
    %%      maxL2 = max(l2)
    %%  end
    %%  if max(l3) > 1.1*grid.c(1,1)*source_vec(nS).step
    %%      maxL3 = max(l3)
    %%  end
end


% Display message
%msg = strcat({'Computed given source with '}, int2str(nPoints), {' steps.'});
%disp(msg{1});
