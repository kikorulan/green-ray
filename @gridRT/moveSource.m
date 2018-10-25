function grid = moveSource(grid, gridOriginal, nS, typeMove)
% MOVESOURCE moves the source nS from gridOriginal to grid
%grid = moveSource(grid, gridOriginal)
% INPUTS
% grid: gridRT object to copy the source
% gridOriginal: gridRT object to find the source to copy
% nS: number of source to copy
% typeMove: choose between 'add', 'replace' and 'append'
%
% OUTPUTS:
% grid: gridRT object 
% 
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke
switch typeMove
    case 'add'
        grid.source = [grid.source(1:nS-1) gridOriginal.source(nS) grid.source(nS:end)];
    case 'replace'
        grid.source(nS) = gridOriginal.source(nS);
    case 'append'
        nSources = length(grid.source);
        grid.source(nSources+1) = gridOriginal.source(nS);
    otherwise
        error('Wrong argument for typeMove');
end
%grid.deleteSource(nS + 1);
