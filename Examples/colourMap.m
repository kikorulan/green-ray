%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Colour map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Normal color map %%%
% Initialize associative array
colourMapV = containers.Map('KeyType', 'double', 'ValueType', 'any');
colourMapV(1) = [1 0 0];
colourMapV(2) = [0 1 0];
colourMapV(3) = [0 0 1];
colourMapV(4) = [0.8 0.3 0.3];
colourMapV(5) = [0.3 0.8 0.3];
colourMapV(6) = [0.3 0.3 0.8];

%%% Rainbow color map %%%
% Initialize associative array
colourMapR = containers.Map('KeyType', 'double', 'ValueType', 'any');
colourMapR(1) = [1 1 0]; % yellow
colourMapR(2) = [1 0 0]; % red
colourMapR(3) = [1 0 1]; % magenta
colourMapR(4) = [0 0 1]; % blue
colourMapR(5) = [0 1 1]; % cyan
colourMapR(6) = [0 1 0]; % green

%%% Blue color map %%%
% Initialize associative array
colourMapBlue = containers.Map('KeyType', 'double', 'ValueType', 'any');
nColours = 15;
vect = 0:1/(nColours-1):1;
for i = 1:nColours
    colourMapBlue(i) = [0 vect(i) 1-vect(i)];
end

%%% Red color map %%%
% Initialize associative array
colourMapRed = containers.Map('KeyType', 'double', 'ValueType', 'any');
nColours = 15;
vect = 0:1/(nColours-1):1;
for i = 1:nColours
    colourMapRed(i) = [1-vect(i) vect(i) 0];
end
clear i

