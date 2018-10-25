% General grid with C values

classdef gridC < handle
    properties (GetAccess = public, SetAccess = private)
        % grid dimensions in grid points
        Nx = 1;
        Ny = 1;
        Nz = 1;
        
        % grid point spacing in metres
      	dx = 0;
        dy = 0;
        dz = 0;
        
        % set dimension
        dim = 0;
        
        % set axis
        xAxis = [];
        yAxis = [];
        zAxis = [];

        % Sound speed
        c = []; % sound speed
        n = []; % eta (inverse of sound speed)
        Gn = []; % gradient of eta

        % Phase
    	phase = [];
        % u - pressure value
        u = [];
    end
    %========================================
    % Constructor function
    %========================================
    methods
        function grid = gridC(varargin)
            % assign the input values to the grid object
            switch nargin
                case 4
                    % 2D uniform grid
                    grid.Nx = varargin{1};
                    grid.dx = varargin{2};
                    grid.Ny = varargin{3};
                    grid.dy = varargin{4};

                    % set the number of dimensions
                    grid.dim = 2;                    

                    % set axis
                    grid.xAxis = 0:grid.dx:(grid.Nx-1)*grid.dx;
                    grid.yAxis = 0:grid.dy:(grid.Ny-1)*grid.dy;
                        
                case 6
                    % 3D uniform grid
                    grid.Nx = varargin{1};
                    grid.dx = varargin{2};
                    grid.Ny = varargin{3};
                    grid.dy = varargin{4};
                    grid.Nz = varargin{5};
                    grid.dz = varargin{6};
                    % set the number of dimensions
                    grid.dim = 3;

               otherwise
                    error(num2str(nargin));
            end
            grid.phase = inf(grid.Nx, grid.Ny, grid.Nz);
            grid.c = ones(grid.Nx, grid.Ny, grid.Nz);
            grid.u = zeros(grid.Nx, grid.Ny, grid.Nz);
        end
    end

    %========================================
    % Methods for modifying attributes
    %========================================
    methods (Access = private)
    	% Sound speed %%%%%%%%%%%%%%%%%%%%%%%
        % Obtain a particular c value
        c = getC(grid, point);
    	% Sound speed inverse %%%%%%%%%%%%%%%%%%%%%%%
        % Obtain a particular n value
        n = getN(grid, point);
        % Get eta and gradient of eta
        [n, Gn] = getNGN(grid, point);
        % Pressure %%%%%%%%%%%%%%%%%%%%%%%%%%
        % Obtain a particular u value
        U = getU(grid, point);
    end
    methods
    	% Sound speed %%%%%%%%%%%%%%%%%%%%%%%
        % Modify a particular c index
        grid = setC(grid, point, cVal);
        % Modify c
        function grid = setCMatrix(grid, v)
            if (~isequal(size(v), size(grid.c)))
                error('Dimensions for given v and grid.c do not match');
            else
                grid.c = v;
                grid.n = 1./grid.c;
            end
            grid.setGNMatrix();
        end

        % Modify a particular n index
        grid = setN(grid, point, nVal);
        % Modify n
        function grid = setNMatrix(grid, v)
            if (~isequal(size(v), size(grid.n)))
                error('Dimensions for given v and grid.n do not match');
            else
                grid.n = v;
                grid.c = 1./grid.n;
            end
            grid.setGNMatrix();
        end
        % Create nGn matrix
        grid = setGNMatrix(grid);

        % Phase %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Obtain a particular phase value
        phase = getPhase(grid, point);
        % Modify a particular phase index
        grid = setPhase(grid, point, phaseVal);
    	% Modify phase
        function grid = setPhaseMatrix(grid, v)
            if (~isequal(size(v), size(grid.phase)))
                error('Dimensions for given v and grid.phase do not match');
            else
                grid.phase = v;
            end
        end

        % Pressure %%%%%%%%%%%%%%%%%%%%%%%%%%
        % Modify a particular u index
        grid = setU(grid, point, uVal);
    	% Modify u
        function grid = setUMatrix(grid, v)
            if (~isequal(size(v), size(grid.u)))
                error('Dimensions for given v and grid.u do not match');
            else
                grid.u = v;
            end
        end
    end

    %========================================
    % Methods for finding coordinates
    %========================================
    methods (Access = private)
        % Find coordinates
        coord = findCoordinates(grid, point);
        % Find neighbours
        neigh = findNeighbours(grid, point);
    end
    %========================================
    % Methods for computing gradients and hessians
    %========================================
    methods (Access = private)
        % Gradient of the sound speed
        grad = gradSpeed(grid, point);
        % Hessian of the sound speed
        hess = hessSpeed(grid, point);
        % Gradient of the sound speed inverse
        grad = gradEta(grid, point);
    	% Gradient of the phase
        grad = gradPhase(grid, point);
        % Hessian of the phase
        hess = hessPhase(grid, point);
    end

end
