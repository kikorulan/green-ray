%==================================================================================================
% GRIDRT is a high frequency approximation method for PAT 
%==================================================================================================
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke
classdef gridRT < gridC
    %======================================================================
    % Properties
    %======================================================================
    properties (GetAccess = public, SetAccess = private)
        % Source
        %source = Source.empty();    % Vector of sources
        dt = 0;                      % Global time step
        tForward = [];               % Global time signal for the forward propagation
        tReverse = [];               % Global time signal for the reverse propagation
        tMaxBeam = 0;                % Maximum time duration for all beams
        % Delta X                   
        deltaX = 0;                  % for computing the amplitude auxiliary rays
        % Delta P                   
        deltaP = 0;                  % for computing the amplitude auxiliary rays

        % Fundamental filter characteristics
        cFilter = 0;        % Reference sound speed for the filter
        delayFilter = 0;    % Delay of non-causal filter
        tFilter = [];       % Time signal
        Filter = [];        % Filter values
        %%% Sound speed estimation
        LMatrix = [];
    end
    properties (GetAccess = public, SetAccess = public)
        % Matrices for pixel-by-pixel pressure
        pixelABeamFFT = []      % FFT signal for the beam 
        pixelFilterFFT = [];    % FFT of the filter at each pixel
        pixelDelayFilter = [];  % Delay for the filter at each pixel
        tMaxReverse = 0;

        % For reconstruction
        tFilterReverse = [];
        FilterReverse = [];
        delayFilterReverse = [];
        pixelAReverse = [];
        pixelTReverse = [];
        pixelAPropagation = [];
    end
    %======================================================================
    % Constructor function
    %======================================================================
    methods
        function grid = gridRT(varargin)
            grid = grid@gridC(varargin{:});
            % Delta Values
            grid.deltaX = grid.dx*1e-4;
            grid.deltaP = 1e-2;
        end
    end

    %======================================================================
    % Modifying attributes - created for this class
    %======================================================================
    methods (Access = public)
        %%%% Grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Set value to time step
        function grid = setTime(grid, dt, tForward); 
            grid.dt = dt; 
            grid.tForward = 0:dt:tForward;
        end;
        % Set value to deltaX & deltaP
        function grid = setDeltaX(grid, deltaX); grid.deltaX = deltaX; end;
        function grid = setDeltaP(grid, deltaP); grid.deltaP = deltaP; end;   
        % Set values to fundamental filter characteristics
        function grid = setCFilter(grid, c); grid.cFilter = c; end;
        function grid = setDelayFilter(grid, delay); grid.delayFilter = delay; end;
        function grid = setTFilter(grid, t); grid.tFilter = t; end;
        function grid = setFilter(grid, Filter); grid.Filter = Filter; end;
        
    end

    %======================================================================
    % Modify attributes - Extended from Grid C & Implemented for this class
    %======================================================================
    methods (Access = public)
        % Sound speed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Obtain a particular c value
        c = getC(grid, point);
        % Sound speed inverse %%%%%%%%%%%%%%%%%%%%%%%
        % Obtain a particular n value
        n = getN(grid, point);
        % Obtain a particular n value and its gradient
        [n, Gn] = getNGN(grid, point);
        % Sound speed inverse %%%%%%%%%%%%%%%%%%%%%%%
        % Obtain a particular u value
        u = getU(grid, point);
    end
    methods  % Extended from class gridC
        % Sound speed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Modify c
        function grid = setCMatrix(grid, v); grid = grid.setCMatrix@gridC(v); end;
        % Modify a particular c index
        function grid = setC(grid, point, cVal); grid = grid.setC@gridC(point, cVal); end;
        % Randomise the sound speed
        function grid = randomiseC(grid, spatialSDV, inc); grid = grid.randomiseC@gridC(spatialSDV, inc); end;

        % Sound speed inverse %%%%%%%%%%%%%%%%%%%%%%%
        % Modify n
        function grid = setNMatrix(grid, v); grid = grid.setNMatrix@gridC(v); end;
        % Modify a particular n index
        function grid = setN(grid, point, nVal); grid = grid.setN@gridC(point, nVal); end;
        
        % Sound speed inverse %%%%%%%%%%%%%%%%%%%%%%%
        % Modify u
        function grid = setUMatrix(grid, v); grid = grid.setUMatrix@gridC(v); end;
        % Modify a particular u index
        function grid = setU(grid, point, uVal); grid = grid.setU@gridC(point, uVal); end;
    end
    
    %======================================================================
    % Find Coordinates in the grid
    %======================================================================
    methods (Access = private)
        % Find Coordinates
        coord = findCoordinates(grid, point);
        % Find Neighbours
        neigh = findNeighbours(grid, point);
    end

    %======================================================================
    % Compute gradients and hessians
    %======================================================================
    methods (Access = private)
        % Gradient of the phase
        function grad = gradPhase(grid, point); grad = grid.gradPhase@gridC(point); end;
        % Hessian of the phase
        function hess = hessPhase(grid, point); hess = grid.hessPhase@gridC(point); end;
        % Hessian of the phase for the ODE method
        hess = hessPhaseODE(grid, nS, typeSource);
        % Gradient of the speed
        grad = gradSpeed(grid, point);
        % Hessian of the speed
        hess = hessSpeed(grid, point);
        % Gradient of eta - based on coordinates
        grad = gradEta(grid, point);
        % Gradient of eta - bilinear interpolation
        n = interpolateEta(grid, point);
        % Gradient of eta - bilinear interpolation
        grad = interpolateGradEta(grid, point);
    end

    %======================================================================
    % Related to the source
    %======================================================================
    methods 
        % Create new source
        source = newSource(grid, x0, angleMin, angleMax, nR, step, tauMax);
        % Compute the beam for the given source
        grid = forward_beam(grid, source); 
        % Compute the time signal for the given source
        grid = forward_timeSignal(grid, source);
    end
    %======================================================================
    % Compute the trajectories
    %======================================================================
    methods
        %%% Main Ray t-Hamiltonian
        grid = forward_trajectories(grid, source);
        % Rays, pressure and time signal
        grid = computeHamil(grid, source, mode, del);
    end
    %======================================================================
    % Parallelisation
    %======================================================================
    methods
        % Move source from the original grid to the objective grid 
        grid = moveSource(grid, gridOriginal, nS, typeMove);
        % Copy the grid
        grid = copyGrid(gridOriginal);
        % Compute the forward problem for all sources parallely
        grid = computeForwardParallel(grid, x0, angleMin, angleMax, nRays, tStep, tMax, mode, deleteRays); 
        % Compute the adjoint problem for all sources parallely
        source = computeAdjointParallel(grid, source);
    end
    %======================================================================
    % Compute amplitude & pressure
    %======================================================================
    methods (Access = private)
        % Cpmpute the amplitude - proximal ray method
        A = forward_amplitude(grid, nS);
        % Cpmpute the amplitude - ODE method
        A = forward_amplitudeODE(grid, nS);
        % Cpmpute the amplitude - Gauss beam method
        A = forward_amplitudeGauss(grid, nS);
        % Cpmpute the amplitude - Gauss beam method
        A = forward_amplitudeGB(grid, nS);
        % Compute the reverse amplitude
        revA = forward_revAmplitude(grid, nS);
        % Compute the pressure
        grid = forward_pressure(grid, nS);
        % Convert the attenuation, distance and time info to grid format
        grid = forward_raysToGrid(grid, nS);
    end

    %======================================================================
    % Runge-Kutta
    %======================================================================
    methods (Access = private)
        %%%% Main Ray t-Hamiltonian
        % Derivative of the ray
        [DX, DP] = deriveHamil(grid, X, P);
        % Runge-Kutta step for the hamiltonian
        [Xn, Pn] = stepRK4Hamil(grid, h, X, P);
        [Xn, Pn] = stepRK2Hamil(grid, h, X, P);

        %%%% Amplitude ODE method
        % Derivative of Dx and Dp
        [D2x, D2p] = deriveGauss(grid, source, index, Dx, Dp);
        % Runge-Kutta step for the amplitude
        [Dxn, Dpn] = stepRK2Gauss(grid, source, step, index, Dx, Dp);
        % Derivative for the amplitude
        [Dxn, Dpn] = deriveAmplitude(grid, source, index, Dx, Dp);
        % Runge-Kutta for the amplitude
        [Dxn, Dpn] = stepRK2Amplitude(grid, source, step, index, Dx, Dp);
        [Dxn, Dpn] = stepRK4Amplitude(grid, source, step, index, Dx, Dp);

        % Derivative of q and u
        [Dq, Du] = deriveGB(grid, source, index, q, u);
        % Runge-Kutta step for the amplitude
        [qn, un] = stepRK2GB(grid, source, step, index, q, u);

        %%  %%%% Gaussian beams
        %%  % Derivative for the Riccati equation
        %%  [DY, DN] = deriveGauss(grid, nS, index, Y, N);
        %%  % Derivative for the amplitude equation
        %%  DA = deriveGaussAmplitude(grid, nS, index, A);
        %%  % Runge-Kutta step for the gaussian beams
        %%  [Yn, Nn] = stepRK4Gauss(grid, nS, index, h);
        %%  % Runge-Kutta step for the amplitude
        %%  An = stepRK4GaussAmplitude(grid, nS, index, h);

    end
    %======================================================================
    % Time propagation
    %======================================================================
    methods
        % Compute the impulse response for the given domain
        grid = impulse_additive(grid, sourceType);
        % Compute the impulse response for a dirichlet source for the given domain
        grid = impulse_dirichlet(grid);
        % Compute the time propagation on the given sensors
        grid = inverse_timeSource(grid, source)
        % Compute the time propagation on the given sensors using the filter at each pixel
        grid = inverse_timeSource_fft(grid, source);
        % Compute the filters
        grid = inverse_filter(grid, nFilters);
        % Compute the filter at each pixel
        grid = inverse_filter_fft(grid, nFilters);
        % Compute the reverse beam
        grid = inverse_beam(grid, source);
        % Compute the reverse beam considering the other sensors
        aReverse = inverse_beam_adjoint(grid, source);
        % Compute the reverse signal
        grid = inverse_signal(grid, source);
    end
    %======================================================================
    % Iterative reconstruction
    %======================================================================
    methods
        % Forward operator
        forward_data = operator_forward(grid, source, index, initial_pressure, norm_factor);
        % Inverse operator
        initial_pressure = operator_inverse(grid, source, index, forward_data);
        % Compute an iteration of the iterative reconstruction
        pixelAReverse = iterative_recon(grid, source, step, nIter);
        % Compute an iteration of the iterative reconstruction for a particular sensor
        pixelAReverse = iterative_recon_sensor(grid, source);
    end

    %======================================================================
    % Other
    %======================================================================
    methods
        % Compute an isotime surface with the number of given curves
        x = findIsoTime(grid, nS, nCurves, tMin, tMax);
    end
    %======================================================================
    % Gaussian Beams
    %======================================================================
    methods
        % Compute subrays
        grid = forward_subrays(grid, source, nRsub, dx);
        % Compute the rays
        grid = gauss_beam(grid, nS);
    end
    %======================================================================
    % Plot
    %======================================================================
    methods
        % Plot the rays of the given vector of sources
        h = plot_soundSpeed(grid);
        % Plot the pressure at a particular time stamp
        h = plot_pressure(grid, timeStamp);
    end
    %======================================================================
    % Sound speed estimation
    %======================================================================
    methods
        % Allocate L matrix
        function ssestimation_allocate(grid, source)
            lSource = length(source);
            grid.LMatrix = zeros(lSource*lSource, grid.Nx*grid.Ny);
        end
        % Select the rays that connect the selected source with the others
        [ray_index, prop_time] = ssestimation_rays(grid, source_vec, nS);
        % Assign lengths in matrix
        ssestimation_length(grid, source_vec, nS);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DEPRECATED                                                                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      %======================================================================
%%      % Compute the trajectories
%%      %======================================================================
%%      methods
%%          %%% Main Ray
%%          % Compute ray step
%%          grid = stepTrajec(grid, nS, tauIndex, step);
%%          % Compute global ray
%%          grid = computeTrajec(grid, nS);
%%          % Compute global ray as well as Q and Amplitude
%%          grid = computeTrajecQ(grid, nS, thresholdR, thresholdQ);
%%          % Compute global ray using Snell's law
%%          grid = computeTrajecSnell(grid, nS);
%%          % Find the end points of a source
%%          [x, angle] = findTrajectoryEnd(grid, nS);
%%      end
%%  
%%      %======================================================================
%%      % Compute amplitude & pressure
%%      %======================================================================
%%      methods
%%          % Proximal Ray Method ==================================
%%          %%% Phi and Theta Rays
%%          % Compute source step
%%          grid = stepTrajecPhi(grid, nS, index);
%%          % Delete phi ray from a given point
%%          function grid = deletePhiRay(grid, nS, stepsPhi); 
%%              grid.source(nS).deletePhiRay(stepsPhi, grid.deltaX, grid.deltaP);
%%          end
%%          % Compute the amplitude using Proximal Rays
%%          A = amplitudeProximal(grid, nS);
%%      
%%          % ODE Method ============================================
%%          % Compute the amplitude solving the ODE
%%          A = amplitudeODE(grid, nS);
%%          % Compute the hessian at the source point
%%          hess = hessPhaseODE(grid, nS, point);
%%                          
%%          % Interface Method ======================================
%%          % Compute the rays as straight lines
%%          grid = computeTrajecStraight(grid, nS);
%%          % Compute the amplitude
%%          A = amplitudeSI(grid, nS);
%%          I = interfaceLoss(grid, nS);
%%      end
%%  
%%      %======================================================================
%%      % Compute Trajectory + Amplitude + Pressure
%%      %======================================================================
%%      methods
%%          grid = computeSource(grid, nS);
%%          %%%% Main Ray
%%          % Derivative of the ray at a given point
%%          DP = deriveTrajec(grid, tauIndex, tau, X, nS);
%%          % Runge-Kutta step
%%          xn = stepRK4(grid, tauIndex, h, nS);
%%  
%%          %%%% Amplitude (Proximal ray method)
%%          % Derivative of the ray at a given point
%%          DP = deriveTrajecPhi(grid, tauIndex, tau, X, nS);
%%          % Runge-Kutta step
%%          xn = stepRK4Phi(grid, tau, X, h, nS);
%%          % Derivative of the ray at a given point
%%          Dx = deriveTheta(grid, tau, X, nS);
%%          % Runge-Kutta step
%%          xn = stepRK4Theta(grid, tau, X, h, nS);
%%  
%%          %%% Amplitude (ODE method)
%%          % Derive
%%          [D2x D2p] = deriveAmplitude(grid, nS, tauIndex, Dx, Dp);
%%          % Runge-Kutta step
%%          [Dxn Dpn] = stepRK2Amplitude(grid, nS, step, tauIndex, Dx, Dp);
%%      end
%%  
%%  
%%      methods
%%          % Flip the trajectories of the given source
%%          function grid = flipTrajectories(grid, nS); grid.source(nS).flipTrajectories(grid.deltaX, grid.deltaP); end;
%%          % Sort the trajectories of the given source
%%          function grid = sortTrajectories(grid, nS, method); grid.source(nS).sortTrajectories(method); end;
%%      
%%  
%%          % Compute the reverse amplitude solving the ODE
%%          A = revAmplitudeODE(grid, nS);
%%          % Compute the reverse amplitude for the given initial pressure
%%          A = revAmplitudeProximal(grid, nS);
%%          % Compute the reverse amplitude using SI
%%          A = revAmplitudeSI(grid, nS);
%%  
%%          % Derive - Reversal
%%          [D2x D2p] = deriveRevAmplitude(grid, nS, tauIndex, Dx, Dp);
%%          % Runge-Kutta step - Reversal
%%          [Dxn Dpn] = stepRK2RevAmplitude(grid, nS, step, tauIndex, Dx, Dp);
%%  
%%      end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
