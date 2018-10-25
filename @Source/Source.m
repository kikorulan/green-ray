%================================================================================
% This class contains a source with multiple rays
%================================================================================
% Copyright (C) Kiko Rul·lan, Marta M. Betcke
classdef Source < handle
    %====================================================================
    % Properties
    %====================================================================
    properties (GetAccess = public, SetAccess = private)
        %%% SOURCE VARIABLES %==================
        x0 = 0;                 % Location of the source
        tau = [];               % Tau values
        tBeam = [];             % Time observed in the source for the beam
        aBeam = [];             % Amplitude observed in the source for the beam
        dBeam = [];             % Length of the front wave
        aForward = [];          % Amplitude signal measured at the source (for update)
        aForward_initial = [];  % Amplitude signal measured at the source (initial data)
        step = 0;               % Step taken for the calculation of rays
        nPoints = 0;            % Number of steps for the source

        %%% MAIN RAYS %=========================
        x = [];         % Hamiltonian x
        p = [];         % Hamiltonian p
        n = [];         % Eta
        Gn = [];        % Grad(Eta)
        nGn = [];       % Eta*grad(Eta)
        phi = [];       % Phase values
        xTD = [];       % Traveled distance
        amplitude = []; % Amplitude
        revAmplitude = []; % Reverse Amplitude
        pressure = [];  % Pressure over the ray

        % Jacobian
        Dx = [];
        Dp = [];
        D2phase = [];
        kIndex = [];
        % q - determinant of the Jacobian
        q = [];

        %%% Gaussian Beams
        xGB = []
        qGB = []
        amplitudeGB = []
        pressureGB = []
        Y = [];
        N = [];
        M = [];



    end

    properties (GetAccess = public, SetAccess = public)
        %%% PIXEL MATRICES %====================
        aReverse = [];   % Amplitude signal for the reverse propagation
        aPropagation = [] % Amplitude signal for the reverse propagation
        pixelDistance = [];
        pixelAttenuation = [];
        pixelPressure = [];
        pixelTime = [];
        pixelAReverse = [];
        pixelAPropagation = [];
        pixelAngleCorrection = [];
    end
    %====================================================================
    % Constructor function
    %====================================================================
    methods
        function s0 = Source(x0, n, nGn, phi, deltaX, deltaP, angleMin, angleMax, nR, step, nPoints)
            if (nargin == 11)
                s0.x0 = x0;
                s0.tau = 0:step:step*(nPoints-1);
                s0.step = step;
                s0.nPoints = nPoints;
                
                % Initialize Rays
                j = (1:nR)';
                if (nR == 1) angle = angleMin;
                else angle = angleMin + (j-1)*(angleMax-angleMin)/(nR-1); end;
                np0(:, 1, 1) = cos(angle);
                np0(:, 1, 2) = sin(angle);
                p0 = np0*n;
                
                % Main Ray
                s0.x = repmat(x0, [nR 1]);
                s0.p = p0;
                s0.n = repmat(n, [nR 1]);
                s0.nGn = repmat(nGn, [nR 1]);
                s0.Gn = s0.nGn(:, 1, :)/n;
                s0.phi = 0;
                s0.xTD = repmat(0, [nR 1]);
                
            end       
        end
    end
    %====================================================================
    % Allocate space for rays
    %====================================================================
    methods
        % Rays
        function s0 = allocateRays(s0, nR, nPoints, dim)
            % Main ray
            s0.x   = [s0.x   nan(nR, nPoints-1, dim)];
            s0.p   = [s0.p   nan(nR, nPoints-1, dim)];
            s0.n   = [s0.n   nan(nR, nPoints-1)];
            s0.nGn = [s0.nGn nan(nR, nPoints-1, dim)];
            s0.Gn  = [s0.Gn  nan(nR, nPoints-1, dim)];
            s0.phi = [s0.phi zeros(1, nPoints-1)];
            s0.xTD = [s0.xTD nan(nR, nPoints-1)];
            % Jacobian
            s0.q = nan(nR, nPoints);            
            s0.Dx = nan(nR, nPoints, dim^2);
            s0.Dp = nan(nR, nPoints, dim^2);
            s0.D2phase = nan(nR, nPoints, dim^2);
            s0.kIndex = zeros(nR, nPoints);
            s0.amplitude = nan(nR, nPoints);      
            s0.pressure = zeros(nR, nPoints);
        end
        % Gaussian beams
        function s0 = allocateGauss(s0, nRsub)
            [nR, nPoints, dim] = size(s0.x);
            s0.xGB = nan(nR, nPoints, dim, 2*nRsub);
            %s0.qGB = nan(nR, nPoints);
            s0.amplitudeGB = zeros(nR, nPoints, 1, 2*nRsub);
            s0.pressureGB = zeros(nR, nPoints, 1, 2*nRsub);
            %s0.Y = nan(nR, nPoints, (dim+1)^2);
            %s0.N = nan(nR, nPoints, (dim+1)^2);
            %s0.M = nan(nR, nPoints, (dim+1)^2);
        end
    end
    %====================================================================
    % Insert values functions
    %====================================================================
    methods
        % Main ray %============================
        function s0 = insertX(s0, x, index); s0.x(:, index, :) = x; end;
        function s0 = insertP(s0, p, index); s0.p(:, index, :) = p; end;
        function s0 = insertTau(s0, tau, index); s0.tau(:, index, :) = tau; end;
        function s0 = insertN(s0, n, index); s0.n(:, index, :) = n; end;
        function s0 = insertNGN(s0, nGn, index); s0.nGn(:, index, :) = nGn; end;
        function s0 = insertGN(s0, Gn, index); s0.Gn(:, index, :) = Gn; end;
        function s0 = insertPhi(s0, phi, index); s0.phi(:, index, :) = phi; end;
        function s0 = insertPhiVector(s0, phi); s0.phi = phi; end;
        function s0 = insertXTD(s0, xTD, index); s0.xTD(:, index, :) = xTD; end;

        % Amplitude and Reverse Amplitude %=====
        function s0 = insertA(s0, amplitude, index); s0.amplitude(:, index, :) = amplitude; end;
        function s0 = insertAVector(s0, amplitude); s0.amplitude = amplitude; end;
        function s0 = insertRevA(s0, revAmplitude); s0.revAmplitude = [s0.revAmplitude revAmplitude]; end;
        function s0 = insertRevAVector(s0, revAmplitude); s0.revAmplitude = revAmplitude; end;
        function s0 = insertQ(s0, q, index); s0.q(:, index, :) = q; end;
        function s0 = insertQVector(s0, q); s0.q = q; end;
        function s0 = insertDX(s0, Dx, index); s0.Dx(:, index, :) = Dx; end;
        function s0 = insertDXVector(s0, Dx); s0.Dx = Dx; end;
        function s0 = insertDP(s0, Dp, index); s0.Dp(:, index, :) = Dp; end;
        function s0 = insertDPVector(s0, Dp); s0.Dp = Dp; end;
        function s0 = insertD2Phase(s0, D2phase, index); s0.D2phase(:, index, :) = D2phase; end;
        function s0 = insertKIndex(s0, kI, index); s0.kIndex(:, index, :) = kI; end;
        function s0 = insertPressure(s0, pressure, index); s0.pressure(:, index, :) = pressure; end;
        function s0 = insertPressureVector(s0, pressure); s0.pressure = pressure; end;

        % Gauss beams %=========================
        function s0 = insertXGB(s0, xGB, indexRay); s0.xGB(:, :, :, indexRay) = xGB; end;
        function s0 = insertQGB(s0, qGB, index); s0.qGB(:, index, :) = qGB; end;
        function s0 = insertAmplitudeGB(s0, amplitudeGB, indexRay); s0.amplitudeGB(:, :, :, indexRay) = amplitudeGB; end;
        function s0 = insertPressureGB(s0, pressureGB, indexRay); s0.pressureGB(:, :, :, indexRay) = pressureGB; end;
        function s0 = insertPressureGBindex(s0, pressureGB, indexRay, index); s0.pressureGB(:, index, :, indexRay) = pressureGB; end;
        function s0 = insertYGauss(s0, Y, index); s0.Y(:, index, :) = Y; end;
        function s0 = insertNGauss(s0, N, index); s0.N(:, index, :) = N; end;
        function s0 = insertMGauss(s0, M, index); s0.M(:, index, :) = M; end;
    end

    %====================================================================
    % Insert rays at the specified index
    %====================================================================
    methods
        % Main ray
        function s0 = insertRayX(s0, x, index); s0.x = cat(1, s0.x(1:index, :, :), x, s0.x(index+1:end, :, :)); end;
        function s0 = insertRayP(s0, p, index); s0.p = cat(1, s0.p(1:index, :, :), p, s0.p(index+1:end, :, :)); end;
        function s0 = insertRayNGN(s0, nGn, index); s0.nGn = cat(1, s0.nGn(1:index, :, :), nGn, s0.nGn(index+1:end, :, :)); end;
        function s0 = insertRaGN(s0, Gn, index); s0.Gn = cat(1, s0.Gn(1:index, :, :), Gn, s0.Gn(index+1:end, :, :)); end;
        function s0 = insertRayPhi(s0, phi, index); s0.phi = cat(1, s0.phi(1:index, :, :), phi, s0.phi(index+1:end, :, :)); end;
        function s0 = insertRayXTD(s0, xTD, index); s0.xTD = cat(1, s0.xTD(1:index, :, :), xTD, s0.xTD(index+1:end, :, :)); end;
        % Amplitude
        function s0 = insertRayDX(s0, Dx, index); s0.Dx = cat(1, s0.Dx(1:index, :, :), Dx, s0.Dx(index+1:end, :, :)); end;
        function s0 = insertRayDP(s0, Dp, index); s0.Dp = cat(1, s0.Dp(1:index, :, :), Dp, s0.Dp(index+1:end, :, :)); end;
        function s0 = insertRayQ(s0, q, index); s0.q = cat(1, s0.q(1:index, :, :), q, s0.q(index+1:end, :, :)); end;
    end

    %====================================================================
    % Measure signals
    %====================================================================
    methods
        % Compute the contributions from all rays
        s0 = computeBeam(s0, dt);
        % Set the values for the beam computation
        function s0 = setBeam(s0, tBeam, aBeam); 
            s0.aBeam = aBeam;
            s0.tBeam = tBeam; 
        end;
        % Set the values for the time signal computation
        function s0 = setForwardSignal(s0, aForward); s0.aForward = aForward; end;
        function s0 = setForwardSignal_initial(s0, aForward_initial); s0.aForward_initial = aForward_initial; end;
    end

    %====================================================================
    % Plot
    %====================================================================
    methods
        % Plot the ray trajectories
        h = plot_rays(source, grid, nColours);
        % Plot the ray trajectories
        h = plot_subrays(source, grid, nColours);
        % Plot Q
        h = plot_q(source, nColours, plotImag);
        % Plot amplitude
        h = plot_amplitude(source, nR, nColours, plotImag);
        % Plot amplitude
        [h1, h2] = plot_amplitude_subrays(source, indexRay, plotImag);
        % Plot YNM matrices
        plot_YNM(source, nR);
        % Plot pressure
        h = plot_pressure(source);
        [h1, h2] = plot_pressure_subrays(source, indexRay);
    end
    %====================================================================
    % Select rays
    %====================================================================
    methods
        function s0 = selectRays(s0, ray_index, t_index)
            s0.x = s0.x(ray_index, :, :);
            for i = 1:length(t_index)
                s0.x(i, t_index(i)+1:end, :) = nan;
            end
        end
    end
    %====================================================================
    % Delete rays
    %====================================================================
    methods
        % Delete phi and theta rays and amplitude vector
        function s0 = deleteAmplitude(s0)
            % Amplitude
            s0.amplitude(:) = nan;
            % q
            s0.q(:) = nan;
            s0.Dx(:) = nan;
            s0.Dp(:) = nan;
        end

        % Delete rays
        function s0 = deleteRays(s0)
            % Main ray
            s0.x = [];
            s0.phi = [];
            % Amplitude
            s0.revAmplitude = [];
            % Pressure
            s0.pressure = [];
        end

        function s0 = deleteAuxInfo(s0)
            % Main ray
            s0.p = [];
            s0.n = [];
            s0.nGn = [];
            s0.Gn = [];
            s0.xTD = [];
            s0.tau = [];
            % Amplitude
            %s0.amplitude = [];
            s0.kIndex = [];
            s0.Dx = [];
            s0.Dp = [];
            s0.D2phase = [];

            % Pixel matrices
            %s0.pixelTime = [];
            s0.pixelDistance = [];
            s0.pixelPressure = [];
            s0.pixelAngleCorrection = [];
        end

    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DEPRECATED                                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      %====================================================================    %%
%%      % Modification of trajectories                                           %%
%%      %====================================================================    %%
%%      methods                                                                  %%
%%          s0 = flipTrajectories(s0, deltaX, deltaP);                           %%
%%          s0 = sortTrajectories(s0, method);                                   %%
%%                                                                               %%
%%          % Delete Phi ray from the given point and compute the next initial point
%%          function s0 = deletePhiRay(s0, index, deltaX, deltaP)
%%              switch length(s0.x0)
%%                  case 1
%%                  case 2
%%                      %%%% Phi
%%                      s0.xPhi = s0.xPhi(:, 1:index-1, :);
%%                      s0.pPhi = s0.pPhi(:, 1:index-1, :);
%%                      s0.nGnPhi = s0.nGnPhi(:, 1:index-1, :);
%%                      % Compute next index
%%                      p0 = s0.pR(:, index, :);
%%                      np0 = sqrt(sum(p0.*p0, 1));
%%                      xS = s0.xR(:, index, :) - deltaX*p0./repmat(np0, [2 1 1]);
%%                      pPhi = [p0(1, 1, :)*cos(deltaP) - p0(2, 1, :)*sin(deltaP); ...
%%                                  p0(2, 1, :)*cos(deltaP) + p0(1, 1, :)*sin(deltaP)];
%%                      s0.insertXPhi(xS + deltaX*pPhi./repmat(np0, [2 1 1]));
%%                      s0.insertPPhi(pPhi);
%%                  case 3 %%%%%%%%%%%% NOT IMPLEMENTED %%%%%%%%%%%%%%%
%%                      error('Not implemented');
%%              end
%%          end
%%                                                                               %%
%%      end                                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
