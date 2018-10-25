function source = inverse_timeSource_fft(grid, source)
% INVERSE_TIMESOURCE_FFT computes the time signal for the whole domain given a source
%source = inverse_timeSource_fft(grid, nS)
%
% INPUTS
% grid: object from gridRT class
% source: source to do the conversion
% 
% OUTPUTS
% grid: object from gridRT class
%   - pixelAmplitudeSignal: amplitude signals for the pixels in the domain
%   - tSignal: time signal
%
% Copyright (C) Kiko Rul·lan, Marta M. Betcke

% Filter time array
cMax = max(grid.c(:)); % maximum sound speed
cMin = min(grid.c(:)); % minumum sound speed
dt = grid.dt; % time signal increment
delayPropMax = length(grid.tForward); % maximum propagation delay
delayFilterMax = floor(grid.delayFilter*grid.cFilter/cMin); % maximum delay for the filter
lengthFilter = length(grid.tFilter);

% Length reverse 
lengthForward = length(grid.tForward);
lengthReverse = lengthForward + lengthFilter + delayPropMax + delayFilterMax;


%==================================================
% Choose particular sensor
%==================================================
if length(source) == 1
    % Display initialisation message
    msg = strcat({'Computing time propagation for the given source...'});
    disp(msg{1});

    % Time signal FFT
    aForward = source.aForward;
    aForwardPad = padarray(aForward, [0 lengthFilter + delayPropMax + delayFilterMax], 0, 'post');
    aForwardFFT = fft(aForwardPad);

    % Generate Signal matrix
    aForwardFFTMatrix = repmat(permute(aForwardFFT, [1 3 2]), [grid.Nx grid.Ny 1]);
    % Shift matrix to adjust delay
    delayMatrix = min(floor(source.pixelTime/dt), delayPropMax) - grid.pixelDelayFilter; 
    pixelDelayMatrix = repmat(delayMatrix, [1 1 lengthReverse]);
    shiftVector = permute(0:1:(lengthReverse-1), [1 3 2]);
    pixelShift = repmat(shiftVector, [grid.Nx grid.Ny 1]);
    pixelExp = exp(-2*pi*i/lengthReverse*pixelShift.*pixelDelayMatrix);

    % Convolve in Frequency domain
    attenuationMatrix = repmat(source.pixelAttenuation, [1 1 lengthReverse]);
    source.pixelAPropagation = attenuationMatrix.*ifft(aForwardFFTMatrix.*grid.pixelFilterFFT.*pixelExp, lengthReverse, 3);
    % Adjust filter Delay
    source.pixelAPropagation(:, :, end-delayFilterMax-lengthFilter+1:end) = [];     
    source.pixelAReverse = source.pixelAPropagation(:, :, delayPropMax);
%==================================================
% Sum all sensors
%==================================================
else
    % Display initialisation message
    msg = strcat({'Computing time propagation for all source...'});
    disp(msg{1});

    grid.pixelAPropagation = zeros(grid.Nx, grid.Ny, lengthForward + delayPropMax);
    grid.pixelAReverse = zeros(grid.Nx, grid.Ny);
    for j = 1:length(source)
        % Display initialisation message
        msg = strcat({'Computing time propagation for source '}, int2str(j), {'...'});
        disp(msg{1});

        % Time signal FFT
        aForward = fliplr(source(j).aForward);
        aForwardPad = padarray(aForward, [0 lengthFilter + delayPropMax + delayFilterMax], 0, 'post');
        aForwardFFT = fft(aForwardPad);
        
        % Generate Signal matrix
        aForwardFFTMatrix = repmat(permute(aForwardFFT, [1 3 2]), [grid.Nx grid.Ny 1]);
        % Shift matrix to adjust delay
        delayMatrix = min(floor(source(j).pixelTime/dt), delayPropMax) - grid.pixelDelayFilter; 
        pixelDelayMatrix = repmat(delayMatrix, [1 1 lengthReverse]);
        shiftVector = permute(0:1:(lengthReverse-1), [1 3 2]);
        pixelShift = repmat(shiftVector, [grid.Nx grid.Ny 1]);
        pixelExp = exp(-2*pi*i/lengthReverse*pixelShift.*pixelDelayMatrix);
        
        % Convolve in Frequency domain
        attenuationMatrix = repmat(source(j).pixelAttenuation, [1 1 lengthReverse]);
        pixelAPropagation = attenuationMatrix.*ifft(aForwardFFTMatrix.*grid.pixelFilterFFT.*pixelExp, lengthReverse, 3);
        % Adjust filter Delay
        pixelAPropagation(:, :, end-delayFilterMax-lengthFilter+1:end) = [];
        grid.pixelAPropagation = grid.pixelAPropagation + pixelAPropagation;
        grid.pixelAReverse = grid.pixelAReverse + pixelAPropagation(:, :, delayPropMax);
    end
end
