function s0 = computeBeam(s0, dt)
% Function for finding the total ray in the given source
    [nR lengthRay dim] = size(s0.x);

    % Time signal
    tMax = max(max(s0.phi(s0.phi<inf)));
    s0.tBeam = 0:dt:tMax+dt;

    % Length of the front wave
    rayI = s0.x(2:nR, :, :);
    rayF = s0.x(1:nR-1, :, :);
    normRay = sqrt(sum((rayI - rayF).^2, 3));
    normRay(isinf(normRay) | isnan(normRay)) = 0;
    dBeam0 = zeros(1, lengthRay); % zero distance for first and last rays
    s0.dBeam = 0.5*(cat(1, dBeam0, normRay) + cat(1, normRay, dBeam0));
    dimD = size(s0.dBeam);
    % Amplitude measured by the sensor
    s0.pressure(isinf(s0.pressure) | isnan(s0.pressure)) = 0;

    % Preallocate distances and pressures
    dBeamSpline = nan(nR, length(s0.tBeam));
    pressureSpline = nan(nR, length(s0.tBeam)); 
    % Interpolate pressure and distance values
    for j = 1:nR
        phiFin = s0.phi(j, s0.phi(j, :, 1) < inf);
        % Distance
        dBeamFin = s0.dBeam(j, s0.phi(j, :, 1) < inf);
        dBeamSpline(j, s0.tBeam <= phiFin(end)) = ...
            spline([-1 phiFin], [0 dBeamFin], s0.tBeam(s0.tBeam <= phiFin(end)));
        dBeamSpline(j, s0.tBeam > phiFin(end)) = 0;
        % Pressure
        pressureFin = s0.pressure(j, s0.phi(j, :, 1) < inf);
        pressureSpline(j, s0.tBeam <= phiFin(end)) = ...
            spline([-1 phiFin], [0 pressureFin], s0.tBeam(s0.tBeam <= phiFin(end)));
        pressureSpline(j, s0.tBeam > phiFin(end)) = 0;
    end 
    %s0.aBeam = sum(pressureSpline.*dBeamSpline, 1); % Distance correction
    s0.aBeam = sum(pressureSpline, 1); % No distance correction
end
