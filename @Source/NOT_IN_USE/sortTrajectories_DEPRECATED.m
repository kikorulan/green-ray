%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEPRECATED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s0 = sortTrajectories(s0, method)
% This function sorts the trajectories x of the given source s0, stored in xR
    [dim length nR] = size(s0.xR);
    for n = 1:nR
        s0.xR(:, :, n) = [s0.xR(:, s0.xR(1, :, n) == inf, n) s0.xR(:, s0.xR(1, :, n) < inf, n)];
        s0.pR(:, :, n) = [s0.pR(:, s0.pR(1, :, n) == inf, n) s0.pR(:, s0.pR(1, :, n) < inf, n)];
        %s0.phi(1, :, n) = [s0.phi(1, s0.phi(1, :, n) == inf, n) s0.phi(1, s0.phi(1, :, n) < inf, n)];
    end

    switch dim
        case 1
        case 2
            % Phi and Theta rays
            if(~isempty(s0.xPhi))
                for n = 1:nR
                    s0.xPhi(:, :, n) = [s0.xPhi(:, s0.xPhi(1, :, n) == inf, n) s0.xPhi(:, s0.xPhi(1, :, n) < inf, n)];
                    s0.pPhi(:, :, n) = [s0.pPhi(:, s0.pPhi(1, :, n) == inf, n) s0.pPhi(:, s0.pPhi(1, :, n) < inf, n)];
                end
            end
            % Q
            if(~isempty(s0.q))
                for n = 1:nR
                    s0.Dx(:, :, n) = [s0.Dx(:, isnan(s0.Dx(1, :, n)) , n) s0.Dx(:, ~isnan(s0.Dx(1, :, n)), n)];
                    s0.q(1, :, n) = [s0.q(1, isnan(s0.q(1, :, n)), n) s0.q(1, ~isnan(s0.q(1, :, n)), n)];
                 end
                if(strcmp(method,'ODE')) % Index correction
                    s0.q = [repmat(NaN, [1 1 nR]) s0.q(:, 1:end-1, :)];
                elseif(strcmp(method,'Prox'));
                    %s0.q = [ s0.q(:, 2:end, :) repmat(NaN, [1 1 nR])];
                else error('Method not defined');
                end
            end
            % Amplitude
            if(~isempty(s0.amplitude))
                for n = 1:nR
                    %disp('sortTrajectories')
                    s0.amplitude(1, :, n) = [s0.amplitude(1, isnan(s0.amplitude(1, :, n)) , n) s0.amplitude(1, ~isnan(s0.amplitude(1, :, n)), n)];
                end
                if(strcmp(method,'ODE')) % Index correction
                     s0.amplitude = [repmat(NaN, [1 1 nR]) s0.amplitude(:, 1:end-1, :)];
                elseif(strcmp(method,'Prox'));
                     %s0.amplitude = [s0.amplitude(:, 2:end, :) repmat(NaN, [1 1 nR])];
                elseif(strcmp(method,'SI'));
                else error('Method not defined');
                end
            end
        case 3 %%%%%%%%%%%%%% NOT IMPLEMENTED %%%%%%%%%%%%%%%%%%%%%%%%%
                error('Not implemented');
    end
end 
        
