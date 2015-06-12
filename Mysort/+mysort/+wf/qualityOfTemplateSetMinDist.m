function [Q, Qq, E, thr, P] = qualityOfTemplateSet(T, priors, debug)
    % Compute the quality of a given template set as the expected error for
    % classification via LDA.
    % THE COVARIANCE MATRIC IS ASSUMED TO BE THE IDENTITY
    % Inputs: 
    %   T      - template set. 3D matrix, [time x electrodes x items]
    %   priors - prior probabilities for all templates
    
    if nargin < 3
        debug = 0;
    end
    
    [L, nE, nT] = size(T);
    assert(length(priors) == 1 || length(priors) == nT, 'Either one prior for all or one prior for each!');
    if length(priors) == 1
        priors = priors * ones(1, nT);
    end
    priors = priors(:)';
    
    noisePrior = 1 - sum(priors);
    
    % add noise template
    T = [zeros(1, nE*L); mysort.wf.t2v(T)];
    priors = [noisePrior priors];
    nT = nT+1;
    
    % Compute optimal discrimination thresholds via LDA
    nD = (nT+1)*nT/2;
    thr = zeros(1,nD);
    I = zeros(2, nD);
    E = zeros(1, nD);
    P = zeros(nD, L*nE);
    combinations = zeros(2, nD);
    count = 0;
    if debug 
        fh = mysort.plot.figure([1200 900]);
    end
    for i=1:nT
        for j=i+1:nT
            count = count+1;
            combinations(:, count) = [i j]';
            [x1 , x2] = minDistTwoTemplates( T(i,:) ,  T(j,:) , L , nE );
            Pj = x1 - x2;
            %            P(count,:) = T(i,:) - T(j,:);
            N = norm(Pj);
            %N = norm(P(count,:));
            % Normalizing the projection has the advantage that the
            % standard deviation on the fisher discriminant stays identical
            %thr_ = -log(priors(i)/priors(j)) + .5*(T(i,:)*T(i,:)' - T(j,:)*T(j,:)');
            thr_ = -log(priors(i)/priors(j)) + .5*(x1*x1' - x2*x2');
            if N < 0.000001
                MU = [0 0]';
                % templates are identical, always decide for the more
                % probable template
                if priors(i) < priors(j)
                    thr(count) = 1;
                    i1 = 1;
                    i2 = 0;
                else
                    thr(count) = -1;
                    i1 = 0;
                    i2 = 1;                    
                end
            else
                thr(count) = thr_/N;
                Pj = Pj / N;
                %P(count,:) = P(count,:)/N;
%                MU  = T([i j],:)*P(count,:)';
                %MU  = [x1 ; x2]*P(count,:)';
                MU  = [x1 ; x2]*Pj';

                i1 = normcdf(thr(count), MU(1), 1);
                i2 = 1-normcdf(thr(count), MU(2), 1);
            end
            I(:, count) = [i1 i2]';  
            E(count) = priors([i j])*[i1 i2]';

            if debug
                %%
                vis();
            end
        end
    end
    Q = sum(E);
    Qq = Q/count;
    
    function vis
        if ~debug
            return
        end
        subplot(nT-1,nT-1,(i-1)*(nT-1) + j-1);
        dMU = abs(diff(MU));
        xr = linspace(min(MU)-3, max(MU)+3, 200);
        Y1 = priors(i)*normpdf(xr, MU(1), 1);
        Y2 = priors(j)*normpdf(xr, MU(2), 1);
       
        mima = [min(xr)-.001 max(xr)+.001];
        maxY = max([Y1 Y2]);
%         [~, idx] = min(abs(xr));
%         ythr = Y1(idx);
        minYval = 0.0000001;
        plot(MU, [minYval minYval], 'rx', 'markersize', 12, 'linewidth', 3);
        hold on
        plot(thr(count), minYval, 'og', 'markersize', 12, 'linewidth', 3);
        plot(xr, Y1,'b', 'linewidth', 2);
        plot(xr, Y2,'c', 'linewidth', 2);
        O = cumsum(Y1)-cumsum(Y2);
        O = (O-min(O)) / (max(O) - min(O));
        plot(xr, O, ':m')
        title(sprintf('(%d-%d), I: %.2f|%.2f, E:%.2f', i-1, j-1, i1, i2, E(count)));

        set(gca, 'xlim', mima, 'ylim', [-.05*maxY maxY]);       
        if i==1
    %        set(gca, 'yscale', 'log')
        end
    end
end