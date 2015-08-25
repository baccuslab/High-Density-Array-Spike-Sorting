
function [fetX pcs T pcs_shift] = dimReductionPCA(X, NDims, Templates, maximalElements, upsample)
    T = [];
    if nargin < 2 || isempty(NDims);
        NDims = 4;
    end
    if nargin < 4 || isempty(maximalElements)
        maximalElements = 1000000;
    end
    if nargin < 5 || isempty(upsample)
        upsample = 1;
    end
    pcs_shift = [];
    meanX = [];
    if numel(X) <= maximalElements && upsample == 1
        [pcs fetX] = princomp(X);
        fetX = fetX(:,1:NDims);
    else
        if numel(X) <= maximalElements
            [pcs] = princomp(X);
        else
            r = randperm(size(X,1));
            nRows = ceil(maximalElements/size(X,2));    
            warning('Too many rows selected for PCA. Only subset is used (%d of %d)!', nRows, size(X,1));
            % 'score' added in order to force matlab to use the svd pca
            % algorithm (might solve a bug with negative eigenvalues)
            [pcs score] = princomp(X(r(1:nRows),:),'econ');            
        end
        % compute mean free X
        meanX = mean(X, 1);
        X_mf = (X-repmat(meanX, size(X,1), 1));
    
        if upsample == 1
            fetX = X_mf*pcs(:,1:NDims);     
        else
            warning('This works but if of no use. Upsample and align waveforms on mean before PCA!');
            % compute the shift == 0 case, this is the normal case.
            F(:,:,2*upsample) = X_mf*pcs(:,1:NDims);            
            upsample = upsample*2;
            % resample pcs
            pcs_up = resample(pcs(:,1:NDims), upsample, 1);
            % org samples  :    x     x     x     x     x     x 
            % upsampled 3x :    x x x x x x x x x x x x x x x x x x
            % shift -3     :0     2     3     4     5     6      
            % shift -1     :  0     2     3     4     5     6     
            % shift  0     :    1     2     3     4     5     6
            % shift  1     :      1     2     3     4     5     6
            % shift  2     :        1     2     3     4     5     6
            pcs_shift(:,:,upsample) = pcs(:,1:NDims);
            for shift = 1:upsample-1
                pcs_shift(:,:,upsample-1+shift+1) = pcs_up(shift+1:upsample:end,:);
                F(:,:,shift) = X_mf*squeeze(pcs_shift(:,:,upsample-1+shift+1));
                 % Dont do this, gives artefact around 0 !
%                 P = X_mf*squeeze(pcs_shift(:,:,upsample-1+shift+1));
%                 FF(:,:,1) = fetX;
%                 FF(:,:,2) = P;
%                 [F IF] = max(abs(FF),[],3);
%                 fetX(IF==2) = P(IF==2);
                
                pcs_shift(:,:,upsample-1-shift+1) = [zeros(1, NDims); pcs_up(shift+1+upsample:upsample:end,:)];
                F(:,:,shift+upsample) = X_mf*squeeze(pcs_shift(:,:,upsample-1-shift+1));
%                 P = X_mf*squeeze(pcs_shift(:,:,upsample-1-shift+1));
%                 FF(:,:,1) = fetX;
%                 FF(:,:,2) = P;                
%                 [F IF] = max(abs(FF),[],3);
%                 fetX(IF==2) = P(IF==2);
            end
            
            % find the projection that procuced the maximal value for each
            % row
            [maxis proj_idx] = max(abs(F),[],3);
            [maxis dim_idx]  = max(maxis,[],2);
            maxProjForEachRow = zeros(size(F,1),1);
            fetX = zeros(size(F,1), size(F,2));
            for i=1:size(F,1)
                maxProjForEachRow(i) = proj_idx(i, dim_idx(i));
                fetX(i,:) = F(i,:,maxProjForEachRow(i));
            end
            
        end
    end
    
    if exist('Templates', 'var') && nargout>2 && ~isempty(Templates)
        if isempty(meanX)
            meanX = mean(X,1);
        end
        T = (Templates - repmat(meanX, size(Templates,1), 1))*pcs(:,1:NDims);
    end
    