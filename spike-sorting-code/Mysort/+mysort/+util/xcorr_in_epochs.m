
function [c totalEpochLength] = xcorr_in_epochs(X,epochs,maxLag,minEpochLength)
    % N -> number of parallel time series (channels)
    % T -> number of data points
    [nC,T]=size(X);

    if nargin < 4
        minEpochLength = maxLag;
    end
    % initialize noise xorr
    %normlag = T+[-maxLag:0 -1:-1:-maxLag];
    
    c = zeros(1,2*maxLag+1);
    epochLength = epochs(:,2)-epochs(:,1)+1;
    epochs = epochs(epochLength>minEpochLength,:);
    epochLength = epochLength(epochLength>minEpochLength);
    
    totalEpochLength = sum(epochLength);
    if nC == 1
        for e=1:size(epochs,1)
            % Compute the block diagonal elements
            c = c + xcorr(X(epochs(e,1):epochs(e,2)), maxLag, 'none');
            % normalize xcorr. Take care, that the crosscorrelation 
            % is calculated with zeropadding. this means for larger lags
            % we have slightly less data to compute the crosscorrelation
            % No, dont do this, this can cause negative eigenvalues.
            %xx = xx./normlag;
        end
    else
        for e=1:size(epochs,1)
            % Compute the block diagonal elements
            c = c + xcorr(X(1,epochs(e,1):epochs(e,2)),X(2,epochs(e,1):epochs(e,2)),maxLag, 'none');
            % normalize xcorr. Take care, that the crosscorrelation 
            % is calculated with zeropadding. this means for larger lags
            % we have slightly less data to compute the crosscorrelation
            % No, dont do this, this can cause negative eigenvalues.
            %xx = xx./normlag;
        end
    end
    c = c/totalEpochLength;
end