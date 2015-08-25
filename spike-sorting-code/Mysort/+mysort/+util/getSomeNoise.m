
function Noise = getSomeNoise(X, epochs, L, N)
    numChannel = size(X,1);
    % cut randomly N samples of length L from epochs from X
    Noise = zeros(N,L*numChannel);
    numEpochs = size(epochs,1);
    % Correct for length of sample, only valid sample positions should be
    % used for probability estimation
    epochLength = epochs(:,2) - epochs(:,1) - L;
    epochLength(epochLength<0) = 0;
    
    totalLength = sum(epochLength);
    cumSumLength = cumsum(epochLength);
    for i=1:N
        % choose one epoch with probability equal to relative length of epoch
        e = ceil(rand*totalLength);
        for k=1:numEpochs
            if cumSumLength(k) >= e
                % we have choosen epoch k, now choose a starting index in
                % the k-th epoch
                idx = ceil(rand*epochLength(k))-1;
                start = round(epochs(k,1)+idx);
                ende = round(start+L-1);
                temp = X(:, start:ende)';
                Noise(i,:) = temp(:);
                break
            end
        end 
    end    
end