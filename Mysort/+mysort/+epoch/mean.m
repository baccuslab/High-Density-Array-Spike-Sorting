
function m = mean(epochs, X)

    % x = []; % Dirty, online mean would be better
    % for i=1:size(epochs,1)
    %     xt = X(:, epochs(i,1):epochs(i,2));
    %     x = [x xt];
    % end
    % m = mean(x,2);

    nrNoiseEpochs = size(epochs,1);
    epochLengths = epochs(:,2) - epochs(:,1) + 1;
    totalLengths = cumsum(epochLengths);
    idx = zeros(1,totalLengths(end));
    startIdx = 1;
    endIdx = totalLengths(1);
    for i=1:nrNoiseEpochs
        idx(startIdx:endIdx) = epochs(i,1):epochs(i,2);
        startIdx = endIdx + 1;
        endIdx = totalLengths(min(nrNoiseEpochs,i+1));
    end

    % disp(' ')
    % display([' --> length noise ' num2str(length(idx))])
    % disp(' ')

    m = mean(X(:,idx),2);

    % disp('done epochsMean')