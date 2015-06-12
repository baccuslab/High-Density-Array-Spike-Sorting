function [Y Hd] = butter_bp_filtfilt(X, order, fc1, fc2, srate, burnin_time)
    if nargin < 5
        srate = 2; % makes nyq = 1;
    end
    if nargin < 6
        burnin_time = 20*order;
    end
    nyq = srate/2;
    fc1 = fc1/nyq;
    fc2 = fc2/nyq;

    d = fdesign.bandpass('n,f3dB1,f3dB2', order, fc1, fc2);
    Hd = design(d,'butter');

    Y = filter(Hd, X);
    Y(1:burnin_time, :) = 0;
    Y = flipud(filter(Hd, flipud(Y)));
    Y(end-burnin_time+1:end,:) = 0;
    
    
    % This is with "chunking" over the channels
%     chunker = mysort.util.Chunker(size(X,2), 'chunkSize', 1,...
%         'progressDisplay', 'console');
% 
%     Y = zeros(size(X));
%     while chunker.hasNextChunk()
%         idx = chunker.getNextChunk(); idx = idx(1);
%         Y(:,idx) = filter(Hd, X(:,idx));
%         Y(1:burnin_time, idx) = 0;
%         Y(:,idx) = flipud(filter(Hd, flipud(Y(:,idx))));
%         Y(end-burnin_time+1:end,idx) = 0;
%     end