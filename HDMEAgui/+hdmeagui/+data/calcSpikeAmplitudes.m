function D = calcSpikeAmplitudes(D)
    disp('Calculating spike amplitudes...');
    T = D.templates;
    nT = length(T.names);
    maxJitter = min(T.Tf-T.cutleft, T.cutleft);
%     maxJitter = 4;
    if ~isfield(D, 'singleChannelDataModified') || isempty(D.singleChannelDataModified)
        D.singleChannelDataModified = D.singleChannelDataSorted;
        % normalize with noise std
        D.singleChannelDataModified(:,3) = ...
            D.singleChannelDataModified(:,3)./D.smad(D.singleChannelDataModified(:,2))';
    end
    % times channels amplitudes
    if isfield(D, 'precutSpikes')
        thr = 0;
    else
        thr = 3.5;
    end
    mainIdx = find(D.singleChannelDataModified(:,3) > thr);
    st1 = D.singleChannelDataModified(mainIdx,1);
    for t=1:nT
        if ~T.accepted(t) || T.amps_modified(t)
            continue
        end
        idx1 = find(T.bIdx(:,t));
        st2 = D.singleChannelDataModified(idx1(T.selIdx{t}),1);
        disp('Matching spike trains');
        tic
        fprintf('St1: %d St2: %d\n', length(st1), length(st2));
        M = mysort.spiketrain.findNearSpikes(st1, st2, 30);
        if isempty(M)
            continue
        end
        D.singleChannelDataModified(mainIdx(M(:,1)), 3) = -1;
%         m_chanIdx = D.singleChannelDataModified(M(:,1), 2);
%         temp = mysort.util.v2m(T.selected_template(t,:), T.nC);
%         temp(temp>0) = 0;
%         toc
%         disp('Updating amplitudes');
%         tic
%         for j = -maxJitter:maxJitter
%             tempj = temp(:, T.cutleft+j+1);
%             jitteridx = M(:,3) == j;
%             D.singleChannelDataModified(M(jitteridx,1), 3) = ...
%                 D.singleChannelDataModified(M(jitteridx,1), 3) + ...
%                 tempj(m_chanIdx(jitteridx));
%         end
        toc
        T.amps_modified(t) = true;
    end    
    D.templates = T;
end