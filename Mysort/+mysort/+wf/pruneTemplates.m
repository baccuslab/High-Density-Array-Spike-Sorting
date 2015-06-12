function T = pruneTemplates(T, varargin)
    P.maxChannels = [];
    P.minChannels = 1;
    P.absThreshold = 0;
    P.setInvalidChannelsTo = 0;
    P = mysort.util.parseInputs(P, varargin, 'error');
    
    [Tf nC nT] = size(T);
    
    if ~isempty(P.maxChannels)
        for temp=1:nT
            [maxis maxTimes] = max(abs(squeeze(T(:,:,temp))),[],1);
            maxis(isnan(maxis)) = -1;
            nGreaterThanTreshold = sum(maxis > P.absThreshold);
            nP = max(P.minChannels, min(P.maxChannels, nGreaterThanTreshold));
            sr = sortrows([maxis' (1:nC)']);
            T(:,sr(1:(nC-nP),2),temp) = P.setInvalidChannelsTo;
        end
    end
   