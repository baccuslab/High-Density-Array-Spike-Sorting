function DATA = cutBrushedSpikes(DATA, INTER, CONFIG)
    T = DATA.templates;
    t = INTER.tIdx;
    fprintf('Cutting brushed spikes...\n');
    tic
 % cut spikes
    DATA.cutSpikesTf   = CONFIG.CutSpikesTf;
    DATA.cutSpikesCutleft = CONFIG.CutSpikesCutleft;
    
    
    
    tf = CONFIG.CutSpikesTf;
    cl = CONFIG.CutSpikesCutleft;
    
    
    elIdx = DATA.useElectrodes;
    nC_ = sum(elIdx);
%     if isa(DATA.X, 'mysort.mea.CMOSMEA')
%         DATA.cutSpikes = DATA.X.getCutWaveforms(sT, cl, tf-cl-1);
%         DATA.cutSpikes = mysort.util.vChannelMul(DATA.cutSpikes, 1./DATA.smad);
%     else
        
        if isfield(DATA, 'precutSpikes')
            new_idx = find(INTER.bIdx);
            new_idx = ceil(new_idx/4);
            T.cutSpikes{t} = DATA.precutSpikes(new_idx,:);
        else
            sT = DATA.singleChannelDataSorted(INTER.bIdx, 1);
            nS = length(sT);
            T.cutSpikes{t} = zeros(nS, nC_ * tf);
            for i=1:nS
                if sT(i)-cl >= 1 && sT(i)+tf-cl-1 <= size(DATA.X,1)
                    wf = DATA.X(sT(i)-cl:sT(i)+tf-cl-1, find(elIdx))';
                    wf = wf./repmat(DATA.smad(find(elIdx))', 1, size(wf,2));
                    T.cutSpikes{t}(i,:) = mysort.wf.m2v(wf);
                end
            end
        end
%     end
    toc
    fprintf('Done.\n');
    DATA.templates = T;