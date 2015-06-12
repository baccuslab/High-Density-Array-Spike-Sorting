function DATA = cutSpikes(DATA, INTER, CONFIG)
    error('This function is obsolete!');
    fprintf('Need to recut spikes...\n');

    % cut spikes
    DATA.cutSpikesTf   = CONFIG.CutSpikesTf;
    DATA.cutSpikesCutleft = CONFIG.CutSpikesCutleft;
    
    sT = DATA.singleChannelSpikeTimes;
    
    tf = CONFIG.CutSpikesTf;
    cl = CONFIG.CutSpikesCutleft;
    
    nS = length(sT);
    elIdx = DATA.useElectrodes;
    nC_ = sum(elIdx);
    if isa(DATA.X, 'mysort.mea.CMOSMEA')
        DATA.cutSpikes = DATA.X.getCutWaveforms(sT, cl, tf-cl-1);
        DATA.cutSpikes = mysort.util.vChannelMul(DATA.cutSpikes, 1./DATA.smad);
    else
        DATA.cutSpikes = zeros(nS, nC_ * tf);
        for i=1:nS;
            wf = DATA.X(sT(i)-cl:sT(i)+tf-cl-1, find(elIdx))';
            wf = wf./repmat(DATA.smad(elIdx), 1, size(wf,2));
            DATA.cutSpikes(i,:) = mysort.util.m2v(wf);
        end      
    end
    fprintf('Done.\n');
