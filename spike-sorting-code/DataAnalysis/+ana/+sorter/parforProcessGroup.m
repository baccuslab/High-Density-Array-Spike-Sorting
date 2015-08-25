function [r spikes_cut] = parforProcessGroup(data, P)
    % merge spikes
    allspikes = ana.sorter.mergeSpikes(data.detUp, data.pksUp, data.detDown, data.pksDown, P.spikeDetection.mergeSpikesMaxDist);
    r.ts = allspikes(:,1); 
    DS = mysort.ds.Matrix(data.X);
    
    % cutspikes
    spikes_cut = DS.getWaveform(r.ts,...
                P.spikeCutting.cutLeft, P.spikeCutting.Tf); 
            
	% estimate noise cov
    if isempty(data.noise)
        s1 = r.ts-P.noiseEstimation.minDistFromSpikes;
        s2 = r.ts+P.noiseEstimation.minDistFromSpikes;
        L  = min(P.noiseEstimation.maxLength, size(DS,1));
        noise.epochs = mysort.epoch.flip(mysort.epoch.merge([s1 s2]), L);

        Cest = mysort.noise.Covest2(DS, 'maxLag', P.featureExtraction.Tf,...
            'maxSamples', L, 'noiseEpochs', noise.epochs, 'forceMethod', 'xcorr');
        noise.C_time = mysort.noise.ccol2Cte(Cest.CCol, P.featureExtraction.Tf);
%         noise.C_time_aligned = mysort.noise.ccol2Cte(Cest.CCol, P.spikeAlignment.Tf);
        noise.meanNoiseStd = sqrt(mean(diag(noise.C_time)));
        noise.CestS = Cest.toStruct();  
        r.noise = noise; clear noise;
    end
    
    % upsample spikes
    spikes_cut = mysort.wf.tResample(mysort.wf.v2t(spikes_cut, size(DS,2)), P.alignment.upsample, 1);
    % remove artefacts not needed for alignment
    spikes_cut([1:P.alignment.upsample_cutEnds*P.alignment.upsample end-P.alignment.upsample_cutEnds*P.alignment.upsample+1:end],:,:) = [];
    spikes_cut = mysort.wf.t2v(spikes_cut);