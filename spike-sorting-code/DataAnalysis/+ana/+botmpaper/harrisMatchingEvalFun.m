function R = harrisMatchingEvalFun(T, cutSpikeWaveforms, correctSpikes, noiseSnippets, C)
    nSp = size(cutSpikeWaveforms,1);
    nSpCor = size(correctSpikes,1);
    nT = size(T,1);
    nNoise = size(noiseSnippets,1);
    
    R.Eucl = [];
    R.EuclCor = [];
    R.EuclNoise = [];

    R.Maha = [];
    R.MahaCor = [];
    R.MahaNoise = [];

    R.Conv = [];
    R.ConvCor = [];
    R.ConvNoise = [];

    R.Matched = [];
    R.MatchedCor = [];
    R.MatchedNoise = [];

    R.botm = [];
    R.botmCor = [];
    R.botmNoise = [];

    R.nEucl = [];
    R.nEuclCor = [];
    R.nEuclNoise = [];

    R.nConv = [];
    R.nConvCor = [];
    R.nConvNoise = [];

    R.nBotm = [];
    R.nBotmCor = [];
    R.nBotmNoise = [];
%     C_ = .5*C + .5*diag(diag(C));
%     R.U_ = chol(R.C_);
    C_ = C;
    edges = linspace(-10e8, 10e8, 20000);
    botm_edges = linspace(-10e4, 10e4, 20000);
    
    for i=1:nT
        F = C_\T(i,:)';

%         R.Eucl(:,i)    = sum((cutSpikeWaveforms-repmat(T(i,:), nSp, 1)).^2,2);
%         R.EuclCor(:,i) = sum((correctSpikes-repmat(T(i,:), nSpCor, 1)).^2,2);
%         R.EuclNoise(:,i) = sum((noiseSnippets-repmat(T(i,:), nNoise, 1)).^2,2);
% 
%         R.Maha(:,i)    = sum( ((cutSpikeWaveforms-repmat(T(i,:), nSp, 1))/C_).*(cutSpikeWaveforms-repmat(T(i,:), nSp, 1)) ,2);
%         R.MahaCor(:,i) = sum( ((correctSpikes-repmat(T(i,:), nSpCor, 1))/C_).*(correctSpikes-repmat(T(i,:), nSpCor, 1)) ,2);
%         R.MahaNoise(:,i) = sum( ((noiseSnippets-repmat(T(i,:), nNoise, 1))/C_).*(noiseSnippets-repmat(T(i,:), nNoise, 1)) ,2);
% 
%         R.Conv(:,i)    = cutSpikeWaveforms*T(i,:)';
%         R.ConvCor(:,i) = correctSpikes*T(i,:)';
%         R.ConvNoise(:,i) = noiseSnippets*T(i,:)';
% 
%         R.Matched(:,i)    = cutSpikeWaveforms*F;
%         R.MatchedCor(:,i) = correctSpikes*F;
%         R.MatchedNoise(:,i) = noiseSnippets*F;    

        R.botm(:,i)    = cutSpikeWaveforms*F - .5*T(i,:)*F;
        R.botmCor(:,i) = correctSpikes*F  - .5*T(i,:)*F;
        R.botmNoise(:,i) = noiseSnippets*F  - .5*T(i,:)*F;

%         R.edges = linspace(-10e8, 10e8, 20000);
%         R.botm_edges = linspace(-10e4, 10e4, 20000);
% 
%         R.nEucl(:,i) = histc(R.Eucl(:,i), edges);
%         R.nEuclCor(:,i) = histc(R.EuclCor(:,i), edges);
%         R.nEuclNoise(:,i) = histc(R.EuclNoise(:,i), edges);
% 
%         R.nMaha(:,i) = histc(R.Maha(:,i), botm_edges);
%         R.nMahaCor(:,i) = histc(R.MahaCor(:,i), botm_edges);
%         R.nMahaNoise(:,i) = histc(R.MahaNoise(:,i), botm_edges);
% 
%         R.nConv(:,i) = histc(R.Conv(:,i), edges);
%         R.nConvCor(:,i) = histc(R.ConvCor(:,i), edges);
%         R.nConvNoise(:,i) = histc(R.ConvNoise(:,i), edges);    
% 
%         R.nMatched(:,i) = histc(R.Matched(:,i), botm_edges);
%         R.nMatchedCor(:,i) = histc(R.MatchedCor(:,i), botm_edges);
%         R.nMatchedNoise(:,i) = histc(R.MatchedNoise(:,i), botm_edges);
% 
%         R.nBotm(:,i) = histc(R.botm(:,i), botm_edges);
%         R.nBotmCor(:,i) = histc(R.botmCor(:,i), botm_edges);
%         R.nBotmNoise(:,i) = histc(R.botmNoise(:,i), botm_edges);
    end
    R.C_ = C;
end
