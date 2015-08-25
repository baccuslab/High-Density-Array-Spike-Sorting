function P = waveforms2D(wfs, electrodePositions, varargin)
% wfs is tensor time x channels x items
% electrodePositions is matrix, first column x, second column y
P.AxesHandle = [];
P.IDs = [];
P.absThreshold = [];
P.maxNumberOfChannels = [];
P.plotArgs = {'-', 'color', [.5 .5 .5], 'linewidth', 1};
P.plotAllTraces = 1;
P.plotMedian = 0;
P.plotMedianPlotArgs = {'-', 'color', [.0 .0 .0], 'linewidth', 2};
P.plotMeanConfidenceWithSigma = [];
P.plotElNumbers = [];
P.maxWaveforms = 3000;
P.channelIdx = [];
P.scaling = [];
P.width = 15;
P = mysort.util.parseInputs(P, varargin, 'error');
%assert(size(wfs,2) == size(electrodePositions,1), 'electrode number does not match!');
assert(isempty(P.IDs) || length(P.IDs) == 1 || length(P.IDs) == size(wfs,3), 'Length of IDs does not match number of items in wfs')

if isempty(P.AxesHandle)
    P.AxesHandle = axes();
end

if ~iscell(P.plotArgs)
    P.plotArgs = {P.plotArgs};
end

if isempty(wfs)
    return
end

% prepare IDs
if isempty(P.IDs)
    P.IDs = ones(1, size(wfs,3));
elseif length(P.IDs) == 1
    P.IDs = ones(1, size(wfs,3))*P.IDs;
end

uIDs = unique(P.IDs);
nU = length(uIDs);

if ~isempty(P.absThreshold)
    maxs = max(abs(wfs),[],1);
    electrodePositions(maxs<P.absThreshold,:) = [];
    wfs(:,maxs<P.absThreshold,:) = [];
end

if ~isempty(P.maxNumberOfChannels) && size(wfs,2) > P.maxNumberOfChannels
    mwfs = max(abs(wfs),[],3);
    maxs = max(abs(mwfs),[],1);
    maxs = sortrows([maxs(:) (1:length(maxs))'], 1);
    idx = maxs(end-P.maxNumberOfChannels+1:end,2);
    electrodePositions = electrodePositions(idx,:);
    wfs = wfs(:,idx,:);
end

if ~isempty(P.maxWaveforms) & P.maxWaveforms < size(wfs,3)
    rp = randperm(size(wfs,3));
    wfs = wfs(:,:, rp(1:P.maxWaveforms));
end

EP = electrodePositions;
[Tf, nC, nWf] = size(wfs);
pIdx = nan((Tf+1)*nC,1);
pIdy = nan((Tf+1)*nC,1);
Y = nan((Tf+1)*nC,nWf);
if isempty(P.scaling)
    P.scaling = 15/max(abs(wfs(:)));
end
for i=1:nC
    s1 = ((Tf+1)*(i-1)+1);
    idx = s1:s1+Tf-1;
    pIdx(idx,1) = EP(i,1) + P.width*(0:Tf-1)/Tf;
    Y(idx,:)   = EP(i,2) - squeeze(wfs(:,i,:))*P.scaling;
end

set(P.AxesHandle, 'NextPlot', 'add');

if P.plotAllTraces
for u = 1:nU
    if nU > 1
        col = 0.75 + 0.25* mysort.plot.vectorColor(uIDs(u));
        P.plotArgs = [P.plotArgs{:}, {'color'}, {col}];
    end
    h = plot(P.AxesHandle, pIdx, Y(:,P.IDs == uIDs(u) ), P.plotArgs{:} );
end
end

if P.plotMedian && ~isempty(P.plotMedianPlotArgs)
    hold on
    for u = 1:nU
        MED = median( Y(:,P.IDs == uIDs(u) ), 2);
        
        if nU > 1
            col = mysort.plot.vectorColor(uIDs(u));
            P.plotMedianPlotArgs = [P.plotMedianPlotArgs{:}, {'color'}, {col}];
        end
        
        h(2) = plot(P.AxesHandle, pIdx, MED, P.plotMedianPlotArgs{:});
        if ~isempty(P.plotMeanConfidenceWithSigma)
            var_mean = P.plotMeanConfidenceWithSigma/sqrt(nWf);
            h(3) = plot(P.AxesHandle, pIdx, MED+3*var_mean, P.plotMedianPlotArgs{:});
            h(3) = plot(P.AxesHandle, pIdx, MED-3*var_mean, P.plotMedianPlotArgs{:});
        end
    end
end
if ~isempty(P.plotElNumbers)
    for i=1:nC
        x = EP(i,1);
        y  = EP(i,2);
        text(x+P.width, y, num2str(P.plotElNumbers(i)), 'parent', P.AxesHandle);
    end
end
if ~isempty(P.channelIdx)
    for i = 1:length(P.channelIdx)
        
        
        %P.referenceElectrodes(i)
        %el = P.referenceElectrodes(i)
        %plot(P.AxesHandle, pIdx(P.referenceElectrodes(i)), -10*ones(10,1) , P.plotMedianPlotArgs{:});
        %x = pIdx(P.referenceElectrodes(i))
        %y = pIdy(P.referenceElectrodes(i))
        
        %plot(P.AxesHandle, pIdx(P.referenceElectrodes(i)), 
        rectangle('Position', [ EP(P.channelIdx(i), 1) EP(P.channelIdx(i), 2)-0.5*P.width P.width P.width ] );% elW elH ]);
    end
end

axis(P.AxesHandle, 'tight');
axis(P.AxesHandle, 'ij');

P.h = h;
P.pIdx = pIdx;