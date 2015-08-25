
function idx = getBurstIndex(spiketrain, maxdist)
    % Calculates for the spikes in the spike train the index of that spike
    % in a possible burst. A single spike will get a one. A spike that is
    % right after another spike, gets a 2. A spike that is very close after
    % a spike with index n gets the index n+1
    if size(spiketrain,1)>size(spiketrain,2)
        spiketrain = spiketrain';
    end
    
    diffs = [0 diff(spiketrain)];
    idx = ones(1, length(spiketrain));
    for i=2:length(spiketrain)
        if diffs(i) <= maxdist
            idx(i) = idx(i-1)+1;
        end
    end
        