
function spikesX = extractWaveform(X, epochs)  
    spikesX =  [];
    if isempty(epochs)
        return
    end
    assert(size(epochs,2) == 2, 'epochs must be a two column matrix!');
    nSpikes = size(epochs,1);
    nC = size(X,1);
    Tf = epochs(1,2)-epochs(1,1)+1;
    spikesX = zeros(nSpikes, nC*Tf);

    
    for i=1:size(epochs,1)
        startSample = epochs(i,1);
        endSample = epochs(i,2);       
        zL = 0;
        if startSample<=0 
            zL = -startSample+1;
            startSample=1;
        end
        zR = 0;
        if endSample>size(X,2)
            zR = endSample-size(X,2);
            endSample = size(X,2);
        end
        for c=1:nC
            spikesX(i,(1:Tf) + (c-1)*Tf) = ...
                [zeros(1,zL) X(c, startSample:endSample) zeros(1,zR)];
        end
    end
          