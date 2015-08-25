%------------------------------------------------------------------
function cp = computeChannelPairs4Channels(channelIdx)
    nC = length(channelIdx);
    nXC = nC*(nC-1)/2+nC;
    cp = zeros(nXC, 2);
    count = 1;
    for i=1:length(channelIdx)
        for j=i:length(channelIdx)
            cp(count,:) = [channelIdx(i) channelIdx(j)];
            count = count+1;
        end
    end
end