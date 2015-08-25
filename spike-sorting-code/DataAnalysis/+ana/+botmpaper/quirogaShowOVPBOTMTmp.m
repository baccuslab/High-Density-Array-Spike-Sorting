%% 
ovpgdf = [];
for i=1:3
    for k=1:size(peaksOvp{i},1)
        mySpikeTime   = peaksOvp{i}(k,1);
        mySpikeHeight = peaksOvp{i}(k,2);
        mySpikeId     = peakId(i,peaksOvp{i}(k,3));
        keepMe = true;
        for i2=1:3
            if i2==i
                continue
            end
            otherPeaks = find(abs(peaksOvp{i2}(:,1) - mySpikeTime) < 7);
            if isempty(otherPeaks)
                % no other spike here
            elseif peaksOvp{i2}(otherPeaks,2) < mySpikeHeight
                % other spike is smaller
            else
                keepMe = false;
            end 
        end
        if keepMe
            ovpgdf(end+1,:) = [mySpikeId mySpikeTime];
        end
    end
end
[gdf_ovp wasResolved] = mysort.spiketrain.resolveOvpIndexInGdf(ovpgdf, DC.DOvpIndex);
gdf_ovp(:,2) = gdf_ovp(:,2)-cutLeft;

R_ovp = mysort.spiketrain.alignGDFs(tgdf, gdf_ovp, 15, 3, 25);
mysort.plot.printEvaluationTable(R_ovp)
