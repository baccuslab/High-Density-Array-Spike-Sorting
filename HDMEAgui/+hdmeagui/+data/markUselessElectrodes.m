function DATA = markUselessElectrodes(DATA, INTER, CONFIG)
    fprintf('Need to mark useless electrodes...\n');
    DATA.bNeedCutSpikes = 1;
    
    nC = DATA.nC;
    DATA.useElectrodes = ones(1,nC);
    
    for i=1:nC
        sp = find(DATA.singleChannelDataSorted(:,2)==i);
        if isempty(sp)
            DATA.useElectrodes(i) = 0;
        elseif 10 > length(sp) && 5 > sum(DATA.singleChannelDataSorted(:,3)-4.5)  
            DATA.useElectrodes(i) = 0;
        end
    end
    fprintf('Done (%d).\n', sum(DATA.useElectrodes));