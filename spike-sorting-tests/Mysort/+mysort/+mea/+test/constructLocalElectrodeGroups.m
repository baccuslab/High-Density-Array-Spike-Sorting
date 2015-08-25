
    P.maxElPerGroup = 7;
    P.minElPerGroup = 1;
    P.addIfNearerThan = 20; % always add direct neighbors
%     P.maxOverlapPerGroup = 6;
    P.maxDistanceWithinGroup = 52;  %keep over 51.59 !#
    
    xy = [100 100  
          119 100
          200 100 
          300 100 
          340 100
          501 100
          502 100
          503 100
          504 100
          505 100
          506 100
          507 100
          508 100
          509 100];
[groups nGroupsPerElectrode] = mysort.mea.constructLocalElectrodeGroups(xy(:,1), xy(:,2), 'otherP', P)

