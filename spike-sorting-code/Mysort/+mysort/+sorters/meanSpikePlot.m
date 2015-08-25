% function meanSpikePlot(SO)

    
  length(unique(classes))
  figure; plot(R.BIC)
  mysort.plot.spikes(SO.templates, 'classes', 1:size(SO.templates,1))  

  mysort.plot.clustering(SO.clustering.features, SO.sorting(:,1))
  
  C = eye(size(SO.templates,2));
  Cfet = eye(size(SO.clustering.features,2));      
  
  mysort.plot.clusterProjection(SO.clustering.features, ...
                                SO.sorting(:,1),...
                                centers, Cfet, 'iU', Cfet,...
                                'plotOnlyDiag', 1);
                           
  mysort.plot.clusterProjection(SO.spikes, ...
                                SO.sorting(:,1),...
                                SO.templates, C, 'iU', C);
                            
  mysort.plot.clusters(SO.spikes, SO.sorting(:,1),...
                       C, 'iU', C);
             
  mysort.plot.clusters(SO.clustering.features, SO.sorting(:,1),...
                       Cfet, 'iU', Cfet);                   