function dataSources = checkDataSources(dataSources)
    if ~iscell(dataSources)
        dataSources = {dataSources};
    end   
    for i=1:length(dataSources)                    
        if ~isempty(dataSources{i}) 
            if isnumeric(dataSources{i})
                dataSources{i} = mysort.ds.Matrix(dataSources{i});
            else
                assert(isa(dataSources{i}, 'mysort.ds.DataSourceInterface'), 'all datasources must implement the DataSourceInterface!');                                
            end
        end    
    end
