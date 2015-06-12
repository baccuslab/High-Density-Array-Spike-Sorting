function X = check(X)
    % Performs a check on X if it is just multi channel data - then a
    % dataSource Object is created - or already a datasource object.
    
    if isa(X, 'mysort.datasource.DataSourceInterface')
        return
    end 
    X = mysort.datasource.DataSource(X);