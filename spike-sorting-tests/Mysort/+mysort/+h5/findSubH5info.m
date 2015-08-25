function ifo = findSubH5info(h5info, h5path)
    % walks the h5info object and tries to find the subbranch specified by
    % h5path. returns empty object if it cant be found
    ifo = [];
    % check if we have a tree to walk through
    if isempty(h5info)
        return
    end
    
    % check if h5path is empty. that will be treated as the root node
    if isempty(h5path)
        ifo = h5info;
        return
    end
    
    % convert h5path to cell if necessary
    if ischar(h5path)
        h5path = mysort.h5.splitH5path(h5path);
    end
    
    % check if we are the root
    if isfield(h5info, 'GroupHierarchy')
        ifo = mysort.h5.findSubH5info(h5info.GroupHierarchy, h5path);
        return
    end
    
    % check if the current h5info is already the searched one and we dont
    % need to walk further
    idx = strfind(h5info.Name, ['/' h5path{1}]);
    if ~isempty(idx) && (idx+length(['/' h5path{1}])-1 == length(h5info.Name))
        if length(h5path) == 1
            ifo = h5info;
            return
        end
        ifo = mysort.h5.findSubH5info(h5info, h5path(2:end));
        if ~isempty(ifo)
            return
        end
    end
    
    % we are not the searched element, check if one of our Groups is
    if isfield(h5info, 'Groups')
        for i=1:length(h5info.Groups)
            ifo = mysort.h5.findSubH5info(h5info.Groups(i), h5path);
            if ~isempty(ifo)
                return
            end
        end
    end
    % check datasets
    if isfield(h5info, 'Datasets')
        for i=1:length(h5info.Datasets)
            ifo = mysort.h5.findSubH5info(h5info.Datasets(i), h5path);
            if ~isempty(ifo)
                return
            end
        end        
    end
    
