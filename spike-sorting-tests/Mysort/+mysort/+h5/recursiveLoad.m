
function D = recursiveLoad(fname, hinfo)
    if nargin == 1
        % No hdf5 path or hinfo object given, load whole file
        hinfo = hdf5info(fname);
        for i=1:length(hinfo.GroupHierarchy)
            D{i} = mysort.h5.recursiveLoad(fname, hinfo.GroupHierarchy(i));
        end
    elseif ischar(hinfo)
        str = hinfo;
        hinfo = hdf5info(fname);
        hinfo = getHinfoFromString(hinfo.GroupHierarchy, str);
        if isempty(hinfo)
            error('Could not find specified h5path in GroupHierarchy: %s (file: %s)', str, fname);
        end
        D = mysort.h5.recursiveLoad(fname, hinfo);
    else
        isVar = true;
        if isfield(hinfo, 'Datasets')
            isVar = false;
            for i=1:length(hinfo.Datasets)
                [tmp attr] = hdf5read(hinfo.Datasets(i), ...
                    'ReadAttributes',true, 'V71Dimensions', true);
                if ndims(tmp)==3
                    tmp = permute(tmp,[3 1 2]);
                end
                if isa(tmp, 'hdf5.h5string') 
                    if length(tmp)>1
                        TMP = {};
                        for k=1:length(tmp)
                            TMP{k} = char(tmp(k).Data);
                        end
                    else
                        TMP = char(tmp.Data);
                    end
                else
                    TMP = tmp;
                end
                slashidx = strfind(hinfo.Datasets(i).Name, '/');
                fieldn = hinfo.Datasets(i).Name(max(slashidx)+1:end);
                if ischar(fieldn) && ~isempty(strfind(fieldn, '____EMPTY____'))
                    fieldn = fieldn(1:strfind(fieldn, '____EMPTY____')-1);
                    TMP = [];
                end
                [x] = str2double(fieldn);
                if ~isnan(x) && ~isempty(x) && isnumeric(x)
                    if isstruct(TMP)
                        D(x) = TMP;
                    else
                        D{x} = TMP;
                    end
                else
                    D.(fieldn) = TMP;
                end  
            end
        end
        if isfield(hinfo, 'Groups')
            isVar = false;
            for i=1:length(hinfo.Groups)
                slashidx = strfind(hinfo.Groups(i).Name, '/');
                fieldn = hinfo.Groups(i).Name(max(slashidx)+1:end);  
                TMP = mysort.h5.recursiveLoad(fname, hinfo.Groups(i));
                [x status] = str2num(fieldn);
                if status
                    if isstruct(TMP)
                        D(x) = TMP;
                    else
                        D{x} = TMP;
                    end
                else
                    classtoken = strfind(fieldn, '____CLASS____');
                    if ~isempty(classtoken)
                        thisclass = fieldn(classtoken+length('____CLASS____'):end);
                        fieldn = fieldn(1:classtoken-1);
                        eval(['TMP = ' thisclass '(''RESTORE_FROM_STRUCT'', TMP);']);
                    end
                    D.(fieldn) = TMP;
                end  
            end        
        end
        if isVar
            slashidx = strfind(hinfo.Name, '/');
            fieldn = hinfo.Name(max(slashidx)+1:end);  
            tmp = hdf5read(hinfo, 'V71Dimensions', true);
            if isa(tmp, 'hdf5.h5string') 
               if length(tmp)>1
                    TMP = {};
                    for k=1:length(tmp)
                        TMP{k} = char(tmp(k).Data);
                    end
                else
                    TMP = char(tmp.Data);
                end
            else
                TMP = tmp;
            end            

            [x status] = str2num(fieldn);
            if status
                if isstruct(TMP)
                    D(x) = TMP;
                else
                    D{x} = TMP;
                end
            else
                classtoken = strfind(fieldn, '____CLASS____');
                if ~isempty(classtoken)
                    thisclass = fieldn(classtoken+length('____CLASS____'):end);
                    fieldn = fieldn(1:classtoken-1);
                    eval(['TMP = ' thisclass '(''RESTORE_FROM_STRUCT'', TMP);']);
                end
                D = TMP;
            end  
        end
    end
    
    %%% ------------------------------------------------------
    function D = storeData()
        [x status] = str2num(fieldn);
        if status
            if isstruct(TMP)
                D(x) = TMP;
            else
                D{x} = TMP;
            end
        else
            D.(fieldn) = TMP;
        end        
    end
    %%% ------------------------------------------------------
    function hinfo = getHinfoFromString(hinfo, str) 
        if isempty(str)
            return
        end
        [searchname rest] = getFirstName(str);
        if isfield(hinfo, 'Name') && (strcmp(hinfo.Name, str) || strcmp(hinfo.Name, searchname))
            return
        end
        if isfield(hinfo, 'Datasets')
            for f=1:length(hinfo.Datasets)
                lastname = getLastName(hinfo.Datasets(f).Name);
                if strcmp(lastname, searchname)
                    hinfo = getHinfoFromString(hinfo.Datasets(f), rest);
                    return
                end
            end
        end
        if isfield(hinfo, 'Groups')
            for f=1:length(hinfo.Groups)
                lastname = getLastName(hinfo.Groups(f).Name);
                if strcmp(lastname, searchname)
                    hinfo = getHinfoFromString(hinfo.Groups(f), rest);
                    return
                end
            end
        end
        hinfo = [];
    end
    %%% ------------------------------------------------------
    function [firstname rest]= getFirstName(str)
        slashidx = strfind(str, '/');
        s1 = slashidx(1)+1;
        if length(slashidx) == 1
            s2 = length(str);
        else
            s2 = slashidx(2)-1;
        end  
        firstname = str(s1:s2);
        rest  = str(s2+1:end);
    end
    %%% ------------------------------------------------------
    function lastname= getLastName(str)
        slashidx = strfind(str, '/');
        lastname = str(max(slashidx)+1:end);
    end
end
