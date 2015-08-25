
function str = printTable(x, varargin)
    P.rowLabel = {};
    P.colLabel = {};
    P.ylabelTop = {};
    P.nDecSpaces = 2;
    P.topLeftLabel = ' ';
    P.markRowMax = [];
    P.markColMax = [];
    P.markRowMin = [];
    P.markColMin = [];
    P.hlineAfterRows = [];
    P.vlineAfterCols = [];
    P.printColMean = 0;
    P.printColSum  = 0;
    P.repreatHeaderAtEnd = 0;
    P.headerSeparator = '|';
    P.removeZeroRows = 0;
    
    P = mysort.util.parseInputs(P, 'printTable', varargin);
    
    str = '';
    
    if P.removeZeroRows
        sx = sum(x~=0,2);
        x(sx==0,:) = [];
        if ~isempty(P.rowLabel)
            P.rowLabel(sx==0) = {};
        end
    end
    
    nRows = size(x,1);
    if ~isempty(P.rowLabel)
        assert(length(P.rowLabel) == nRows, sprintf('Number of rows and row labels must be identical! (label:%d, rows:%d)', length(P.rowLabel), nRows));
    end
    nCols = size(x,2);
    if ~isempty(P.colLabel)
        assert(length(P.colLabel) == nCols, 'Number of columns and column labels must be identical!');
    end

    

    rowLabelLen = 0; rowLabFormat = [];
    if ~isempty(P.rowLabel)
        rowLabelLen = max(max(cellfun(@length,P.rowLabel), length(P.topLeftLabel)));
        rowLabelFormat = ['%-' num2str(rowLabelLen) 's|'];
    else
        rowLabelFormat = '';
    end
    
    [colDataFormat colLabelFormat maxLens] = getPrintFormats();
    
    headerStr = buildHeaderStr();
    lineStr   = buildLineStr();

    str = [str sprintf([headerStr '\n'])];
    str = [str sprintf([lineStr '\n'])];
    
    maxColDF = [use_str P.markColMax];
    maxRowDF = [use_str P.markRowMax];
    minColDF = [use_str P.markColMin];
    minRowDF = [use_str P.markRowMin];
    
    if markerLen == 0
        dataSpacer = '';
    else
        dataSpacer = sprintf(['%' num2str(markerLen) 's' ],' ');
    end
    
%     [m colMaxis] = max(x,[],1);
%     [m colMinis] = min(x,[],1);
    for i=1:nRows
        rowStr = '';
        if ~isempty(P.rowLabel)
            rowStr = sprintf(rowLabelFormat, P.rowLabel{i});
        end
        str = [str printRow(x(i,:))];
        if any(P.hlineAfterRows==i)
            str = [str sprintf([lineStr '\n'])];
        end
    end
    if P.printColMean
       m = mean(x,1);
       rowStr = sprintf(rowLabelFormat, 'Mean');   
       str = [str sprintf([lineStr '\n'])];
       str = [str printRow(m)];
    end
    if P.printColSum
       m = sum(x,1);
       rowStr = sprintf(rowLabelFormat, 'Total');   
       str = [str sprintf([lineStr '\n'])];
       str = [str printRow(m)];
    end
    if P.repreatHeaderAtEnd
        str = [str headerStr];
    end
    if nargout > 0
        return
    else
        disp(str);
    end
    
    %----------------------------------------------------------------------
    function [colDataFormat colLabelFormat maxLens] = getPrintFormats()
        % estimate length of fields for every column
        maxLens = zeros(1,nCols);
        if ~isempty(P.colLabel)
            maxLens = cellfun(@length,P.colLabel)+1;
        end
        f_str = ['%.' sprintf('%d',P.nDecSpaces) 'f '];
        i_str = ['%' sprintf('%d',P.nDecSpaces) 'd '];
        s_str = ['%s '];

        markerLen = length(P.markRowMax)+length(P.markRowMin);
        colLabelFormat = {}; colDataFormat = {};
        for ii=1:nCols
            for j=1:nRows
                use_str = f_str;
                if iscell(x) 
                    if ischar(x{j,ii})
                        use_str = s_str;
                    elseif int32(x{j,ii}) == x{j,ii}
                        use_str = i_str;
                    end
                    tmp = sprintf(use_str, x{j,ii});
                else
                    if int32(x(j,ii)) == x(j,ii)
                        use_str = i_str;
                    end
                    tmp = sprintf(use_str, x(j,ii));
                end
                
                maxLens(ii) = max(maxLens(ii), length(tmp)+markerLen);
            end  

            if P.printColSum
                use_str = f_str;
                if int32(sum(x(:,ii))) == sum(x(:,ii))
                    use_str = i_str;
                end
                tmp = sprintf(use_str, sum(x(:,ii)));
                maxLens(ii) = max(maxLens(ii), length(tmp)+markerLen);            
            end
            if any(P.vlineAfterCols==ii)
                spacer = '|';
            else
                spacer = ' ';
            end
            if iscell(x)
                if ischar(x{1,ii})
                    colDataFormat{ii} = ['%' sprintf('%d',maxLens(ii)-1) 's' spacer];
                elseif int32(x{1,ii}) == x{1,ii}
                    colDataFormat{ii} = ['%' sprintf('%d',maxLens(ii)-1) 'd' spacer];
                else
                    colDataFormat{ii} = ['%' sprintf('%d',maxLens(ii)-1) '.' sprintf('%d',P.nDecSpaces) 'f' spacer];
                end                
            else
                if int32(x(1,ii)) == x(1,ii)
                    colDataFormat{ii} = ['%' sprintf('%d',maxLens(ii)-1) 'd' spacer];
                else
                    colDataFormat{ii} = ['%' sprintf('%d',maxLens(ii)-1) '.' sprintf('%d',P.nDecSpaces) 'f' spacer];
                end
            end
            colLabelFormat{ii} = ['%' num2str(maxLens(ii)-1) 's'];
        end        
    end

    %----------------------------------------------------------------------        
    function S = printRow(r)
        S = '';

%         m = min(r); 
%         rowMin = find(r==m);
        for k=1:nCols
            if ~isempty(P.markRowMax)
                m = max(r);
                rowMax = find(r==m);
                if any(k==rowMax) 
                    rowStr = [rowStr sprintf(maxRowDF, r(k))];
                end
%             elseif colMaxis(k)==i && ~isempty(P.markColMax)
%                 rowStr = [rowStr sprintf(maxColDF, str_x{r,k})];                            
%             elseif any(k==rowMin) && ~isempty(P.markRowMin)
%                 rowStr = [rowStr sprintf(minRowDF, str_x{r,k})];
%             elseif colMinis(k)==i && ~isempty(P.markColMin)
%                 rowStr = [rowStr sprintf(minColDF, str_x{r,k})];  
            elseif k<nCols
                if iscell(r)
                    rowStr = [rowStr sprintf(colDataFormat{k}, r{k}) dataSpacer];
                else
                    rowStr = [rowStr sprintf(colDataFormat{k}, r(k)) dataSpacer];
                end
            else
                if iscell(r)
                    rowStr = [rowStr sprintf(colDataFormat{k}, r{k})]; 
                else
                    rowStr = [rowStr sprintf(colDataFormat{k}, r(k))]; 
                end
            end
        end        
        S = sprintf([rowStr '|\n']);
    end

    %----------------------------------------------------------------------
    function headerStr = buildHeaderStr()
        headerStr = '';
        if ~isempty(P.colLabel)
            if ~isempty(P.rowLabel)
                headerStr = sprintf(rowLabelFormat, P.topLeftLabel);
            end
            for ii=1:nCols
                if ii<nCols
                    headerStr = [headerStr sprintf(colLabelFormat{ii}, strrep(P.colLabel{ii}, '%', '%%')) P.headerSeparator];
                else
                    headerStr = [headerStr sprintf(colLabelFormat{ii}, strrep(P.colLabel{ii}, '%', '%%')) ' ' P.headerSeparator];
                end
            end
    %         headerStr = [headerStr '|'];
        end
    end

    %----------------------------------------------------------------------
    function lineStr = buildLineStr()
        lineStr = '';
        if ~isempty(P.rowLabel)
            for ii=1:length(sprintf(rowLabelFormat, P.topLeftLabel))-1
                lineStr = [lineStr '-'];
            end
            lineStr = [lineStr '|'];
        end    
        for ii=1:sum(maxLens)
            lineStr = [lineStr '-'];
        end
        lineStr = [lineStr '|'];
    end
end