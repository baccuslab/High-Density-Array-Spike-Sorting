
function str = printEvaluationTable(R, varargin)
    P.ignoreEmptyGT = false;
    P.ignoreEmptySorted = true;
    P.printPercTpDet = true;
    P.printPercTpGt  = true;
    [P uP] = mysort.util.parseInputs(P, varargin, 'split');
    uP = mysort.util.deflateP(uP);

    x = R.table(:,2:end);
    
    sx = sum(x(:,2:end)~=0,2);
    x(sx==0,:) = [];
    rlab = [];
    if ~isempty(R.tableRowLabel)
        rlab = R.tableRowLabel(sx~=0);
    end
    
    colLabels = R.tableLabelShort(2:end);
    % remove gt units with 0 spikes and with 0 detected spikes
    if P.ignoreEmptyGT
        keepidx = x(:,2)>0;
        x = x(keepidx,:);
        rlab = rlab(keepidx);
    end
    if P.ignoreEmptySorted
        keepidx = x(:,5)>0;
        x = x(keepidx,:);
        rlab = rlab(keepidx);
    end
    if P.printPercTpDet
        pdet = round(100*x(:,11)./x(:,5));
        pdet(isnan(pdet)) = 0;
        x = [x(:, 1:11) pdet x(:,12:end)];
        colLabels = [colLabels(1:11) 'Det%' colLabels(12:end)];
    end
    if P.printPercTpGt
        pgt = round(100*x(:,11)./x(:,2));
        pgt(isnan(pgt))=0;
        x = [x(:, 1:11) pgt x(:,12:end)];
        colLabels = [colLabels(1:11) 'Gt%' colLabels(12:end)];
    end
    
    if nargout > 0 
        str = mysort.plot.printTable(x, 'colLabel', colLabels, ...
                                            'rowLabel', rlab, ...
                                            'printColSum',1, ...
                                            'vlineAfterCols', 10,...
                                            'removeZeroRows', 0,...
                                            uP{:});
    else
        mysort.plot.printTable(x, 'colLabel', colLabels, ...
                                            'rowLabel', rlab, ...
                                            'printColSum',1, ...
                                            'vlineAfterCols', 10,...
                                            'removeZeroRows', 0,...
                                            uP{:});        
    end