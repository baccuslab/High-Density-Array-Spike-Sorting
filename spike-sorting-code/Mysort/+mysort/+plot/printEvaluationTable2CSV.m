
function printEvaluationTable2CSV(fname, R)

    x = R.table(:,2:end);
    sx = sum(x(:,2:end)~=0,2);
    x(sx==0,:) = [];
    rlab = [];
    if ~isempty(R.tableRowLabel)
        rlab = R.tableRowLabel(sx~=0);
    end

    % remove gt units with 0 spikes
    keepidx = x(:,2)>0;
    x = x(keepidx,:);
    rlab = rlab(keepidx);
    
    dlmwrite(fname, x, 'delimiter', '\t');