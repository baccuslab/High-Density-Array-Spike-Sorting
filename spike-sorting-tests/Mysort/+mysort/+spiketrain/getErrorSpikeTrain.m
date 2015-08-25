function errorSt = getErrorSpikeTrain(evaluation, errorType)
    def = mysort.util.defs();

    if strcmp(errorType, 'FP')
        TYPE = def.fp;
        st = evaluation.St2;
        label = evaluation.spikeLabel2;
    elseif strcmp(errorType, 'FN')
        TYPE = def.fn;
        st = evaluation.St1;
        label = evaluation.spikeLabel1;
    elseif strcmp(errorType, 'FNO')
        TYPE = def.fno;
        st = evaluation.St1;
        label = evaluation.spikeLabel1;           
    elseif strcmp(errorType, 'CL')
        TYPE = def.cl;
        st = evaluation.St1;
        label = evaluation.spikeLabel1; 
    elseif strcmp(errorType, 'CLO')
        TYPE = def.clo;
        st = evaluation.St1;
        label = evaluation.spikeLabel1;                 
    else
        error('Unknown Error Type!');
    end
    errorSt = {};
    for i=1:length(st)
        idx = find(label{i}==TYPE);
        errorSt{i} = st{i}(idx);
    end            
end