function DATA = calcCurrentTemplate(DATA, INTER, CONFIG)
    T = DATA.templates;
    t = INTER.tIdx;
    if isempty(t)
        % no template selected
        fprintf('no template selected\n');
        return
    end
    
    nC = DATA.nUsedElectrodes;    
    
    if isempty(INTER.bIdx) || isempty(T.selIdx{t})
       return 
    end
    
    if isempty(T.bIdx)
        T.bIdx = false(length(INTER.bIdx), t);    
    end

    T.bIdx(:, t) = INTER.bIdx;

    if ~any(T.bIdx(:,t))
        return
%         T.selected_template(t,:) = [];
%         T.selected_template_cleaned(t,:) = 
    end
    
    temp = median(T.cutSpikes{t}(T.selIdx{t},:),1);
    T.selected_template(t,:) = temp;
    temp = cleanTemplate(temp);
    T.selected_template_cleaned(t,:) = temp;
    
    temp = mysort.wf.vSubsel(temp, DATA.nUsedElectrodes, (1:20));    
    tempNorm = norm(temp);
    assert(tempNorm>0, 'template norm is not positive!');
    temp = temp./tempNorm;  
    
    T.spikeProjections{t} = mysort.wf.vSubsel(T.cutSpikes{t}, DATA.nUsedElectrodes, (1:20))*temp';               

    if any(T.unselIdx{t})
        temp = median(T.cutSpikes{t},1);
        T.brushed_template(t,:) = temp;
        temp = cleanTemplate(temp);
        T.brushed_template_cleaned(t,:) = temp;              
        
        temp = median(T.cutSpikes{t}(~T.selIdx{t},:),1);
        T.unselected_template(t,:) = temp;
        temp = cleanTemplate(temp);
        T.unselected_template_cleaned(t,:) = temp;                      
    else
        T.brushed_template(t,:) = T.selected_template(t,:);
        T.brushed_template_cleaned(t,:) = T.selected_template_cleaned(t,:);   
        T.unselected_template(t,:) = T.selected_template(t,:);     
        T.unselected_template_cleaned(t,:) = T.selected_template(t,:);    
    end
    
    DATA.templates = T;
    
    
    %----------------------------------------------------------------------
    function t = cleanTemplate(t)
        medM = mysort.wf.v2m(t, nC);
        maxPeak = max(abs(medM), [], 2);
        medMCleaned = medM;
        if isfield(DATA, 'precutSpikes')
            % do nothing
        else
            medMCleaned(maxPeak<.75,:) = 0;
        end
        
        medVCleaned = mysort.wf.m2v(medMCleaned);
        t = medVCleaned;
    end
end