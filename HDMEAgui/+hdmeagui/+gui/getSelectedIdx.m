function DATA = getSelectedIdx(DATA, INTER, CONFIG)
    T = DATA.templates;
    t = INTER.tIdx;
    DATA.bSelectedIdxChanged = 0;
    
    if INTER.bBrushingChanged
        T.selIdx{t} = true(sum(T.bIdx(:,t)),1);
        T.unselIdx{t} = ~T.selIdx{t};
    elseif isfield(INTER, 'NormHistThresh') && ...
       ~isempty(INTER.NormHistThresh) && ~isempty(T.spikeProjections{t})
   
        INTER.NormHistThresh = min(INTER.NormHistThresh, max(T.spikeProjections{t}));
        newSelIdx = (T.spikeProjections{t} >= INTER.NormHistThresh);
        
        if isempty(T.selIdx{t}) || length(newSelIdx) ~= length(T.selIdx{t})...
            || any(newSelIdx ~= T.selIdx{t}) 
            fprintf('Selection change.\n');
            T.selIdx{t} = newSelIdx;
            T.unselIdx{t} = ~newSelIdx;
            if ~isempty(INTER.NormHistThresh)
                T.normHistThresh(t) = INTER.NormHistThresh;
            else
                T.normHistThresh(t) = -1;
            end
            DATA.bSelectedIdxChanged = 1;
        end
    end
    DATA.templates = T;
%     fprintf('No selection change.\n');