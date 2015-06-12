function H = initInput(H, UD, C, GC)
    % check if the CONFIG exists
    if isfield(UD, 'CONFIG')
        H.CONFIG = mysort.util.parseInputs(C, UD.CONFIG, 'merge');
    else
        H.CONFIG = C;
    end

    % check if the CONFIG exists
    if isfield(UD, 'GUI_CONFIG')
        H.GUI_CONFIG = mysort.util.parseInputs(GC, UD.GUI_CONFIG, 'merge');
    else
        H.GUI_CONFIG = GC;
    end

    if ischar(UD)
        if strcmp(UD(end-2:end), 'mat')
            disp('Loading from File!');
            error('not implemented');
            %L = load(UD);            
        elseif strcmp(UD(end-1:end), 'h5')
            h5file = UD;
            
            % Set data source
            DS = mysort.mea.CMOSMEA(h5file, ...
                'prefilter', 1, ...
                'preprocess', 1, ...
                'hpf', 350, 'lpf', 8000, 'filterOrder', 10);
            
            D = DS.getPreprocessedData();
            TCA = D.timesChansAmpsSorted;
            
            % Create Waveform manager
            H.WF = mysort.wf.WfManagerBuffCut(DS, TCA(:,1), TCA(:,2), ...
                H.GUI_CONFIG.CutSpikesCutleft, H.GUI_CONFIG.CutSpikesTf);
            
            % Register features
            mima = mysort.wf.feature.MinMax(H.WF);
            mima.needsUpdate = false;
            mima.F = -D.timesChansAmpsSorted(:,3:3);
            H.WF.registerFeature(mima);
            
            % Create Template Manager
            H.T = mysort.wf.TemplateManager(H.WF);
            
        else
            error('not supported');
        end
    elseif isa(UD, 'mysort.wf.WfManager')
        H.WF = UD;
        H.T = mysort.wf.TemplateManager(H.WF);
    else
        error('not supported');
    end