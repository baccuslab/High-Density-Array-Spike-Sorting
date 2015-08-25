function VIEW = plotFischerProj(handles)    
    %warning('this function needs a proper estimate of C!');

    CONF = handles.CONFIG;
    
    ax = handles.FischerProjAxes;
    G = handles.GUI;
    T = G.T;
    C = handles.CONFIG;
    
    T = G.T;
    t = G.selTemplateIdx;

    cla(ax);
    if isempty(t) || isempty(T.templates(t).sourceWfIdx)
        return
    end
    
    set(ax, 'Nextplot', 'replace');

    classes = ones(1, size(T.cutSpikes{t},1));
    classes(T.selIdx{t}) = 2;
    templates = [T.brushed_template_cleaned(t,:); T.selected_template_cleaned(t,:)];
    
%     mcorrSpikes = [D.cutSpikes(~D.selectedIdx,:) - repmat(D.brushedTemplate,  sum(~D.selectedIdx), 1);
%                    D.cutSpikes( D.selectedIdx,:) - repmat(D.selectedTemplate, sum( D.selectedIdx), 1)];
%     
    C = eye(size(templates,2)); iU = C;
    %C = cov(mcorrSpikes); iU = [];

    mysort.plot.clusterProjection(T.cutSpikes{t}, classes, templates, C, ...
        'iU', iU, 'axesHandles', ax, 'binning', [-200:200], 'xbinning', [-200:200], ...
        'colors', {CONF.color_unselected_spikes, CONF.color_selected_spikes});
    xlabel('fischer proj');
    ylabel('hist');
    
    