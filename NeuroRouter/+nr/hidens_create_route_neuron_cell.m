function npos = hidens_create_route_neuron_cell(electrodeMatlabIdx, stimMatlabIdx, electrodePriorityCosts)
    if nargin < 3
        electrodePriorityCosts = -ones(length(electrodeMatlabIdx), 1);
    end
    els = nr.hidens_get_all_electrodes(2);
    npos = {};
    for i=1:length(electrodeMatlabIdx)
        npos{i}.el_idx = els.el_idx(electrodeMatlabIdx(i));
        npos{i}.x = els.x(electrodeMatlabIdx(i));
        npos{i}.y = els.y(electrodeMatlabIdx(i));
        npos{i}.label = sprintf('neuron%d', i);
        if electrodePriorityCosts(i) ~= -1
            npos{i}.cost = electrodePriorityCosts(i);
        end
        if any(electrodeMatlabIdx(i)==stimMatlabIdx)
            npos{i}.isStim = true;
        else
            npos{i}.isStim = false;
        end
    end