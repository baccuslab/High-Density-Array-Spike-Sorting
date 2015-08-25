function savefigs(fhs, dpath, varargin)
    % Saves multiple figures (handles in fhs) to the folder in dpath
    for i=1:length(fhs)
        name = get(fhs(i), 'name');
        mysort.plot.savefig(fhs(i), fullfile(dpath, name), varargin{:});
    end