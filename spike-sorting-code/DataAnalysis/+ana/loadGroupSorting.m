function S = loadGroupSorting(exp, fname, run_name, groupidx)
spath = '/net/bs-filesvr01/export/group/hierlemann/recordings/HiDens/SpikeSorting/Roska/';

sort_path = fullfile(spath, exp, fname, 'sortings', run_name, sprintf('group%03d', groupidx));

tmp = load(fullfile(sort_path, [run_name '.P.mat']));
S = tmp.S;

filelist = fieldnames(S.files);
for i=1:length(filelist)
    if exist(S.files.(filelist{i}), 'file')
        tmp = load(S.files.(filelist{i}));
        fname = fieldnames(tmp);
        for k=1:length(fname)
            S.(fname{k}) = tmp.(fname{k});
        end
    end
end
