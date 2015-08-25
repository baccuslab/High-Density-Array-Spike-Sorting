function [ntk2 ntk]=ntk_load(ntk_in, load_chunk_size, varargin)
%
% use this function to load a chunck from an ntk file (measured and
% simulated)
%
%


ntkstreamproc_args={};
nofilters=0;
noconfig=0;
verbose=0;

i=1;
while i<=length(varargin)
    if not(isempty(varargin{i}))
        if strcmp(varargin{i}, 'keep_all')
            ntkstreamproc_args{end+1}=varargin{i};  %#ok<AGROW>
        elseif strcmp(varargin{i}, 'keep_disconnected')
            ntkstreamproc_args{end+1}=varargin{i};  %#ok<AGROW>
        elseif strcmp(varargin{i}, 'keep_damaged')
            ntkstreamproc_args{end+1}=varargin{i};  %#ok<AGROW>
        elseif strcmp(varargin{i}, 'keep_dummy')
            ntkstreamproc_args{end+1}=varargin{i};  %#ok<AGROW>
        elseif strcmp(varargin{i}, 'keep_only')
            ntkstreamproc_args{end+1}=varargin{i};  %#ok<AGROW>
            i=i+1;
            ntkstreamproc_args{end+1}=varargin{i};  %#ok<AGROW>
        elseif strcmp(varargin{i}, 'ignore_oversampling')
            ntkstreamproc_args{end+1}=varargin{i};  %#ok<AGROW>
        elseif strcmp(varargin{i}, 'time_align')
            ntkstreamproc_args{end+1}=varargin{i};  %#ok<AGROW>
        elseif strcmp(varargin{i}, 'time_align_precision')
            ntkstreamproc_args{end+1}=varargin{i};  %#ok<AGROW>
            i=i+1;
            ntkstreamproc_args{end+1}=varargin{i}; %#ok<AGROW>
        elseif strcmp(varargin{i}, 'nofiltering')
            ntkstreamproc_args{end+1}=varargin{i};  %#ok<AGROW>
            nofilters=1;
        elseif strcmp(varargin{i}, 'digibits')
            ntkstreamproc_args{end+1}=varargin{i};  %#ok<AGROW>
        elseif strcmp(varargin{i}, 'images_v1')
            ntkstreamproc_args{end+1}=varargin{i};  %#ok<AGROW>
        elseif strcmp(varargin{i}, 'noconfig') %don't generate filters
            noconfig=1;
        elseif strcmp(varargin{i}, 'calc_mean') %calculate prehp mean
            ntkstreamproc_args{end+1}=varargin{i};  %#ok<AGROW>
        elseif strcmp(varargin{i}, 'calc_at_rail') %calculate at rail percentage
            ntkstreamproc_args{end+1}=varargin{i};  %#ok<AGROW>
        elseif strcmp(varargin{i}, 'verbose') %don't generate filters
            verbose=1;
        else
            fprintf('unknown argument at pos %d\n', 2+i);
        end
    end
    i=i+1;
end



% load data
ntkstream_args={};
if noconfig
    ntkstream_args{end+1}='noconfig';
end
if verbose
    ntkstream_args{end+1}='verbose';
end

MAX_CHUNK=3*60*20000;
if load_chunk_size>MAX_CHUNK && isfield(ntk_in, 'nrcmds')
    fprintf('loading data in loop mode, please don''t use frame data in ntk\n');



else
    ntk=ntk_stream(ntk_in,load_chunk_size, ntkstream_args{:});
    ntk.images = ntk_in.images;
    %rotate channels if needed
    ntk=ntk_stream_rotate(ntk);
    ntk2=ntk_stream_process(ntk, ntkstreamproc_args{:});
    ntk.images.last_frame = ntk2.images.last_frame;
end
    





