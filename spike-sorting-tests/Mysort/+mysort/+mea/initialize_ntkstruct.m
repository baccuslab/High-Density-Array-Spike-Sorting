function y=initialize_ntkstruct(fname_or_gntk, varargin)
%
% use this function to initilize an ntk struct
%
%
%
%

filter=[];
hhp=[];
llp=[];
bbp=[];
ntkavg=[];
lpf=4400;
hpf=5;
nofilters=0;
noconfig=0;
verbose=0;
keep_discarded=0;
use_local_filters=0;

i=1;
while i<=length(varargin)
    if strcmp(varargin{i}, 'filter')
        i=i+1;
        filter=varargin{i};
    elseif strcmp(varargin{i}, 'hhp')  %dfilt highpass
        i=i+1;
        hhp=varargin{i};
    elseif strcmp(varargin{i}, 'llp')  %dfilt lowpass
        i=i+1;
        llp=varargin{i};
    elseif strcmp(varargin{i}, 'bbp')  %dfilt bandpass
        i=i+1;
        llp=varargin{i};
    elseif strcmp(varargin{i}, 'lpf') %1st order lowpass cutoff
        i=i+1;
        lpf=varargin{i};
    elseif strcmp(varargin{i}, 'hpf') %1st order highpass cutoff
        i=i+1;
        hpf=varargin{i};
    elseif strcmp(varargin{i}, 'nofilters') %don't generate filters
        nofilters=1;
	elseif strcmp(varargin{i}, 'use_local_filters') %use locally stored ones
        use_local_filters=1;    
    elseif strcmp(varargin{i}, 'noconfig')
        noconfig=1;
    elseif strcmp(varargin{i}, 'verbose')
        verbose=1;
    elseif strcmp(varargin{i}, 'keep_discarded') %keep file even though it has been discarded
        keep_discarded=1;
    else
        fprintf('unknown argument at pos %d\n', 2+i);
    end
    i=i+1;
end

if ischar(fname_or_gntk)
    y.fname=fname_or_gntk;
    if not(exist(y.fname, 'file'))
        y.fname=[ '../proc/' fname_or_gntk];
        if not(exist(y.fname, 'file'))
            error('file %s does not exist\n',  y.fname);
            return
        end
    end
    
    [pathstr, name, ext, versn] = fileparts(y.fname); %#ok<NASGU>
    if strcmp(ext, '.ntk')
        %load ntk file
        ntkstream_args={};
        if noconfig
            ntkstream_args{end+1}='noconfig';
        end
        if verbose
            ntkstream_args{end+1}='verbose';
        end
        y=ntk_stream(y,0, ntkstream_args{:});
        y.timestamp(11)=' ';    % convert into MATLAB-compatible standard
        y=load_params(y, 'all');
        if isfield(y.recfile_param, 'discard') && ~keep_discarded
            if y.recfile_param.discard
                fprintf('******************************\nskip file %s, discarded\n******************************\n', y.fname);
                y=[];
                y.eof=true;
                return
            end
        end
    elseif strcmp(ext, '.mat')   % load testbench 'mat' file
        
        z=load(y.fname);
        y=z.gntk_data1;
        y.fname=fname_or_gntk;
        %runde1=1;
        y.simul=1;
        y.recfile_param=[];
        y.channel_nr=1:length(y.el_idx);
        y.sr=20000;
        y.pos=1;
        y.eof=0;
    else
        error('unknown neurodata file extention')
    end
else
    y=fname_or_gntk;
    if not(isfield(y,'fname'))
        fprintf('Warning: Inexistent field ''fname'' in structure gntk\n');
        y.fname='';
    end
    y.simul=1;
    y.recfile_param=[];
    y.channel_nr=1:length(y.el_idx);
    y.sr=20000;
    y.pos=1;
    y.eof=false;
end

%set the filers
if isempty(filter) && not(nofilters)
    if exist('fdesign') && ~use_local_filters
        if isempty(hhp) && isempty(llp) && isempty(bbp)
            bp=fdesign.bandpass('n,f3dB1,f3dB2', 2, hpf, lpf, y.sr);
            %y.filters.bbp=butter(bp);
            y.filters.bbp=design(bp,'butter','sosscalenorm','l1');
            % [b, a] = sos2tf(y.filters.bbp.sosMatrix, prod(y.filters.bbp.ScaleValues));
            % y.filters.bp.b = b;
            % y.filters.bp.a = a;
            % To filter data do output = filter(y.filters.bbp.b, y.filters.bbp.a, data)
        else
            if not(isempty(bbp))    % is this needed...?
                y.filters.bbp=bbp;
            end
            if isempty(hhp)
                hp = fdesign.highpass('n,fst,ast',1,hpf,3,y.sr);
                %y.filters.hhp=cheby2(hp);
                y.filters.hhp=design(hp,'cheby2','sosscalenorm','l1');
            else
                y.filters.hhp=hhp;
            end
            if isempty(llp)
                lp = fdesign.lowpass('n,fst,ast',1,lpf,3,y.sr);
                %y.filters.llp=cheby2(lp);
                y.filters.llp=design(lp,'cheby2','sosscalenorm','l1');
            else
                y.filters.llp=llp;
            end
            
        end
        if not(isempty(ntkavg))
            y.filters.ntkavg.iirFactor=ntkavg(1);
            y.filters.ntkavg.step=ntkavg(2);
        end
    else
        [b a]=get_local_filter(hpf,lpf);    % if fdesign not available...
        y.filters.coeffs.b=b;
        y.filters.coeffs.a=a;
    end
    
else
    if not(nofilters)
        y.filters=filter;
    end
end

% set the images structure to some reasonable defauls...
y.images.last_frame = 0;




