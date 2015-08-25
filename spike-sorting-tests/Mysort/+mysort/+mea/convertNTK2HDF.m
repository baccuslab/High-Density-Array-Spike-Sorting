function P = convertNTK2HDF(ffile, varargin)
P.outFile = [];
P.prefilter = 0;
P.hpf = 300;
P.lpf = 7000;
P.filterOrder = 4;
P.fir_filterOrder = 110;
P.filtfilt = true;
P.deflation = 1;
P.chunkLen = 1000;
P.downSample = 1;
P.chunkIndividualChannels = 0;
P.sessionName = '/Sessions/Session0';
P.save_as_binary = true;
P = mysort.util.parseInputs(P, varargin);

if iscell(ffile)
    assert(isempty(P.outFile)|| (iscell(P.outFile)&&length(P.outFile)==length(ffile)),'if multiple ntk files are specified, outFile must be empty or contain one file per ntk file.');
    P_ = P;
    for i=1:length(ffile)
        if iscell(P.outFile)
            P_.outFile = P.outFile{i};
        end
        mysort.mea.convertNTK2HDF(ffile{i},P_);
    end
    return
end
assert(exist(ffile, 'file')>0, ['NTK File does not exist! (' ffile ')']);
[pathstr, name, ext] = fileparts(ffile);
if isempty(P.outFile)
    P.outFile = fullfile(pathstr, [name '.h5']);
end
%     delete(P.outFile);
assert(~exist(P.outFile, 'file'), ['Output File does already exist! ' P.outFile]);
assert(int32(P.downSample)==P.downSample, 'downSample must be an integer!');
% NOT NECESSARY ANYMORE! assert(P.prefilter || P.downSample==1, 'You must prefilter before downsampling to avoid aliasing!');
assert(P.downSample==1 || P.lpf >0, 'You must high pass filter before dowsampling, otherwise you will get artefacts at the endpoints!');

siz = 10000;
ntk=initialize_ntkstruct(ffile, 'nofilters');
[high_density_data ntk] = ntk_load(ntk, siz);
nC = size(ntk.data, 1);

sr = ntk.sr;
gainmultiplier = 1;
filterType = 'None      ';

if P.downSample > 1
    sr = ntk.sr/P.downSample;
    if P.prefilter
        assert(P.lpf <= .5*sr, 'The lowpass filter must be below the Nyquist Frequency (AFTER downsampling, if enabled)!');
    end
end

% Prefiltered
if P.prefilter
    gainmultiplier = 256;
    h5Type = 'H5T_NATIVE_SHORT';
    dfun = @(x) int16(gainmultiplier*x);
    
    % create filter object
    %         hd = mysort.mea.filter_design(P.hpf, P.lpf, ntk.sr, P.filterOrder);
    %         filterType = 'IIR butter';
    b  = mysort.mea.filter_design_fir(P.hpf, P.lpf, sr, P.fir_filterOrder);
    filterType = 'FIR fircls';
    d = P.downSample;
    h = P.hpf;
    l = P.lpf;
    fo = P.filterOrder;
elseif P.downSample > 1
    % do we want to have a different gain after only downsampling?
    h = 0;
    l = 0;
    fo = 0;
    d = P.downSample;
    gainmultiplier = 256;
    h5Type = 'H5T_NATIVE_SHORT';
    dfun = @(x) int16(gainmultiplier*x);
else
    h5Type = 'H5T_NATIVE_USHORT';
    dfun = @(x) uint16(x);
    gainmultiplier = 1;
    h = 0;
    l = 0;
    fo = 0;
    d = 1;
end
% Set File as being in process
proc = mysort.h5.createVariableAndOrFile(P.outFile, '/bFileIsInProcess', [1 1], [1 1], 'H5T_NATIVE_INT');
proc(1,1) = int32(1);

pref = mysort.h5.createVariableAndOrFile(P.outFile, '/Sessions/Session0/filter/prefiltered', [1 1], [1 1], 'H5T_NATIVE_INT');
pref(1,1) = int32(P.prefilter);
clear pref
high = mysort.h5.createVariableAndOrFile(P.outFile, '/Sessions/Session0/filter/highpass', [1 1], [1 1], 'H5T_NATIVE_INT');
high(1,1) = int32(h);
clear high
low = mysort.h5.createVariableAndOrFile(P.outFile, '/Sessions/Session0/filter/lowpass', [1 1], [1 1], 'H5T_NATIVE_INT');
low(1,1) = int32(l);
clear low
down = mysort.h5.createVariableAndOrFile(P.outFile, '/Sessions/Session0/filter/downsamplefactor', [1 1], [1 1], 'H5T_NATIVE_INT');
down(1,1) = int32(d);
clear down
ord = mysort.h5.createVariableAndOrFile(P.outFile, '/Sessions/Session0/filter/order', [1 1], [1 1], 'H5T_NATIVE_INT');
ord(1,1) = int32(fo);
clear type
ord = mysort.h5.createVariableAndOrFile(P.outFile, '/Sessions/Session0/filter/type', [1 20], [1 20], 'H5T_C_S1');
ord(1,1:length(filterType)) = filterType;
clear ord
gd = mysort.h5.createVariableAndOrFile(P.outFile, '/Sessions/Session0/filter/gainmultiplier', [1 1], [1 1], 'H5T_NATIVE_INT');
gd(1,1) = int32(gainmultiplier);
clear gd



% CHIP ID
chipid = mysort.h5.createVariableAndOrFile(P.outFile, '/Sessions/Session0/chipid', [1 1], [1 1], 'H5T_NATIVE_INT');
chipid(1,1) = int32(ntk.chipid);
clear chipid

% GAIN
gain = mysort.h5.createVariableAndOrFile(P.outFile, '/Sessions/Session0/gain', [1 4], [1 4], 'H5T_NATIVE_DOUBLE');
gain(1,2:4) = [ntk.gain1 ntk.gain2 ntk.gain3];
gain(1,1) = prod(gain(1,2:4)); % total gain
clear gain

% ADC range and resolution
gain = mysort.h5.createVariableAndOrFile(P.outFile, '/Sessions/Session0/adc_resolution', [1 1], [1 1], 'H5T_NATIVE_DOUBLE');
gain(1,1) = ntk.adc_resolution;
gain = mysort.h5.createVariableAndOrFile(P.outFile, '/Sessions/Session0/adc_range', [1 1], [1 1], 'H5T_NATIVE_DOUBLE');
gain(1,1) = ntk.adc_range;
clear gain;

% SR
sr_ = mysort.h5.createVariableAndOrFile(P.outFile, '/Sessions/Session0/sr', [1 1], [1 1], 'H5T_NATIVE_INT');
sr_(1,1) = int32(sr);
clear sr_

% VERSION
version = mysort.h5.createVariableAndOrFile(P.outFile, '/Sessions/Session0/version', [1 1], [1 1], 'H5T_NATIVE_INT');
version(1,1) = int32(ntk.version);
clear version

% CHANNEL LIST
names = {'channel_nr', 'connected', 'x', 'y', 'idx', 'dummy', 'damaged'};
type_id = H5T.create('H5T_COMPOUND',length(names)*32);
for i=1:length(names)
    H5T.insert(type_id, names{i}, (i-1)*32, 'H5T_NATIVE_INT');
end
cl = mysort.h5.createVariableAndOrFile(P.outFile, '/Sessions/Session0/channel_list', [size(ntk.channel_list,2)], [size(ntk.channel_list,2)], type_id);
% get one object
a = cl(1);
isConnected = zeros(size(ntk.channel_list,2),1);
elx = [];
ely = [];
for i=1:size(ntk.channel_list,2)
    % set values into object
    a.channel_nr = int32(ntk.channel_list(i));
    if ~isempty(ntk.channels{i,1}.els)
        a.connected = int32(1);
        isConnected(i) = 1;
        a.x = int32(ntk.channels{i,1}.els{1,1}.x*1000);
        a.y = int32(ntk.channels{i,1}.els{1,1}.y*1000);
        a.idx = int32(ntk.channels{i,1}.els{1,1}.idx);
        a.dummy = int32(ntk.channels{i,1}.els{1,1}.dummy);
        a.damaged = int32(ntk.channels{i,1}.els{1,1}.damaged);
    else
        a.connected = int32(0);
        a.x = int32(0);
        a.y = int32(0);
        a.idx = int32(0);
        a.dummy = int32(0);
        a.damaged = int32(0);
    end
    % set correct index in file to values of object
    cl(i)=a;
    elx(1,i) = a.x;
    ely(1,i) = a.y;
end
clear a cl
isConnected = isConnected == 1;
nC_effective = length(ntk.channel_list);
% assert(sum(isConnected) == nC_effective, 'Must be identical');

% NEW CL FORMAT
x = mysort.h5.createVariableAndOrFile(P.outFile, '/Sessions/Session0/channel_nr', [1 nC_effective], [1 nC_effective], 'H5T_NATIVE_DOUBLE');
x(1,1:nC_effective) = ntk.channel_list;
clear x
x = mysort.h5.createVariableAndOrFile(P.outFile, '/Sessions/Session0/channel_posx', [1 nC_effective], [1 nC_effective], 'H5T_NATIVE_DOUBLE');
x(1,1:nC_effective) = elx;
clear x
x = mysort.h5.createVariableAndOrFile(P.outFile, '/Sessions/Session0/channel_posy', [1 nC_effective], [1 nC_effective], 'H5T_NATIVE_DOUBLE');
x(1,1:nC_effective) = ely;
clear x
x = mysort.h5.createVariableAndOrFile(P.outFile, '/Sessions/Session0/channel_connected', [1 nC_effective], [1 nC_effective], 'H5T_NATIVE_DOUBLE');
x(1,1:nC_effective) = double(isConnected(:)');
clear x

% SIG
X = ntk.data';
if P.downSample>1
    X = double(X);
    X = bsxfun(@minus, X, mean(X,1));
    X = resample(X, 1, P.downSample);
end
if P.prefilter
    X = double(X);
    %         dummy = filtfilthd(hd, X);
    %         clear dummy
    %         X = filtfilthd(hd, X); % do this twice for burn-in!
    X(:,isConnected) = bsxfun(@minus, X(:,isConnected), mean(X(:,isConnected),1));
    X(:,isConnected) = conv2(X(:,isConnected), b(:), 'same');
end
X(:,isConnected) = dfun(X(:,isConnected));

if P.chunkIndividualChannels
    chunkDims = [P.chunkLen 1];
else
    chunkDims = [P.chunkLen nC];
end

% reading 1mio samples of
% a) one channel
% b) 102 channels
% with
% 1) chunkDims (over channels) = 1 or
% 2) chunkDims = 102
% takes:
% 1) a) 0.1 sec
%    b) 3.2 sec
% 2) a) 1.1 sec
%    b) 1.2 sec
% for 1) the files are around 30% bigger

dims = [size(X,1) nC];
maxDims = [-1 nC];
if ~P.save_as_binary
    sig = mysort.h5.createVariableAndOrFile(P.outFile, [P.sessionName '/sig'], ...
        dims, maxDims, h5Type, chunkDims, P.deflation);
else    
    [pathstr,name,ext] = fileparts(P.outFile);
    P.binFile = fullfile(pathstr, [name, '.dat']);
    if exist(P.binFile, 'file')
        warning('Binary File exists already: %s', P.binFile);
%         s = input('Overwrite? (y/n)', 's');
%         if strcmp(s, 'y')
            delete(P.binFile);
%         else
%             error('Binary File exists already, converting aborted');
%         end
    end
    % create a binary file where the data is stored
    %sig = mysort.ds.binaryFileMatrix(P.binFile, [1 dims(2)], 'writable', true);
    sig = mysort.ds.binaryFileMatrix(P.binFile, [1 nC_effective], 'writable', true);
     
    % Save a link to the binary file into /sig:
    binFileName = [name, '.dat'];
    sig_link = mysort.h5.createVariableAndOrFile(P.outFile, [P.sessionName '/sig'], [1 length(binFileName) ], [1 length(binFileName) ], 'H5T_C_S1');
    sig_link(1,1:length(binFileName)) = binFileName;
    clear sig_link;
end
sig(1:size(X,1),isConnected) = X;

chunkSize = 500000;
nextIdx = size(sig,1)+1;
%     fprintf('Writing Data...');
while ~ntk.eof
    [high_density_data ntk] = ntk_load(ntk, chunkSize);
    X = ntk.data';
    if P.downSample > 1
        X = double(X);
        X = bsxfun(@minus, X, mean(X,1));
        X = resample(X, 1, P.downSample);
    end
    if P.prefilter
        %             X = filtfilthd(hd, double(X)); % do this twice for burn-in!
        X = double(X);
        %         dummy = filtfilthd(hd, X);
        %         clear dummy
        %         X = filtfilthd(hd, X); % do this twice for burn-in!
        X(:,isConnected) = bsxfun(@minus, X(:,isConnected), mean(X(:,isConnected),1));
        X(:,isConnected) = conv2(X(:,isConnected), b(:), 'same');
    end
    X(:,isConnected) = dfun(X(:,isConnected));
    sig(nextIdx:nextIdx+size(X,1)-1,isConnected) = X;
    nextIdx = nextIdx+size(X,1);
end
clear sig
bin_dims = mysort.h5.createVariableAndOrFile(P.outFile, [P.sessionName '/bin_dims'], [1 2], [1 2], 'H5T_NATIVE_LONG');
bin_dims(1, :) = [nextIdx-1 nC_effective];
clear bin_dims
% Set File as being done
proc(1,1) = int32(0);
clear proc