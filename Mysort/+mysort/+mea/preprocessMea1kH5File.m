function P = preprocessMea1kH5File(ffile, map, varargin)
P.outFile = [];
P.prefilter = 0;
P.hpf = 300;
P.lpf = 7000;
P.filterOrder = 4;
P.fir_filterOrder = 110;
P.filtfilt = true;
P.bUseFPGA_IIR_Filter = false;
P.deflation = 1;
P.chunkLen = [];
P.downSample = 1;
P.chunkIndividualChannels = 0;
P.subtractMeanOverAllChannels = true;
P.restrictToTimePeriod = [];
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
        mysort.mea.preprocessMea1kH5File(ffile{i},P_);
    end
    return
end

%%
% Im H5 file ist spalte 0 channel 0.
% 1023 ist der letzte Kanal dann DAC und !command counter"
% DAC ist nicht im nrk file 
% 

%% FROM JELENA FOR THE FPGA PREP FILTER
% X_fixP_10b = round(X_unfilt*2)/2;
% Xin = X_fixP_10b';
%  
% Flp=3000;
% Fhp=500;
%  
% f_bp = mysort.mea.filter_design(Fhp,Flp,20000);
% mysort.mea.filter_init(f_bp, Xin(1:nS,:));
% X_BP=filter(f_bp,Xin)';
%%

assert(exist(ffile, 'file')>0, ['H5 File does not exist! (' ffile ')']);
[pathstr, name, ext] = fileparts(ffile);
if isempty(P.outFile)
    P.outFile = fullfile(pathstr, [name '.h5']);
end
%     delete(P.outFile);
assert(~exist(P.outFile, 'file'), ['Output File does already exist! ' P.outFile]);
assert(int32(P.downSample)==P.downSample, 'downSample must be an integer!');
% NOT NECESSARY ANYMORE! assert(P.prefilter || P.downSample==1, 'You must prefilter before downsampling to avoid aliasing!');
assert(P.downSample==1 || P.lpf >0, 'You must high pass filter before dowsampling, otherwise you will get artefacts at the endpoints!');

ntk = struct();
ntk.data = mysort.h5.matrix(ffile, '/sig', true);
ntk.chipid = -1;
ntk.gain1 = 1;
ntk.gain2 = 1;
ntk.gain3 = 1;
ntk.adc_resolution = 13.8;
ntk.adc_range = 1;
[nSamples, nC_] = size(ntk.data);

if ~isempty(P.restrictToTimePeriod)
    error('Not Implemented At the moment');
    assert(length(P.restrictToTimePeriod) == 2, 'P.restrictToTimePeriod must be a time period in samples [start end]!');
    assert(~any(round(P.restrictToTimePeriod) ~= P.restrictToTimePeriod), 'P.restrictToTimePeriod must be in intergers!');
    assert(P.restrictToTimePeriod(1) > 0 && P.restrictToTimePeriod(2) < nSamples, 'P.restrictToTimePeriod out of bounds!');
    assert(P.restrictToTimePeriod(1) < P.restrictToTimePeriod(2), 'P.restrictToTimePeriod must be increasing');
end

nC = map.noChans;
assert(length(nC) == 1, 'map.noChans must be a number!')
assert(nC>0, 'map.noChans must be greater 0!')

% Make sure the map vectors are all in the right format
map.elNo  = map.elNo(:)';
map.chans = map.chans(:)';
map.mposx = map.mposx(:)';
map.mposy = map.mposy(:)';


ntk.sr = 20000;
sr = ntk.sr;
ntk.version = -1;
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
    d = P.downSample;
    if P.bUseFPGA_IIR_Filter
        % INIT IIR Filter
        filterType = 'IIR butterworth FPGA'; 
        h = 500;
        l = 3000; 
        fo = 2;
        f_bp = mysort.mea.filter_design_fpga(h, l, sr, fo);      
    else
        % INIT FIR Filter
        b  = mysort.mea.filter_design_fir(P.hpf, P.lpf, sr, P.fir_filterOrder);
        filterType = 'FIR fircls';
        
        h = P.hpf;
        l = P.lpf;
        fo = P.fir_filterOrder;
    end
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

pref = mysort.h5.createVariableAndOrFile(P.outFile, [P.sessionName '/filter/prefiltered'], [1 1], [1 1], 'H5T_NATIVE_INT');
pref(1,1) = int32(P.prefilter);
clear pref
high = mysort.h5.createVariableAndOrFile(P.outFile, [P.sessionName '/filter/highpass'], [1 1], [1 1], 'H5T_NATIVE_INT');
high(1,1) = int32(h);
clear high
low = mysort.h5.createVariableAndOrFile(P.outFile, [P.sessionName '/filter/lowpass'], [1 1], [1 1], 'H5T_NATIVE_INT');
low(1,1) = int32(l);
clear low
down = mysort.h5.createVariableAndOrFile(P.outFile, [P.sessionName '/filter/downsamplefactor'], [1 1], [1 1], 'H5T_NATIVE_INT');
down(1,1) = int32(d);
clear down
ord = mysort.h5.createVariableAndOrFile(P.outFile, [P.sessionName '/filter/order'], [1 1], [1 1], 'H5T_NATIVE_INT');
ord(1,1) = int32(fo);
clear type
ord = mysort.h5.createVariableAndOrFile(P.outFile, [P.sessionName '/filter/type'], [1 20], [1 20], 'H5T_C_S1');
ord(1,1:length(filterType)) = filterType;
clear ord
gd = mysort.h5.createVariableAndOrFile(P.outFile, [P.sessionName '/filter/gainmultiplier'], [1 1], [1 1], 'H5T_NATIVE_INT');
gd(1,1) = int32(gainmultiplier);
clear gd



% CHIP ID
chipid = mysort.h5.createVariableAndOrFile(P.outFile, [P.sessionName '/chipid'], [1 1], [1 1], 'H5T_NATIVE_INT');
chipid(1,1) = int32(ntk.chipid);
clear chipid

% GAIN
gain = mysort.h5.createVariableAndOrFile(P.outFile, [P.sessionName '/gain'], [1 4], [1 4], 'H5T_NATIVE_DOUBLE');
gain(1,2:4) = [ntk.gain1 ntk.gain2 ntk.gain3];
gain(1,1) = prod(gain(1,2:4)); % total gain
clear gain

% ADC range and resolution
gain = mysort.h5.createVariableAndOrFile(P.outFile, [P.sessionName '/adc_resolution'], [1 1], [1 1], 'H5T_NATIVE_DOUBLE');
gain(1,1) = ntk.adc_resolution;
gain = mysort.h5.createVariableAndOrFile(P.outFile, [P.sessionName '/adc_range'], [1 1], [1 1], 'H5T_NATIVE_DOUBLE');
gain(1,1) = ntk.adc_range;
clear gain;

% SR
sr_ = mysort.h5.createVariableAndOrFile(P.outFile, [P.sessionName '/sr'], [1 1], [1 1], 'H5T_NATIVE_INT');
sr_(1,1) = int32(sr);
clear sr_

% VERSION
version = mysort.h5.createVariableAndOrFile(P.outFile, [P.sessionName '/version'], [1 1], [1 1], 'H5T_NATIVE_INT');
version(1,1) = int32(ntk.version);
clear version

% CHANNEL LIST
names = {'channel_nr', 'connected', 'x', 'y', 'idx', 'dummy', 'damaged'};
type_id = H5T.create('H5T_COMPOUND',length(names)*32);
for i=1:length(names)
    H5T.insert(type_id, names{i}, (i-1)*32, 'H5T_NATIVE_INT');
end
cl = mysort.h5.createVariableAndOrFile(P.outFile, [P.sessionName '/channel_list'], nC, nC, type_id);
% get one object
a = cl(1);
isConnected = zeros(nC,1);
count = 0;
for i=1:nC
    % set values into object
    a.channel_nr = int32(map.elNo(i));
    if a.channel_nr > 1
        a.connected = int32(1);
        isConnected(i) = 1;
        a.x = int32(map.mposx(i)*1000)*17.5;
        a.y = int32(map.mposy(i)*1000)*17.5;
        a.idx = int32(map.chans(i));
        a.dummy = int32(0);
        a.damaged = int32(0);
        
        % set correct index in file to values of object
        count = count+1;
        cl(count)=a;
        
    else
        a.connected = int32(0);
        a.x = int32(0);
        a.y = int32(0);
        a.idx = int32(0);
        a.dummy = int32(0);
        a.damaged = int32(0);
    end

end
clear a cl
nC_effective = sum(isConnected);
assert(nC_effective == count, 'Must be identical');
isConnected = isConnected == 1;

% New CL Format
idx = isConnected==1;
idx = idx(:);
x = mysort.h5.createVariableAndOrFile(P.outFile, [P.sessionName '/channel_nr'], [1 nC_effective], [1 nC_effective], 'H5T_NATIVE_DOUBLE');
x(1,1:nC_effective) = map.elNo(idx);
clear x
x = mysort.h5.createVariableAndOrFile(P.outFile, [P.sessionName '/channel_posx'], [1 nC_effective], [1 nC_effective], 'H5T_NATIVE_DOUBLE');
x(1,1:nC_effective) = map.mposx(idx)*17.5;
clear x
x = mysort.h5.createVariableAndOrFile(P.outFile, [P.sessionName '/channel_posy'], [1 nC_effective], [1 nC_effective], 'H5T_NATIVE_DOUBLE');
x(1,1:nC_effective) = map.mposy(idx)*17.5;
clear x
x = mysort.h5.createVariableAndOrFile(P.outFile, [P.sessionName '/channel_connected'], [1 nC_effective], [1 nC_effective], 'H5T_NATIVE_DOUBLE');
x(1,1:nC_effective) = ones(1,nC_effective);
clear x


if isempty(P.chunkLen) || (~isempty(P.chunkLen) && P.chunkLen == 0)
    % This disables chunking and deflation
    chunkDims = [];
else
    if P.chunkIndividualChannels
        chunkDims = [P.chunkLen 1];
    else
        chunkDims = [P.chunkLen nC_effective];
    end
end

isConnectedIdx = find(isConnected);

if isempty(P.restrictToTimePeriod)
    lastSample = 0;
else
    nSamples = P.restrictToTimePeriod(2) - P.restrictToTimePeriod(1) +1;
    lastSample = P.restrictToTimePeriod(1)-1;
end

dims = [nSamples nC_effective];
maxDims = [nSamples nC_effective];

if ~P.save_as_binary
    sig = mysort.h5.createVariableAndOrFile(P.outFile, [P.sessionName '/sig'], ...
        dims, maxDims, h5Type, chunkDims, P.deflation);
else
    
    [pathstr,name,ext] = fileparts(P.outFile);
    P.binFile = fullfile(pathstr, [name, '.dat']);
    
    % create a binary file where the data is stored
    %sig = mysort.ds.binaryFileMatrix(P.binFile, [1 dims(2)], 'writable', true);
    sig = mysort.ds.binaryFileMatrix(P.binFile, [1 nC_effective], 'writable', true);
     
    % Save a link to the binary file into /sig:
    binFileName = [name, '.dat'];
    sig_link = mysort.h5.createVariableAndOrFile(P.outFile, [P.sessionName '/sig'], [1 length(binFileName) ], [1 length(binFileName) ], 'H5T_C_S1');
    sig_link(1,1:length(binFileName)) = binFileName;
    clear sig_link;

    bin_dims = mysort.h5.createVariableAndOrFile(P.outFile, [P.sessionName '/bin_dims'], [1 2], [1 2], 'H5T_NATIVE_LONG');
    bin_dims(1, :) = [nSamples nC_effective];
    clear bin_dims;
end
    
if P.prefilter && P.bUseFPGA_IIR_Filter
    % If we use an IIR filter, we need to burn it in
    X = ntk.data(1:min(end, 5000),isConnectedIdx);
    X = double(X);
    X = bsxfun(@minus, X, mean(X,1));
    if P.subtractMeanOverAllChannels
        X = X-repmat(mean(X,2), 1, size(X,2));
    end    
    mysort.mea.filter_init(f_bp, X);
end

chunkSize = 250000;
nextIdx = 1;
%     fprintf('Writing Data...');

timer = mysort.util.ProcessTimer(ceil(nSamples/chunkSize));
while lastSample < nSamples
    timer.next();
    timer.showProgress();

    thisLastSample = min(nSamples, lastSample+chunkSize);
    X = ntk.data(lastSample+1:thisLastSample, isConnectedIdx);
    
    if P.downSample > 1
        X = double(X);
        X = bsxfun(@minus, X, mean(X,1));
        if P.subtractMeanOverAllChannels
            X = X-repmat(mean(X,2), 1, size(X,2));
        end
        X = resample(X, 1, P.downSample);
    end
    if P.prefilter
        %             X = filtfilthd(hd, double(X)); % do this twice for burn-in!
        X = double(X);
        %         dummy = filtfilthd(hd, X);
        %         clear dummy
        %         X = filtfilthd(hd, X); % do this twice for burn-in!
        X = bsxfun(@minus, X, mean(X,1));
        if P.subtractMeanOverAllChannels
            X = X-repmat(mean(X,2), 1, size(X,2));
        end
        if P.bUseFPGA_IIR_Filter
            X = filter(f_bp, X);
        else
            X = conv2(X, b(:), 'same');
        end
    end
    X = dfun(X);
    sig(nextIdx:nextIdx+size(X,1)-1,:) = X;
    nextIdx = nextIdx+size(X,1);
    lastSample = thisLastSample;  
end

% i = 1;
% while i < nSamples
%     timer.next();
%     timer.showProgress();
%     
%     j = min(nSamples, i+P.chunkSize-1);
%     
%     %sig(i:j, :) = data(i:j, :);
%     %sig(i:j, :) = data.sessionList.X.sessionList.getRawData(i:j, :);
%     sig(i:j, :) = data.getData(i:j, :);
%     
%     %assert( ~any(any( sig(i:j, :) ~= data.getData(i:j, :) )), 'error!') 
%     
%     i = i + P.chunkSize;
% end
% disp('Done copying.')



disp('Done preprocessing.')
clear sig

% Set File as being done
proc(1,1) = int32(0);
clear proc