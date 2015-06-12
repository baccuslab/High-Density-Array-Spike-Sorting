function P = copyH5toBinary(inFile, outFile, binFile, varargin)

P.chunkSize = 250000;
P.sessionName = '/Sessions/Session0';
P = mysort.util.parseInputs(P, varargin);

assert(exist(inFile, 'file')>0, ['H5 File does not exist! (' inFile ')']);
[pathstr, name, ext] = fileparts(inFile);
if isempty(outFile)
    outFile = fullfile(pathstr, [name '.h5']);
end
assert(~exist(outFile, 'file'), ['Output File does already exist! ' outFile]);

data = mysort.h5.matrix(inFile, [P.sessionName '/sig'], true);
[nSamples, nC] = size(data);

proc = mysort.h5.createVariableAndOrFile(outFile, '/bFileIsInProcess', [1 1], [1 1], 'H5T_NATIVE_INT');
proc(1,1) = int32(1);

pref = mysort.h5.createVariableAndOrFile(outFile, [P.sessionName '/filter/prefiltered'], [1 1], [1 1], 'H5T_NATIVE_INT');
pref(1,1) = h5read(inFile, [P.sessionName '/filter/prefiltered']);

clear pref
high = mysort.h5.createVariableAndOrFile(outFile, [P.sessionName '/filter/highpass'], [1 1], [1 1], 'H5T_NATIVE_INT');
high(1,1) = h5read(inFile, [P.sessionName '/filter/highpass']);
clear high
low = mysort.h5.createVariableAndOrFile(outFile, [P.sessionName '/filter/lowpass'], [1 1], [1 1], 'H5T_NATIVE_INT');
low(1,1) = h5read(inFile, [P.sessionName '/filter/lowpass']);
clear low
down = mysort.h5.createVariableAndOrFile(outFile, [P.sessionName '/filter/downsamplefactor'], [1 1], [1 1], 'H5T_NATIVE_INT');
down(1,1) = h5read(inFile, [P.sessionName '/filter/downsamplefactor']);
clear down
ord = mysort.h5.createVariableAndOrFile(outFile, [P.sessionName '/filter/order'], [1 1], [1 1], 'H5T_NATIVE_INT');
ord(1,1) = h5read(inFile, [P.sessionName '/filter/order']);
clear type
filterType = h5read(inFile, [P.sessionName '/filter/type']); filterType = cell2mat(filterType);
l = length(filterType);
ord = mysort.h5.createVariableAndOrFile(outFile, [P.sessionName '/filter/type'], [1 l], [1 l], 'H5T_C_S1');
ord(1,1:l) = filterType;
clear ord l filterType
gd = mysort.h5.createVariableAndOrFile(outFile, [P.sessionName '/filter/gainmultiplier'], [1 1], [1 1], 'H5T_NATIVE_INT');
gd(1,1) = h5read(inFile, [P.sessionName '/filter/gainmultiplier']);
clear gd

% CHIP ID
chipid = mysort.h5.createVariableAndOrFile(outFile, [P.sessionName '/chipid'], [1 1], [1 1], 'H5T_NATIVE_INT');
chipid(1,1) = h5read(inFile, [P.sessionName '/chipid']);
clear chipid

% GAIN
gain = mysort.h5.createVariableAndOrFile(outFile, [P.sessionName '/gain'], [1 4], [1 4], 'H5T_NATIVE_DOUBLE');
gain(1,:) = h5read(inFile, [P.sessionName '/gain']);
clear gain

% ADC range and resolution
gain = mysort.h5.createVariableAndOrFile(outFile, [P.sessionName '/adc_resolution'], [1 1], [1 1], 'H5T_NATIVE_DOUBLE');
gain(1,:) = h5read(inFile, [P.sessionName '/adc_resolution']);
clear gain
gain = mysort.h5.createVariableAndOrFile(outFile, [P.sessionName '/adc_range'], [1 1], [1 1], 'H5T_NATIVE_DOUBLE');
gain(1,:) = h5read(inFile, [P.sessionName '/adc_range']);
clear gain;

% SR
sr_ = mysort.h5.createVariableAndOrFile(outFile, [P.sessionName '/sr'], [1 1], [1 1], 'H5T_NATIVE_INT');
sr_(1,1) = h5read(inFile, [P.sessionName '/sr']);
clear sr_

% VERSION
version = mysort.h5.createVariableAndOrFile(outFile, [P.sessionName '/version'], [1 1], [1 1], 'H5T_NATIVE_INT');
version(1,1) = h5read(inFile, [P.sessionName '/version']);
clear version

% CHANNEL LIST
channel_list =  h5read(inFile, [P.sessionName '/channel_list']);
names = fieldnames(channel_list);
type_id = H5T.create('H5T_COMPOUND',length(names)*32);
for i=1:length(names)
    H5T.insert(type_id, names{i}, (i-1)*32, 'H5T_NATIVE_INT');
end
cl = mysort.h5.createVariableAndOrFile(outFile, [P.sessionName '/channel_list'], nC, nC, type_id);
% get one object
a = cl(1);
isConnected = zeros(nC,1);
count = 0;
for i=1:nC
    for j = 1:length(names)
        a.(names{j}) = channel_list.(names{j})(i);
    end
    cl(i) = a;
end
clear a cl

x = mysort.h5.createVariableAndOrFile(outFile, [P.sessionName '/channel_nr'], [1 nC], [1 nC], 'H5T_NATIVE_DOUBLE');
x(1,:) = h5read(inFile, [P.sessionName '/channel_nr']);
clear x
x = mysort.h5.createVariableAndOrFile(outFile, [P.sessionName '/channel_posx'], [1 nC], [1 nC], 'H5T_NATIVE_DOUBLE');
x(1,:) = h5read(inFile, [P.sessionName '/channel_posx']);
clear x
x = mysort.h5.createVariableAndOrFile(outFile, [P.sessionName '/channel_posy'], [1 nC], [1 nC], 'H5T_NATIVE_DOUBLE');
x(1,:) = h5read(inFile, [P.sessionName '/channel_posy']);
clear x
x = mysort.h5.createVariableAndOrFile(outFile, [P.sessionName '/channel_connected'], [1 nC], [1 nC], 'H5T_NATIVE_DOUBLE');
x(1,:) = h5read(inFile, [P.sessionName '/channel_connected']);
isConnected = x;
clear x

%dims = [nSamples nC];
%maxDims = [nSamples nC];

% create a binary file where the data is stored
sig = mysort.ds.binaryFileMatrix(binFile, [1 nC], 'writable', true);

% Save a link to the binary file into /sig:
[pathstr,name,ext] = fileparts(binFile);
binFileName = [name, '.dat'];
sig_link = mysort.h5.createVariableAndOrFile(outFile, [P.sessionName '/sig'], [1 length(binFileName) ], [1 length(binFileName) ], 'H5T_C_S1');
sig_link(1,1:length(binFileName)) = binFileName;
clear sig_link;

bin_dims = mysort.h5.createVariableAndOrFile(outFile, [P.sessionName '/bin_dims'], [1 2], [1 2], 'H5T_NATIVE_LONG');
bin_dims(1, :) = [nSamples nC];
clear bin_dims;

fprintf('Writing Data...');

% data is in double, however, we want it in uint:
%dfun = @(x) uint16(x);

timer = mysort.util.ProcessTimer(ceil(nSamples/P.chunkSize));
i = 1;
while i < nSamples
    timer.next();
    timer.showProgress();
    
    j = min(nSamples, i+P.chunkSize-1);
    
    %sig(i:j, :) = data(i:j, :);
    %sig(i:j, :) = data.sessionList.X.sessionList.getRawData(i:j, :);
    sig(i:j, :) = data.getData(i:j, :);
    
    %assert( ~any(any( sig(i:j, :) ~= data.getData(i:j, :) )), 'error!') 
    
    i = i + P.chunkSize;
end
disp('Done copying.')
clear sig

% Set File as being done
proc(1,1) = int32(0);
clear proc