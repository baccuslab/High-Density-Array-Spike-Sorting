%%
fname = 'Trace_id843_2011-12-06T11_27_53_5.stream.h5';
var_str = '/Sessions/Session0/sig';
plist = 'H5P_DEFAULT';
rmode = 'H5F_ACC_RDONLY';
channelSet = [1:2:92];
nC = length(channelSet);
T = 100000;

%% Read from mea object
fprintf('\nUse mea object\n');
tic
mea = mysort.mea.CMOSMEA(fname);
X0 = mea(1:T, channelSet);
toc
%% Read individual rows

% Init hyperslab with first selected channel (row) and add (OR)
% the other rows. Use only one H5D.read operation! (much faster
% than reading every row individually)
mode = 'H5S_SELECT_SET';
fprintf('\nRead individual rows\n');
tic 
fid = H5F.open(fname,rmode,plist); 
dset_id = H5D.open(fid, var_str);
dims =  [1 T];  
h5_dims = dims; 
mem_space_id = H5S.create_simple(2, h5_dims, h5_dims);
file_space_id = H5D.get_space(dset_id);
block = [1 T];
h5_block = block; 
X1 = uint8(zeros(fliplr(dims)));
for c = 1:nC
    myC = channelSet(c);
    offset = [myC-1 0];      
    h5_offset = offset; 
    H5S.select_hyperslab(file_space_id,mode,h5_offset,[],[1 1],h5_block);
    X1(:,c) = H5D.read(dset_id, 'H5ML_DEFAULT', mem_space_id, file_space_id, plist);
end
H5D.close(dset_id);
H5F.close(fid);
toc

%% Build hyperslab from individual rows
block = [1 T];
h5_block = block; 
% Init hyperslab with first selected channel (row) and add (OR)
% the other rows. Use only one H5D.read operation! (much faster
% than reading every row individually)
mode = 'H5S_SELECT_SET';
fprintf('\nBuild hyperslab from individual rows\n');
tic
fid = H5F.open(fname,rmode,plist); 
dset_id = H5D.open(fid, var_str);
dims =  [nC T];  
h5_dims = dims; 
mem_space_id = H5S.create_simple(2, h5_dims, h5_dims);
file_space_id = H5D.get_space(dset_id);
for c = 1:nC
    myC = channelSet(c);
    offset = [myC-1 0];      
    h5_offset = offset; 
    H5S.select_hyperslab(file_space_id,mode,h5_offset,[],[1 1],h5_block);
    if c == 1
        mode = 'H5S_SELECT_OR';
    end
end
X2 = H5D.read(dset_id, 'H5ML_DEFAULT', mem_space_id, file_space_id, plist);
H5D.close(dset_id);
H5F.close(fid);
toc

%% Build big hyperslab over all rows and remove those that are too much
c1 = min(channelSet);
c2 = max(channelSet);
channelSet_ = c1:c2;
nC_ = c2-c1+1;
block = [nC_ T];
h5_block = block; 
% Init hyperslab with first selected channel (row) and add (OR)
% the other rows. Use only one H5D.read operation! (much faster
% than reading every row individually)
mode = 'H5S_SELECT_SET';
fprintf('\nRead one big hyperslab\n');
tic
fid = H5F.open(fname,rmode,plist); 
dset_id = H5D.open(fid, var_str);
dims =  [nC_ T];  
h5_dims = dims; 
mem_space_id = H5S.create_simple(2, h5_dims, h5_dims);
file_space_id = H5D.get_space(dset_id);
offset = [c1-1 0];      
h5_offset = offset; 
H5S.select_hyperslab(file_space_id,mode,h5_offset,[],[1 1],h5_block);
X3 = H5D.read(dset_id, 'H5ML_DEFAULT', mem_space_id, file_space_id, plist);
H5D.close(dset_id);
H5F.close(fid);
toc

%% test function
fprintf('\nRead with mea.h5read (one big hyperslab)\n');
tic
X4 = mysort.mea.h5read(fname, var_str, channelSet, 1:T);
toc

%%
sum(sum(abs(X4-X1)))
sum(sum(abs(X4-X2)))

%%
fprintf('\nRead whole data with matlab highend function hdf5read\n');
tic
X = hdf5read(fstr, '/Sessions/Session0/sig');
toc                 

%%
fprintf('\nRead half the data with mea.h5read\n');
tic
X = mysort.mea.h5read(fname, '/Sessions/Session0/sig', 1:131, 1:round(mea.size(1)/2));
toc                 
