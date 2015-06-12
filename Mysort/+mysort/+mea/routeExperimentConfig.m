function routeExperimentConfig(expPath, ELs4Michele)
currentPath = pwd;
% create experimental folders
matlabFolder = fullfile(expPath, 'Matlab');
bMatFolderExisted = exist(matlabFolder, 'file');
if ~bMatFolderExisted
    mkdir(matlabFolder);
end
mkdir(fullfile(expPath, 'Configs'));
% load electrode positions, sorted from best electrode to worst elecotrdes

% THIS IS YOUR OUTPUT
M = ELs4Michele;
% THIS IS YOUR OUTPUT


req_electrodes=4; % MAX NUMBER OF ELECTRODES YOU WANT TO ROUTE
spec_electrodes=10; % LIST LENGTH FROM WHICH YOU WANT TO GET ELETRODES

% STUPID LOOP FROM APST CODE, IGNORE IT
n_rout = cell(1,size(M,1));
for i = 1:size(M,1)

    n_rout{i}.el_idx=M(i,:,3);
    n_rout{i}.x=M(i,:,1);
    n_rout{i}.y=M(i,:,2);
    
end


% make text file

npos_multiloc={};
figure('color','w');hold on

for n=1:length(n_rout)
    
    plot(n_rout{n}.x,n_rout{n}.y,'go');
   
    for ii=1:length(n_rout{n}.el_idx)
        npos_multiloc{end+1}.label=sprintf('neuron%d', n);
        npos_multiloc{end}.x=n_rout{n}.x(ii);
        npos_multiloc{end}.y=n_rout{n}.y(ii);
        npos_multiloc{end}.cost=ii;
        npos_multiloc{end}.elcnt=req_electrodes;
        npos_multiloc{end}.multiloc=1;
    end
    
  
    
end
axis equal
axis ij
 
% check that you have the function below

cd(fullfile(expPath, 'Matlab'))
hidens_write_neuroplacement('multiloc.neuropos.nrk', 'npos', npos_multiloc, 'size', 20);


% execute NeuroDishRouter

fnames={'multiloc'};
nr_exe='NeuroDishRouter';
for fn=1:length(fnames)
    fprintf('running: %s -n -v 2 -l %s -s %s\n', nr_exe, [pwd '/../Configs/matlab_specs/' fnames{fn} '.neuropos.nrk'], [pwd '/../Configs/' fnames{fn}]);
    unix(sprintf('%s -n -v 2 -l %s -s %s\n', nr_exe, [pwd '/../Configs/matlab_specs/' fnames{fn} '.neuropos.nrk'], [pwd '/../Configs/' fnames{fn}]));
end


% visualize routed els

fn=1;
fname=['../Configs/' fnames{fn} '.el2fi.nrk2'];
fid=fopen(fname);
elidx=[];
tline = fgetl(fid);
while ischar(tline)
    [tokens] = regexp(tline, 'el\((\d+)\)', 'tokens');
    elidx(end+1)=str2double(tokens{1});
    tline = fgetl(fid);
end
fclose(fid);

els=hidens_get_all_electrodes(2);

figure('color','w'); % plot selected vs routed
hold on
for n=1:length(n_rout)
   l{1}=plot(n_rout{n}.x,n_rout{n}.y,'go');
end
l{2}=plot(els.x(elidx+1), els.y(elidx+1), 'rx');
axis equal
axis ij

res_avgrank=nan(1,length(fnames));
res_total_els=nan(1,length(fnames));
res_ng2=nan(1,length(fnames));
res_ng3=nan(1,length(fnames));

selected_els=zeros(1, length(els.x));
selected_els(elidx)=1;
    
nr_connected=zeros(1,length(n_rout));
avg_rank=zeros(1,length(n_rout));
for n=1:length(n_rout)
    con=selected_els(n_rout{n}.el_idx);
    nr_connected(n)=sum(con);
    if sum(con)==0
        avg_rank(n)=spec_electrodes;
    else
        avg_rank(n)=sum((1:spec_electrodes).*con)/sum(con);
    end
end

res_avgrank(fn)=mean(avg_rank);
res_total_els(fn)=length(elidx);
res_ng2(fn)=sum(nr_connected>=2);
res_ng3(fn)=sum(nr_connected>=3);

figure('color','w')
subplot(2,1,1)
hist(nr_connected, 0:spec_electrodes);
title([fnames{fn} ' Routing: Connected Electrodes per Neuron'])
ylabel('Cnt')
xlabel('# of Electrodes')

subplot(2,1,2)
hist(avg_rank, 0:spec_electrodes);
title([fnames{fn} ' Routing: Avg. Rank per Neuron'])
ylabel('Cnt')
xlabel('Avg Rank')

fprintf('Routing     Total-Els    N with >=2 els    N with >=3 els   Avg. Rank\n');
fprintf('---------------------------------------------------------------------\n');
for fn=1:length(fnames)
    fprintf('%-17s %3d               %3d               %3d         %1.1f\n', fnames{fn}, res_total_els(fn), res_ng2(fn), res_ng3(fn), res_avgrank(fn));
end

if ~bMatFolderExisted
    rmdir(matlabFolder);
end

cd(currentPath);