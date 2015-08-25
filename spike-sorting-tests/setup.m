% rmpath('C:\SVN\Matlab\');
% rmpath('C:\SVN\Matlab\mysortpackage\');
% rmpath('C:\SVN\Matlab\dataviewerpackage');

this_path = pwd;
fprintf('Folder located in %s\n', this_path);
cmosmea_external_matlab = fullfile(this_path, '..', '..', '..', '..', 'cmosmea_external', 'matlab', 'trunk');
cmosmea_external_meabench = fullfile(this_path, '..', '..', '..', '..', 'cmosmea_external', 'meabench', 'trunk', 'matlab');
isLinux = isempty(strfind(computer, 'WIN'));
% offlineFile = fullfile(this_path, 'DataAnalysis', 'offline_status_indicator.txt');
% try
%     delete(offlineFile)
% end
   
if isLinux
    nlxPath = 'nlx2mat-unix-master';
else
    nlxPath = 'NeuralynxMatlabImportExport_v501';
end
paths = {%this_path
         cmosmea_external_matlab
         fullfile(cmosmea_external_matlab, 'SpikeSorter')
         fullfile(cmosmea_external_matlab, 'misc')
%          fullfile(cmosmea_external_matlab, 'tools', 'm2html')
         genpath(cmosmea_external_meabench)
         
         
         fullfile(this_path, 'BOTM')
         fullfile(this_path, 'DataAnalysis')
         fullfile(this_path, 'DataViewer')
         fullfile(this_path, 'HDMEAgui')
         fullfile(this_path, 'GridComputing')
         fullfile(this_path, 'GuiUtils')
         fullfile(this_path, 'Code')
         fullfile(this_path, 'NeuroRouter')
         fullfile(this_path, 'Circular')
         fullfile(this_path, 'Sandbox', 'Vrushali')
         fullfile(this_path, '..', 'setup', 'software', 'matlab')
         fullfile(this_path, 'External')         
             fullfile(this_path, 'External', 'UIToolbox')
             fullfile(this_path, 'External', 'UIToolbox', 'Patch')
             fullfile(this_path, 'External', 'UIToolbox', 'layoutHelp')
           %fullfile(this_path, 'External', 'discrim')  % LDA package
             fullfile(this_path, 'External', 'arnoOnken')
             fullfile(this_path, 'External', 'profile_history')
             fullfile(this_path, 'External', 'FastICA')
%              fullfile(this_path, 'External', 'SLMtools')
             fullfile(this_path, 'External', 'MeanShift')
             fullfile(this_path, 'External', 'Optometrika')
             fullfile(this_path, 'External', 'boundedline_20140701')
             fullfile(this_path, 'External', nlxPath)   
             fullfile(this_path, 'External', 'matlab_bgl-4.0.1', 'matlab_bgl')
             fullfile(this_path, 'External', 'CircStat2011f')   
             
         fullfile(this_path, 'External', 'uiinspect')
         fullfile(this_path, 'Misc')
         fullfile(this_path, 'Mysort')
         fullfile(this_path, 'FileUtils')};
t1 = tic;
for k=1:size(paths)
    fprintf('Adding to path: %s\n', paths{k});
    addpath(paths{k});
end
mysort.addUMS2000ToPath(this_path);
t2 = toc(t1);
fprintf('This took %.2f sec.\n', t2);

fprintf('Checking setup... ');
try
    mysort.wf.m2v([]);
    fprintf('success!\n');
catch
    fprintf('failed!\n');
end

clear t1 t2
clear k cmosmea_external this_path paths ans
[colvec C marker] = mysort.plot.vectorColor(1);
set(0,'DefaultAxesColorOrder',C);

