%%
postProcessingOnly = false; %true;

pd = pdefs();
basePath   = pd.Mea1kShared;

dataFiles = {fullfile(basePath, '140821','data_analysed', 'Roland', 'preprocessed_binary', '0001.h5') }
destFolder = fullfile(basePath, '140821', 'data_analysed', 'Roland', 'testrun');
if ismac
    dataFiles = {'/Users/rolandd/tmp/0001.h5'};
    destFolder = '/Users/rolandd/tmp/destinationfolder';
end

jobName = 'test_grid_sorting'
nG = 8

%%
sj = grid.SortJob(jobName, dataFiles, destFolder, 'nGroups', nG, 'postProcOnly', postProcessingOnly);
sj.createBOTMGroups();
sj.copyDataFiles();
sj.setTaskParameters();
sj.prepareTasks();

%%
if ismac
    if ~postProcessingOnly
        for i = 1:nG
            grid.SortJob.run(['/Users/rolandd/tmp/' jobName '/taskFiles/taskFile' num2str(i) '.mat']);
        end
    end
end

%%
all_tasks_completed = sj.waitForTasksToFinish();

if all_tasks_completed

    sj.copyBackResults();
    
    tic
    sj.postprocessBOTMSorting();
    toc
    
    sj.runQC()
end
