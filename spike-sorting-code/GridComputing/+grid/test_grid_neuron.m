%%
pd = pdefs();
basePath   = pd.Mea1kShared;

sortingFolder = fullfile(  basePath, '140821', 'data_analysed', 'Roland', 'botm_sorted', 'files_0to0');
templateFile = fullfile( sortingFolder, 'hamsterDS_binary_results.mat');
if ismac
    sortingFolder = '/Users/rolandd/tmp/destinationfolder';
    templateFile = fullfile( sortingFolder, 'mytest_binary_results.mat');
end

jobName = 'mytest_neuronvisualisation'

nN = 4;
sj = grid.NeuronJob(jobName, templateFile, 'nNeurons', nN );
sj.loadTemplateFile();
sj.setParameters('createGraphs', 'test')

sj.copyDataFilesToScratch();

%sj.setResume(true);
sj.prepareTasks();

%%
if ismac
    %parfor i = 1:nN
    for i = 1:sj.nTasks
        grid.NeuronJob.run(['/Users/rolandd/tmp/' jobName '/taskFiles/taskFile' num2str(i) '.mat'])
    end
end

%%
all_tasks_completed = sj.waitForTasksToFinish();

if all_tasks_completed
    
    if ismac
        destFolder = '/Users/rolandd/tmp/destinationfolder';
    else
        destFolder = fullfile(  basePath, '140821', 'data_analysed', 'Roland', 'botm_sorted', 'analyzed_neurons_test')
    end
    
    sj.copyBackResults(destFolder);
end

