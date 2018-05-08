%% if statement to call parcluster or regular analysis function
function [outputs] = CycIF_execute(jobs, parallel, imageDirectory, saveDirectory,...
    experiment, mag, maxCycle, FOVlimits, channelNames,...
    bugGFP, bugmCherry, punctaChannels)

if parallel == true
    c=parcluster;
    c.AdditionalProperties.WallTime = '6:00:00';
    c.AdditionalProperties.QueueName = 'short';
    c.AdditionalProperties.AdditionalSubmitArgs = '--mem-per-cpu=2G'
    batchJob = createJob(c);
end

    
for j = 1:size(jobs,1)
    timepoint = jobs{j,1};
    row = jobs{j,2};
    column = jobs{j,3};
    if parallel == false
        [morphOut, fluorescenceOut] = CycIF_head(imageDirectory, saveDirectory,...
            experiment, timepoint, row, column, mag, maxCycle, FOVlimits, channelNames,...
            bugGFP, bugmCherry, punctaChannels);
    else
        createTask(batchJob, @CycIF_head, 2, {imageDirectory, saveDirectory,...
            experiment, timepoint, row, column, mag, maxCycle, FOVlimits, channelNames,...
            bugGFP, bugmCherry, punctaChannels})
    end
end
if parallel == true
    submit(batchJob);
    wait(batchJob);
%     for i = 1:length(batchJob.Tasks)
%         task = batchJob(i);
%         if strcmp(task.State, 'error')
%             fprintf(1, 'Task %d: %s', i, task.ErrorMessage);
%         end
%     end
    %batchJob.delete to delete log files
end

    
