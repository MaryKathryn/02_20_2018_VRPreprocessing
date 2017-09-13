% Example function
% 
%%% Goal: Plot rasters for every trial for one neuron

%% Load all necessary data
s = filesep ;

dirResultsSess = 'E:\ALVR\Results\Theo\20170413\';
dirResultsNeur = [dirResultsSess 'NeuralData' s];
dirResultsBehav = regexprep(dirResultsNeur,'Neural','Behav');

load([dirResultsBehav ls([dirResultsBehav '*ALNovel_Block1*'])])
% load([dirResultsNeur ls([dirResultsNeur '*NeuralTrials*'])])
spikeFile = [dirResultsNeur ls([dirResultsNeur '*hdf5'])]


%% Ensure only valid & full trials are kept
% Removes 

Options.defaultToZero           = 0;  % Sets all other options to 0
Options.RemoveIncompleteTrials  = 1;  % Removes invalid trials (user stopped, timed out, crash)
Options.RemoveShortTrials       = 1; 
Options.TrialConditions         = 0;  % XMaze only trial summary table
Options.SelectTrials            = 0;  % Select trials based on specific trial criteria
Options.TrialCriteria           = {}; 
Options.Save                    = 0;
[SelectTrials,~] = SelectTrials(taskStruct, Options);

trlIDs = fieldnames(taskStruct.Trials);
trls2remove = trlIDs(~ismember(trlIDs,SelectTrials.Trials));
taskStruct.Trials = rmfield(taskStruct.Trials,trls2remove);


%%


