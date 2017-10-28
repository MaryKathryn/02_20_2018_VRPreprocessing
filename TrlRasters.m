function TrlRasters(spikeFile,MLStruct,varargin)

% This function is built to parse the neural data in to ML trials, which
% can then be used to run subsequent time- or event-locked analyses.
%
% This script will plot rasters for all trials contained in
% MLStruct.(session).Trials. You can choose to alter the
% plotted trials within this script by calling
% [trlIDs2plot, ~] = SelectTrials(MLStruct,Options) with options containing
% whatever trial criteria you would like met. See SelectTrials help for
% more info.
%
% Inputs
%   spikeFile   - String. Full filepath of your spike kdf5 file
%
%   MLStruct    - Structure. Behavioural trial structure you would like to
%                 plot. Because trial epochs are hard-coded in the function
%                 right now, this version of this structure only works for
%                 Unreal 3 associative memory tasks (XMaze, ALFixed,
%                 ALNovel)
%                 Note: Use SelectTrials() before running this to only
%                 include trials you want to plot in MLStruct.
%
%   neurons     - Cell. Contains all unitNames you would like to include
%                 in the plot.
%                 Default: {'all'} includes all neurons found in the hdf5
%                 file
% 
%   mMultiunit  - Boolean. Combine all units on a single channel into
%                 multi-unit activity
% 
%   bPlot       - Boolean. Plot trial rasters for all trials
% 
%   bSave       - Boolean. Save plotted trial rasters
%                 
% 
%%% Usage example %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
% dirResultsSess = 'E:\ALVR\Results\Theo\20170413\';
% dirResultsNeur = [dirResultsSess 'NeuralData' filesep];
% dirResultsBehav = regexprep(dirResultsNeur,'Neural','Behav');
% 
% load([dirResultsBehav ls([dirResultsBehav '*ALNovel_Block1*'])])
% spikeFile = [dirResultsNeur ls([dirResultsNeur '*hdf5'])]
% 
% % Keep only valid trials
% Options.defaultToZero           = 0;  % Sets all other options to 0
% Options.RemoveIncompleteTrials  = 1;  % Removes invalid trials (user stopped, timed out, crash)
% Options.RemoveShortTrials       = 1; 
% Options.TrialConditions         = 0;  % XMaze only trial summary table
% Options.SelectTrials            = 1;  % Select trials based on specific trial criteria
% Options.TrialCriteria           = {'OutcomeWord' {'Correct' 'Incorrect'}}; 
% Options.Save                    = 0;
% [selectedTrials,~] = SelectTrials(taskStruct, Options);
% % Keep only selectedTrials in taskStruct.Trials
% trlIDs = fieldnames(taskStruct.Trials);
% trls2remove = trlIDs(~ismember(trlIDs,selectedTrials.Trials));
% taskStruct.Trials = rmfield(taskStruct.Trials,trls2remove);
% 
% TrlRasters(spikeFile,taskStruct,'bMultiUnit',true,'bPlot',true,'bSave',true)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parse inputs
p = inputParser;
%default parameters
defaultNeurons = {'all'};
defaultBMultiUnit = false;
defaultBPlot = true;
defaultBSave = false;
% Parse inputs
addParameter(p,'neurons',defaultNeurons,@iscell);
addParameter(p,'bMultiUnit',defaultBMultiUnit,@islogical);
addParameter(p,'bPlot',defaultBPlot,@islogical);
addParameter(p,'bSave',defaultBSave,@islogical);
parse(p,varargin{:});
% Define parsed variables
neurons     = p.Results.neurons;
bMultiUnit  = p.Results.bMultiUnit;
bPlot       = p.Results.bPlot;
bSave       = p.Results.bSave;

%% Set up
s = filesep;

% Name the session
session =  datestr(MLStruct.SessionInfo.Date,'yyyymmdd');

%% Retrieve timestamps for all units


%%% Read session and recording info
display(cell2mat(h5read(spikeFile,'/Info/SessionInfo/Date/')));
display(cell2mat(h5read(spikeFile,'/Info/SessionInfo/MonkeyName/')));
display(cell2mat(h5read(spikeFile,'/Info/SessionInfo/FileName/')));
display(cell2mat(h5read(spikeFile,'/Info/RecordingInfo/FileComment/')));

%%% To retrieve all units that belong to a channel
chs = cell2mat(h5read(spikeFile,'/Info/UnitInfo/Channels/'));
chs = cellstr(chs);
chs = cellfun(@strtrim,chs,'uni',0);
unitNum = h5read(spikeFile,'/Info/UnitInfo/UnitNum/');

%%% Timestamps
tsInfo = h5info(spikeFile,'/Data/TimeStamps');
if strcmp(neurons,'all')
    unitNames = {tsInfo.Datasets.Name}';
else % Need to change this to ensure neurons are the same format as those used in the hdf5 file
    unitNames = neurons;
end
% Retrieve timestamps
TimeStamps = cell(length(unitNames),1);
for k = 1:length(unitNames)
    TimeStamps{k} = h5read(spikeFile,['/Data/TimeStamps/' unitNames{k}]);
end
if bMultiUnit
    uniqueChs = unique(chs);
    multiUnits = cell(length(uniqueChs),1);
    for ch = 1:length(uniqueChs)
        
        channelInds = strcmp(uniqueChs{ch},chs);
        multiUnits{ch} = sort(vertcat(TimeStamps{channelInds}));
        
    end
    TimeStamps = multiUnits;
    chs = uniqueChs;
    unitNames = cellfun(@(x)([x '_multi']),chs,'uni',0);
end

%%% Trial IDs
mlTrlIDs = fieldnames(MLStruct.Trials);
cerTrlIDs = cellstr(cell2mat(...
    h5read(spikeFile,'/Info/TrialInfo/TrialIDs/') ) );
cerSOTs = h5read(spikeFile,'/Info/TrialInfo/TrialStartTs/');
% cerEOTs = h5read(spikeFile,'/Info/TrialInfo/TrialEndTs/'); % Not using this.
% % Instead use cerSOTs + (mlEOTs-mlSOTs)

%%% Remove trials from cerebus trlIDs and SOTs that are not inlcuded in
%%% MLStruct
cerTrls2remove = ~ismember(cerTrlIDs,mlTrlIDs);
cerTrlIDs(cerTrls2remove) = [];
cerSOTs(cerTrls2remove) = [];
if ~isequal(mlTrlIDs,cerTrlIDs)
    mlTrls2remove = mlTrlIDs(~ismember(mlTrlIDs,cerTrlIDs));
    MLStruct.Trials = rmfield(MLStruct.Trials,mlTrls2remove);
    mlTrlIDs = fieldnames(MLStruct.Trials);
end


%% Create raster data and plot multi-neuron rasters for each trial of interest

for k=1:length(mlTrlIDs) % Trial counter
    
    
    % Temporarily store trial behavioural data
    temp_trl = MLStruct.Trials.(mlTrlIDs{k});
    trlLength = temp_trl.EOT_Time-temp_trl.SOT_Time;
    ctxOnInd = floor((temp_trl.ContextOnset-temp_trl.SOT_Time)*1000);
    goalsOnInd = floor((temp_trl.GoalsOnset-temp_trl.SOT_Time)*1000);
    decisInd = floor((temp_trl.DecisionOnset-temp_trl.SOT_Time)*1000);
    tempRaster = false(length(TimeStamps),floor((trlLength)*1000));
    
    
    % For each channel
    for unit=1:length(TimeStamps)
        
        % Retrieve spike times relative to the start of the trial
        spikeTimes = TimeStamps{unit}(...
            TimeStamps{unit} > cerSOTs(k) & ...
            TimeStamps{unit} <(cerSOTs(k)+trlLength)) - cerSOTs(k);
        spikeInds = floor(spikeTimes*1000); % Convert to raster with 1ms resolution
        spikeInds(spikeInds==0) = 1; % Ensure there is no spike at raster index 0
        
        tempRaster(unit,spikeInds) = true;
        
    end
    
    if bPlot
        close all
        hRaster = figure('Units','Normalized','Position',[0 0 1 1]);
        hold on;
        % Using imagesc is not suitable for creating plots you intend to work
        % with later in illustrator, as it creates a bitmapped image rather
        % than a vector graphic. However, it is suitable and fast for data
        % vizualization.
        hRax = imagesc(tempRaster);
        colormap([1 1 1; 0 0 0])
        plot([ctxOnInd ctxOnInd],ylim,'-b')
        plot([goalsOnInd goalsOnInd],ylim,'-r')
        plot([decisInd decisInd],ylim,'-g')
        ylim([0.5 length(TimeStamps)+.5])
        xlim([0 size(tempRaster,2)])
        xlabel('Time (ms)')
        if bMultiUnit
            ylabel('Channel number (multi-unit)')
        else
            ylabel('Neuron number')
        end
        title(regexprep(...
            {strtok(MLStruct.SessionInfo.FileName,'.') ; ...
            mlTrlIDs{k} },...
            '_','-') )
        legend({'Context Onset' 'Goals Onset' 'Decision'},...
            'Location','NorthOutside','Orientation','horizontal',...
            'box','off')
    end
    
    if bSave
        
        fileNameParts = strsplit(MLStruct.SessionInfo.FileName,'_');
        taskName = strjoin(fileNameParts(3:end),'_');
        taskName = strtok(taskName,'.');
        
        dirSave = [strtok(spikeFile,'N') 'TrialRasters' ...
            filesep taskName filesep];
        if ~isdir(dirSave)
            mkdir(dirSave)
        end
        saveStr = ['Raster_' mlTrlIDs{k} '.png'];
        
        saveas(hRaster,[dirSave saveStr])
    end
    
end




end % Channel counter






