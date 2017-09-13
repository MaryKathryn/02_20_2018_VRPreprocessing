function ProcessNeuralData_HighLevel(Options,varargin)

% This file pre-processes recorded neural data recorded in HcTask and saves
% the resultant LFP data as a *.mat file, and/or the spike data as an HDF5
% file. It was made by Roberto Gulli and Guillaume Doucet, built and
% adapted from Neuroshare's Example.m file. This function is meant for
% the import and analysis of neural file information, spike times, spike
% waveforms, LFP data, and events from a sorted Plexon or Cerebus neural
% file. This is intended to be the globally used import function for neural
% data collected in the Martinez-Trujillo lab in or after 2014.
%
% This was built on the architecture of the RG_nsImportSpike file used in
% the analysis of Sergio's optogenetics data (November 2013). The Sync
% Pulse portion of this script was written originally by GD, and
% modified and incorporated into this script by RG.
%
% Note, HDF5 files can have specific subfields accessed directly without
% loading the entire file in to RAM, saving a great deal of time when this
% data needs to be accessed for subsequent analyses. LFP import was done
% prior to the lab switching to HDF5 as a preferred data format.
%
%
% USAGE NOTES:
% Include all Neuroshare scripts in their own subfolder. This can be
% located ~/Analyses/Low-level functions/. The Neuroshare subfolder folder
% should be named 'Neuroshare_library'. Ensure this folder also contains
% the DLLs provided by the recording system manufacturer. Note, NeuroShare
% functions on PCs ONLY.
%
%
% CRITICAL RECORDING INSTRUCTIONS AND NOTES:
% 1) Ensure that the Cerebus neural recording file is continuously
% recorded, not incremented with pauses occurring between trials. If these
% pauses are included in the recording, the continuous data (LFP) is not
% included in the .nev file, and this function must be altered to read the
% continuous data from the .ns2  ile instead.
%
% To ensure Continuous recording is enabled, so that Monkeylab does not
% send words to stop/start the neural recording in the intertrial interval.
% Settings > DAQ > Click Continuous Recording
%
% 2) Ensure the SyncPulse is on in your behavioural task on the MonkeyLab
% computer:
% 
% In Monkeylab > Settings > DAQ > A/D O:
%       Event Name: SyncPulseRead
%       Field Name: SyncPulse
%       Sample Time: 1
%       Check A/DO
%
% 3) Ensure that the analog input ainp9 in the Central settings on the
% Cerebus computer or Channel 20 on the plexon system are
% being sampled at 1kHz
% 
% 4) In the hardware configuration, the digital inputs need to be set to
% read 16 bits on word strobe, otherwise the events (words/trial markers)
% sent from Monkeylab will not be read properly.
% Digital input > Function > 16-bit on Word Strobe
%
% 5) If you are using the lever, ensure the lever is enabled for your
% subject. Settings > DAQ > A/D 4:
%       Event Name: LeverStateUpdate > Field
%       Name: Lever
%       Sample Time: 1
%
% 6) The .nev file needs to be ***_sorted.nev or ***.nev file for this
% script to automatically detect them.


% DATA RETRIEVAL 
%%% To see the architecture of your HDF5 file
% h5disp(fileName)
% h5info(fileName)

%%% To read the data just written in to the file
% cell2mat(h5read(fileName,'/Info/SessionInfo/Date/'))
% cell2mat(h5read(fileName,'/Info/SessionInfo/MonkeyName/'))
% cell2mat(h5read(fileName,'/Info/SessionInfo/FileName/'))
% cell2mat(h5read(fileName,'/Info/RecordingInfo/FileComment/'))   

%%% To retrieve all units that belong to a channel
% fileName = 'D:\HcTask\Results\Woody\20141217\W_20141217_SpikeData.hdf5';
% chs      = cell2mat(h5read(fileName,'/Info/UnitInfo/Channels/'));
% unnums   = h5read(fileName,'/Info/UnitInfo/UnitNum/');
% unitIDs  = [chs repmat('_',length(unnums),1) num2str(unnums)];
% unitIDs  = cellstr(unitIDs);
% unitIDs  = regexprep(unitIDs,' ','')

%%% Rasters
% h5readatt(fileName,'/Data/TrialRasters','PadTime')
% Rasters = h5read(fileName,'/Data/TrialRasters');
% Note: Remove padding of 2's before moving on with this data

%%% Timestamps
% tsInfo = h5info(fileName,'/Data/TimeStamps');
% hdfUnitNames = {tsInfo.Datasets.Name}';
% for k = 1:length(hdfUnitNames)
%     TimeStamps{k} = h5read(fileName,['/Data/TimeStamps/' hdfUnitNames{k}]);
% end

%%% Waveforms
% wfInfo = h5info(fileName,'/Data/WaveForms');
% hdfUnitNames = {wfInfo.Datasets.Name}';
% for k = 1:length(hdfUnitNames)
%     WaveForms{k} = h5read(fileName,['/Data/WaveForms/' hdfUnitNames{k}]);
% end
% hdfUnitNames = regexprep({wfInfo.Datasets.Name}',' ','');


%% EDITS
% July 21, 2016
% - Removed 'ResultsFolder' input to the function (was unused anyway)
% - Re-created VGS task rasters so that they span the ITC-18 start time to
%   EOT time. The SOT is when the visual stimulus was presented. The ITC-18
%   start time is the start of the inter-trial interval.


%% Set up

% % DEFAULT FUNCTION OPTIONS
% 1 to output the selected data, 0 to skip
if isempty(Options)
    Options.Spikes = 1;
    Options.LFPs = 1;
    Options.Trials = 1;
else
    if ~isfield(Options, 'Spikes')
        Options.Spikes = 1;
    end
    
    if ~isfield(Options, 'LFPs')
        Options.LFPs = 1;
    end
    
    if ~isfield(Options, 'Trials')
        Options.Trials = 1;
    end
    
end


% Defines Slash for separating folders
sl = filesep; % File folder seperator (slash) different on Mac and PC
dirAnalyses = mfilename('fullpath');
pathParts = strsplit(dirAnalyses,sl);
dirRoot = find(~cellfun(@isempty,(strfind(pathParts,':'))));
projectName = pathParts(dirRoot+1);
dirAnalyses = [strjoin(pathParts(1:end-1),sl) sl];
dirTask = dirAnalyses(1:max(strfind(dirAnalyses,[sl 'Analyses'])));
dirData = [dirTask 'Data' sl];

% Parse inputs
p = inputParser;
%default parameters
defaultMonkeyName = 'Theo';
defaultSessions = {'all'};
defaultBSave = false;
% Parse inputs
addParameter(p,'monkeyName',defaultMonkeyName,@ischar);
addParameter(p,'sessions',defaultSessions,@iscell);
addParameter(p,'bSave',defaultBSave,@islogical);
parse(p,varargin{:});
% Define parsed variables
monkeyName  = p.Results.monkeyName;
sessions    = p.Results.sessions;
if isrow(sessions)
    sessions = sessions';
end
if p.Results.bSave
    bSave = 'true';
else
    bSave = 'false';
end




% % LOAD NEUROSHARE DLL
% Gets current file directory (removes the filename from the mfilename function
dirAnalyses = mfilename('fullpath');
dirAnalyses = dirAnalyses(1:max(strfind(dirAnalyses,sl)));
% Checks that a location was identified, and prompts user for a location if necessary
if isempty(dirAnalyses)
    dirAnalyses = uigetdir(cd, 'Select directory containing Analysis scripts');
end
% Ensures that there is a slash at the end of the path
if ~strcmp(dirAnalyses(end), sl)
    dirAnalyses = [dirAnalyses sl];
end
%Checks for Neuroshare_library subfolder folder and generates the path
[~,dirNSLibrary] = system('dir /S/B Neuroshare_library');
dirNSLibrary = [strtrim(dirNSLibrary) sl];
if isempty(dirNSLibrary)
    dirNSLibrary = uigetdir('C;\', 'Select NeuroShare Library Folder');
    dirNSLibrary = [dirNSLibrary sl];
end
%Adds the NeuroShare directory path to search paths
addpath(dirNSLibrary);
% Determine whether the computer is 32 or 64 bits, to load the correct DLL
CompInfo = computer;
if strcmp(CompInfo(end-1:end), '64')
    LibraryFile = 'nsNEVLibrary64.dll';
else
    LibraryFile = 'nsNEVLibrary.dll';
end
% Load the corresponding DLL
% Note: DLLs are produced by the vendors (e.g., Plexon, Blackrock) and
% allow Neuroshare to read many proprietary filetypes.
[nsresult] = ns_SetLibrary([dirNSLibrary LibraryFile]);
if (nsresult ~= 0)  % Check that the DLL was loaded properly
    disp('DLL was not found!');
    return
end



%% CHOOSE WHICH SESSIONS TO RUN
dirResults = [dirTask 'Results' sl monkeyName sl ];
% Get sessions
if strcmp('all',sessions)
    dirDataFolders = dir(...
        [dirResults '2*']);
    sessions = {dirDataFolders.name}';
end
SelectedDirectories = cellstr([ repmat(dirResults,length(sessions),1) cell2mat(sessions) ]);

%% LOOP THROUGH SELECTED SESSIONS
for sessionNum = 1: numel(SelectedDirectories)
    
    %% INITIALIZE STRUCTURES & SET UP
    display(['Processing neural data for ' SelectedDirectories{sessionNum}]);
    sessionTimer = tic;
    SpikeData = struct();
    LFPData = struct();
    NeuralTrials = struct();
    
    
    dirResultsSess = [SelectedDirectories{sessionNum} sl];
    dirDataSess    = strrep(dirResultsSess,'Results','Data');
    
    dirResultsNeur  = [dirResultsSess 'NeuralData' sl];
    dirResultsBehav = [dirResultsSess 'BehavData' sl];
    
    
    %%% FIND NEV FILE
    NevFile = dir([dirDataSess 'NSP1_001.nev']); % RG 20170503: Changed to process unsorted NSP1_001
    if isempty(NevFile)
        cd(dirAnalyses);
        warning('No nev file was found, aborting.');
        SpikeData = struct();
        LFPData = struct();
        NeuralTrials = struct();
        return
    end
    % If multiple nev files are contained in folder the user needs to specify
    % which one to use
    if size(NevFile,1)>1
        [FileIndex, Ok] = listdlg('liststring', {NevFile.name}, 'selectionmode', 'single', 'name', 'Select correct Nev file', 'listsize', [300 300]);
        
        if ~Ok
            FileIndex=1;
        end
        
        NevFile = NevFile(FileIndex);
    end
    File = [dirDataSess NevFile.name];
    
    
    
    %%% FIND NS6 FILE
    FileFound = dir([dirDataSess 'NSP1_001.ns6']); % RG 20170503: Changed to process unsorted NSP1_001
    if size(FileFound,1) ~=1
        [LFPFilename, LFPPathname] = uigetfile({'*.ns6'},...
            'Select analog recording file',dirDataSess);
        LFPFile = [LFPPathname LFPFilename];
    else
        LFPFile = [dirDataSess FileFound.name];
    end
    if isempty(LFPFile)
        warning('No ns6 file was found, aborting.');
        SpikeData = struct();
        LFPData = struct();
        NeuralTrials = struct();
        return
    end
    
    
    
    %% LOAD BEHAVIOURAL DATA
    
    display(['Loading behavioural data for ' SelectedDirectories{sessionNum}]);
    behavData = tic;
    behavFiles = cellstr(ls([dirResultsBehav '*.mat']));
    
    allTasksStruct = struct();
    for task = 1:length(behavFiles)
        taskStruct = struct();
        load([dirResultsSess 'BehavData' sl behavFiles{task}])
        if task ==1
            allTasksStruct.SessionInfo = taskStruct.SessionInfo;
            allTasksStruct.Trials = taskStruct.Trials;
        else
            allTasksStruct.Trials = catstruct(allTasksStruct.Trials,taskStruct.Trials);
        end
    end
    [allTasksStruct.Trials,perm] = orderfields(allTasksStruct.Trials);
    toc(behavData)
    
    %% %% PROCESS SYNC PULSE AND COMPUTE LAG TIMES
    if Options.Trials
        
        display(['Processing SyncPulse for ' ...
            SelectedDirectories{sessionNum}]);
        syncPulseTimer = tic;
        [NeuralTrials] = HcTask_ProcessSyncPulse(...
            allTasksStruct,LFPFile);
        display(['Done processing SyncPulse for ' SelectedDirectories{sessionNum}]);
        toc(syncPulseTimer)
        
    end % Options.Trials
    
    
    
    %% PROCESS UNIT DATA
    % Spike data is contained in the '.nev' file for the Cerebus system, the
    % '.plx' file from the Plexon system. At this point the data should be
    % manually sorted using Offline Sorter, and saved as a '.nev' file.
    
    if Options.Spikes
        
        
        % Ensure the workspace contains the NeuralTrials structure
        if ~Options.Trials
            if ~isempty(dir([dirResultsNeur '*NeuralTrials*mat']))
                NeuralTrialsFile = dir([dirResultsNeur '*NeuralTrials*.mat']);
                load([dirResultsNeur NeuralTrialsFile.name]);
                
            else
                [NeuralTrials] = HcTask_ProcessSyncPulse(...
                                    allTasksStruct,LFPFile);
            end
        end
        
        display(['Processing unit data for ' SelectedDirectories{sessionNum}]);
        unitTimer = tic;
        HcTask_ProcessUnitData(allTasksStruct,NeuralTrials,File,dirResultsNeur)
        display(['Done processing unit data for ' SelectedDirectories{sessionNum}]);
        toc(unitTimer)
    end
    
    
    
    %% PROCESS LFP DATA
    
    if Options.LFPs
        % Ensure the workspace contains the NeuralTrials structure
        if ~Options.Trials
            if ~isempty(dir([dirResults '*NeuralTrials*mat']))
                NeuralTrialsFile = dir([dirResults '*NeuralTrials*.mat']);
                load([dirResults NeuralTrialsFile.name]);
            else
                [NeuralTrials] = HcTask_ProcessSyncPulse(EyeCalStruct,XMazeStruct,FreeRoamStruct,EyeEndStruct,File,dirResults);
            end
        end
        
        display(['Processing LFP data for ' SelectedDirectories{sessionNum}]);
        lfpTimer = tic;
        HcTask_ProcessLFPData(XMazeStruct,LFPFile,dirResults);
        display(['Done processing LFP data for ' SelectedDirectories{sessionNum}]);
        toc(lfpTimer)
    end
    
    
    
    
    
    
    clear XMazeStruct EyeCalStruct EyeEndStruct FreeRoamStruct
    
    
    
    %% End session loop
    display(['Done processing neural data for '  SelectedDirectories{sessionNum}]);
    toc(sessionTimer)
    display(char(10));
    
    
end % SESSION LOOP


end % End of function



