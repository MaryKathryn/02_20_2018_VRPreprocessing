function ProcessBehavData_HighLevel(varargin)

% Pre-processes behavioural data recorded in any task and saves the
% resultant *.mat file.
% 
% This function relies on a *unique* low-level function for each task type
% you have included in your project. Include 
%
% 
% =========================================================================
% INPUTS
% =========================================================================
% 
% taskNames:       A cell containing strings for each task you would like
%                  to analyze
%       'EyeCal' - Pre-processes data from the visually guided
%                  saccade task (aka calibration task, aka 
%                  cued saccade task)
%                  Written by RG/BC/GD
% 
%       'EyeEnd' - Pre-processes data from the visually guided
%                  saccade task (aka calibration task, aka 
%                  cued saccade task) done at the end of the
%                  session
%                  Written by RG/BC/GD
% 
%        'XMaze' - (Old; may not work) 
%                  Pre-processes data from Rob's Associative
%                  memory task used in the hippocampal recordings
%                  Written by RG/GD
% 
%     'FreeRoam' - (Old; may not work) 
%                  Pre-processes data from Rob's Foraging task
%                  Written by RG/GD
% 
%           'AL' - Pre-processes data from Rob's Associative learning task.
%                  Includes the fixed (Orange/Steel-Purple/Wood)
%                  associatiation that was done at the beginning and end of
%                  every recording session, and the novel association(s)
%                  that were also done each day. % 
%                  Written by RG
% 
% monkeyNames:      A cell containing the names of which monkeys' data to
%                   analyze
%           'all' - Default. Runs through all folders stored in ~/<taskName>/Data/
%                   Note: Will also create a corresponding folder in
%                   ~/HcTask/Results/ if there is not one
% 
% 
% sessions:         A cell containing strings of session dates, if you 
%                   choose to specify which sessions to analyze
%           'all' - Default. Runs through each date stored in each monkey's
%                   data folder. 
%                   Note: Will also create a corresponding folder in
%                   ~/HcTask/Results/ if there is not one
%                   
%                   
% bSave:            Boolean identifier, whether to save taskStruct output or not. 
%                   false (default) or true
% 
% 
% 
% =========================================================================
%      Usage notes 
% =========================================================================
% 
% 1) FOLDER ORGANIZATION
% This function should be stored in the directory ~/<taskName>/Analyses
% It is important for for the pipeline that the folder hierarchy is
% preserved across machines. In all cases, the folder separation notation
% is determined upon function call, but this pipeline does assume a
% specific folder hierarchy. 
% 
% Folder hierarchy
% <taskName>
%   / Analyses
%       Stores high-level functions that are directly called within the 
%       analysis Master function
%       /Low-level Functions
%        Stores functions that are not directly called within your analysis
%        pipeline Master file. For functions that are analysis-specific,
%        they should be in a subfolder named after the higher-level 
%        function that calls it. 
%   / Data
%       / <monkeyNames>
%            / <sessions>
%              Each session folder should have all of the recorded data,
%              recording notes, and subfolders for each of the recorded
%              tasks. Named according to date in YYYYMMDD format. 
%                   / <taskNames>
%                     Each task should have all of its corresponding *.mat
%                     files (one for each trial run that day)
%   / Results
%       /SummaryResults
%       / <monkeyNames>
%            / <sessions>
% 
% This will run on all sessions that are in your ~/<ProjectName>/Data/ folder. 
% Each session should be a subfolder named in the style YYYYMMDD
% By default, the monkey names are Woody and Raul.
%  
% 
% *********************** Version History ********************************* 
% The original functions were written in part by RG with a deep
% contribution by GD and BC. GD is the expert for low-level Monkeylab and
% Unreal questions. BC is the expert for deep detail of any aspect of eye
% movements and eye data processing.
% 
% All of these functions were reformatted, cleaned and rewritten for
% distribution by RG in April/May 2017. Any further changes from this point
% forward should be dated and documented here. 
% 
% Change log: 
% - 20170921, by Rogelio Luna. "ODRTask_PrepBehavData" function added to prep data
%   from ODR task.
% - 20170914, by Rogelio Luna. A subfunction to look for the session folder
%   that will be processed was added. This is called from the Task loop
%   section. If the folder is not found, the loop skip it and continues to
%   next task in line.
% - 20170913, by Rogelio Luna. Minor style and typo corrections done.
%   "KMTasks_PrepBehavData" function added to prep data from "KeyMap" and
%   "KeyMapWM" tasks.
%  -20170915 - addition to task KeyMapWM to process multiple files of the same task type
% *************************************************************************

%% Initialization

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
defaultTaskNames = {'AL'};
defaultMonkeyNames = {'all'};
defaultSessions = {'all'};
defaultBSave = false;
% Parse inputs
addParameter(p,'taskNames',defaultTaskNames,@iscell);
addParameter(p,'monkeyNames',defaultMonkeyNames,@iscell);
addParameter(p,'sessions',defaultSessions,@iscell);
addParameter(p,'bSave',defaultBSave,@islogical);
parse(p,varargin{:});
% Define parsed variables
taskNames   = p.Results.taskNames;
monkeyNames = p.Results.monkeyNames;
sessions    = p.Results.sessions;
if p.Results.bSave
    bSave = 'true';
else
    bSave = 'false';
end

% Retrieve monkeyNames if they are not already specified. 
if ismember(monkeyNames,'all')
    dirDataFiles = dir(dirData);
    dirDataFiles(ismember({dirDataFiles.name},{'.','..'})) = [];
    dirDataFolders = {dirDataFiles([dirDataFiles.isdir]).name};
    monkeyNames = dirDataFolders;
    if ~iscell(monkeyNames)
        monkeyNames = {monkeyNames};
    end
end


% Monkey loop
for monkey = 1:length(monkeyNames)
    
    monkeyName = monkeyNames{monkey};
    
    % Get sessions
    if strcmp('all',sessions)
        dirDataFolders = dir(...
            [dirData monkeyName sl '2*']);
        sessions = {dirDataFolders.name}';
    end
    
    
    % Session loop
    for sess = 1:length(sessions)
        
        % Display monkey and session date
        session = sessions{sess};
        display(['Loading behavioural data for ' monkeyName ' ' session]);
        dirDataSess = [dirData monkeyName sl sessions{sess} sl];

        
        % Task loop.
        % CALLS TO TASK-SPECIFIC LOW-LEVEL FUNCTIONS
        for task = 1:length(taskNames)
            
            taskStruct = struct();
            taskName = taskNames{task};
            go = thereisfolder;             % Check that folder exists...
            if ~go, continue, end
            
            switch taskName
                case 'AL'
                    ALVR_PrepALData(dirDataSess,bSave)
                case 'EyeCal'
                    CSTask_PrepEyeCalEndData(dirDataSess, 'EyeCal',bSave)
                case 'EyeEnd'
                    CSTask_PrepEyeCalEndData(dirDataSess, 'EyeEnd',bSave)
                case 'KeyMap'
                    KMTasks_PrepBehavData(dirDataSess, 'KeyMap', bSave)
                case 'KeyMapWM'
                    KMTasks_PrepBehavData(dirDataSess, 'KeyMapWM', bSave)
                       files = dir(dirDataSess);
                       files = {files.name};
                       folders_of_interest = find(cell2mat(cellfun(@(x)(any(regexpi(x,taskName))),files,'uni',0)));
                       for folder_number = 1:length(folders_of_interest)
                           KMTasks_PrepBehavData(dirDataSess, files{folders_of_interest(folder_number)}, bSave)
                       end
                case {'ODR', 'ODRTask'}
                    ODRTask_PrepBehavData(dirDataSess, taskName, bSave)
            end

        end % task loop

        display(['Done processing ' num2str(sess) ' of ' num2str(length(sessions)) ' sessions.'])
    end % session loop
end % monkey loop

%% ===== Subfunctions =====
    function [go] = thereisfolder
        if exist([dirDataSess taskName],'dir')
            disp(['Processing folder ' dirDataSess taskName])
            go = true;
        else
            disp(['Folder ' dirDataSess taskName ' was not found!'])
            go = false;
            pause(2)
            disp(' ')
        end
    end % thereisfolder subfunction

end % function
