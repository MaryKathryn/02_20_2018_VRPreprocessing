function [SelectedTrials] = HcTask_SelectTrials(MLStruct, Options)

% Function built by RG to return trialIDs that meet user-specified
% criteria. 
% 
% 
%%% Trial selection
% trialCriteria is userinput, and contains any trial selection criteria for
% desired for processing in the next steps. This MUST be the first cell
% input into the function.
% 
% Options:
% 
%  Options.RemoveInvalidTrials
%       Boolean. Removes any trials with the outcome other than 'Correct'
%       or 'Incorrect' (e.g. {'Time Run Out', 'userStoppedTrial',
%       'BadTrial', 'badTrial'}).
% 
%  Options.RemoveShortTrials;
%       Boolean. Only in VR tasks. This must be modified if there are any
%       tasks where the trial start pathnode is moved.
%       Remove all trials that started with the animal in the black box or
%       at the TrialStart path node location in the maze. This happens on
%       any time the task resumes after being paused (e.g. the first trial
%       of the day), and on trials and trials after a time out. (e.g.
%       subject completes a North trial and then just sits there; the
%       subsequent South trial times out, and the subject is transported
%       back to the trial start path node in the middle of the maze).
% 
%  Options.SelectTrials
%       Boolean. Select trials according to Options.TrialCriteria       
% 
%  Options.TrialCritera is a cell (x by 2):
%       The first column is the criterion, and the second the value. Everything
%       on a single line is interpreted as an OR, different lines are
%       interpreted as an AND. For example:
%       {
%       'OutcomeWord',                           {'Correct', 'Incorrect'};
%       {'GoalWestColor' 'GoalEastColor'},       Red
%       'Direction',                             'North'
%       }
%       This would return any trial where the outcome is EITHER Correct or
%       Incorrect AND where ANY of the two presented goals is Red AND where the
%       direction is North.
% 
%  Options.Save;
%       Save output the session's results folder (or prompts user for save
%       folder)
% 
% 
%%% Example Usage %%%
% 
% 


%% Initialize Variables
Slash = filesep;
AllTrials = fieldnames(MLStruct.Trials);
NbrTrials = numel(AllTrials);


% Remove invalid/incomplete trials
if ~isfield(Options,'RemoveInvalidTrials')
    Options.RemoveInvalidTrials = 0;
end
% RemoveShortTrials
if ~isfield(Options,'RemoveShortTrials')
    Options.RemoveShortTrials = 0;
end
% Select trials based on parameters
if ~isfield(Options,'SelectTrials')
    Options.SelectTrials = 0;
elseif ~isfield(Options.TrialCriteria)
    help HcTask_SelectTrials
    tempIDs = fieldNames(MLStruct.Trials);
    display(fieldnames(MLStruct.Trials.(tempIDs{1})));
    Options.TrialCriteria = input(...
        'Input cell for SelectTrialOptions.TrialCriteria now.');
end
% Save Selected trials
if ~isfield(Options,'Save')
    Options.Save = 0;
end



%% Remove trials that are shorter than 150ms
% This only typically happens in FreeRoam, when a trial starts and the
% monkey already happens to be in spot where the goal appeared. These
% trials don't have enough usable eye data or position data.

thresh = 75;
fewSamples = cell2mat(struct2cell(structfun(@(x)...
    (length(x.EyeSamples)),MLStruct.Trials,'uni',0)...
    )) < thresh+1;

if any(fewSamples)
    MLStruct.Trials = rmfield(MLStruct.Trials,AllTrials(fewSamples));
    AllTrials = fieldnames(MLStruct.Trials);
    NbrTrials = numel(AllTrials);
end



%% Option: Remove invalid trials (i.e., not Correct or Incorrect)
% Note: Trials that are completed correctly are either coded as 'Correct'
% or 'Incorrect' in the taskStruct.Trials.(trlID).OutcomeWord in EyeCal and
% most VR tasks. However, in the XMazeStruct, the reliable trial outcome
% field is 'OutcomeCorr'.
if Options.RemoveInvalidTrials == 1
    
    isInvalidTrial = false(NbrTrials,2);
    isInvalidTrial = ~sum(cell2mat(struct2cell(structfun(@(x) ...
        strcmp(x.OutcomeWord, {'Correct' 'Incorrect'}), MLStruct.Trials, ... Look for all trials with x.OutcomeWord matching any outcomesToLookFor
        'uni',false))),2);
    
    if any(strfind(MLStruct.SessionInfo.FileName,'XMaze'))
        
        isInvalidTrial(:,2) = ~sum(cell2mat(struct2cell(structfun(@(x) ...
        strcmp(x.OutcomeCorr, {'Correct' 'Incorrect'}), MLStruct.Trials, ... Look for all trials with x.OutcomeWord matching any outcomesToLookFor
        'uni',false))),2);
        
        isInvalidTrial = any(isInvalidTrial,2);
        
    end
    
    MLStruct.Trials = rmfield(MLStruct.Trials,AllTrials(isInvalidTrial));
    AllTrials = fieldnames(MLStruct.Trials);
    NbrTrials = numel(AllTrials);
        
end % Remove Incomplete Trials
    


%% Option: Remove short trials
if Options.RemoveShortTrials
    
    isShortTrial = [];
    % Find trials that started in the black box
    isShortTrial = [ AllTrials(...
                        structfun(@(x)(...
                        x.PositionMatrix_deg(1,2)>700),MLStruct.Trials)) ; ...
                    isShortTrial];
    % Find trials that started at the trial start pathnode 
    isShortTrial = [ AllTrials(...
                        structfun(@(x)(...
                        all(x.PositionMatrix_deg(1,1:2) == [175.89 0])),MLStruct.Trials)) ; ...
                    isShortTrial];            
    % Remove trials
    MLStruct.Trials = rmfield(MLStruct.Trials,isShortTrial);
    AllTrials = fieldnames(MLStruct.Trials);
    NbrTrials = numel(AllTrials);
    
end % Remove Short Trials
    
 

%% Trial selection
if Options.SelectTrials == 1 
    
    trialCriteria = Options.TrialCriteria;
    
    
    % Select trials according to trialCriteria
    
    
    NbrCriteria = size(trialCriteria,1);
    AllMatch = zeros(NbrTrials, NbrCriteria);
    
    % If trialCriteria were specified, return trials where the condition is
    % met for a given field of MLStruct.Trials.(Trial)
    for CptrCriteria = 1:NbrCriteria % CptrCriteria = criteria counter in french.
        
        if iscell(trialCriteria{CptrCriteria,1})
            NbrFields = numel(trialCriteria{CptrCriteria,1});
        else
            NbrFields = 1;
        end
        
        Match = zeros(NbrTrials, NbrFields);
        
        Conditions = trialCriteria{CptrCriteria, 2};
        
        for CptrField = 1:NbrFields
            
            %Read the condition and variable for each row
            % Defines the field in which certain conditions will be looked for
            if NbrFields >= 2
                Fields = char(trialCriteria{CptrCriteria,1}(CptrField));
            else
                Fields = char(trialCriteria{CptrCriteria,1});
            end
            
            %Find the all of the conditions for that field
            tempConditions = struct2cell(structfun(@(x) (x.(Fields)), MLStruct.Trials, 'uniformoutput', false));
            
            if ischar(tempConditions{1})
                tempMatch = sum(cell2mat(cellfun(@(x) (strcmp(x, Conditions)), tempConditions, 'uni',false)),2);
            else
                tempMatch = cell2mat(cellfun(@(x) sum(ismember(x, Conditions)), tempConditions, 'uni',false));
            end
            
            Match(:,CptrField) = tempMatch;
        end
        
        AllMatch(:,CptrCriteria) = sum(Match,2);
    end
    
    %Only keep trials that match ALL contidions
    GoodTrials = sum(AllMatch, 2) == size(AllMatch,2);
    
    % Create the output varible and store the trial ID headers that meet
    % the crieteria.
    % If more than one condition is specified, only trials that meet all of
    % the conditions are included.
    SelectedTrials = AllTrials(GoodTrials);

% If no trial critera have been selected, keep all of the trials for
% subsequent processing
else 
    SelectedTrials = fieldnames(MLStruct.Trials);
    trialCriteria = [];
    
end


%% Prepare output variable
% NOTE: Changed name of output variable and folder to SelectTrials from
% SelectedTrials to match the function name

SelectTrials.SessionInfo = MLStruct.SessionInfo;
SelectTrials.Trials      = SelectedTrials;
SelectTrials.Criteria    = trialCriteria;
SelectTrials.Options     = Options;


%% Option: Save output variable
if Options.Save ==1
    
    FileName = MLStruct.SessionInfo.FileName;
    
    %Get the session date from the filename
    % Looks for 4-8 consecutive digits in the session name i.e. :
    %       20150101 or 150101 or 0101
    % We aren't using the 4 and 6 digits dates anymore but just to make sure
    % the script works.
    [startIndex,endIndex] = regexp(FileName, '\d{4,8}');
    SessionDate = FileName(startIndex:endIndex);
    
    %Gets directory of executed script, should be ..\HcTask\Analyses\
    AnalysesFolder = mfilename('fullpath');
    AnalysesFolder = strsplit(AnalysesFolder, 'Analyses');
    AnalysesFolder = AnalysesFolder{1};
    
    
    if ~isempty(AnalysesFolder) && ~isempty(SessionDate)
        OutputFolder = [AnalysesFolder 'Results' Slash MLStruct.SessionInfo.MonkeyName Slash SessionDate Slash];
    else
        OutputFolder = uigetdir(AnalysesFolder, 'Select Output Dir');
    end
    
     OutputFolder = [OutputFolder Slash 'SelectTrials'];
    if ~exist(OutputFolder, 'dir')
        mkdir(OutputFolder)
    end
%     save(OutputFile, 'SelectTrials', 'trialCriteria', 'Options')
    save([OutputFolder Slash 'SelectTrials_' FileName], 'SelectTrials')
    
    if Options.TrialConditions == 1
        save([OutputFolder Slash 'TrialConditions_' FileName], 'TrialConditions')
    end
 
end



end
