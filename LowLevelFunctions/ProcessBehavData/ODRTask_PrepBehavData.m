function [taskStruct] = ODRTask_PrepBehavData(dirDataSess, taskName, bSave)
% Extracts session and trials information from the oculomotor delayed
% response (ODR) task folders. This function was adapted to integrate it in
% the VR_preprocessing toolbox.
 
% Gets the string containing the session's folder and selects all files to
% open
dirTask = [dirDataSess taskName filesep];
pathParts = strsplit(dirTask, filesep);
session = pathParts{~cellfun('isempty', strfind(pathParts, '201'))};

Contents = dir(dirTask);

% We need to sort the files in their date order because by default the
% algorithm sorts them like: Trial1, Trial10, Trial11,... Trial2, Trial
% 20, and so on. Next, the valid trials are collected.

[~ , Order] = sort([Contents.datenum]);
Contents = Contents(Order);

ValidTrials = {Contents(cell2mat(cellfun(@(x) (~isempty(x)), ...
    (strfind({Contents.name}, 'Trial')), 'uni', false))).name};

% Loads the EDF file that has the eye samples for the whole session
EDFFile = {Contents(cell2mat(cellfun(@(x) (~isempty(x)), ...
    (strfind({Contents.name}, 'edf')), 'uni', false))).name};

if ~strcmp('el.edf',EDFFile) && ~isempty(EDFFile)
    disp('Eyelink Data File found')
    edfEyeSamplefile = edfmex(char([dirTask EDFFile{1}]));
    disp(' ')
else
    disp('Eyelink Data File WAS NOT found. "Bad" samples will be used.')
    pause(2)
    edfEyeSamplefile = [];
end

%% Define static variables
% Define column labels for the Summary table created at the end
SummaryLabels = {'TrialNumber' ...
                 'WMTrial' ...
                 'FPColor' ...
                 'TargetColor' ...
                 'XYPointsPosition' ...
                 'TrialOutcome' ...
                 'CorrectCounter' ...
                 'FailCounter' ...
                 'StopStimulusTime' ...
                 'HideFPTime' ...
                 'Trial SOT' ...
                 'Trial EOT'};
        
% Define the states of interest to be retrieved. These strings will be
% later used as labels for the Summary table as well.
statesOfInterest = {'interTrialState'; ...
                    'PreFixationStartState'; ...
                    'WaitFixation'; ...
                    'PreTargetHold'; ...
                    'TargetOnset'; ...
                    'PostTargetHold'; ...
                    'FixPointOffset'; ...
                    'endTrialState'};


%% Extract session and trial information

% Initiate the waitbar
hWaitBar = waitbar(0, ['Going through ' taskName ' trials...']);
                
% Defines TrialSummary headers
taskStruct.TrialSummary = cell(1,10);
    %{'Trial Number','Direction','Training Phase',...
    %     'Context','Test','Test','Test','Test','Test','Goal West Color',
    %     'Goal East Color','Outcome'};

% define Ses to hold the trial params to compare across the session
% Ses = [];
for trl = 1:numel(ValidTrials)
%     ValidTrials{trl}
    tic
    
    % Sets some counters
    caseCounterA = 1;
    
    % Initialize dynamic variables
    Sys_SetupID         = 00;  %Need setup id for eye calibration data
    Ses_MonkeyName      = 'Default Monkey';
    Ses_Date            = '01-Jan-0000';
    Trl_TrialID         = '000000000000000';
    Trl_SOT             = 0;
    Trl_EOT             = 0;
    %     Trl_Behav           = '';
    Trl_Outcome         = '';
    Fail_Counter        = 0;
    Correct_Counter     = 0;
    Sync_ITC18StartTime = 0;
    Sync_Pulse          = []; % Columns 1 - 2: Pulse Time
    Eye_Samples         = []; % Columns 1 - 3: EyeX EyeY MonkeyLab_Time but taken from EDF
    Eye_SamplesBad      = []; % Columns 1 - 3: EyeX EyeY MonkeyLab_Time
    Eye_Calibration     = []; % Columns 1 - 3: GainX GainY Offset; Rows 1 - 2: EyeX EyeY
    Eye_ScreenHeight    = [];
    Eye_ScreenWidth     = [];
    Eye_ScreenDistance  = [];
    Eye_EDFtimeDiff     = []; % Holds the difference in time between MLtime and eyelink time
    Eye_LinkTime        = [];
    EDFTrialStartTime   = [];
    Task_States         = {};
    Eye_ScreenFrames    = 0;
    Task_Params         = struct;
    StateTimes          = struct;

    
    % Load the event data for the trial.
    events = cell(1,1);  %Need to predefine the variable to remove conflict with events function
    
    % Concatenates file path with file name, no need to change cd.
    load(char([dirTask ValidTrials{trl}] ));
    eventsNumber = length(events);
    
    %% Separates data according to type of event
    for k = 1:eventsNumber
        
        switch events{k}.name
            
            case 'SyncPulseRead'                % (9) ---------------------
                
                % =========================================================
                % Changed the event name to SyncPulseRead and the field
                % name to SyncPulse because it was not reading it properly.
                % =========================================================
                
                % Get Pulse points
                tempPulse = double(events{k}.SyncPulse);
                tempPulseTime = zeros(size(tempPulse,1), size(tempPulse,2));
                
                % Validate that we have data
                if ~isempty(tempPulse)
                    tempPulseTime(1, size(tempPulse,2)) = events{k}.time;
                    
                    % Add Data to variable
                    Sync_Pulse = double([Sync_Pulse; tempPulse' tempPulseTime']);
                end
                
            case 'eyeSample'                    % (4) ---------------------
                
                % Since eye samples are sent one at a time, we do not need
                % to convert time to ML, we already have it.
                
                % File should have a sampleTimeSecs field that time stamps
                % the time at which the sample is taken and not when the
                % event is sent (i. e., "time"). However some older files
                % might not have it, so need to check.
                
                if isfield(events{k}, 'sampleTimeSecs')
                    Eye_SamplesBad = [Eye_SamplesBad; ...
                        double(events{k}.eyeUnits) ...
                        double(events{k}.sampleTimeSecs)];
                    Eye_EDFtimeDiff = [Eye_EDFtimeDiff; ...
                        double(events{k}.sampleTimeSecs) - ...
                        double(events{k}.eyeLinkSampleTime)...
                        *.001];
                    Eye_LinkTime = [Eye_LinkTime; ...
                        double(events{k}.eyeLinkSampleTime)...
                        *.001];
                    
                else %if not, then use time. CAREFUL: there might be a
                    %delay between "sampleTimeSecs" and "time".
                    
                    Eye_SamplesBad = [Eye_SamplesBad; ...
                        double(events{k}.eyeUnits) ...
                        double(events{k}.time)];
                    Eye_EDFtimeDiff = [Eye_EDFtimeDiff; ...
                        double(events{k}.time) - ...
                        double(events{k}.eyeLinkSampleTime)...
                        *.001];
                end
                
            case 'frame'
                Eye_ScreenFrames(events{k}.frame) = events{k}.flipTime;
                
                
            case 'savingController::taskDialog' % (1) ---------------------
                
                % Create "taskParams" variable, which contains the
                % different delay times, inter-trial intervals, Fixation
                % point and Target sizes, FP and Target windows sizes, and
                % some other. All of these parameters are stored in a
                % single structure, which will be attach as a field in the
                % output structure. Please note that this function is not
                % tracking changes over any of these parameters from trial
                % to trial.
                
                % Checks that parameters are coming from the right GUI
                if ~strfind(events{k}.FileName,'ODRTask')
                    error('Information contained within this file may not correspond to the ODR task')
                end
                
                % Retrieves task's parameters:
                Task_Params.FailDelay       = events{k}.txt_FailDelay;
                Task_Params.WMPercent       = events{k}.txt_HidePercent;
                Task_Params.RT              = events{k}.txt_RT;
                Task_Params.WMenabled       = events{k}.chk_HideTarget;
                Task_Params.Reward          = events{k}.txt_Reward;
                Task_Params.TargetFix       = events{k}.txt_TargetFix;
                Task_Params.PostTargDelay   = [events{k}.txt_maxPostTarget ...
                    events{k}.txt_minPostTarget];
                Task_Params.TargDelay       = [events{k}.txt_maxTarget ...
                    events{k}.txt_minTarget];
                Task_Params.PreTargDelay    = [events{k}.txt_maxPreTarget ...
                    events{k}.txt_minPreTarget];
                Task_Params.PreTrialEyeDisp = events{k}.txt_PreTrialEyeDisplay;
                Task_Params.PreTrialFixTime = events{k}.txt_PreTrialFix;
                Task_Params.ITI             = [events{k}.txt_maxITI ...
                    events{k}.txt_minITI];
                Task_Params.TargetSize      = events{k}.txt_TargetSize;
                Task_Params.TargetWinSize   = events{k}.txt_TargetWin;
                Task_Params.FixWinSize      = events{k}.txt_FixWin;
                Task_Params.RandFixPoint    = events{k}.chk_RandFixPoint;
                Task_Params.TargetColor     = events{k}.pop_TargetColor;
                Task_Params.FixColor        = events{k}.pop_FixColor;
                Task_Params.BGIntensity     = events{k}.txt_BGIntensity;
                Task_Params.FixPtSize       = events{k}.txt_FixPtSize;
                Task_Params.NbrV            = events{k}.txt_NbrV;
                Task_Params.NbrH            = events{k}.txt_NbrH;
                Task_Params.VSep            = events{k}.txt_VSeparation;
                Task_Params.HSep            = events{k}.txt_HSeparation;
                
                
            case 'savingController::eyeDialog'  % (6) ---------------------
                
                % The eye calibration data might also be found in the
                % saving of the eyeDialog window.
                Eye_Calibration(1,1) = events{k}.a;
                Eye_Calibration(1,2) = events{k}.b;
                Eye_Calibration(1,3) = events{k}.c;
                Eye_Calibration(1,4) = events{k}.d;
                Eye_Calibration(1,5) = events{k}.e;
                
                Eye_Calibration(2,1) = events{k}.f;
                Eye_Calibration(2,2) = events{k}.g;
                Eye_Calibration(2,3) = events{k}.h;
                Eye_Calibration(2,4) = events{k}.i;
                Eye_Calibration(2,5) = events{k}.j;
                
                QuadrantCorrect(1,1) = events{k}.m_1;
                QuadrantCorrect(1,2) = events{k}.m_2;
                QuadrantCorrect(1,3) = events{k}.m_3;
                QuadrantCorrect(1,4) = events{k}.m_4;
                
                QuadrantCorrect(2,1) = events{k}.n_1;
                QuadrantCorrect(2,2) = events{k}.n_2;
                QuadrantCorrect(2,3) = events{k}.n_3;
                QuadrantCorrect(2,4) = events{k}.n_4;
                
            case 'savingController::screenDialog'  % (7) ---------------------
                
                Eye_ScreenHeight    = events{k}.screenHeightCM*2;
                Eye_ScreenWidth     = events{k}.screenWidthCM*2;
                Eye_ScreenDistance  = events{k}.screenDistanceCM;
                
                
            case 'stateAction'
                % Get the trial states, from which different stimulus onset
                % and offset times can be retrieved
                
                Task_States{caseCounterA,1} = events{k}.stateName;
                Task_States{caseCounterA,2} = events{k}.time;
                
                
                caseCounterA = caseCounterA + 1;

                
                % 6) endOfTrial -----------------------------------------------
            case 'endOfTrial'                   % (7) ---------------------
                
                % Get Trial outcome
                Trl_Outcome = events{k}.eot;
                
                % Get End of trial time
                Trl_EOT = events{k}.time;
                

            case 'startOfTrial'                 % (8) ---------------------
                
                % Get Setup ID from trial ID;
                % Trial ID  = DDD DSS SSS ssT TTT
                %   D = Days since Jan 1 2012
                %   S = Seconds of day
                %   s = setup ID
                %   T = Trial number
                Trl_TrialID = events{k}.trialID;
                Sys_SetupID = Trl_TrialID(10:11);  % 01 = Cerebus, 02 = Plexon
                
                %Get Start of trial Time
                Trl_SOT = events{k}.time;
                
                % Get Monkey Name
                Ses_MonkeyName = events{k}.MonkeyName;
                
                % Get Session data
                %                 Ses_Task = events{k}.taskName;
                Ses_Date = events{k}.date;
                
                % get start of Trial for EDF calculations
                if isempty(EDFTrialStartTime)
                    EDFTrialStartTime = events{k}.time;
                end

                
            case 'stimulusMovieReady'
                
                FixPointColor = events{k}.FixPointColor;
                TargetColor = events{k}.TargetColor;
                PointsPosXY = [events{k}.FixPointX ...
                    events{k}.FixPointY; ...
                    events{k}.TargetX ...
                    events{k}.TargetY];
                WMTrial = events{k}.IsODR;
                
                % PENDING: recover both fixation and target point size (or
                % diameter)
                
%             %% ERROR: this time actually corresponds to the EOT time                
%             case 'stopStimulus'
%                 % Get the stop stimulus time
%                 stopStimulusTime = events{k}.time;
%                 
%             %% ERROR: this time actually corresponds to the EOT time                 
%             case 'hideFixationPoint'
%                 % Get the hide fixation point time
%                 hideFixPointTime = events{k}.time;
                
                
            case 'fixationPointGrid'
                % Get the grid of possible X and Y positions for the
                % fixation point
                fpGridX = events{k}.fpgX;
                fpGridY = events{k}.fpgY;
                
                
            case 'requestStimulusMovie'
                HideTarget = events{k}.HideTargetFlag;
                
                
            case 'rtpCounters'
                % Get the trial outcome.
                Fail_Counter = events{k}.Fail;
                Correct_Counter = events{k}.Correct;
                
                
            case 'ITC18TimeZero'                % (10)---------------------
                % Get ITC 18 start collection time
                Sync_ITC18StartTime = events{k}.time;
                
        end
        
    end
    
    %% Pre-processes eye data
    
    if ~strcmp(Trl_TrialID,'000000000000000')
        
        if ~isempty(Eye_EDFtimeDiff)
            %             edfDiff = Eye_EDFtimeDiff-Eye_EDFtimeDiff(1);
            
            % find the index in eyesamplesBad(:,3) that is greater than the
            % start of the trial, and the index in eyesamplesBad(:,3) that
            % is less that EOT.
            % edfDiff = EDFtimeDiff
            SOTIndex = find(Eye_SamplesBad(:,3) > EDFTrialStartTime,1);
            EOTindex = find(Eye_SamplesBad(:,3) < Trl_EOT,1,'last');
            %             edfDiff = Eye_EDFtimeDiff(SOTIndex:EOTindex)-Eye_EDFtimeDiff(SOTIndex);
            %          figure
            %          plot(edfDiff)
            
            % Find the eyelink time that aligns with the start of the trial
            % (the end of the last trial) and the end of this trial. get the
            % indices of these times, and then get the eye samples for these
            % periods
            if ~isempty(edfEyeSamplefile) && ~isempty(SOTIndex)
                Trl_SOTEDFindex = find(edfEyeSamplefile.FSAMPLE.time > ...
                    (EDFTrialStartTime-Eye_EDFtimeDiff(SOTIndex))*1000,1);
                
                Trl_EOTEDFindex = find(edfEyeSamplefile.FSAMPLE.time < ...
                    (Trl_EOT-Eye_EDFtimeDiff(SOTIndex))*1000,1,'last');
                
                Eye_Samples = [double(edfEyeSamplefile.FSAMPLE.hx(1,Trl_SOTEDFindex:Trl_EOTEDFindex))' ...
                    double(edfEyeSamplefile.FSAMPLE.hy(1,Trl_SOTEDFindex:Trl_EOTEDFindex))'...
                    double(edfEyeSamplefile.FSAMPLE.time(Trl_SOTEDFindex:Trl_EOTEDFindex))'*.001+Eye_EDFtimeDiff(1)...
                    double(edfEyeSamplefile.FSAMPLE.pa(1,Trl_SOTEDFindex:Trl_EOTEDFindex))'];
            else
                disp('    Eye Samples bad were used for this trial!')
                pause(1);
                Eye_Samples = Eye_SamplesBad;
            end
        end
        
        %% Collects times from the states of interest
        if ~isempty(Task_States)
            for k = 1:numel(statesOfInterest)
                found = strcmpi(Task_States(:,1),statesOfInterest{k});
                if numel(find(found)) > 1
                    Itemp = find(found,1,'last');
                    found = false(numel(found),1);
                    found(Itemp) = true;
                end
                if any(found)
                   eval(['StateTimes.' statesOfInterest{k} ' = Task_States{found,2};'])
                else
                    eval(['StateTimes.' statesOfInterest{k} ' = [];'])
                end
            end
        else
            error('States time values were not found')
        end
        
        
        %% Create the output MLStruct variable
        
        % Session-specific information
        taskStruct.SessionInfo.MonkeyName      = Ses_MonkeyName;      % (8)
        taskStruct.SessionInfo.Date            = Ses_Date;            % (8)
        taskStruct.SessionInfo.SetupID         = Sys_SetupID;         % (8)
        
        %If no eye calibration in last trial, use the last valid values
        if ~isempty(Eye_Calibration)
            taskStruct.SessionInfo.EyeCalibration  = Eye_Calibration;
            taskStruct.SessionInfo.QuadrantCorrect = QuadrantCorrect;
        end
        
        taskStruct.SessionInfo.ScreenHeight    = Eye_ScreenHeight;
        taskStruct.SessionInfo.ScreenWidth     = Eye_ScreenWidth;
        taskStruct.SessionInfo.ScreenDistance  = Eye_ScreenDistance;
        taskStruct.SessionInfo.FileName        =  ...
                                                [Ses_MonkeyName(1) '_' ...
                                                session '_' ...
                                                taskName '.mat'];
        % Trial-specific information
        taskStruct.Trials.(['ID_' Trl_TrialID]).Task_Params    = Task_Params;         % (2)
        taskStruct.Trials.(['ID_' Trl_TrialID]).FixPointColor  = FixPointColor;       % (2)
        taskStruct.Trials.(['ID_' Trl_TrialID]).TargetColor    = TargetColor;         % (2)
        taskStruct.Trials.(['ID_' Trl_TrialID]).PointsPosXY    = PointsPosXY;         % (2)
        taskStruct.Trials.(['ID_' Trl_TrialID]).WMTrial        = WMTrial;             % (1)
%         taskStruct.Trials.(['ID_' Trl_TrialID]).StopStimTime   = stopStimulusTime;    % ERROR: see above and fix it if needed
%         taskStruct.Trials.(['ID_' Trl_TrialID]).HideFPTime     = hideFixPointTime;    % ERROR: see above and fix it if needed
        taskStruct.Trials.(['ID_' Trl_TrialID]).FPGridX        = fpGridX;             % ()
        taskStruct.Trials.(['ID_' Trl_TrialID]).FPGridY        = fpGridY;             % ()
        taskStruct.Trials.(['ID_' Trl_TrialID]).HideTarget     = HideTarget;          % ()
        taskStruct.Trials.(['ID_' Trl_TrialID]).Task_States    = Task_States;         % (2)
        taskStruct.Trials.(['ID_' Trl_TrialID]).ScreenFrames   = Eye_ScreenFrames;    % (2)
        taskStruct.Trials.(['ID_' Trl_TrialID]).CorrectCount   = Correct_Counter;     % (2)
        taskStruct.Trials.(['ID_' Trl_TrialID]).FailCount      = Fail_Counter;        % (2)
        taskStruct.Trials.(['ID_' Trl_TrialID]).EyeSamplesBad  = Eye_SamplesBad;      % ()
        taskStruct.Trials.(['ID_' Trl_TrialID]).EyeSamples     = Eye_Samples;         % ()
        taskStruct.Trials.(['ID_' Trl_TrialID]).ITC18StartTime = Sync_ITC18StartTime; % (10)
        taskStruct.Trials.(['ID_' Trl_TrialID]).SOT_Time       = Trl_SOT;             % (6)
        taskStruct.Trials.(['ID_' Trl_TrialID]).EOT_Time       = Trl_EOT;             % (7)
        taskStruct.Trials.(['ID_' Trl_TrialID]).Outcome        = Trl_Outcome;         % (7)
        % ODRStruct.Trials.(['ID_' Trl_TrialID]).Behavior       = Trl_Behav;           % (7)
        taskStruct.Trials.(['ID_' Trl_TrialID]).SyncPulse      = Sync_Pulse;          % (7)
        taskStruct.Trials.(['ID_' Trl_TrialID]).States         = StateTimes;
        
        % Trials summary information
        taskStruct.TrialSummary(trl,:) = ...
            {['ID_' Trl_TrialID] ...    % Trial ID
            WMTrial ...                 % WM or Perception trial
            FixPointColor ...           % Fixation Point Colour
            TargetColor ...             % Target Colour
            PointsPosXY ...             % XY position of grid points
            Trl_Outcome ...             % Trial outcome
            Correct_Counter ...         % Counter for correct trials
            Fail_Counter ...            % Counter for fail trials
            Trl_SOT ...                 % Start of Trial time
            Trl_EOT};                   % End of trial time
%             stopStimulusTime ...        % Stimulation stopping time     % ERROR, see above.
%             hideFixPointTime ...        % Time of Fixation Point hiding % ERROR, see above.
        
         StateFields = fieldnames(StateTimes);
%          lastcol = numel(ODRStruct.TrialSummary(trl,:));
         for kk = 1:numel(StateFields)
             eval(['ODRStruct.TrialSummary{trl,12+kk}' ...
                   ' = StateTimes.(StateFields{kk});'])
         end

        display(['ODR Trial ' Trl_TrialID(end-2:end) ...
            ' processing time : ' num2str(toc, '%3.3f') ' seconds']);
    else
        display(['ODR Trial contained in file ' ValidTrials{trl} ...
            ' was not processed.'])
    end
    waitbar(trl/numel(ValidTrials));
    
end
close(hWaitBar);

% Sometimes the order gets screwed up when copying files over
taskStruct.Trials = orderfields(taskStruct.Trials);

disp('Data prep of ODR folder is done')

%% Processes the Eye signals recovered from this session

% taskStruct = HcTask_FormatEyeDataOffScreen(taskStruct, 0);
% taskStruct = ML_FormatEyeDataOnScreen(taskStruct, 0);
bPreSmoothed = false;

%   Cleaning and smoothing eye data ====================================
timerFormatEye = tic;
display('Cleaning, formatting and smoothing eye data');

% Modify this to handle Eye_SamplesBad as an input rather than looking
% for it within taskStruct
taskStruct = HcTask_FormatEyeDataOffScreen(taskStruct, bPreSmoothed);

% Classify eye movements
bPreSmoothed = true;    

display('Classifying eye movements');
hWaitBar3 = waitbar(0, 'Classifying eye movements...');
trlIDs = fieldnames(taskStruct.Trials);
for t = 1:length(trlIDs)
    if isfield(taskStruct.Trials.(trlIDs{t}), 'EyeDegrees')
    taskStruct.Trials.(trlIDs{t}).SaccadeData = ...
        HcTask_SaccadeProcessing(taskStruct.Trials.(trlIDs{t}).EyeDegrees, bPreSmoothed);
    else
        taskStruct.Trials.(trlIDs{t}).SaccadeData = [];
        disp([' --- Warning:  EyeDegrees field not found in ' trlIDs{t}])
        pause(2)
    end
    waitbar(t/length(trlIDs));
end
close(hWaitBar3);

toc(timerFormatEye);
%   ======================================================================

%   Save taskStruct ====================================================
if bSave
    display('Saving taskStruct')
    timerSave = tic;
    
    % taskStruct.BlockInfo.SkippedTrials = SkippedTrials;
    
    dirResultsSession = strrep(dirDataSess,'Data','Results');
    if ~isdir(dirResultsSession)
        mkdir(dirResultsSession)
    end
    
    %Update FileName Info in structure
    if exist([dirResultsSession filesep 'BehavData'],'dir') ~= 7
        disp(['Creating folder ' dirResultsSession filesep 'BehavData'])
        mkdir([dirResultsSession filesep 'BehavData'])
    end
    save([dirResultsSession 'BehavData' filesep taskStruct.SessionInfo.FileName], 'taskStruct');
    
    clear taskStruct
    toc(timerSave);
end
%   ======================================================================

end % End of function
