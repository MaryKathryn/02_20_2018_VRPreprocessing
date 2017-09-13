function [XMazeStruct] = HcTask_PrepXMazeData(SessionDir)

% Create a string contianing the file location
if ispc
    Slash = '\';
else
    Slash = '/';
end

%Selects all files to open
Contents = dir(SessionDir);

%We need to sort the files in their date order because by default the
%algorithm sorts them like: Trial1, Trial10, Trial11,... Trial2, Trial
%20,....

[~ , Order] = sort([Contents.datenum]);
Contents = Contents(Order);

ValidTrials = {Contents(cell2mat(cellfun(@(x) (~isempty(x)), (strfind({Contents.name}, 'Trial')), 'uni', false))).name};
EDFFile = {Contents(cell2mat(cellfun(@(x) (~isempty(x)), (strfind({Contents.name}, 'edf')), 'uni', false))).name};


%% Define static variables
Task_ColorList      = {'Invisible', 'Red', 'Green', 'Blue', 'Yellow', 'Orange', 'Purple', 'Cyan', 'Gray'}; %#ok<*NASGU>
Task_TextureList    = {'ML_Material_1', 'ML_Material_2', 'ML_Material_3', 'ML_Material_4', 'ML_Invisible'};
Task_TextureName    = {'Grass', 'Wood', 'Steel', 'Swan', 'Invisible'};
tempDxn             = {'N' 'S'};
Goals               = {'A' 'B' 'C' 'D'};
goalBaseStr         = 'Task_ColorList{events{k}.pop_GoalCtx';
rewBaseStr          = 'events{k}.RwdCtx';
chkBaseStr          = 'events{k}.chk_GoalCtx';
XMazeStruct.TrialSummary = {'Trial Number','Direction','Training Phase','Context','Goal West Color','Goal East Color','Outcome'};

%% Extract session and trial information

hWaitBar = waitbar(0, 'Going through XMaze trials.');

% load the EDF file that has the eye samples for the whole session
if ~strcmp('el.edf',EDFFile)&~isempty(EDFFile)
    edfEyeSamplefile = edfmex(char([SessionDir Slash EDFFile{1}] ));
else
    edfEyeSamplefile = [];
end

% define Ses to hold the trial params to compare across the session
  Ses                 = [];
for trl = 1:numel(ValidTrials)
    tic
    
    % Initialize dynamic variables
    eventsNumber        = [];
    Sys_SetupID         = 00;  %Need setup id for eye calibration data
    Ses_MonkeyName      = 'Default Monkey';
    Ses_Date            = '01-Jan-0000';
    Trl_TrialID         = '000000000000000';
    Trl_SOT             = 0;
    Trl_EOT             = 0;
    Trl_Outcome         = '';
    Sync_ITC18StartTime = 0;
    Sync_Pulse          = []; % Columns 1 - 2: Pulse Time
    Unreal_States       = {};
    Unreal_PosRots      = []; % Columns 1 - 3: PosX PosY Rot
    Unreal_Times        = []; % ML Time of UDK samples
    Task_Cue            = '';
    Task_Context        = 0;
    Task_GoalWest       = '';
    Task_GoalEast       = '';
    Task_TrainingPhase  = [];
    Eye_Samples         = []; % Columns 1 - 3: EyeX EyeY MonkeyLab_Time but taken from EDF
    Eye_SamplesBad      = []; % Columns 1 - 3: EyeX EyeY MonkeyLab_Time
    Eye_Calibration     = []; % Columns 1 - 3: GainX GainY Offset; Rows 1 - 2: EyeX EyeY
    Eye_ScreenHeight    = [];
    Eye_ScreenWidth     = [];
    Eye_ScreenDistance  = [];
    Eye_EDFtimeDiff     = []; % Holds the difference in time between MLtime and eyelink time
    EyeLinkTime         = [];
    EDFTrialStartTime   = [];
    % Load the event data for the trial.
    events = cell(1,1);  %Need to predefine the variable to remove conflict with events function
    
    % Concatenates file path with file name, no need to change cd.
    load(char([SessionDir Slash ValidTrials{trl}] ));
    eventsNumber = length(events);
    
    for k = 1:eventsNumber
        
        % Separates data according to event type
        switch events{k}.name
            
            case 'savingController::taskDialog' % (1) ---------------------
                
                % Create trlParams variable, which contains the goal
                % colour, goal value and wall material for the two contexts
                % and two directions as read from the HcTask GUI. Later on,
                % this information will be stored for the first trial as an
                % array in the Ses variable, and referenced on all
                % subsequent trials to ensure that none of these parameters
                % have changed. If they have changed, a warning is
                % reported.
                % This information is also later sorted, invisible goals
                % are taken out, and the colour hierarchy is stored for use
                % in later analyses.
                for ctxNum = 1:2
                    for dxnNum = 1:2
                        tempCtx = ['Context' num2str(ctxNum) tempDxn{dxnNum}];
                        trlParams.(tempCtx)      = {0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0};
                        for goalNum = 1:4
                            goalStr = [goalBaseStr num2str(ctxNum) tempDxn{dxnNum} Goals{goalNum} '}'];
                            rewStr  = [rewBaseStr num2str(ctxNum) tempDxn{dxnNum} Goals{goalNum}];
                            chkStr  = [chkBaseStr num2str(ctxNum) tempDxn{dxnNum} Goals{goalNum}];
                            % Create trlParams variable
                            trlParams.(tempCtx){1,goalNum} = Goals{goalNum}; % Goal position in ML GUI (A, B, C or D)
                            trlParams.(tempCtx){2,goalNum} = eval(goalStr);  % Goal colour in ML GUI (7 possible)
                            trlParams.(tempCtx){3,goalNum} = eval(rewStr);   % Reward duration in ML GUI (# of seconds)
                            trlParams.(tempCtx){4,goalNum} = eval(chkStr);   % Whether this reward was checked/activated in ML GUI
                        end
                        % Remove any columns with 'Invisible' specified as
                        % the reward colour.
                        if any(cell2mat(trlParams.(tempCtx)(4,:))==0)
                            rmvColumn = find(cell2mat(trlParams.(tempCtx)(4,:))==0);
                            trlParams.(tempCtx)(:,rmvColumn) = [];
                        end
                        % Append reward value hierarchy, and reorder columns according to hierarchy
                        [~,tempHier] = sort([trlParams.(tempCtx){3,:}],'descend');       % Returns the reward duration hierarchy for this context
                        trlParams.(tempCtx) = [trlParams.(tempCtx); num2cell(tempHier)]; % Appends the hierarchy to the trlParams.Context** variable
                        trlParams.(tempCtx) = trlParams.(tempCtx)(:,tempHier);           % Reorders the columns, such that the highest reward val is column 1
                    end
                end
                
                % Retrieve trial-specific Goal information using the Goals
                % specified in the PreStartTrial event.
                findWestGoalInd   = 'find(strcmp(Task_GoalWest,trlParams.(tempContext)(2,:)))';
                Task_GoalWestRew  = cell2mat(trlParams.(tempContext)(3,eval(findWestGoalInd)));
                Task_GoalWestHier = cell2mat(trlParams.(tempContext)(5,eval(findWestGoalInd)));
                findEastGoalInd   = 'find(strcmp(Task_GoalEast,trlParams.(tempContext)(2,:)))';
                Task_GoalEastRew  = cell2mat(trlParams.(tempContext)(3,eval(findEastGoalInd)));
                Task_GoalEastHier = cell2mat(trlParams.(tempContext)(5,eval(findEastGoalInd)));
                
                % Save Training Phase as recorded in the GUI during the
                % training session
                Task_TrainingPhase = events{k}.pop_trainingPhase;
                
                % Check if Goals or Reward durations have changed between
                % the first analyzed trial and the current trial.
                if ~isequal(trlParams.Context1N, trlParams.Context1S)
                    warning(['Context 1 North and South parameters not matching on Trial' Trl_TrialID(end-3:end) '!'])
                end
                if ~isequal(trlParams.Context2N, trlParams.Context2S)
                    warning(['Context 2 North and South parameters not matching on Trial' Trl_TrialID(end-3:end) '!'])
                end
                if isempty(Ses) % if Ses has not yet be defined, define it as TrlParams
                    % Store the context and goal information to be
                    % referenced in subsequent trials
                    Ses = trlParams;
                elseif ~isequal(trlParams, Ses)
                    % Check that the values in this trial match the
                    % reference cell values. If they do not, send an error
                    % message stating what does not match, and what the
                    % non-matching values are.
                    warning(['Task parameters not consistent between Trial' Trl_TrialID(end-3:end) ' and the first analyzed trial!']);
                end
                %
            case 'PreStartTrialEvent'           % (2) ---------------------
                
                % Get trial-specific Context, Direction, Cue, and Training Phase information
                Task_Context  = events{k}.Context;                          % Context output is a number, not a string
                Task_Cue      = events{k}.Textures;                         % Output is a string with the Material name as seen in the Content Browswer in UDK
                Task_CueMat   = cell2mat(Task_TextureName(find(strcmp(Task_Cue,Task_TextureList)))); %#ok<*FNDSB>
                if events{k}.bIsNorth == '1'
                    Task_Direction = 'North';
                elseif events{k}.bIsNorth == '0'
                    Task_Direction = 'South';
                else
                    warning('Unknown trial direction! (North vs South)');
                end
                
                % Get trial-specific Goal informations:
                tempGoals        = regexp(events{k}.GoalsColor,'/','split'); % The first item in this list is always the goal colour presented in the West side of the maze
                Task_GoalWest    = tempGoals{1};
                Task_GoalEast    = tempGoals{2};
                tempContext      = ['Context' num2str(Task_Context) 'N'];
                
                
                % ========================================================
                % Replaced Task_Direction == 'North' with strcmp which is
                % more stable and efficient. The == will compare each
                % letter individually and return 1 or 0 for each.
                if strcmp(Task_Direction, 'North')
                    Task_GoalLeft  = Task_GoalWest;
                    Task_GoalRight = Task_GoalEast;
                elseif strcmp(Task_Direction, 'South');
                    Task_GoalLeft  = Task_GoalEast;
                    Task_GoalRight = Task_GoalWest;
                else
                    warning('Could not determine GoalLeft/GoalRight identity!')
                end
                %
                %                 if Task_Direction == 'North'
                %                     Task_GoalLeft  = Task_GoalWest;
                %                     Task_GoalRight = Task_GoalEast;
                %                 elseif Task_Direction == 'South'
                %                     Task_GoalLeft  = Task_GoalEast;
                %                     Task_GoalRight = Task_GoalWest;
                %                 else
                %                     warning('Could not determine GoalLeft/GoalRight identity!')
                %                 end
                % ========================================================
                
                %
            case 'UnrealPlayerState'            % (3) ---------------------
                if ~isempty(events{k}.SampleTimes)
                    %Get Player State informations:
                    Unreal_States = [Unreal_States; regexp(events{k}.PlayerStates,'/','split')'];
                    
                    % Get Player Position and Rotation informations:
                    Unreal_PosRots = [Unreal_PosRots; events{k}.PlayerPositions(:,1:2) events{k}.PlayerRotations];
                    
                    % Get Unreal Time and convert it in MonkeyLab time to sync
                    Unreal_TempTime = regexp(events{k}.SampleTimes,'/','split')';
                    
                    % Convert UDK times from string to numerical values
                    % Times are now in vector: [YYYY MM DD hh mm ss.fff]
                    Unreal_TempTime = datevec(Unreal_TempTime, 'HH:MM:SS.FFF');
                    
                    if isfield(events{k},'MLQueryTimeVect')

                    Unreal_TempTime = Unreal_TempTime(:,4:6);
                    
                    %Now since the two computers are synched to ~1ms
                    %precision we just need to compute the way to translate
                    %the "vector" times (i.e. [2014 09 23 16 28 10.444] or
                    %[YYYY MM DD HH MM SS.FFF] into "GetSecs" times (i.e.
                    %4.555704997871990e+05). To do this we saved the
                    %variables : MLQueryTimeGetSecs and MLQueryTimeVect.
                    %This operation takes < 0.1 ms so we can assume both of
                    %these times to be equivalent. To get the UDK times in
                    %ML times we simply need to compute the difference in
                    %SECONDS between the UDK times and the MLQueryTimeVect
                    %and add this to the MLQueryTimeGetSecs.
                    
                    MLVectTime = repmat(events{k}.MLQueryTimeVect, size(Unreal_TempTime,1),1);
                    %Keep only the time (HH:MM:SS.FFF)
                    MLVectTime = MLVectTime(:,4:6);
                    
                    %Get the time difference in seconds
                    TimeDiff = Unreal_TempTime - MLVectTime;
                    TimeDiff = 3600*(TimeDiff(:,1)) + 60 * (TimeDiff(:,2)) + (TimeDiff(:,3));
                    
                    
                    Unreal_TempTime = events{k}.MLQueryTimeGetSecs + TimeDiff;
                    
                    % if this is a really old file, then use the old method
                    % of aligning MLtime and UDKtime
                    else
                    %----------------------------------------------------------
                    % HOW: We save the time the query has been received in UDK
                    % and when it was sent in ML. By equalling the
                    % UDKQueryTimes with the MLQueryTime we can rearrange the
                    % times. However we first need to convert the UDK times
                    % from string to time values.
                    %----------------------------------------------------------
                    
                    % Calculate delta time between query time
                    Unreal_QueryTime = datevec(events{k}.UDKQueryTimes);
                    
                    % Replicate line vector in matrix for substraction
                    Unreal_QueryTime = repmat(Unreal_QueryTime, size(Unreal_TempTime,1), 1);
                    
                    %Substract Query Time from sample times
                    Unreal_TempTime = Unreal_TempTime - Unreal_QueryTime;
                    
                    % Converts Times in seconds : HH:MM:SS.FFF = HH * 3600 + MM * 60 + SS.FFF
                    Unreal_TempTime = (Unreal_TempTime(:,4)*3600) + (Unreal_TempTime(:,5)*60) + (Unreal_TempTime(:,6));
                    
                    % Add the Query times in secs to the MonkeyLab Query Time
                    Unreal_TempTime = Unreal_TempTime + events{k}.MLQueryTime;
                    end
                    % Puts the values in the PosRotTime matrix plus
                    Unreal_Times = [Unreal_Times; Unreal_TempTime]; %#ok<*AGROW>
                end
                %
            case 'eyeSample'                    % (4) ---------------------
                
                % Since eye samples are sent one at a time, we do not need
                % to convert time to ML, we already have it.
                
                
                % File should have a sampleTimeSecs field that time stamps
                % the time at which the sample is taken and not when the
                % event is sent (time). However some older files might not
                % have it, so need to check.
                if isfield(events{k}, 'sampleTimeSecs')
                    Eye_SamplesBad = [Eye_SamplesBad; double(events{k}.eyeUnits) double(events{k}.sampleTimeSecs)];
                    Eye_EDFtimeDiff = [Eye_EDFtimeDiff; double(events{k}.sampleTimeSecs) - ...
                        double(events{k}.eyeLinkSampleTime)*.001];
                    EyeLinkTime = [EyeLinkTime; double(events{k}.eyeLinkSampleTime)*.001];
                else %if not, then use time. CAREFUL: there might be a delay between sampleTimeSecs and time.
                    Eye_SamplesBad = [Eye_SamplesBad; double(events{k}.eyeUnits) double(events{k}.time)];
                    Eye_EDFtimeDiff = [Eye_EDFtimeDiff; double(events{k}.time) - ...
                        double(events{k}.eyeLinkSampleTime)*.001];
                end
                
                
                %Removed because of conflicting information between
                %eyeCalibration and savingController::eyeDialog. Probably
                %because the eyeCalibration event doesn't take into account
                %last eye calibration trial.
                % 5) eyeCalibration -------------------------------------------
                %             case 'eyeCalibration'               % (5) ---------------------
                
                % Gets the info from the EyeCalibration Event
                % We only get how to convert the units to degrees, we will
                % need further processing to convert to pixels.
                %                 Eye_Calibration = events{k}.units2Deg;
                %                 QuadrantCorrect = events{k}.QuadrantCorrect;
                
                
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
                
                % 6) endOfTrial -----------------------------------------------
            case 'endOfTrial'                   % (7) ---------------------
                
                % Get Trial outcome
                Trl_Outcome = events{k}.eot;
                
                % Get End of trial time
                Trl_EOT = events{k}.time;
                
                %
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
                
                % Get Session date
                Ses_Date = events{k}.date;
                
                % get start of Trial for EDF calculations
                if isempty(EDFTrialStartTime)
                    EDFTrialStartTime = events{k}.time;
                end
                
                % 8) SyncPulseRead --------------------------------------------
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
                
                
                % 9) ITC18TimeZero --------------------------------------------
            case 'ITC18TimeZero'                % (10)---------------------
                
                % Get ITC 18 start collection time
                Sync_ITC18StartTime = events{k}.time;
                
            case 'sendTrialID'
                EDFTrialStartTime = events{k}.time;
        end
        
    end
    
    
    
    % Retrieve which goal the subject selected from the last PlayerState.
    if isempty(Unreal_States)
        continue
    end
    tempOutcome = cell2mat(Unreal_States(length(Unreal_States)));
    % Ensure that this result includes a Goal.
    % If true, determine which goal was selected (colour & hierarchy).
    % Otherwise, warn that this trial was not completed.
    if ~isempty(strfind(tempOutcome, 'Goal'))
        if     tempOutcome(end) == 'W';
            Task_OutcomeColor = Task_GoalWest;
            Task_OutcomeHier  = Task_GoalWestHier;
            Task_OutcomeDirec = 'West';
            if strcmp(Task_Direction,'North')
                Task_OutcomeSide = 'Left';
            elseif strcmp(Task_Direction,'South')
                Task_OutcomeSide = 'Right';
            else
                warning('Could not determine Outcome side! (Left vs Right)');
            end              % Determine if Left or Right
            if ~isempty(Task_GoalWestHier) && ~isempty(Task_GoalEastHier)
                if Task_GoalWestHier(1) < Task_GoalEastHier(1);
                    Task_OutcomeCorr = 'Correct';
                else
                    Task_OutcomeCorr = 'Incorrect';
                end % Determine if Correct or Incorrect
            else
                Task_OutcomeCorr = 'BadTrial';
            end
        elseif tempOutcome(end) == 'E';
            Task_OutcomeColor = Task_GoalEast;
            Task_OutcomeHier  = Task_GoalEastHier;
            Task_OutcomeDirec = 'East';
            if strcmp(Task_Direction,'North')
                Task_OutcomeSide = 'Right';
            elseif strcmp(Task_Direction,'South')
                Task_OutcomeSide = 'Left';
            else
                warning('Could not determine Outcome side! (Left vs Right)');
            end               % Determine if Left or Right
            if ~isempty(Task_GoalWestHier) && ~isempty(Task_GoalEastHier)
                if Task_GoalEastHier(1) < Task_GoalWestHier(1);
                    Task_OutcomeCorr = 'Correct';
                else
                    Task_OutcomeCorr = 'Incorrect';
                end % Determine if Correct or Incorrect
            else
                Task_OutcomeCorr = 'BadTrial';
            end
        end
    else
        warning(['Last PlayerState in Trial ' Trl_TrialID(end-2:end) ' is not a goal location!']);
        if strcmp(Trl_Outcome,'userStoppedTrial');
            Task_OutcomeCorr = 'UserStopped';
            warning('Trial terminated by the user.');
        else
            Task_OutcomeCorr = 'Unknown';
            warning('Trial NOT terminated by user.');
        end
        Task_OutcomeColor = [];
        Task_OutcomeHier  = [];
        Task_OutcomeDirec = [];
        Task_OutcomeSide = [];
    end
    
    
    if ~strcmp(Trl_TrialID,'000000000000000')
        
        if ~isempty(Eye_EDFtimeDiff)
            edfDiff = Eye_EDFtimeDiff-Eye_EDFtimeDiff(1);
            
            %find the index in eyesamplesBad(:,3) that is greater than the start of the trial,
            %and the index in eyesamplesBad(:,3) that is less that EOT.
            % edfDiff = EDFtimeDiff
            SOTIndex = find(Eye_SamplesBad(:,3)>EDFTrialStartTime,1);
            EOTindex = find(Eye_SamplesBad(:,3)<Trl_EOT,1,'last');
            edfDiff = Eye_EDFtimeDiff(SOTIndex:EOTindex)-Eye_EDFtimeDiff(SOTIndex);
            %          figure
            %          plot(edfDiff)
            
            % Find the time eyelink time that aligns with the start of the trial
            % (the end of the last trial) and the end of this trial. get the
            % indices of these times, and then get the eye samples for these
            % periods
            if ~isempty(edfEyeSamplefile)
            Trl_SOTEDFindex = find(edfEyeSamplefile.FSAMPLE.time>...
                (EDFTrialStartTime-Eye_EDFtimeDiff(SOTIndex))*1000,1);
            Trl_EOTEDFindex = find(edfEyeSamplefile.FSAMPLE.time<...
                (Trl_EOT-Eye_EDFtimeDiff(SOTIndex))*1000,1,'last');
            Eye_Samples = [double(edfEyeSamplefile.FSAMPLE.hx(1,Trl_SOTEDFindex:Trl_EOTEDFindex))' ...
                double(edfEyeSamplefile.FSAMPLE.hy(1,Trl_SOTEDFindex:Trl_EOTEDFindex))'...
                double(edfEyeSamplefile.FSAMPLE.time(Trl_SOTEDFindex:Trl_EOTEDFindex))'*.001+Eye_EDFtimeDiff(1)...
                double(edfEyeSamplefile.FSAMPLE.pa(1,Trl_SOTEDFindex:Trl_EOTEDFindex))'];
            else
                Eye_Samples = Eye_SamplesBad;
            end
        end
        
        % Create the output MLStruct variable
        
        % Session-specific information
        XMazeStruct.SessionInfo.MonkeyName      = Ses_MonkeyName;      % (8)
        XMazeStruct.SessionInfo.Date            = Ses_Date;            % (8)
        XMazeStruct.SessionInfo.SetupID         = Sys_SetupID;         % (8)
        XMazeStruct.SessionInfo.Context1NParams = trlParams.Context1N; % (1)
        XMazeStruct.SessionInfo.Context1SParams = trlParams.Context1S; % (1)
        XMazeStruct.SessionInfo.Context2NParams = trlParams.Context2N; % (1)
        XMazeStruct.SessionInfo.Context2SParams = trlParams.Context2S; % (1)
        
        %If no eye calibration in last trial, use the last valid values
        if ~isempty(Eye_Calibration)
            XMazeStruct.SessionInfo.EyeCalibration  = Eye_Calibration;     %
            XMazeStruct.SessionInfo.QuadrantCorrect = QuadrantCorrect;     %
        end
        
        XMazeStruct.SessionInfo.ScreenHeight    = Eye_ScreenHeight;
        XMazeStruct.SessionInfo.ScreenWidth     = Eye_ScreenWidth;
        XMazeStruct.SessionInfo.ScreenDistance  = Eye_ScreenDistance;
        
        
        % Trial-specific information
        XMazeStruct.Trials.(['ID_' Trl_TrialID]).Context        = Task_Context;        % (2)
        XMazeStruct.Trials.(['ID_' Trl_TrialID]).ContextName    = Task_Cue;            % (2)
        XMazeStruct.Trials.(['ID_' Trl_TrialID]).ContextMat     = Task_CueMat;         % (2)
        XMazeStruct.Trials.(['ID_' Trl_TrialID]).Direction      = Task_Direction;      % (2)
        XMazeStruct.Trials.(['ID_' Trl_TrialID]).TrainingPhase  = Task_TrainingPhase;  % (1)
        XMazeStruct.Trials.(['ID_' Trl_TrialID]).GoalWestColor  = Task_GoalWest;       % (2)
        XMazeStruct.Trials.(['ID_' Trl_TrialID]).GoalWestRew    = Task_GoalWestRew;    % (2)
        XMazeStruct.Trials.(['ID_' Trl_TrialID]).GoalWestHier   = Task_GoalWestHier;   % (2)
        XMazeStruct.Trials.(['ID_' Trl_TrialID]).GoalEastColor  = Task_GoalEast;       % (2)
        XMazeStruct.Trials.(['ID_' Trl_TrialID]).GoalEastRew    = Task_GoalEastRew;    % (2)
        XMazeStruct.Trials.(['ID_' Trl_TrialID]).GoalEastHier   = Task_GoalEastHier;   % (2)
        XMazeStruct.Trials.(['ID_' Trl_TrialID]).GoalLeft       = Task_GoalLeft;       % (2)
        XMazeStruct.Trials.(['ID_' Trl_TrialID]).GoalRight      = Task_GoalRight;      % (2)
        XMazeStruct.Trials.(['ID_' Trl_TrialID]).OutcomeCorr    = Task_OutcomeCorr;    % (3)
        XMazeStruct.Trials.(['ID_' Trl_TrialID]).OutcomeColor   = Task_OutcomeColor;   % (3)
        XMazeStruct.Trials.(['ID_' Trl_TrialID]).OutcomeHier    = Task_OutcomeHier;    % (3)
        XMazeStruct.Trials.(['ID_' Trl_TrialID]).OutcomeDirec   = Task_OutcomeDirec;   % (3)
        XMazeStruct.Trials.(['ID_' Trl_TrialID]).OutcomeSide    = Task_OutcomeSide;    % (3)
        XMazeStruct.Trials.(['ID_' Trl_TrialID]).OutcomeWord    = Trl_Outcome;         % (7) Change this to be 'Completed' rather than 'Correct'
        XMazeStruct.Trials.(['ID_' Trl_TrialID]).Unreal_States  = Unreal_States;       % (3)
        XMazeStruct.Trials.(['ID_' Trl_TrialID]).Unreal_PosXPosYRot = Unreal_PosRots;  % (3)
        XMazeStruct.Trials.(['ID_' Trl_TrialID]).Unreal_Times   = Unreal_Times;        % (3)
        XMazeStruct.Trials.(['ID_' Trl_TrialID]).EyeSamplesBad  = Eye_SamplesBad;         %
        XMazeStruct.Trials.(['ID_' Trl_TrialID]).EyeSamples     = Eye_Samples;         %
        XMazeStruct.Trials.(['ID_' Trl_TrialID]).ITC18StartTime = Sync_ITC18StartTime; % (10)
        XMazeStruct.Trials.(['ID_' Trl_TrialID]).SOT_Time       = Trl_SOT;             % (6)
        XMazeStruct.Trials.(['ID_' Trl_TrialID]).EOT_Time       = Trl_EOT;             % (7)
        XMazeStruct.Trials.(['ID_' Trl_TrialID]).SyncPulse      = Sync_Pulse;          %
        
        
        XMazeStruct.TrialSummary(trl+1,:) = {['ID_' Trl_TrialID] Task_Direction Task_TrainingPhase...
        Task_CueMat Task_GoalWest Task_GoalEast Trl_Outcome};
    end
    
    
    
    
    waitbar(trl/numel(ValidTrials));
    display(['XMaze Trial ' Trl_TrialID(end-2:end) ' processing time : ' num2str(toc, '%3.3f') ' seconds']);
    
end

close(hWaitBar);

XMazeStruct = HcTask_FormatEyeDataOffScreen(XMazeStruct, 1); % Modified by Ben Corrigan summer 2017


%% Retrieve the moment of goal appearance and first rotation for all trials
% Get index of the moment the goals appear for this trial

trlIDs = fieldnames(XMazeStruct.Trials);
for trl = 1:length(trlIDs)
    goalsAppearInd = [];
    %     if strcmp(XMazeStruct.Trials.(trlIDs{trl}).Direction, 'North')
    %         if any(XMazeStruct.Trials.(trlIDs{trl}).PositionMatrix_deg(:,1) >= 768)
    %             goalsAppearInd  = ...
    %                 find(XMazeStruct.Trials.(trlIDs{trl}).PositionMatrix_deg(:,1) ...
    %                 >= 768, 1, 'first');
    %         end
    %     elseif strcmp(XMazeStruct.Trials.(trlIDs{trl}).Direction, 'South')
    %         if any(XMazeStruct.Trials.(trlIDs{trl}).PositionMatrix_deg(:,1) <= -786)
    %             goalsAppearInd  = ...
    %                 find(XMazeStruct.Trials.(trlIDs{trl}).PositionMatrix_deg(:,1)...
    %                 <= -786, 1, 'first');
    %         end
    %     end
    
    %GoalsOnset Time
    goalsAppearInd = find(strcmp(XMazeStruct.Trials.(trlIDs{trl}).Unreal_States, ['Goals' XMazeStruct.Trials.(trlIDs{trl}).Direction 'On']),1, 'first');
    XMazeStruct.Trials.(trlIDs{trl}).GoalsOnset = ...
        XMazeStruct.Trials.(trlIDs{trl}).Unreal_Times(goalsAppearInd);
    
    %Context Onset Time: since the change in player state doesn't work, we
    %will rely on position and by taking into acount the collision radius
    %of the Pawn, which is: ~34 units
    if strcmp(XMazeStruct.Trials.(trlIDs{trl}).Direction, 'North')
        if any(XMazeStruct.Trials.(trlIDs{trl}).Unreal_PosXPosYRot(:,1) >= (-703-34))
            
            XMazeStruct.Trials.(trlIDs{trl}).ContextOnset = XMazeStruct.Trials.(trlIDs{trl}).Unreal_Times(find(XMazeStruct.Trials.(trlIDs{trl}).Unreal_PosXPosYRot(:,1) >= (-703-34), 1, 'first'));
            %
            %             goalsAppearInd  = ...
            %                 find(XMazeStruct.Trials.(trlIDs{trl}).PositionMatrix_deg(:,1) ...
            %                 >= 768, 1, 'first');
        else
            XMazeStruct.Trials.(trlIDs{trl}).ContextOnset = [];
        end
    elseif strcmp(XMazeStruct.Trials.(trlIDs{trl}).Direction, 'South')
        if any(XMazeStruct.Trials.(trlIDs{trl}).Unreal_PosXPosYRot(:,1) <= (704+48))
            
            XMazeStruct.Trials.(trlIDs{trl}).ContextOnset = XMazeStruct.Trials.(trlIDs{trl}).Unreal_Times(find(XMazeStruct.Trials.(trlIDs{trl}).Unreal_PosXPosYRot(:,1) <= (704+48), 1, 'first'));
            
            %             goalsAppearInd  = ...
            %                 find(XMazeStruct.Trials.(trlIDs{trl}).PositionMatrix_deg(:,1)...
            %                 <= -786, 1, 'first');
        else
            XMazeStruct.Trials.(trlIDs{trl}).ContextOnset = [];
        end
    end
    
    %     XMazeStruct.Trials.(trlIDs{trl}).PositionMatrix_deg(goalsAppearInd,6);
    
    % Determine at which index the monkey turned away from looking straight
    % North or straight south by thresholding the Unreal_PosRot value. A
    % difference in rotation of one degree in either direction from the
    % subject's rotation when the goals appeared is used as the threshold.
    goalsRotTime = [];
    if ~isempty(goalsAppearInd) && (strcmp(XMazeStruct.Trials.(trlIDs{trl}).OutcomeCorr, 'Correct')||...
            strcmp(XMazeStruct.Trials.(trlIDs{trl}).OutcomeCorr, 'Incorrect'))
        
        %Find the first sample where the rotation crosses the threshold
        RotThreshInd = find(abs(...
            XMazeStruct.Trials.(trlIDs{trl}).Unreal_PosXPosYRot(goalsAppearInd:end,3)-...
            XMazeStruct.Trials.(trlIDs{trl}).Unreal_PosXPosYRot(goalsAppearInd,3)) > 10,...
            1, 'first')+goalsAppearInd;
        
        %Go back to the first sample where the rotation was initiated
        RotStartInd = find(...
            XMazeStruct.Trials.(trlIDs{trl}).Unreal_PosXPosYRot(goalsAppearInd:RotThreshInd,3) == ...
            XMazeStruct.Trials.(trlIDs{trl}).Unreal_PosXPosYRot(goalsAppearInd,3),1,'last')+goalsAppearInd;
        
        goalsRotTime = XMazeStruct.Trials.(trlIDs{trl}).Unreal_Times(RotStartInd) ;
    end
    XMazeStruct.Trials.(trlIDs{trl}).DecisionOnset = goalsRotTime;
end



end
