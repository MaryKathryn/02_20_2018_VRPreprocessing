function [] = KMTasks_PrepBehavData(dirDataSess, taskName, bSave)
% Extracts session and trials information from the "KeyMap" and/or
% "KeyMapWM" task folders.
% 
% This function was first created by RG, BC and GD to prepare data from the
% freeroam task. RL and MR modified the original function to prepare data
% from KeyMap and KeyMapWM tasks. RL further adapted this function to
% integrate it in the VR_preprocessing toolbox. 

% Create a string contaning the file location
dirTask = [dirDataSess taskName filesep];
pathParts = strsplit(dirTask, filesep);
session = pathParts{~cellfun(@isempty,strfind(pathParts,'2'))};

%Selects all files to open
Contents = dir(dirTask);

%We need to sort the files in their date order because by default the
%algorithm sorts them like: Trial1, Trial10, Trial11,... Trial2, Trial
%20,....

[~ , Order] = sort([Contents.datenum]);
Contents = Contents(Order);

ValidTrials = {Contents(cell2mat(cellfun(@(x) (~isempty(x)), (strfind({Contents.name}, 'Trial')), 'uni', false))).name};
EDFFile = {Contents(cell2mat(cellfun(@(x) (~isempty(x)), (strfind({Contents.name}, 'edf')), 'uni', false))).name};

%% Extract session and trial information

hWaitBar = waitbar(0, ['Going through ' taskName ' trials...']);

% load the EDF file that has the eye samples for the whole session
if ~isempty(EDFFile) && ~strcmp('el.edf', EDFFile)
    disp(['Loading data from EDF file called ' dirDataSess taskName filesep EDFFile{1}])
    edfEyeSamplefile = edfmex(char([dirDataSess taskName filesep EDFFile{1}]));
else
    disp('No EDF file was found. "Bad" samples will be used.')
    pause(3)
    edfEyeSamplefile = [];
end

%Static var
Task_Goals = [];

for trl = 1:numel(ValidTrials)
    tic
    
    % Initialize dynamic variables
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
    Task_GoalID         = [];
    Eye_Samples         = []; % Columns 1 - 3: EyeX EyeY MonkeyLab_Time pulled from EDF
    Eye_SamplesBad      = []; % Columns 1 - 3: EyeX EyeY MonkeyLab_Time
    Eye_Calibration     = []; % Columns 1 - 3: GainX GainY Offset; Rows 1 - 2: EyeX EyeY
    Eye_ScreenHeight    = [];
    Eye_ScreenWidth     = [];
    Eye_ScreenDistance  = [];
    Eye_EDFtimeDiff     = []; % Holds the difference in time between MLtime and eyelink time
    EyeLinkTime         = [];
    
    % Load the event data for the trial.
    events = cell(1,1);  %Need to predefine the variable to remove conflict with events function
    
    % Concatenates file path with file name, no need to change cd.
    load(char([dirDataSess taskName filesep ValidTrials{trl}]));
    eventsNumber = length(events);
    
    for k = 1:eventsNumber
        
        % Separates data according to event type
        switch events{k}.name
            case 'PathNodesInfos'
                
                %Gets the location of all goals in the environment
                temp = events{k};
                fields = fieldnames(temp);
                GoalsIDs = find(cell2mat(cellfun(@(x) (isempty(strfind('Goal', x(1)))), fields, 'uni', false)));
                
                for m=1:numel(GoalsIDs)
                    temp = rmfield(temp, fields{GoalsIDs(m)});
                end
                
                Task_Goals = temp;
                
            case 'PreStartTrialEvent'           % (2) ---------------------
                
                %Get which goal was selected
                Task_GoalID  = events{k}.GoalID;
                EDFTrialStartTime = events{k}.time;
                
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
        end
        
    end
    
     % Retrieve which goal the subject selected from the last PlayerState.
    if isempty(Unreal_States)
        continue
    end
    
    if ~strcmp(Trl_TrialID,'000000000000000')
        
        if ~isempty(Eye_EDFtimeDiff)
         edfDiff = Eye_EDFtimeDiff-Eye_EDFtimeDiff(1);
   
         %find the index in eyesamplesBad(:,3) that is greater than SOT,
         %and the index in eyesamplesBad(:,3) that is less that EOT.
         % edfDiff = EDFtimeDiff
         SOTIndex = find(Eye_SamplesBad(:,3)>EDFTrialStartTime,1);
         EOTindex = find(Eye_SamplesBad(:,3)<Trl_EOT,1,'last');
         edfDiff = Eye_EDFtimeDiff(SOTIndex:EOTindex)-Eye_EDFtimeDiff(SOTIndex);
         if max(abs(edfDiff))>.002
             
             figure
             plot(edfDiff)
         end
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
        taskStruct.SessionInfo.MonkeyName      = Ses_MonkeyName;      % (8)
        taskStruct.SessionInfo.Date            = Ses_Date;            % (8)
        taskStruct.SessionInfo.SetupID         = Sys_SetupID;         % (8)
        taskStruct.SessionInfo.Goals           = Task_Goals;          % (1)
        
        %If no eye calibration in last trial, use the last valid values
        if ~isempty(Eye_Calibration)
            taskStruct.SessionInfo.EyeCalibration  = Eye_Calibration;     %
            taskStruct.SessionInfo.QuadrantCorrect = QuadrantCorrect;     %
        end
        
        taskStruct.SessionInfo.ScreenHeight    = Eye_ScreenHeight;
        taskStruct.SessionInfo.ScreenWidth     = Eye_ScreenWidth;
        taskStruct.SessionInfo.ScreenDistance  = Eye_ScreenDistance;
        taskStruct.SessionInfo.FileName        =  ...
                                                [Ses_MonkeyName(1) '_' ...
                                                session '_' ...
                                                taskName '.mat'];
        % Trial-specific information
        taskStruct.Trials.(['ID_' Trl_TrialID]).GoalID         = Task_GoalID;         % (2)
        taskStruct.Trials.(['ID_' Trl_TrialID]).Unreal_States  = Unreal_States;       % (3)
        taskStruct.Trials.(['ID_' Trl_TrialID]).Unreal_PosXPosYRot = Unreal_PosRots;  % (3)
        taskStruct.Trials.(['ID_' Trl_TrialID]).Unreal_Times   = Unreal_Times;        % (3)
        taskStruct.Trials.(['ID_' Trl_TrialID]).EyeSamplesBad  = Eye_SamplesBad;      %
        taskStruct.Trials.(['ID_' Trl_TrialID]).EyeSamples     = Eye_Samples;         %
        taskStruct.Trials.(['ID_' Trl_TrialID]).ITC18StartTime = Sync_ITC18StartTime; % (10)
        taskStruct.Trials.(['ID_' Trl_TrialID]).SOT_Time       = Trl_SOT;             % (6)
        taskStruct.Trials.(['ID_' Trl_TrialID]).EOT_Time       = Trl_EOT;             % (7)
        taskStruct.Trials.(['ID_' Trl_TrialID]).SyncPulse      = Sync_Pulse;          %
        taskStruct.Trials.(['ID_' Trl_TrialID]).OutcomeWord    = Trl_Outcome;
     
 taskStruct.TrialSummary(trl+1,:) = {['ID_' Trl_TrialID] Trl_Outcome,Task_GoalID};
   
    end
    
    waitbar(trl/numel(ValidTrials));
    display([taskName ' Trial ' Trl_TrialID(end-2:end) ' processing time : ' num2str(toc, '%3.3f') ' seconds']);
    
end

close(hWaitBar);
% Sometimes the order gets screwed up when copying files over
taskStruct.Trials = orderfields(taskStruct.Trials);

% taskStruct = ML_FormatEyeDataOffScreen(taskStruct, 1);
% taskStruct = ML_FormatEyeDataOnScreen(taskStruct, 1);

% taskStruct = HcTask_FormatEyeDataOffScreen(taskStruct, 0); % Modified by Ben Corrigan summer 2016
bPreSmoothed = false;
ProcessEyesAndSave


%% ===== SUBFUNCTIONS =====
    function ProcessEyesAndSave
        
        % Cleaning and smoothing eye data ====================================
        timerFormatEye = tic;
        display('Cleaning and smoothing eye data');
        % Modify this to handle Eye_SamplesBad as an input rather than looking
        % for it within taskStruct
        taskStruct = HcTask_FormatEyeDataOffScreen(taskStruct, bPreSmoothed);
        bPreSmoothed = 1;
        
        % Classify eye movements
        hWaitBar3 = waitbar(0, 'Classifying eye movements.');
        trlIDs = fieldnames(taskStruct.Trials);
        for t = 1:length(trlIDs)
            taskStruct.Trials.(trlIDs{t}).SaccadeData = ...
                HcTask_SaccadeProcessing(taskStruct.Trials.(trlIDs{t}).EyeDegrees,bPreSmoothed);
            waitbar(t/length(trlIDs));
        end
        close(hWaitBar3);
        
        toc(timerFormatEye);
        %======================================================================
        
        % Save taskStruct ====================================================
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
            save([dirResultsSession filesep 'BehavData' filesep taskStruct.SessionInfo.FileName], 'taskStruct');
            
            clear taskStruct
            toc(timerSave);
        end
        %======================================================================
        
    end

end


