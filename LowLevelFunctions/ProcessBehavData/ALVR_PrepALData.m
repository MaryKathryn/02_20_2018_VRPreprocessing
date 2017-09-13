function [taskStruct] = ALVR_PrepALData(dirDataSess,bSave)

% Function written by RG to pre-process behavioural data from the virtual
% reality associative learning task. Outputs taskStruct containing all
% behavioural data.
%
%
% =========================================================================
% Output architecture:
% =========================================================================
% %
% % taskStruct.SessionInfo.
% %          MonkeyName
% %                Date
% %             SetupID
% %     Context1NParams
% %     Context1SParams
% %     Context2NParams
% %     Context2SParams
% %      EyeCalibration
% %     QuadrantCorrect
% %        ScreenHeight
% %         ScreenWidth
% %      ScreenDistance
% %            FileName
% %
% % taskStruct.Trials.(trlID). ->
% %     TrlParams.
% %         TrainingPhase
% %         TrlDirection
% %         ContextMat
% %         GoalWestColor
% %         GoalEastColor
% %         GoalRightColor
% %         GoalLeftColor
% %         CorrectGoalDirec
% %         CorrectGoalColor
% %         CorrectGoalSide
% %
% %     Events_Outcomes.
% %         SyncPulse
% %         ITC18Start_Time
% %         SOT_Time
% %         ContextOnset
% %         GoalsOnset
% %         DecisionOnset
% %         EOT_Time
% %         DecisionChange
% %         ChosenDirec
% %         ChosenSide
% %         ChosenColor
% %         OutcomeCorr
% %         OutcomeWord
% %         RewardDuration
% %
% %     Place_EyeData. ->
% %         PositionMatrix_deg
% %         EyeSamples
% %         PupilSize
% %         SaccadeData
% %             StartTime
% %             EndTime
% %             PeakTime
% %             PeakVelocity
% %             Duration
% %             Amplitude
% %             StartPointX
% %             StartPointY
% %             EndPointX
% %             EndPointY
% %             Direction
% %             PostSaccadicOscillationEnd
% %             isSaccade
% %             isPSO
% %             isFixation
% %             isSmoothPursuit
% %             isOffScreen
% %         GazeTargets
% %         GazePosition
% %
% %     RawData. ->
% %         Unreal_PosXPosYRot
% %         Unreal_States
% %         Unreal_Times
% %
% =========================================================================
%
%
% =========================================================================
% Event types
% =========================================================================
% 'EyelinkTimeZero'
% 'ITC18TimeZero'
%       Time of the ITI; first timestamp of the trial; precedes SOT_Time
% 'PreStartTrialEvent'
%       Contains trial-specific parameters
% 'SpeedAndFixPointUnreal'
%       Not used
% 'SyncPulseRead'
%       Contains SyncPulse sent to the Cerebus computer that is read on
%       channel ainp9 on both systems
% 'ToggleFixPointUnreal'
%       Not used
% 'UnrealPlayerState'
%       Contains Unreal Player states; defined in Kismet. One per frame
% 'endOfTrial'
% 'eyeCalibration'
%       Not used. Because of conflicting information between eyeCalibration
%       and savingController::eyeDialog. Probably because the
%       eyeCalibration event doesn't take into account last eye calibration
%       trial.
% 'eyeSample'
%       File should have a sampleTimeSecs field that time stamps the time
%       at which the sample is taken and not when the event is sent (time)
% 'playSound'
%       Not used. Correct or incorrect sound at the end of the trial.
% 'requestNewEyeCalibration'
%       Not used. At beginning of trial. Unsure of intended use (RG)
% 'rtpCounters'
%       Not used. For real time plots.
% 'saveDataToFile'
% 'savingController::DAQDialogValues'
% 'savingController::eyeDialog'
% 'savingController::fakeMonkeyDialog'
%       Not used
% 'savingController::screenDialog'
%       Not used
% 'savingController::taskDialog'
%       Contains all critical user-defined task parameters. Used to create
%       BlockInfo
% 'savingInfo'
%       Not used
% 'sendTrialID'
% 'startCollectingData'
% 'startOfTrial'
% 'startReward'
%       Contains reward length. Should be consistent with taskDialog.
% 'startStateSystem'
% 'stateAction'
%       Not used
% 'stopCollectingData'
% 'taskDialogValues'
%       Not used. Redundant with 'savingController::taskDialog'
%
%
% =========================================================================
% USAGE NOTES
% =========================================================================
% All AL task folders should be in a seperate subfolder of the Data folder
% for this session
% e.g. dir('E:\8a_AssociativeLearning\Data\Theo\20170414\')
%
% .                        DL_imp_end.png           NSP0_001.ns6             VL_imp_end.txt
% ..                       DL_spikePanel_end.png    NSP1_001.ccf             VL_spikePanel_end.png
% ALFixedEnd               DL_spikePanel_start.png  NSP1_001.nev             VL_spikePanel_start.png
% ALFixedStart             EyeCal                   NSP1_001.ns2             recording_notes.txt
% ALNovel                  EyeEnd                   NSP1_001.ns6
% DL_ct_start.png          NSP0_001.ccf             VL_ct_start.png
% DL_ct_start.txt          NSP0_001.nev             VL_ct_start.txt
% DL_imp.txt               NSP0_001.ns2             VL_imp_end.png
%
%
% =========================================================================
% Version History
% =========================================================================
%
% 20170501:
% Function Created by Roberto Gulli. The oldest predecessor to this
% function was built by Guillaume Doucet and Roberto Gulli. The functions
% called here to clean and process eye movements were written by Ben
% Corrigan. Details of these functions can be found in Corrigan, Gulli,
% Doucet & Martinez-Trujillo, 2017.
%
% Please document any modification to this function here:
%
%
%
%
%
%
% Example session: Theo April 13, 2017



%% Set up
sl = filesep;
pathParts = strsplit(strtrim(dirDataSess),filesep);
monkeyName = pathParts{find(strcmp(pathParts,'Data'))+1};
session = pathParts{end-1};

% Find all directories with trials from the associative learning task.
dirDataSessFiles = dir(dirDataSess);
dirDataSessFiles(ismember({dirDataSessFiles.name},{'.','..'})) = [];
dirDataSessFolders = {dirDataSessFiles([dirDataSessFiles.isdir]).name};
taskFolders = dirDataSessFolders(~cellfun(@isempty,strfind(dirDataSessFolders,'AL')));


% Task color and materials list
Task_ColorList      = {'Invisible', 'Red', 'Green', 'Blue', 'Yellow', 'Orange', 'Purple', 'Cyan', 'Gray'};
Task_TextureList    = {'ML_Material_1', 'ML_Material_2', 'ML_Material_3', 'ML_Material_4', 'ML_Invisible'};
Task_TextureName    = {'Swan', 'Wood', 'Steel', 'Grass', 'Invisible'};


%% Run through this entire function once for each task block
for task = 1:length(taskFolders)
    
    taskName = taskFolders{task};
    dirTask = [dirDataSess taskName sl];
    
    
    % Selects all files to open. Then we need to sort the files in their date
    % order because by default the algorithm sorts them like: Trial1, Trial10,
    % Trial11,... Trial2, Trial20,...
    taskFiles  = dir(dirTask);
    [~,dirFilesOrder] = sort([taskFiles.datenum]);
    taskFiles  = taskFiles(dirFilesOrder);
    
    
    % LOAD THE EDF FILE
    % Contains eye samples for the whole session
    EDFFile = {taskFiles(cell2mat(cellfun(@(x) (~isempty(x)), ...
        (strfind({taskFiles.name}, 'edf')), 'uni', false))).name};
    % Load
    if ~strcmp('el.edf',EDFFile)&&~isempty(EDFFile)
        edfEyeSamplefile = edfmex(char([dirTask EDFFile{1}] ));
    else
        edfEyeSamplefile = [];
    end
    
    
    % Retrieve trial filenames
    trlNames = {taskFiles(cell2mat(cellfun(@(x)(~isempty(x)), ...
        (strfind({taskFiles.name}, 'Trial')), 'uni', false))).name};
    
    hWaitBar = waitbar(0, 'Going through trials...'); % Starting waiting bar...
    blockCounter = 1;
    blockTrlCounter = 1;
    % TRIAL LOOP
    for trl = 1:numel(trlNames)
        timerTrl = tic;
        disp(['  Processing file ' trlNames{trl} '...'])
        
        
        % Initialize dynamic variables:
        eventsNumber        = [];
        Trl_TrialID         = '000000000000000';
        Trl_TrainingPhase   = uint8([]);
        Trl_Context         = '';
        Trl_GoalWest        = '';
        Trl_GoalEast        = '';
        Trl_OutcomeWord     = '';
        Trl_GoalWestValue   = [];
        Trl_GoalEastValue   = [];
        Trl_GoalLeftValue   = [];
        Trl_GoalRightValue  = [];
        Rew_Duration        = [];
        Trl_OutcomeDirec    = '';
        Trl_OutcomeColor    = '';
        Trl_OutcomeSide     = '';
        Sync_ITC18StartTime = 0;
        Trl_SOT             = 0;
        Trl_EOT             = 0;
        Sync_Pulse          = []; % Columns 1 - 2: Pulse, Time
        Unreal_States       = {};
        Unreal_PosXPosYRot  = []; % Columns 1 - 3: PosX PosY Rot
        Unreal_Times        = []; % ML Time of UDK samples
        EyeSamples          = []; % Columns 1 - 3: EyeX EyeY MonkeyLab_Time but taken from EDF
        Eye_SamplesBad      = []; % Columns 1 - 3: EyeX EyeY MonkeyLab_Time
        Eye_Calibration     = []; % Columns 1 - 3: GainX GainY Offset; Rows 1 - 2: EyeX EyeY
        Eye_EDFtimeDiff     = []; % Holds the difference in time between MLtime and eyelink time
        EyeLinkTime         = [];
        InvalidTrials       = [];
        
        
        % Load events
        events = cell(1,1);  % Must initialize variable to avoid conflicts with the native function 'events'
        % Concatenates file path with file name, no need to change cd.
        load(char([dirTask trlNames{trl}]));
        eventNames = cellfun(@(x)(x.name),events,'uni',0)';
        display([num2str(length(events)) ' events in this trial.']);
        
        %% SessionInfo ====================================================
        sotInd = find(strcmp(eventNames,'startOfTrial'),1,'last');
        eyeDlgInd = find(strcmp(eventNames,'savingController::eyeDialog'),1,'last');
        screenDlgInd = find(strcmp(eventNames,'savingController::screenDialog'),1,'last');
        unrealStatesInd = find(strcmp(eventNames,'UnrealPlayerState'),1,'first');
        if any(cellfun(@isempty,{sotInd,eyeDlgInd,screenDlgInd,unrealStatesInd}))
            warning(['Skipped trial ' trlNames{trl} ' due to missing event types!'])
            InvalidTrials = [InvalidTrials; trlNames(trl)];
            continue
        elseif length(eventNames)<200
            warning(['Skipped trial ' trlNames{trl} '; too few events!'])
            InvalidTrials = [InvalidTrials; trlNames(trl)];
            continue
        end
        
        
        if blockTrlCounter == 1
            % MonkeyName
            taskStruct.SessionInfo.MonkeyName = monkeyName;
            % Date
            taskStruct.SessionInfo.Date = events{sotInd}.date;
            % SetupID
            taskStruct.SessionInfo.SetupID = events{sotInd}.trialID(10:11); % 02 = Recording Room 35B
            % EyeCalibration
            Eye_Calibration(1,1) = events{eyeDlgInd}.a;
            Eye_Calibration(1,2) = events{eyeDlgInd}.b;
            Eye_Calibration(1,3) = events{eyeDlgInd}.c;
            Eye_Calibration(1,4) = events{eyeDlgInd}.d;
            Eye_Calibration(1,5) = events{eyeDlgInd}.e;
            Eye_Calibration(2,1) = events{eyeDlgInd}.f;
            Eye_Calibration(2,2) = events{eyeDlgInd}.g;
            Eye_Calibration(2,3) = events{eyeDlgInd}.h;
            Eye_Calibration(2,4) = events{eyeDlgInd}.i;
            Eye_Calibration(2,5) = events{eyeDlgInd}.j;
            taskStruct.SessionInfo.EyeCalibration = Eye_Calibration;
            % QuadrantCorrect
            QuadrantCorrect(1,1) = events{eyeDlgInd}.m_1;
            QuadrantCorrect(1,2) = events{eyeDlgInd}.m_2;
            QuadrantCorrect(1,3) = events{eyeDlgInd}.m_3;
            QuadrantCorrect(1,4) = events{eyeDlgInd}.m_4;
            QuadrantCorrect(2,1) = events{eyeDlgInd}.n_1;
            QuadrantCorrect(2,2) = events{eyeDlgInd}.n_2;
            QuadrantCorrect(2,3) = events{eyeDlgInd}.n_3;
            QuadrantCorrect(2,4) = events{eyeDlgInd}.n_4;
            taskStruct.SessionInfo.QuadrantCorrect = QuadrantCorrect;
            % ScreenHeight
            taskStruct.SessionInfo.ScreenHeight = events{screenDlgInd}.screenHeightCM*2;
            % ScreenWidth
            taskStruct.SessionInfo.ScreenWidth = events{screenDlgInd }.screenWidthCM*2;
            % ScreenDistance
            taskStruct.SessionInfo.ScreenDistance = events{screenDlgInd }.screenDistanceCM;
            % FileName
            taskStruct.SessionInfo.FileName = ...
                [monkeyName(1) '_' ...
                session '_' ...
                taskName '_'...
                'Block' num2str(blockCounter) '.mat'];
        end
        % SessionInfo =====================================================
        
        
        
        %% Block info =====================================================
        % From the task dialog box
        % Task blocks are defined by the combination of contexts and
        % colours. If it is the first trial in a task block, compile all of
        % the task dialog parameters in taskStruct.SessionInfo.
        %
        % If irregularities in the task dialog box parameters are
        % detected, skip this trial.
        %
        % If this is not the first trial of the task block, and the context
        % and/or colors used have changed, save taskStruct before
        % proceeding to preprocess this trial, and reset the block counter
        % to 1.
        taskDlg = events{find(strcmp('savingController::taskDialog',eventNames))};
        
        % Retrieve active colors for each context and trial direction
        activeCtx1ClrsN = [taskDlg.chk_GoalCtx1NA taskDlg.chk_GoalCtx1NB taskDlg.chk_GoalCtx1NC taskDlg.chk_GoalCtx1ND] .* ...
            [taskDlg.pop_GoalCtx1NA taskDlg.pop_GoalCtx1NB taskDlg.pop_GoalCtx1NC taskDlg.pop_GoalCtx1ND];
        activeCtx1ClrsS = [taskDlg.chk_GoalCtx1SA taskDlg.chk_GoalCtx1SB taskDlg.chk_GoalCtx1SC taskDlg.chk_GoalCtx1SD] .* ...
            [taskDlg.pop_GoalCtx1SA taskDlg.pop_GoalCtx1SB taskDlg.pop_GoalCtx1SC taskDlg.pop_GoalCtx1SD];
        activeCtx2ClrsN = [taskDlg.chk_GoalCtx2NA taskDlg.chk_GoalCtx2NB taskDlg.chk_GoalCtx2NC taskDlg.chk_GoalCtx2ND] .* ...
            [taskDlg.pop_GoalCtx2NA taskDlg.pop_GoalCtx2NB taskDlg.pop_GoalCtx2NC taskDlg.pop_GoalCtx2ND];
        activeCtx2ClrsS = [taskDlg.chk_GoalCtx2SA taskDlg.chk_GoalCtx2SB taskDlg.chk_GoalCtx2SC taskDlg.chk_GoalCtx2SD] .* ...
            [taskDlg.pop_GoalCtx2SA taskDlg.pop_GoalCtx2SB taskDlg.pop_GoalCtx2SC taskDlg.pop_GoalCtx2SD];
        
        % Skip this trial if any of the task dialog parameters were not
        % correctly set (see PDF of behavioural task parameters for
        % details; contact RG or BC)
        if  any(...
                taskDlg.lst_TextureCtx1N ~= taskDlg.lst_TextureCtx1S || ... Ctx1 textures are equal
                taskDlg.lst_TextureCtx2N ~= taskDlg.lst_TextureCtx2S || ... Ctx2 textures are equal
                ~isequal(activeCtx1ClrsN,activeCtx1ClrsS) || ... Ctx 1 North and South colors are identical
                ~isequal(activeCtx2ClrsN,activeCtx2ClrsS) || ... Ctx 2 North and South colors are identical
                ~isequal(activeCtx1ClrsN(activeCtx1ClrsN~=0),fliplr(activeCtx2ClrsN(activeCtx2ClrsN~=0))) || ... Active colors are inverted across contexts
                (taskDlg.chk_AutoTraining && taskDlg.pop_trainingPhase>5) || ... Training phase is 5 or below (if active)
                taskDlg.chk_CuedMaze || ... Not active for this task
                taskDlg.chk_FreeRoam || ... Not active for this task
                ~taskDlg.chk_RandGoals || ... Goal locations should be randomized
                ~taskDlg.chk_CurrentPosStart || ... Trials are continuous (player not transported during ITI)
                taskDlg.chk_UseFogs || ... Not used for this task
                ~taskDlg.chk_ContTrials) %Trials are continuous (player not transported during ITI)
            warning(['Skipped trial ' trlNames{trl} ' due to nonsensical trial parameters (experimenter error)!'])
            InvalidTrials = [InvalidTrials; trlNames(trl)];
            continue
        end
        % Context Textures
        ctx1 = Task_TextureName{taskDlg.lst_TextureCtx1N};
        ctx2 = Task_TextureName{taskDlg.lst_TextureCtx2N};
        % Context Objects
        BlockInfo.([ctx1 'Colors']) = Task_ColorList(activeCtx1ClrsN(activeCtx1ClrsN~=0));
        BlockInfo.([ctx2 'Colors']) = Task_ColorList(activeCtx2ClrsN(activeCtx2ClrsN~=0));
        
        % Store or compare BlockInfo
        if blockTrlCounter == 1
            
            % If this is the first trial of the block, retain BlockInfo
            taskStruct.BlockInfo = BlockInfo;
            
        else
            if ~isequal(BlockInfo,taskStruct.BlockInfo) % Check to see if current trial dialog values have changed from the previous trial
                % If critical block information (contexts and objects) has
                % changed, save the taskStruct as a complete block, and start a
                % new taskStruct using these block parameters
               
                if ~isfield(taskStruct,'Trials')    
                    % If this is the first good trial of the block, retain BlockInfo
                    taskStruct.BlockInfo = BlockInfo;
            
                else
                    % Save taskStruct with old block info
                    ProcessEyesAndSave
                     
                    %Re-initialize taskStruct with this block info
                    taskStruct = struct();
                    taskStruct.BlockInfo = BlockInfo;
                    blockCounter = blockCounter + 1;
                    blockTrlCounter = 0;
                end
            end
        end
        
        %% Other task dialog parameters
        % TrainingPhase
        TrainingPhase = taskDlg.pop_trainingPhase * taskDlg.chk_AutoTraining;
        % BlockInfo =======================================================
        
        
        %% TrialInfo ======================================================
        eventsNumber = length(events);
        for k = 1:eventsNumber
            
            %% Retrieves data according to event type:
            switch events{k}.name
                
                %% PreStartTrialEvent
                case 'PreStartTrialEvent'
                    
                    % Trial context
                    tempCtx = str2num(events{k}.Textures(end));
                    Trl_Context  = Task_TextureName{...
                        tempCtx};   % Context output is a number, not a string
                    
                    % Trial direction
                    if events{k}.bIsNorth == '1'
                        Trl_Direction = 'North';
                    elseif events{k}.bIsNorth == '0'
                        Trl_Direction = 'South';
                    end
                    
                    % Goal Colors
                    % Allocentric
                    tempGoals        = regexp(events{k}.GoalsColor,'/','split'); % Always 'WESTcolor/EASTcolor'
                    Trl_GoalWest    = tempGoals{1};
                    Trl_GoalEast    = tempGoals{2};
                    % Egocentric
                    if strcmp(Trl_Direction, 'North')
                        Trl_GoalLeft  = Trl_GoalWest;
                        Trl_GoalRight = Trl_GoalEast;
                    elseif strcmp(Trl_Direction, 'South')
                        Trl_GoalLeft  = Trl_GoalEast;
                        Trl_GoalRight = Trl_GoalWest;
                    end
                    
                case 'UnrealPlayerState'
                    if ~isempty(events{k}.SampleTimes)
                        %Get Player State informations:
                        Unreal_States = [Unreal_States; regexp(events{k}.PlayerStates,'/','split')'];
                        
                        % Get Player Position and Rotation informations:
                        Unreal_PosXPosYRot = [Unreal_PosXPosYRot; events{k}.PlayerPositions(:,1:2) events{k}.PlayerRotations];
                        
                        % Get Unreal Time and convert it in MonkeyLab time to sync
                        Unreal_TempTime = regexp(events{k}.SampleTimes,'/','split')';
                        
                        % Convert UDK times from string to numerical values
                        % Times are now in vector: [YYYY MM DD hh mm ss.fff]
                        Unreal_TempTime = datevec(Unreal_TempTime, 'HH:MM:SS.FFF');
                        Unreal_TempTime = Unreal_TempTime(:,4:6);
                        
                        % Now since the two computers are synched to ~1ms
                        % precision we just need to compute the way to translate
                        % the "vector" times (i.e. [2014 09 23 16 28 10.444] or
                        % [YYYY MM DD HH MM SS.FFF] into "GetSecs" times (i.e.
                        % 4.555704997871990e+05). To do this we saved the
                        % variables : MLQueryTimeGetSecs and MLQueryTimeVect.
                        % This operation takes < 0.1 ms so we can assume both of
                        % these times to be equivalent. To get the UDK times in
                        % ML times we simply need to compute the difference in
                        % SECONDS between the UDK times and the MLQueryTimeVect
                        % and add this to the MLQueryTimeGetSecs.
                        
                        MLVectTime = repmat(events{k}.MLQueryTimeVect, size(Unreal_TempTime,1),1);
                        %Keep only the time (HH:MM:SS.FFF)
                        MLVectTime = MLVectTime(:,4:6);
                        
                        %Get the time difference in seconds
                        TimeDiff = Unreal_TempTime - MLVectTime;
                        TimeDiff = 3600*(TimeDiff(:,1)) + 60 * (TimeDiff(:,2)) + (TimeDiff(:,3));
                        
                        Unreal_TempTime = events{k}.MLQueryTimeGetSecs + TimeDiff;
                        
                        % Puts the values in the PosRotTime matrix plus
                        Unreal_Times = [Unreal_Times; Unreal_TempTime]; %#ok<*AGROW>
                    end
                    
                case 'eyeSample'
                    
                    % File should have a sampleTimeSecs field that time stamps
                    % the time at which the sample is taken and not when the
                    % event is sent (time). However some older files might not
                    % have it, so need to check.
                    % Eye_SamplesBad are timestamped by Monkeylab. It
                    % contains bits of data that are oversampled (eye
                    % offscreen) or undersampled (accumulated sub-ms
                    % differences between computer clocks). All of these
                    % times are corrected below.
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
                    
                case 'endOfTrial'
                    
                    % Get Trial outcome
                    Trl_OutcomeWord = events{k}.eot;
                    
                    % Get End of trial time
                    Trl_EOT = events{k}.time;
                    
                    %% startOfTrial
                    
                case 'startOfTrial'
                    
                    % Get Setup ID from trial ID;
                    % Trial ID  = DDD DSS SSS ssT TTT
                    %   D = Days since Jan 1 2012
                    %   S = Seconds of day
                    %   s = setup ID
                    %   T = Trial number
                    Trl_TrialID = events{k}.trialID;
                    Sys_SetupID = Trl_TrialID(10:11); % 01 = Cerebus, 02 = Plexon
                    
                    %Get Start of trial Time
                    Trl_SOT = events{k}.time;
                    
                    % Get Session date
                    Ses_Date = events{k}.date;
                    
                    %% SyncPulseRead
                    
                case 'SyncPulseRead'
                    
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
                    
                    %% ITC18TimeZero
                    
                case 'ITC18TimeZero'
                    % Get ITC 18 start collection time
                    Sync_ITC18StartTime = events{k}.time;
                    
                    %% sendTrialID
                    
                case 'sendTrialID'
                    EDFTrialStartTime = events{k}.time;
                    
                    %% startReward
                    
                case 'startReward'
                    Rew_Duration = events{k}.RewardDuration;
            end
        end
        
        %% Secondary trail parameter calculations
        % Check Trl_OutcomeWord
        if ~ismember(Trl_OutcomeWord,{'Correct' 'Incorrect'})
            
            InvalidTrials = [InvalidTrials; trlNames(trl)];

        else
            
            % Trl_GoalWestValue
            Trl_GoalWestValue = find(strcmp(...
                taskStruct.BlockInfo.([Trl_Context 'Colors']),...
                Trl_GoalWest)) == 1 ; % Remove == 1 to keep relative ranking rather than isHighest
            % Trl_GoalEastValue
            Trl_GoalEastValue = find(strcmp(...
                taskStruct.BlockInfo.([Trl_Context 'Colors']),...
                Trl_GoalEast)) == 1 ; % Remove == 1 to keep relative ranking rather than isHighest
            % Trl_GoalLeftValue
            Trl_GoalLeftValue = find(strcmp(...
                taskStruct.BlockInfo.([Trl_Context 'Colors']),...
                Trl_GoalLeft)) == 1 ; % Remove == 1 to keep relative ranking rather than isHighest
            % Trl_GoalRightValue
            Trl_GoalRightValue = find(strcmp(...
                taskStruct.BlockInfo.([Trl_Context 'Colors']),...
                Trl_GoalRight)) == 1 ; % Remove == 1 to keep relative ranking rather than isHighest
            
            % Trl_OutcomeColor, Trl_OutcomeDirec
            if ~isempty(strfind(Unreal_States{end},'Goal'))
                % West
                if strfind(Unreal_States{end},'W')
                    Trl_OutcomeDirec = 'West';
                    Trl_OutcomeColor = Trl_GoalWest;
                    % East
                elseif strfind(Unreal_States{end},'E')
                    Trl_OutcomeDirec = 'East';
                    Trl_OutcomeColor = Trl_GoalEast;
                    %Other
                else
                    error('Unrecognized goal location')
                end
            end
            % Trl_OutcomeSide
            if any(strcmp(Unreal_States{end}(end-1:end),{'NW','SE'}))
                Trl_OutcomeSide = 'Left';
            elseif any(strcmp(Unreal_States{end}(end-1:end),{'NE','SW'}))
                Trl_OutcomeSide = 'Right';
            else
                error('Unrecognized goal location')
            end
            
        end % TrlOutcome conditional
        
        % Correct abberant eye sampling====================================
        if ~isempty(Eye_EDFtimeDiff)
            edfDiff = Eye_EDFtimeDiff-Eye_EDFtimeDiff(1);
            
            % find the index in eyesamplesBad(:,3) that is greater than the
            % start of other trial, and the index in eyesamplesBad(:,3)
            % that is less that EOT.
            % edfDiff = EDFtimeDiff
            SOTIndex = find(Eye_SamplesBad(:,3)>EDFTrialStartTime,1);
            EOTindex = find(Eye_SamplesBad(:,3)<Trl_EOT,1,'last');
            edfDiff = Eye_EDFtimeDiff(SOTIndex:EOTindex)-Eye_EDFtimeDiff(SOTIndex);
            %          figure
            %          plot(edfDiff)
            
            % Find the time eyelink time that aligns with the start of the trial
            % (the end of the last trial) and the end of this trial. Get the
            % indices of these times, and then get the eye samples for these
            % periods
            Trl_SOTEDFindex = find(edfEyeSamplefile.FSAMPLE.time>...
                (EDFTrialStartTime-Eye_EDFtimeDiff(SOTIndex))*1000,1);
            Trl_EOTEDFindex = find(edfEyeSamplefile.FSAMPLE.time<...
                (Trl_EOT-Eye_EDFtimeDiff(SOTIndex))*1000,1,'last');
            EyeSamples = [double(edfEyeSamplefile.FSAMPLE.hx(1,Trl_SOTEDFindex:Trl_EOTEDFindex))' ...
                double(edfEyeSamplefile.FSAMPLE.hy(1,Trl_SOTEDFindex:Trl_EOTEDFindex))'...
                double(edfEyeSamplefile.FSAMPLE.time(Trl_SOTEDFindex:Trl_EOTEDFindex))'*.001+Eye_EDFtimeDiff(1)...
                double(edfEyeSamplefile.FSAMPLE.pa(1,Trl_SOTEDFindex:Trl_EOTEDFindex))'];
        end
        % =================================================================
        
        
        %% ContextOnset, GoalsOnset, DecisionOnset=========================
        % Get index of the moment the goals appear for this trial:
        
        % ContextOnset
        % Since the change in player state doesn't work, we will rely on
        % position and by taking into acount the collision radius of the
        % Pawn, which is: ~34 units
        if strcmp(Trl_Direction, 'North')
            if any(Unreal_PosXPosYRot(:,1) >= (-703-34))
                ContextOnset = ...
                    Unreal_Times(...
                    find(Unreal_PosXPosYRot(:,1) >= (-703-34), 1, 'first')...
                    );
            else
                ContextOnset = [];
            end
        elseif strcmp(Trl_Direction, 'South')
            if any(Unreal_PosXPosYRot(:,1) <= (704+48))
                ContextOnset = ...
                    Unreal_Times(...
                    find(Unreal_PosXPosYRot(:,1) <= (704+48), 1, 'first')...
                    );
            else
                ContextOnset = [];
            end
        end
        
        % GoalsOnset Time:
        goalsAppearInd = [];
        goalsAppearInd = find(strcmp(Unreal_States,...
            ['Goals' Trl_Direction 'On']),1, 'first');
        GoalsOnset = ...
            Unreal_Times(goalsAppearInd);
        
        
        % Determine at which index the monkey turned away from looking straight
        % North or straight south by thresholding the Unreal_PosRot value. A
        % difference in rotation of one degree in either direction from the
        % subject's rotation when the goals appeared is used as the threshold.
        goalsRotTime = [];
        if ~isempty(goalsAppearInd) && any(strcmp(Trl_OutcomeWord, {'Correct','Incorrect'}))
            
            %Find the first sample where the rotation crosses the threshold
            RotThreshInd = find(abs(...
                Unreal_PosXPosYRot(goalsAppearInd:end,3)-...
                Unreal_PosXPosYRot(goalsAppearInd,3)) > 10,...
                1, 'first')+goalsAppearInd;
            
            %Go back to the first sample where the rotation was initiated
            RotStartInd = find(...
                Unreal_PosXPosYRot(goalsAppearInd:RotThreshInd,3) == ...
                Unreal_PosXPosYRot(goalsAppearInd,3),1,'last')+goalsAppearInd;
            
            goalsRotTime = Unreal_Times(RotStartInd) ;
        end
        DecisionOnset = goalsRotTime;
        
        
        
        %% Create the output MLStruct variable=============================
        % Trial-specific information
        taskStruct.Trials.(['ID_' Trl_TrialID]).TrainingPhase  = TrainingPhase;      % (1)
        taskStruct.Trials.(['ID_' Trl_TrialID]).Direction      = Trl_Direction;      % (2)
        taskStruct.Trials.(['ID_' Trl_TrialID]).Context        = Trl_Context;        % (2)
        taskStruct.Trials.(['ID_' Trl_TrialID]).GoalWestColor  = Trl_GoalWest;       % (2)
        taskStruct.Trials.(['ID_' Trl_TrialID]).GoalWestValue  = Trl_GoalWestValue;  % (2)
        taskStruct.Trials.(['ID_' Trl_TrialID]).GoalEastColor  = Trl_GoalEast;       % (2)
        taskStruct.Trials.(['ID_' Trl_TrialID]).GoalEastValue  = Trl_GoalEastValue;  % (2)
        taskStruct.Trials.(['ID_' Trl_TrialID]).GoalLeft       = Trl_GoalLeft;       % (2)
        taskStruct.Trials.(['ID_' Trl_TrialID]).GoalLeftValue  = Trl_GoalLeftValue;  % (2)
        taskStruct.Trials.(['ID_' Trl_TrialID]).GoalRight      = Trl_GoalRight;      % (2)
        taskStruct.Trials.(['ID_' Trl_TrialID]).GoalRightValue = Trl_GoalRightValue; % (2)
        taskStruct.Trials.(['ID_' Trl_TrialID]).OutcomeColor   = Trl_OutcomeColor;   % (3)
        taskStruct.Trials.(['ID_' Trl_TrialID]).OutcomeDirec   = Trl_OutcomeDirec;   % (3)
        taskStruct.Trials.(['ID_' Trl_TrialID]).OutcomeSide    = Trl_OutcomeSide;    % (3)
        taskStruct.Trials.(['ID_' Trl_TrialID]).OutcomeWord    = Trl_OutcomeWord;     % (7)
        taskStruct.Trials.(['ID_' Trl_TrialID]).RewardDuration = Rew_Duration;        % (6)
        taskStruct.Trials.(['ID_' Trl_TrialID]).Unreal_States  = Unreal_States;       % (3)
        taskStruct.Trials.(['ID_' Trl_TrialID]).Unreal_PosXPosYRot = Unreal_PosXPosYRot;  % (3)
        taskStruct.Trials.(['ID_' Trl_TrialID]).Unreal_Times   = Unreal_Times;        % (3)
        taskStruct.Trials.(['ID_' Trl_TrialID]).EyeSamples     = EyeSamples;         % (3)
        taskStruct.Trials.(['ID_' Trl_TrialID]).ITC18StartTime = Sync_ITC18StartTime; % (10)
        taskStruct.Trials.(['ID_' Trl_TrialID]).SOT_Time       = Trl_SOT;             % (6)
        taskStruct.Trials.(['ID_' Trl_TrialID]).ContextOnset   = ContextOnset;        % (6)
        taskStruct.Trials.(['ID_' Trl_TrialID]).GoalsOnset     = GoalsOnset;          % (6)
        taskStruct.Trials.(['ID_' Trl_TrialID]).DecisionOnset  = DecisionOnset;       % (6)
        taskStruct.Trials.(['ID_' Trl_TrialID]).EOT_Time       = Trl_EOT;             % (7)
        taskStruct.Trials.(['ID_' Trl_TrialID]).SyncPulse      = Sync_Pulse;          %
        % =================================================================
        
        
        blockTrlCounter = blockTrlCounter + 1;
        waitbar(trl/numel(trlNames));
        display(['     Trial ' Trl_TrialID(end-2:end) ...       %**
            ' processing time : ' num2str(toc(timerTrl), '%3.3f') ' seconds']);
        
    end % Trial loop
    
    
    
    
    
    close(hWaitBar);
    % TrialInfo============================================================
    
    
    
    if isfield(taskStruct,'Trials')
        
        % Save taskStruct with old block info
        ProcessEyesAndSave
        
    end
    
    
    
    
end % task loop




%% Nested function
    function ProcessEyesAndSave
        
        %%Cleaning and smoothing eye data======================================
        timerFormatEye = tic;
        display('Cleaning and smoothing eye data');
        % Modify this to handle Eye_SamplesBad as an input rather than looking
        % for it within taskStruct
        taskStruct = HcTask_FormatEyeDataOffScreen(taskStruct, 1);
        
        % Classify eye movements
        hWaitBar3 = waitbar(0, 'Classifying eye movements.');
        trlIDs = fieldnames(taskStruct.Trials);
        for t = 1:length(trlIDs)
            
            if isfield(taskStruct.Trials.(trlIDs{t}),'PositionMatrix_deg')
                
                taskStruct.Trials.(trlIDs{t}).SaccadeData = ...
                    HcTask_SaccadeProcessing(taskStruct.Trials.(trlIDs{t}).PositionMatrix_deg(:,4:6),1);
                waitbar(t/length(trlIDs));
                
            else
                continue
            end
            
        end
        close(hWaitBar3);
        
        toc(timerFormatEye);
        %======================================================================
        
        %% Save taskStruct=====================================================
        if bSave
            display('Saving taskStruct')
            timerSave = tic;
            
            taskStruct.BlockInfo.InvalidTrials = InvalidTrials;
            
            dirResultsSession = strrep(dirDataSess,'Data','Results');
            if ~isdir(dirResultsSession)
                mkdir(dirResultsSession)
            end
            
            dirSave = [dirResultsSession 'BehavData' filesep];
            if ~isdir(dirSave)
                mkdir(dirSave)
            end
            
            %Update FileName Info in structure
            save([dirSave taskStruct.SessionInfo.FileName], 'taskStruct');
            
            clear taskStruct
            toc(timerSave);
        end
        % Save taskStruct======================================================
        
    end

end % function end