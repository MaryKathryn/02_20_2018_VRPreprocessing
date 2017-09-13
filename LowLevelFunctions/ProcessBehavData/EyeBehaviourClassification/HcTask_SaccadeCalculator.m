%% UltimateSaccadeCalculator
% open up a results folder - go through each of the XMaze, EyeCal, EyeEnd, and FreeRoam
% create a structure 'SaccadeData' that contains SessionInfo, and both
% smoothed and unsmoothed Saccade structures organized by TrialID (for all
% session types. Also generate a structure that holds the Offscreen
% saccades offset and Onset. These are not differentiated from blinks or
% loss of signal, as it cannot be guatanteed that there is not a saccade
% during this period.
%% Saccade Structures include:
% SaccadeStartTime,
% SaccadeEndTime,
% SaccadePeakTime,
% SaccadePeakVelocity,
% SaccadeDuration,
% SaccadeAmplitude,
% SaccadeStartPointX,
% SaccadeEndPointX, 
% SaccadeStartPointY, 
% SaccadeEndPointY,
% SaccadeDirection
%
% This program will ask for the desired Directory, and then load all
% available session types, and then process the saccades using
% HcTask_SaccadeProcessing. This calculates the acceleration threshhold for
% each trial, and then processes each putative saccade based on direction
% of movement and velocity of eye trace.
% One file is created for each folder, and that can be processed for any of
% the available trials, regardless of session, to analyze the correct
% trials, simply run HcTask_SelectTrials, which will remove unneeded
% trails.
%% Get folders to run analysis
if ispc
    Slash = '\';
else
    Slash = '/';
end
% ParentDirectory = mfilename('fullpath');
% Select directories
[SelectedSessions,ResultsDirectory] = HcTask_chooseDirectories();

%% Cycle through the sessions
for sessionNum = 1:length(SelectedSessions)
    %% Set up files
    display(['Calculating SaccadeInfo for ' SelectedSessions{sessionNum}]);
    sessionTimer = tic;
    sessionDir = [ResultsDirectory SelectedSessions{sessionNum}];
    dirFiles = dir(sessionDir);
    %get the files in the folder
    dirFiles = {dirFiles.name};
    %Cycle through the SessionTypes - fitrst load EyeCal
    EyeCalFile = cell2mat(dirFiles(~cellfun(@isempty,(strfind(dirFiles,'EyeCal')))));
    EyeEndFile = cell2mat(dirFiles(~cellfun(@isempty,(strfind(dirFiles,'EyeEnd')))));
    XMazeFile = cell2mat(dirFiles(~cellfun(@isempty,(strfind(dirFiles,'XMaze')))));
    FreeRoamFile = cell2mat(dirFiles(~cellfun(@isempty,(strfind(dirFiles,'Free')))));
    SaccadeDataFile = cell2mat(dirFiles(~cellfun(@isempty,(strfind(dirFiles,'SaccadeData')))));

    SessInfo = 0;
    
    
    cd(sessionDir)
    
    if ~isempty(XMazeFile)
        load([sessionDir Slash XMazeFile]);
%         XMazeStruct = HcTask_FormatEyeDataOnScreen(XMazeStruct,1);
%         XMazeStruct = HcTask_FormatEyeDataOffScreen(XMazeStruct,1);
%         save(XMazeFile,'XMazeStruct')
        SaccadeData.SessionInfo = XMazeStruct.SessionInfo;
        SessInfo = 1;
    else display('no XMazeData');XMazeStruct = [];
    end
    if ~isempty(EyeEndFile)
        load([sessionDir Slash EyeEndFile]);
%         EyeEndStruct = HcTask_FormatEyeDataOnScreen(EyeEndStruct,0);
%         EyeEndStruct = HcTask_FormatEyeDataOffScreen(EyeEndStruct,0);
%         save(EyeEndFile,'EyeEndStruct')
        if SessInfo~= 1
            SaccadeData.SessionInfo =EyeEndStruct.SessionInfo;
            SessInfo = 1;
        end
    else display('no EyeEnd data'); EyeEndStruct = [];
    end
    if ~isempty(EyeCalFile);
        load([sessionDir Slash EyeCalFile]);
%         EyeCalStruct = HcTask_FormatEyeDataOnScreen(EyeCalStruct,0);
%         EyeCalStruct = HcTask_FormatEyeDataOffScreen(EyeCalStruct,0);
%         save(EyeCalFile,'EyeCalStruct')
        if SessInfo~= 1
            SaccadeData.SessionInfo =EyeCalStruct.SessionInfo;
            SessInfo = 1;
        end
    else display('no EyeCal data'); EyeCalStruct = [];
    end
    if ~isempty(FreeRoamFile)
        load([sessionDir Slash FreeRoamFile]);
%         FreeRoamStruct = HcTask_FormatEyeDataOnScreen(FreeRoamStruct,1);
%         FreeRoamStruct = HcTask_FormatEyeDataOffScreen(FreeRoamStruct,1);
%         save(FreeRoamFile,'FreeRoamStruct')
        if SessInfo~= 1
            SaccadeData.SessionInfo =FreeRoamStruct.SessionInfo;
            SessInfo = 1;
        end
    else display('no FreeRoam data'); FreeRoamStruct = [];
    end
    if ~isempty(SaccadeDataFile)
        load([sessionDir Slash SaccadeDataFile]);
    end
    Types = {'EyeCalStruct'; 'EyeEndStruct';'XMazeStruct';'FreeRoamStruct'};
    BigStruct = struct('XMazeStruct',XMazeStruct,'FreeRoamStruct',FreeRoamStruct...
        ,'EyeCalStruct',EyeCalStruct,'EyeEndStruct',EyeEndStruct); 
        SaccadeData = [];
    
    %% Go Through Files
    %cycle through the loaded files, and compute the saccade data
    if ~isempty(EyeCalStruct) && ~isempty(EyeEndStruct) && ~isempty(XMazeStruct) &&...
            ~isempty(FreeRoamStruct)
        for sessionType = 1:4
            
            if ~isempty(BigStruct.(Types{sessionType}))
                display(['Calculating SaccadeInfo for ' Types{sessionType}]);
                
                %for each trial, compute the
                trialIDs = fieldnames(BigStruct.(Types{sessionType}).Trials);
                for trl = 1:length(trialIDs)
                    % [Types{sessionType} trl]
                    
                    if strncmp(Types{sessionType},'Eye',3)
                        if ~isempty(BigStruct.(Types{sessionType}).Trials.(trialIDs{trl}).EyeSamples)
                            % Get the saccade information for Onscreen
                            % Saccades, using the nonsmoothed info, because it
                            % is better at ending
%                             SaccadeData.NonSmoothTrials.(trialIDs{trl}) = ...
%                                 HcTask_SaccadeProcessing(BigStruct.(Types{sessionType})...
%                                 .Trials.(trialIDs{trl}).NonSmoothedEyes,0);
                            SaccadeData.SmoothTrials.(trialIDs{trl}) = ...
                                HcTask_SaccadeProcessing(BigStruct.(Types{sessionType})...
                                .Trials.(trialIDs{trl}).EyeDegreesOnScreen,1);
                            %the Offscreen saccades
%                             SaccadeData.OffScreen.(trialIDs{trl}) = HcTask_SaccadeProcessingOffScreen(...
%                                 BigStruct.(Types{sessionType}).Trials.(trialIDs{trl}).NonSmoothedEyes,0);
                        end
                        
                    else
                        if ~isempty(BigStruct.(Types{sessionType}).Trials.(trialIDs{trl}).NonSmoothedEyes)
                            
%                             SaccadeData.NonSmoothTrials.(trialIDs{trl}) = HcTask_SaccadeProcessing(BigStruct.(Types{sessionType})...
%                                 .Trials.(trialIDs{trl}).NonSmoothedEyes,0);
                            SaccadeData.SmoothTrials.(trialIDs{trl}) = ...
                                HcTask_SaccadeProcessing(BigStruct.(Types{sessionType})...
                                .Trials.(trialIDs{trl}).PositionMatrix_degOnScreen(:,4:6),1);
                            %the Offscreen saccades
                            SaccadeData.OffScreen.(trialIDs{trl}) = HcTask_SaccadeProcessingOffScreen(...
                                BigStruct.(Types{sessionType}).Trials.(trialIDs{trl}).NonSmoothedEyes,0);
                           
                        end
                    end
                end
            end
        end
    end
%% run analyses on the movements
    if ~isempty(EyeCalStruct) && ~isempty(EyeEndStruct) && ~isempty(XMazeStruct) &&...
            ~isempty(FreeRoamStruct)
        %run on the non smoothed data
%         SaccadeData = HcTask_EyeCalFixationCalculator(EyeCalStruct,EyeEndStruct,SaccadeData,1);
%         SaccadeData = HcTask_XMazeFixationCalculator(XMazeStruct,SaccadeData,1);
%         SaccadeData = HcTask_FreeRoamFixationCalculator(FreeRoamStruct,SaccadeData,1);
%         %run on the smoothed data
        SaccadeData = HcTask_EyeCalFixationCalculator(EyeCalStruct,EyeEndStruct,SaccadeData,2);
        SaccadeData = HcTask_XMazeFixationCalculator(XMazeStruct,SaccadeData,2);
        SaccadeData = HcTask_FreeRoamFixationCalculator(FreeRoamStruct,SaccadeData,2);
%        % MS on nonsmoothed data
%         SaccadeData.MainSequenceDegrees = HcTask_MainSequencePlot(XMazeStruct,...
%             FreeRoamStruct,EyeCalStruct,EyeEndStruct,SaccadeData);
%         
% MS on smoothed Data

        SaccadeData.MainSequenceDegreesSmooth = HcTask_MainSequencePlotSmooth(XMazeStruct,...
            FreeRoamStruct,EyeCalStruct,EyeEndStruct,SaccadeData);
        SaccadeData.SessionInfo.FullSession = 1;
    else
        SaccadeData.SessionInfo.FullSession = 0;
    end
      save([SelectedSessions{sessionNum} '_SaccadeData'],'SaccadeData')
% WArndingDidNotSAve = 1
     SaccadeData = [];
end
