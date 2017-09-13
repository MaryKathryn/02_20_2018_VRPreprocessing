function [] = HcTask_appendEyeMovements()
% 
% Add 'SaccadeData' structure to each trial in all tasks. 
%
% Saccade Structures include:
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
% One structure is added to each trial, and that can be processed for any of
% the available trials, regardless of session, to analyze the correct
% trials, simply run HcTask_SelectTrials, which will remove unneeded
% trails.
% 
% Usage note:
% External functions:
% Mean Angle
%   https://www.mathworks.com/matlabcentral/fileexchange/12631-mean-angle/content/meanangle.m


%% Get folders to run analysis
sl = filesep;
monkeyNames = {'Woody' 'Raul'};
taskTypes = {'EyeCal' 'XMaze' 'FreeRoam' 'EyeEnd'};
dirAnalyses = mfilename('fullpath');
dirTask = dirAnalyses(1:max(strfind(dirAnalyses,[sl 'Analyses'])));
dirRes = [dirTask 'Results' sl];


for mky = 1:length(monkeyNames)
    
    monkeyName = monkeyNames{mky};
    % Return all sessList with a results folder for this monkey
    dirMnkyRes = [dirRes monkeyName sl];
    sessList = cellstr(ls([dirMnkyRes '20*']));
    
    %% Cycle through the sessList
    for sessNum = 1:length(sessList)
        
        
        %% Set up files
        timerSess = tic;
        
        % Session date
        session = sessList{sessNum};
        % Session ID
        sessID = [monkeyName(1) session];
        % Session specific results folder
        dirSess = [dirMnkyRes session sl];
        display([char(10) 'Calculating SaccadeInfo for ' sessID]);
        
        
        % Load behavioural data
        EyeCalStruct = struct();
        XMazeStruct = struct();
        FreeRoamStruct = struct();
        EyeEndStruct = struct();
        
        display (['Loading behavioural data for ' sessID])
        timerLoad = tic;
        eyecalFile = [dirSess ls([dirSess '*EyeCal.mat'])];
        load(eyecalFile)
        xmazeFile = [dirSess ls([dirSess '*XMaze.mat'])];
        load(xmazeFile)
        freeroamFile = [dirSess ls([dirSess '*FreeRoam.mat'])];
        load(freeroamFile)
        eyeendFile = [dirSess ls([dirSess '*EyeEnd.mat'])];
        load(eyeendFile)
        toc(timerLoad)
        
        
        for task = 1:length(taskTypes)
            timerTask = tic;
            taskName = taskTypes{task};
            display(['Calculating SaccadeInfo for ' taskName]);
            trlIDs = [];
            
            switch taskName
                case 'EyeCal'
                    if ~isempty(EyeCalStruct)
                        trlIDs = fieldnames(EyeCalStruct.Trials);
                        for trl = 1:length(trlIDs)
                            EyeCalStruct.Trials.(trlIDs{trl}).SaccadeData = ...
                                HcTask_SaccadeProcessing(EyeCalStruct.Trials.(trlIDs{trl}).EyeDegrees,1);
                        end
                        save(eyecalFile,'EyeCalStruct');
                        display([taskName ' : ' num2str(toc(timerTask)) 's'])
                    end
                    
                case 'XMaze'
                    if ~isempty(XMazeStruct)
                        trlIDs = fieldnames(XMazeStruct.Trials);
                        for trl = 1:length(trlIDs)
                            XMazeStruct.Trials.(trlIDs{trl}).SaccadeData = ...
                                HcTask_SaccadeProcessing(XMazeStruct.Trials.(trlIDs{trl}).PositionMatrix_deg(:,4:6),1);
                        end
                        save(xmazeFile,'XMazeStruct');
                        display([taskName ' : ' num2str(toc(timerTask)) 's'])
                    end
                    
                case 'FreeRoam'
                    if ~isempty(FreeRoamStruct)
                        trlIDs = fieldnames(FreeRoamStruct.Trials);
                        for trl = 1:length(trlIDs)
                            if isfield(FreeRoamStruct.Trials.(trlIDs{trl}),'PositionMatrix_deg')
                                FreeRoamStruct.Trials.(trlIDs{trl}).SaccadeData = ...
                                    HcTask_SaccadeProcessing(FreeRoamStruct.Trials.(trlIDs{trl}).PositionMatrix_deg(:,4:6),1);
                            else
                                FreeRoamStruct.Trials = rmfield(FreeRoamStruct.Trials,trlIDs{trl});
                            end
                            % Note: Had to check for PositionMatrix_deg
                            % field, because there is a trial in W20141216
                            % (number 54/57) that is only 8ms long. It is
                            % labelled as correct because the trial started
                            % with the monkey in the red fog (Goal 60)
                        end
                        save(freeroamFile,'FreeRoamStruct');
                        display([taskName ' : ' num2str(toc(timerTask)) 's'])
                    end
                    
                case 'EyeEnd'
                    if ~isempty(EyeEndStruct)
                        trlIDs = fieldnames(EyeEndStruct.Trials);
                        for trl = 1:length(trlIDs)
                            EyeEndStruct.Trials.(trlIDs{trl}).SaccadeData = ...
                                HcTask_SaccadeProcessing(EyeEndStruct.Trials.(trlIDs{trl}).EyeDegrees,1);
                        end
                        save(eyeendFile,'EyeEndStruct');
                        display([taskName ' : ' num2str(toc(timerTask)) 's'])
                    end
                    
            end
            toc(timerSess)
            
            
        end % task loop
        
    end % session loop
    
end % monkey loop

end % function
