function [NeuralTrials] = HcTask_ProcessSyncPulse(allTasksStruct,LFPFile)

sl = filesep;
pathParts = strsplit(LFPFile,filesep);
sessInd = find(~cellfun(@isempty,strfind(pathParts,'2017')));
monkeyName = pathParts{sessInd-1};
session = pathParts{sessInd};
dirDataSess    = cell2mat(strcat(pathParts(1:end-1),sl));
dirResultsSess = regexprep(dirDataSess,'Data','Results');
dirSave = [dirResultsSess 'NeuralData' sl];

MLtoCer_CommunicationTime = 0.008; % Seconds. This is stable across tasks, monkeys, and days. 


%%% OPEN NS6 FILE
[~, hfile] = ns_OpenFile(LFPFile);
% Get file information
[~, FileInfo] = ns_GetFileInfo(hfile); % Gives you EntityCount, TimeStampResolution and TimeSpan and the time at which the file was created
NeuralTrials.SessionInfo = FileInfo;
NeuralTrials.SessionInfo.MonkeyName = monkeyName;
NeuralTrials.SessionInfo.FileName = [monkeyName(1) '_' session '_NeuralTrials.mat'];


%%% RETRIEVE EVENTS
timerNTWords = tic;
display('Retrieving event data (words)');

% Retreive entity (channel) information
[~, LFPEntityInfo] = ns_GetEntityInfo(hfile, 1:FileInfo.EntityCount);
LFPAnalogList = find([LFPEntityInfo.EntityType] == 2);

% Find the Entity index of the digital words sent from Monkeylab to the
% recording system. On the Plexon, this channel is named 'Strobed'. On the
% Cerebus, it is named 'digin'.
evInd = find(strcmp('digin',{LFPEntityInfo.EntityLabel}), 1);
if isempty(evInd)
    evInd = find(strcmp('Strobed',{LFPEntityInfo.EntityLabel}));
end

if isempty(evInd)
    errordlg('Could not find analog channel containing the events!');
end

% Get basic information from the event channel
[~, nsEventInfo] = ns_GetEventInfo(hfile, evInd);

% Retrieve event words and associated timestamps
[~, evTs, evWords] = ...
    ns_GetEventData(hfile, evInd, 1:LFPEntityInfo(evInd).ItemCount);

NeuralTrials.Events.Info = nsEventInfo;
NeuralTrials.Events.Words = evWords;
NeuralTrials.Events.Timestamps = evTs;

display(['Time retrieving event data: ' num2str(toc(timerNTWords))]);
% display(char(10));



%% Start SyncPulse: Alignment between Monkeylab and Neural recording comps
% The following sections of the script are only relevant for data recorded
% in Monkeylab, 2014 and later. Previous versions of Monkeylab and McLab
% did not incorporate a Sync Pulse to measure lags between the two
% computers.
% The following following Sync Pulse sections are set up to read the
% Strobed analog data sent from Monkeylab to the analog input 9 on the
% Cerebus computer.
% NOTE: Sync pulse channel will be named differently on the Plexon system.
% Named ainp9 (analog input 9) on the Cerebus system.




% Find the Sync Pulse channel, name ainp9 on Cerebus.
if strcmp(LFPFile(end-2:end),'ns6')
    SyncPulseInd = ...
        LFPAnalogList( find ( ...
        strcmp({LFPEntityInfo(LFPAnalogList).EntityLabel}, 'ainp9')) ); %#ok<*FNDSB>
end
if isempty(SyncPulseInd)
    errordlg('Sync Pulse Channel Not Found!');
    return
end

% Retreive Sync Pulse data and timestamps
display('Retrieving Sync Pulse from the neural recording file.');
timerReadNS6 = tic;

[~, ~, tempPulse] = ...
    ns_GetAnalogData(hfile,SyncPulseInd,1,LFPEntityInfo(SyncPulseInd).ItemCount);
% TimeStamps
[~, tempTimeStamp_Pulse] = ...
    ns_GetTimeByIndex(hfile,SyncPulseInd,1:LFPEntityInfo(SyncPulseInd).ItemCount);


% Gets all Trial defining words
NT_EOTInds = find(NeuralTrials.Events.Words ==1017);
NT_SOTInds = find(NeuralTrials.Events.Words==1016);
NT_TrlIDStartInds = find(NeuralTrials.Events.Words==1018);
NT_TrlIDEndInds = find(NeuralTrials.Events.Words==1019);
% Get timestamps associated with these words
NT_EOTTs = NeuralTrials.Events.Timestamps(NT_EOTInds);
NT_SOTTs = NeuralTrials.Events.Timestamps(NT_SOTInds);

% Starts with all of the Cerebus EOTInds, then finds the closest previous
% SOTInd, then the closest previous TrialIDEndInd, then the closest
% previous TrialIDStartInd
% for k = 1:length(NT_EOTInds)
%     
%     NT_ValidTrlIDStartInds(k)   = NT_TrlIDStartInds(...
%         find(NT_TrlIDStartInds < NT_EOTInds(k),1,'last') );
%     
%     NT_ValidTrlIDEndInds(k)     = NT_TrlIDEndInds(...
%         find(NT_TrlIDEndInds > NT_ValidTrlIDStartInds(k),1,'first') );
%     
%     NT_ValidTrlSOTInds(k)       = NT_SOTInds(...
%         find(NT_SOTInds > NT_ValidTrlIDEndInds(k),1,'first') );
%     
% %     NT_ValidTrlEOTInds(k)       = NT_EOTInds(...
% %         find(NT_EOTInds > NT_SOTInds(k),1,'first') );
%     
%     NT_TrialID{k,1}             = ['ID_' num2str(NeuralTrials.Events.Words (...
%         NT_ValidTrlIDStartInds(k)+1:NT_ValidTrlIDEndInds(k)-1)', '%03d')];
%     
%     NeuralTrials.Trials.(NT_TrialID{k,1}).SOT_Time = NeuralTrials.Events.Timestamps(...
%                                                         NT_ValidTrlSOTInds(k) );
% end
% NeuralTrials2.Events.Info = nsEventInfo;
% NeuralTrials2.Events.Words = evWords;
% NeuralTrials2.Events.Timestamps = evTs;
% Find TrialIDs that are not included in both NT and taskStruct




% Gets all TrialIDs and creates data containing structure
% Define variables
NeuralTrials.Trials = struct();

for k=1:length(NT_EOTInds)
    
    
    % Retrieve the Trial ID
    % Word 1018 marks the start of the trial ID. 1019 marks the end of
    % the trial ID.
    tempStartTrialID = ...
        find(NeuralTrials.Events.Words(1:NT_EOTInds(k))==1018);
    tempEndTrialID = ...
        find(NeuralTrials.Events.Words(1:NT_EOTInds(k))==1019);
    % If the recording is started or stopped in the middle of a trial, one
    % of the above variables will be empty; in which case, skip this trial,
    % and move on to the next iteration of the loop.
    if isempty(tempStartTrialID) || isempty(tempEndTrialID)
        continue
    end
    
    TrialID = num2str(NeuralTrials.Events.Words (...
        tempStartTrialID(end)+1:tempEndTrialID(end)-1)', '%03d');
    FieldNameID = ['ID_' TrialID];
    
    %Gets Pulse Sequence for this trial
    PulseIndices = find(tempTimeStamp_Pulse < NT_EOTTs(k));
    
    PulseSequence = tempPulse(PulseIndices);
    TimeStamp_Sequence = tempTimeStamp_Pulse(PulseIndices);
    
    %Removes used points from temp Variables
    tempPulse(PulseIndices)=[];
    tempTimeStamp_Pulse(PulseIndices)=[];
    
    %Adds sequence to structure
    NeuralTrials.Trials.(FieldNameID).Pulse = PulseSequence ;
    NeuralTrials.Trials.(FieldNameID).TimeStamp_Pulse = TimeStamp_Sequence;

    %Gets the StartTrial word Time
    NeuralTrials.Trials.(FieldNameID).Raw_SOT_Time = max( NT_SOTTs((NT_SOTTs < NT_EOTTs(k))) );
    
    
    
end


%% Remove all trials that are not common to taskStruct and NeuralTrials
AllMLTrials = fieldnames(allTasksStruct.Trials);
AllNeuralTrials = fieldnames(NeuralTrials.Trials);
% Remove NeuralTrials that are not included in any of the allTasksStruct
NeuralTrials.Trials = rmfield(NeuralTrials.Trials,...
                              AllNeuralTrials(...
                                find(~ismember(AllNeuralTrials,AllMLTrials))...
                                ) );
AllNeuralTrials = fieldnames(NeuralTrials.Trials);
% Remove allTasksStruct trials that are not included NeuralTrials
allTasksStruct.Trials = rmfield(allTasksStruct.Trials,...
                              AllNeuralTrials(...
                                find(~ismember(AllMLTrials,AllNeuralTrials))...
                                ) );
AllMLTrials = fieldnames(allTasksStruct.Trials);
% Check
if ~isequal(fieldnames(allTasksStruct.Trials),fieldnames(NeuralTrials.Trials))
    error
end


%% Fit the linear equation between ML SOT and Cerebus SOT
ML_SOT = structfun(@(x)(x.SOT_Time),allTasksStruct.Trials);
Cr_SOT = structfun(@(x)(x.Raw_SOT_Time-MLtoCer_CommunicationTime),NeuralTrials.Trials);
% plot(ML_SOT,(ML_SOT-Cr_SOT))
Time_Diff_eqn = polyfit(ML_SOT, (ML_SOT-Cr_SOT), 1);



% Compute TimeDiff and Cerebus SOT Time
for trl = 1:length(AllNeuralTrials)
    
    % Trial ID
    trlID = AllNeuralTrials{trl};
    
    % Find the time difference at this time point
    NeuralTrials.Trials.(trlID).TimeDiff = ...
        (Time_Diff_eqn(1) * allTasksStruct.Trials.(trlID).SOT_Time) + ...
                      Time_Diff_eqn(2);
    
    % Find the Cerebus Start of trial time at this trials corresponding
    % Monkeylab SOT
    NeuralTrials.Trials.(trlID).SOT_Time = ...
        NeuralTrials.Trials.(trlID).Raw_SOT_Time + ...
        NeuralTrials.Trials.(trlID).TimeDiff;
    
    
end

 
% %% Find lag between the start times from each trial 
% % Start by trying to use the SyncPulse. For any trials with a poor
% % correlation between SyncPulse from Cerebus and ML, use the mean
% % difference in times for the SOT and EOT words. 
% 
% % Remove NeuralTrials that are not included in any of the allTasksStruct
% NeuralTrials.Trials = rmfield(NeuralTrials.Trials,...
%                               AllNeuralTrials(...
%                                 find(~ismember(AllNeuralTrials,AllMLTrials))...
%                                 ) );
% AllNeuralTrials = fieldnames(NeuralTrials.Trials);
% 
% 
% trlwaitbar = waitbar(0,'Computing Lag time for all trials');
% display('Computing lag times for all trials');
% lagTimer = tic;
% 
% for k=1:length(AllNeuralTrials)
%     waitbar(k/length(AllNeuralTrials),trlwaitbar);
%     
%     CurrentTrial = AllNeuralTrials{k};
%     
%     if isfield(allTasksStruct.Trials.(CurrentTrial), 'SyncPulse')
%         if ~isempty(allTasksStruct.Trials.(CurrentTrial).SyncPulse)
%             tempCorr = allTasksStruct.Trials.(CurrentTrial).SyncPulse;
%             MLTimeStampsIndices = find(allTasksStruct.Trials.(CurrentTrial).SyncPulse(:,2)~=0);
%             MLTimes = allTasksStruct.Trials.(CurrentTrial).SyncPulse(:,2);
%             ML_SOT = allTasksStruct.Trials.(CurrentTrial).SOT_Time;
%             ML_EOT = allTasksStruct.Trials.(CurrentTrial).EOT_Time;
%         else
%             tempCorr = NaN;
%             MLTimeStampsIndices = NaN;
%             MLTimes = NaN;
%             ML_SOT  = NaN;
%             ML_EOT  = NaN;
%         end
%     else
%         tempCorr = NaN;
%         MLTimeStampsIndices = NaN;
%         MLTimes = NaN;
%         ML_SOT  = NaN;
%         ML_EOT  = NaN;
%     end
%         
% 
%     
%     if ~isnan(tempCorr)
%         tempXcorr = xcorr(tempCorr(:,1), NeuralTrials.Trials.(CurrentTrial).Pulse);
%         MaxIndex = find(tempXcorr==max(tempXcorr));
%         
%         %The lag is defined as the number of samples difference between ML
%         %and Cerebus:
%         % If the lag < 0: Cerebus is earlier than ML
%         % If lag > 0 : ML earlier
%         % The SyncPuse channel is sampled at 1kHz, so the lag time, though
%         % derived from differences in sample number, is equivalent to the lag
%         % time in milliseconds.
%         % To convert Cerebus time to ML time:
%         % Time of ML sample(x) = Time of Cerebus sample(x + lag)
%         SampleLag = length(NeuralTrials.Trials.(CurrentTrial).Pulse) - MaxIndex;
%         
%         
%         %Get all timestamps for ML Sync Pulse and round to the msec
%         CerebTimeStampIndices = MLTimeStampsIndices + SampleLag;
%         
%         %Removes indices that are over the number of samples or smaller than 0
%         IndicesToRemove = CerebTimeStampIndices > numel(NeuralTrials.Trials.(CurrentTrial).Pulse);
%         
%         CerebTimeStampIndices(IndicesToRemove) = [];
%         MLTimeStampsIndices(IndicesToRemove) = [];
%         
%         IndicesToRemove = CerebTimeStampIndices <= 0;
%         
%         CerebTimeStampIndices(IndicesToRemove) = [];
%         MLTimeStampsIndices(IndicesToRemove) = [];
%         
%         %Get time stamps for both datasets
%         
%         CerebTimes = NeuralTrials.Trials.(CurrentTrial).TimeStamp_Pulse(CerebTimeStampIndices);
%         MLTimes = MLTimes(MLTimeStampsIndices);
%         
%         %This is the time difference between ML and Cerebus. SO to get the ML
%         %time from the Cerebus time : MLTime = CerebTime + TimeDiff;
%         NeuralTrials.Trials.(CurrentTrial).TimeDiff = mean(MLTimes - CerebTimes);
%         
%         % Will save the lag in Time difference and sample difference from
%         % now on.
%         NeuralTrials.Trials.(CurrentTrial).Lag = SampleLag;
%         
%         % Save Start and End of Trial times for each trial from the recording file
% %         MUST DEFINE ML_SOT AND ML_EOT ABOVE SO THAT THEY CAN BE USED HERE
%         NeuralTrials.Trials.(CurrentTrial).SOT = ML_SOT-NeuralTrials.Trials.(CurrentTrial).TimeDiff;
%         NeuralTrials.Trials.(CurrentTrial).EOT = ML_EOT-NeuralTrials.Trials.(CurrentTrial).TimeDiff;
%         
%     else
%         % If no TimeDiff computed, remove this trial from the NeuralTrials
%         % structure
%         NeuralTrials.Trials = rmfield(NeuralTrials.Trials,CurrentTrial);
%     end
%     
% end
% 
% close(trlwaitbar)
% 
% toc(lagTimer)


%% SAVE NEURALTRIALS STRUCTURE
display('Saving NeuralTrials structure.');
saveNT = tic;
if ~isempty(fieldnames(NeuralTrials))

            if ~isdir(dirSave)
                mkdir(dirSave)
            end
            savestr = [dirSave NeuralTrials.SessionInfo.FileName];
            save(savestr,'NeuralTrials','-v7.3')
    
end
toc(saveNT)

% clearvars -except SpikeData NeuralTrials Options Slash FileName MonkeyName SessionDate InputFolder...
%     OutputFolder