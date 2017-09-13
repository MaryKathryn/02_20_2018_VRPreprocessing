function HcTask_ProcessUnitData(MLStruct,NeuralTrials,File,OutputFolder)

 % Import, format and save sorted neural data from an HcTask session.
% 
% File Architecture:
% FileName
% /
% /Info
%   /SessionInfo
%       /Date
%       /Monkey
%       /FileName
%   /RecordingInfo
%       /TimeStampResolution
%       /FileComment
%   /TrialInfo
%       /TrialIDs
%       /TrialStartTs
%       /TrialEndTs
%   /UnitInfo
%       /Channels
%       /UnitNum
% /Data
%       /TrialRasters
%   /TimeStamps
%       /Units
%   /Waveforms
%       /Units
%
%
%
%
% To read Data:
%
% h5disp(fileName,'/Data')
% % % % TrialRasters
% h5readatt(fileName,'/Data/TrialRasters','PadTime')
% Rasters = h5read(fileName,'/Data/TrialRasters');
% Rasters = logical(Rasters);
% % % % Timestamps
% tsInfo = h5info(fileName,'/Data/TimeStamps');
% hdfUnitNames = {tsInfo.Datasets.Name}';
% for k = 1:length(hdfUnitNames)
%     TimeStamps{k} = h5read(fileName,['/Data/TimeStamps/' hdfUnitNames{k}]);
% end
% % % Waveforms
% wfInfo = h5info(fileName,'/Data/WaveForms');
% hdfUnitNames = {wfInfo.Datasets.Name}';
% for k = 1:length(hdfUnitNames)
%     WaveForms{k} = h5read(fileName,['/Data/WaveForms/' hdfUnitNames{k}]);
% end
% clear wfInfo k hdfUnitNames

monkeyName = MLStruct.SessionInfo.MonkeyName;
session = datestr(MLStruct.SessionInfo.Date,'yyyymmdd');

%% Create MLStruct from each behavioural data structure
MLTrlIDs = fieldnames(MLStruct.Trials);
NTTrlIDs = fieldnames(NeuralTrials.Trials);
% Remove MLStruct trials that are not included in Neural Trials. It cannot
% be the case that there are NeuralTrials files that are not in MLStruct,
% but I've included this check anyways
if ~isequal(MLTrlIDs,NTTrlIDs)
    trls2remove = MLTrlIDs(~ismember(MLTrlIDs,NTTrlIDs));
    MLStruct.Trials = rmfield(MLStruct.Trials,trls2remove);
    MLTrlIDs = fieldnames(MLStruct.Trials);
    trls2remove = NTTrlIDs(~ismember(NTTrlIDs,MLTrlIDs));
    NeuralTrials.Trials = rmfield(NeuralTrials.Trials,trls2remove);
    NTTrlIDs = fieldnames(NeuralTrials.Trials);
end
% Ensure the trials are consistent between NeuralTrials and MLStruct 
% That is, since only the same trials are now included, and they are in the
% same order
MLStruct.Trials = orderfields(MLStruct.Trials);
MLTrlIDs = fieldnames(MLStruct.Trials);
NeuralTrials.Trials = orderfields(NeuralTrials.Trials);
NTTrlIDs = fieldnames(NeuralTrials.Trials);



% TrialInfo
trlIDs      = MLTrlIDs;
numTrls     = length(trlIDs);

% MLStruct SOTs and EOTs( willl be used for Rasters)
% Get all trial start and end times from the Behavioural structures. Then
% convert into Cerebus time using TimeDiff
StructSOTs = cell2mat(struct2cell(structfun(@(x)(x.SOT_Time),...
    MLStruct.Trials,'uni',0)));
StructEOTs = cell2mat(struct2cell(structfun(@(x)(x.EOT_Time),...
    MLStruct.Trials,'uni',0)));
TimeDiffs = cell2mat(struct2cell(structfun(@(x)(x.TimeDiff),...
    NeuralTrials.Trials,'uni',0)));
CerebusSOTs = structfun(@(x)(x.SOT_Time),NeuralTrials.Trials);
CerebusEOTs = StructEOTs-TimeDiffs;


% % NeuralTrials SOTs and EOTs (retrieved, but should not be used for Rasters)
% try
%     CerebusSOTs     = struct2cell(structfun(@(x) (x.SOT), NeuralTrials.Trials,'uni',0));
%     CerebusEOTs     = struct2cell(structfun(@(x) (x.EOT), NeuralTrials.Trials,'uni',0));
% catch % If above does not work, it is because NeuralTrials was computed 
%     % with the wrong SOT and EOT from an older version of the SyncPulse
%     % function. Must delete the old NeuralTrials file, and re-run the
%     % HcTask_ProcessNeuralData with Options.Trials = 1
%     clear NeuralTrials
%     errordlg([monkeyName ' ' sessDate char(10) ... 
%         'Decimated NeuralTrials file found!' char(10)...
%         'Re-run with Options.Trials=1'])
%     return
% end
% if length(CerebusEOTs) == length(CerebusSOTs) && length(CerebusEOTs) == length(trlIDs) && length(trlIDs) == length(CerebusSOTs)
%     numTrls = length(trlIDs);
% else
%     errordlg([monkeyName ' ' sessDate char(10) ... 
%     'Uneven number of NT SOTs, EOTs and TrialIDs!'])
% end

%% Open sorted neural data file and retrieve basic file & unit information


[~, hfile] = ns_OpenFile(File);
% RecordingInfo
[~, FileInfo] = ns_GetFileInfo(hfile); % Gives you EntityCount, TimeStampResolution and TimeSpan and the time at which the file was created
% Build catalogue of entities RG: Entities are data entries (channels)
% within the file. These include your discretized (discontinuous) spike
% data, continuous analog data, serial and digital inputs, described below.
[~, EntityInfo] = ns_GetEntityInfo(hfile, 1:FileInfo.EntityCount);

% List of EntityIDs needed to retrieve the information and data
NeuralList = find([EntityInfo.EntityType] == 4);
% Neural entities are the timestamps of the segment entities. Each Neural
% entity corresponds to one defined single- or multi-unit cluster.
SegmentList = find([EntityInfo.EntityType] == 3);
% Segment entities are short, time-stamped snippets of data extracted from
% a high-sample data stream. In our case, this is the waveform information
% (time, voltage) for every spike. Each segment corresponds to one channel
% of the recording system.
AnalogList = find([EntityInfo.EntityType] == 2);
% Analog entities are continuous, sampled data that represent digitized
% analog signals such as position of eyes, LFPs, or TTL signals.
EventList = find([EntityInfo.EntityType] == 1); % Event entities are values from digital and serial inputs (time-stamped
% text or binary data packets). These are used to represent data such as
% trial markers, experimental events, digital input values, and embedded
% user comments.
% How many of a particular entity do we have
cNeural = length(NeuralList);   % Number of sorted units in the NEV file (artifact, multi- and single-units from all channels)
cSegment = length(SegmentList); % Number of input channels on the recording system (electrode, digital and analog channels)
cAnalog = length(AnalogList);   % Number of recorded continuous LFP channels
cEvent = length(EventList);     % Number of channels receiving digital inputs (including digital words channel)

% Spike channels are EntityType 3
% Find spike channels with some activity on them
channels = find([EntityInfo.EntityType]==3 & [EntityInfo.ItemCount] ~= 0);
SegmentLabels = {EntityInfo(SegmentList).EntityLabel}';



%% Prepare datasets for /Info branch of the HDF5 file


% SessionInfo
sessDate = MLStruct.SessionInfo.Date;
monkeyName = MLStruct.SessionInfo.MonkeyName;
spikeDataFileName = [monkeyName(1) '_' session '_SpikeData.hdf5'];

% RecordingInfo
recordingisi = FileInfo.TimeStampResolution;
filecomment  = FileInfo.FileComment;
if isempty(filecomment)
    filecomment = 'empty';
end


% UnitInfo
channelNames = {EntityInfo(NeuralList).EntityLabel};
    % Pad unit names to ensure they are the same length for the HDF5 array
    maxchar = max(cellfun(@length,channelNames));
    channelNames = cellfun(@(x) (sprintf(['% ' num2str(maxchar) 's'], x)),channelNames,'uni',0);
    channelNames = reshape(cell2mat(channelNames'),[length(channelNames) maxchar]);
unNums = cellfun(@(x) (strcmp(x,cellstr(channelNames))), unique(cellstr(channelNames)),'uni',0) ;
    unNums = cell2mat(cellfun(@(x) (cumsum(x).*x),unNums,'uni',0)');
    unNums = uint8(sum(unNums,2)-1);

    

%% Save /Info groups (/SessionInfo, /RecordingInfo, /UnitInfo, /TrialInfo) 

fileName = [OutputFolder spikeDataFileName];

% Write data to HDF5 file
plist='H5P_DEFAULT'; %default properties
fileID = H5F.create(fileName,'H5F_ACC_TRUNC','H5P_DEFAULT','H5P_DEFAULT'); %create new file
    infoGID = H5G.create(fileID,'Info',plist);
            sesInfoGID = H5G.create(infoGID,'SessionInfo',plist); % Basic session information from ML file  
                HDF5InsertData(sesInfoGID,'Date',sessDate);
                HDF5InsertData(sesInfoGID,'MonkeyName',monkeyName);
                HDF5InsertData(sesInfoGID,'FileName',spikeDataFileName);
            H5G.close(sesInfoGID) 
            recInfoGID = H5G.create(infoGID,'RecordingInfo',plist); % Recording file information
%                 HDF5InsertData(recInfoGID,'FileType',FileInfo.FileType)
%                 HDF5InsertData(recInfoGID,'EntityCount',FileInfo.EntityCount)
                HDF5InsertData(recInfoGID,'TimeStampResolution',recordingisi)
%                 HDF5InsertData(recInfoGID,'TimeSpan',FileInfo.TimeSpan)
%                 HDF5InsertData(recInfoGID,'AppName',FileInfo.AppName)
%                 HDF5InsertData(recInfoGID,'Time_Year',FileInfo.Time_Year)
%                 HDF5InsertData(recInfoGID,'Time_Month',FileInfo.Time_Month)
%                 HDF5InsertData(recInfoGID,'Time_Day',FileInfo.Time_Day)
%                 HDF5InsertData(recInfoGID,'Time_Hour',FileInfo.Time_Hour)
%                 HDF5InsertData(recInfoGID,'Time_Min',FileInfo.Time_Min)
%                 HDF5InsertData(recInfoGID,'Time_Sec',FileInfo.Time_Sec)
%                 HDF5InsertData(recInfoGID,'Time_MilliSec',FileInfo.Time_MilliSec)
                HDF5InsertData(recInfoGID,'FileComment',filecomment)
            H5G.close(recInfoGID)
            trlInfoGID = H5G.create(infoGID,'TrialInfo',plist); % Trial information
                HDF5InsertData(trlInfoGID,'TrialIDs',cell2mat(trlIDs))
                HDF5InsertData(trlInfoGID,'TrialStartTs',CerebusSOTs)
                HDF5InsertData(trlInfoGID,'TrialEndTs',CerebusEOTs)
            H5G.close(trlInfoGID)
            unitInfoGID = H5G.create(infoGID,'UnitInfo',plist); % Trial information
                HDF5InsertData(unitInfoGID,'Channels',channelNames)
                HDF5InsertData(unitInfoGID,'UnitNum',unNums)
            H5G.close(unitInfoGID)
        H5G.close(infoGID)
% H5F.close(fileID)

% Display Comments from the recording file:
disp([char(10) 'Recording File Comments: ' char(10) FileInfo.FileComment char(10)]);


% % TO EXPLORE THE HDF5 FILE:
% % To see the architecture of your HDF5 file
% h5disp(fileName)
% h5info(fileName)
% % To read the data just written in to the file
% cell2mat(h5read(fileName,'/Info/SessionInfo/Date/'))
% cell2mat(h5read(fileName,'/Info/SessionInfo/MonkeyName/'))
% cell2mat(h5read(fileName,'/Info/SessionInfo/FileName/'))
% cell2mat(h5read(fileName,'/Info/RecordingInfo/FileComment/'))   
% % To see whether your file has open groups
% H5F.get_obj_ids(fileID,'H5F_OBJ_GROUP',100) % Returns an error if your file is not open
% % To retrieve all units that belong to a channel
% chs = cell2mat(h5read(fileName,'/Info/UnitInfo/Channels/'))
% chs = cellstr(a)
% chnames = unique(b)
% chunits = cellfun(@(x) (strcmp(x,chs)), chnames,'uni',0) 
% h5read(fileName,'/Info/UnitInfo/UnitNum/')



%% Extract unit data
% NOTE: The first unit of any channel is always the multi-unit/noise
% cluster. Isolated units will always be unit numbers > 1, unless the
% multi-unit/noise cluster was deleted entirely during spike sorting.

unitcounter = 1;

for ch = 1:length(channels)
    
    chname = cell2mat(SegmentLabels(channels(ch)));
%     display(['Retrieving spike data for channel ' chname]);
    
    % These functions extract the segment data (spike times, waveforms) for
    % each threshold crossing on your selected channel. Note, this extracts
    % data for all threshold crossings on the selected spike channel,
    % regardless of which unit those spikes actually belong to.
    % Properties of  the filter combine all analog and digital filters for
    % spike data
    % RG question: Unsure of how this would affect the import of tetrode data,
    % as all spikes may not be visible on all channels.
    
    % Returns sample rate and number of data points per waveform.
    [~, nsSegmentInfo] = ns_GetSegmentInfo(hfile, SegmentList(channels(ch)));
    
    % Returns information about the filter settings, recorded voltage and
    % resolution. nsSegmentSourceInfo.Resolution refers to the resolution of
    % continuous analog signal to discretized digital conversion in uV.
    [~, nsSegmentSourceInfo] = ...
        ns_GetSegmentSourceInfo(hfile, SegmentList(channels(ch)), 1);
    
    % Retrive all neural data from this channel (timestamps, waveforms, number
    % of samples per waveform, units)
    [~, chan_ts, chan_wf, ~, chan_unitIDs] = ...
        ns_GetSegmentData ( hfile, SegmentList(channels(ch)), ...
        1:EntityInfo(SegmentList(channels(ch))).ItemCount ) ;
    
    units = unique(chan_unitIDs);

    for unit = 1:length(units)
%         tic
        
        unitInds = find(chan_unitIDs==units(unit));
        timestamps{1,unitcounter}    = chan_ts(unitInds);
        waveforms{1,unitcounter}     = chan_wf(:,unitInds);
        unitID(1,unitcounter)        = uint8(chan_unitIDs(unitInds(1))+1);
%         uname = ['Unit_' num2str(unitID{unit})];
%         SpikeData.Channels.(chname).Units.(uname).Timestamps = ...
%             timestamps{unit};
%         SpikeData.Channels.(chname).Units.(uname).Waveforms = ...
%             waveforms{unit};
        
        display(['Retrieving spike data for ' chname ' Unit' num2str(unitID(unit))]);
%         toc
        
        if unitcounter < cNeural        
            unitcounter = unitcounter+1;
        end
        
    end % Unit counter

end % Channel counter
% display(char(10));


% return



%% Prepare datasets for /Data

display(['Prepping unit data for ' spikeDataFileName]);
dataTimer = tic;

% % % Timestamps & Waveforms
% whos('timestamps');
% whos('waveforms');
divider = repmat('_',unitcounter,1);
hdfUnitNames = cellstr([channelNames divider num2str(unNums)]);


%%% RG 20170505:
%%%%% Rasters from 96-channel arrays can be larger than the 4Gb max chunk
%%%%% size for HDF5 files. Therefore, I'm commenting out this section of
%%%%% the function, and opting to write a helper function to read
%%%%% timestamps and create rasters. This function will have the option to
%%%%% select which units, as well as which trials are retrieved. 
% % % % Trial Rasters
% rasterPadTime = 0.5; % seconds; will pad each trial with additional time before and after each trial
% % Find the longest trial; will pad all trial rasters to this length+(2*rasterPadTime)
% maxTrlLength = max(StructEOTs-StructSOTs) + 2*rasterPadTime;
% maxTrlInds = floor(maxTrlLength*1000); % Want this to have 1ms resolution
% % Initialize the matrix of for all trial rasters.
% % [number of trials (rows) * longest trial length (column) * number of units (pages)]
% Rasters = repmat(uint8(2),numTrls,maxTrlInds,cNeural);
% for trl=1:numTrls % Trial counter
%     
%     % Find this trial raster's start and end time
%     % MUST CHANGE THIS TO USE THE MLSTRUCT SOTS AND EOTS, NOT THE NEURAL
%     % TRIALS INFORMATION
%     rastStart  = CerebusSOTs(trl)-rasterPadTime; 
%     rastEnd    = CerebusEOTs(trl)+rasterPadTime;
%     rastEndInd = floor((rastEnd-rastStart)*1000);
%     % Convert all indices corresponding to trial time to 0
%     Rasters(trl,1:rastEndInd,:) = uint8(0);
% 
%      
%     % Find spike times for each unit
%     % Spike times here are relative to the start of the recording session,
%     % not the trial of interest.
%     % CHECK THIS RETREIVES THE TIME OF EACH SPIKE, NOT AN INDEX
%     trlSpikeTimes = cellfun(@(x) (x(rastStart < x & x < rastEnd)),timestamps,'uni',0);
%     % Convert spike times to spike times relative to the start of the
%     % raster period (Trial Start - raster pad time)
%     rastSpikeTimes = cellfun(@(x) (x - rastStart),trlSpikeTimes,'uni',0);
%     % Convert this time to an index for the trial raster
%     rastSpikeInds = cellfun(@(x) (floor( x*1000)),rastSpikeTimes,'uni',0);
%     
%     % Convert 0's to 1's for all trlSpikeInds
%     for u = 1:cNeural
%         rastSpikeInds{u}(rastSpikeInds{u}==0) = 1;
%         Rasters(trl,rastSpikeInds{u},u) = uint8(1);
%     end
%     
% %     % CHECKS
% %     colors = jet(cNeural)
% %     close all; 
% %     figure('units','normalized','position',[0 0 1 1]); 
% %     hold on; 
% %     for a = 1:cNeural 
% %         hold on
% %         plot(Rasters(trl,:,a),'.','color',colors(a,:)); 
% %         title(num2str(a))
% %         pause; 
% %     end
% end
% 
% toc(dataTimer)
% % whos('Rasters')



%% Save /Data groups/datasets (/TrialRasters, /Timestamps, /Waveforms)

display(['Saving unit data to ' spikeDataFileName]);
saveTimer = tic; 

% Write data to HDF5 file
dataGID = H5G.create(fileID,'Data',plist);
%     %%% Trial rasters 
%     % for each unit
%     HDF5InsertData(dataGID,'TrialRasters',Rasters);
%     h5writeatt(fileName,'/Data/TrialRasters','PadTime',['0.5 seconds before trial start, 0.5s after;' char(10) 'i.e. index 500 of every trial = trial start'])
    %%% Timestamps 
    % group containing one dataset for each unit in the file
    tsGID = H5G.create(dataGID,'TimeStamps',plist);
        for unit = 1:cNeural
            HDF5InsertData(tsGID,hdfUnitNames{unit},timestamps{unit});
        end  
    H5G.close(tsGID) 
    %%% Waveforms 
    % group containing one dataset with each spike waveform for each unit
    % in the file
    wfGID = H5G.create(dataGID,'WaveForms',plist); 
        for unit = 1:length(unNums)
            HDF5InsertData(wfGID,hdfUnitNames{unit},waveforms{unit});
        end
    H5G.close(wfGID)
H5G.close(dataGID)
H5F.close(fileID)

toc(saveTimer)


end