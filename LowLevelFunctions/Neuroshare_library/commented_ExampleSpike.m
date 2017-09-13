function Example()

% Prompt for the correct DLL 
% RG: DLL's are produced by the vendors (e.g., Plexon, Blackrock) and allow
% Neuroshare to read many proprietary filetypes.
disp(' ');  % Blank line
DLLName = input('DLL Name (i.e. nsNEVLibrary64.dll) : ', 's');

% Load the appropriate DLL
[nsresult] = ns_SetLibrary(DLLName);
if (nsresult ~= 0)
    disp('DLL was not found!');
    return
end

% Find out the data file from user
disp(' ');  % Blank line
filename = input('Data file: ', 's');


% Load data file and display some info about the file
% Open data file
[nsresult, hfile] = ns_OpenFile(filename);
if (nsresult ~= 0)
    disp('Data file did not open!');
    return
end


% Get file information
[nsresult, FileInfo] = ns_GetFileInfo(hfile);
% Gives you EntityCount, TimeStampResolution and TimeSpan and the time at
% which the file was created
if (nsresult ~= 0)
    disp('Data file information did not load!');
    return
end

% Define some variables needed for firing rates
stepsize = 0.001; % seconds

if FileInfo.TimeSpan > 150  % Limit the timespan shown in the graphs to 150 seconds
    totaltime = 150;
else
    totaltime = FileInfo.TimeSpan; % seconds
end
time = 0 : stepsize : totaltime;   % Initialize time axis for gaussian plot

% Build catalogue of entities SdT : Entities contains one or more indexed
% data entries that are ordered by increasing time. There are four types of
% entities, listed below. Each entity type contains specific information
% about your neural data.
[nsresult, EntityInfo] = ns_GetEntityInfo(hfile, [1 : 1 : FileInfo.EntityCount]);

% List of EntityIDs needed to retrieve the information and data
NeuralList = find([EntityInfo.EntityType] == 4); %Neural event entity timestamps of event and segment entitities that
%are known to represent neural action potential firing times.  For example, if a segment entity contains sorted neural
%spike waveforms, each sorted unit is also exported as a neural entity.
SegmentList = find([EntityInfo.EntityType] == 3); %Segment entities are short, time-stamped segments of digitized analog signals in which
%the  segments are separated by variable amounts of time. They represent discontinuous analog signals such as
%extracellular spike waveforms from electrodes or groups of electrodes
AnalogList = find([EntityInfo.EntityType] == 2); %Analog entities are continuous, sampled data that represent digitized
%analog signals such as position of eyes, or LFPs.
EventList = find([EntityInfo.EntityType] == 1); %Event entities are discrete events that consist of small time-stamped
%text or binary data packets.  These are used to represent data such as trial markers, experimental events, digital input
%values, and embedded user comments. They do include digital inputs and
%serail inputs to Cerebus

% How many of a particular entity do we have
cNeural = length(NeuralList);
cSegment = length(SegmentList);
cAnalog = length(AnalogList);
cEvent = length(EventList);

%clear FileInfo;

if (cNeural == 0)
    disp('No neural events available!');
end

if (cSegment == 0)
    disp('No segment entities available!');
    return;     % It does not make sense to continue in this particular analysis
    % if there are no segment entities.
end

if (cAnalog == 0)
    disp('No analog entities available!'); %SdT : I beleive you should only expect analog entities in your .ns files, not your .nev
    
end

if (cEvent == 0)
    disp('No event entities available!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Look at Event entities
%Get some event identity info (i.e. Digital input, Serial input/output
%Define Entity of interest
EntityID = input('Which Entity ID to open : ');
%Get info
[nsresult, nsEventInfo] = ns_GetEventInfo(hfile, EntityID);

%Define index to look at
Index = input('Which Index to open (e.g. 1 or [1 2 3]) : ');
%Retrieve event data by index
[nsresult, TimeStamp_event, Data, DataSize] = ns_GetEventData(hfile, EntityID, Index);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Have user PICK A CHANNEL or channels for further analysis
% Show the user how many channels are available
disp(' ');
disp(['There are ' num2str(cSegment) ' channels.']);
disp(' ');
%SdT : In .nev files, Segments include all electrodes channels (128 total, even if only 32 are used), and all 16 analog input channels
channel = input('Which data channels would you like to display? (e.g. 1 or [1 2 3])');
if (channel > cSegment)  % Have to check that the selected channel actually exists
    disp('Channel does not exist');
    return
end
%clear cNeural cSegment;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Load the waveform data and do the analysis on segments
%

% These three functions are for demonstration purposes
%SdT : extracts information from segment entities for the selected channels
%out of the 144 (in case of Cerebus) (sampling rate of source, filter type, etc.)
%Properties of the filter combine all analog and digital filters for spike
%data
[nsresult, nsSegmentInfo] = ns_GetSegmentInfo(hfile, SegmentList(channel));
[nsresult, nsSegmentSourceInfo] = ns_GetSegmentSourceInfo(hfile, SegmentList(channel), 1);

% Load the first N waveforms on each selected channel
%SdT : timestamps_wf is the time stamp of each waveform. Samplecount is the 
%number of samples per waveform (48 for
%all of them. Waveform is the actual y values of each waveforms.
%unitIDs is the identifier of a sorted unit, 0 being noise or unsorted, 1
%to 5 being sorted ID (works with Guillaumes Spike Sorter)
disp(' ');
N = input('How many segments do you want to see, starting from the beginning ?  '); %Set the number of segments you want to analyze
disp(' ');
[nsresult, timestamps_wf, waveforms, sampleCount, unitIDs] = ns_GetSegmentData(hfile, SegmentList(channel), [1 : 1 : N]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Dig in Neural event entities
% Reduce the data set to only the first 150 seconds of data
%
% SdT : Creates a list of labels of all the electrodes from which a Unit as
% been identified (sorted perhaps). (YOull get the same number of labels
% and electrodes if no units are sorted. 1 unit per electrode)
NeuralLabels = strvcat(EntityInfo(NeuralList).EntityLabel);

%cChannel is useful when you slected multiple channels to analyze
%simultaniously
for cChannel = 1  : length(channel)
    % You Have to figure out which Neural entities correspond to the selected segment entities
    %Find the entity label that corresponds to the current channel being
    %analyzed (i.e. elec 1) and match with the list of Neural labels to
    %find the index of this label in the NeuralLabels list.
    list = strmatch(EntityInfo(SegmentList(channel(cChannel))).EntityLabel, NeuralLabels, 'exact');
    
    % Retrieve the data
    %SdT : Gives you the basic info about the units on this channel
    [nsresult, NeuralInfo] = ns_GetNeuralInfo(hfile, NeuralList(list));
    %SdT : gives you the timestamps of spikes of each Unit, padded with the
    %timestamps of the Unit with the biggest Index number.
    %SDT [nsresult, NeuralData] = ns_GetNeuralData(hfile, NeuralList(list), 1, max([EntityInfo(NeuralList(list)).ItemCount]));
    % Error (-7) if you ask for more timestasmps than available for a
    %unit ) The min here will mak sure you
    %have an equal number of spikes per units (cutting down to the lowest number)
    %%%%%%%%%%%%%%%FIX THIS ISSUE!
    %This problem does not exist on the 32-bit mexprog version of
    %Neuroshare. It seems to be a bug specific to the 64-bit version
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %FIX TO QUICKLY OUTPUT NeuralData structure
    clear timestamps_units
    clear NeuralData
    
    for i = 1:length(list)
        [nsresult, NeuralData] = ns_GetNeuralData(hfile, NeuralList(list(i)), 1, [EntityInfo(NeuralList(list(i))).ItemCount]);
        timestamps_units{i} = NeuralData;
    end
    
    SpikeData.(['Channel_' num2str(cChannel)]) = timestamps_units;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %THIS ORIGINAL SECTION OF NEUROSHARE COMMENTED BY SDT BECAUSE OF FIX
%     %Find the spike indices for the first 150 sec of the experiment
%     %SdT : This find function reads first column first till the end, then second
%     %column, retaining the same count. Gives you a single vector with
%     %"indices"
%     if (totaltime == 150)
%         ind = find(NeuralData <= 150);
%     else
%         ind = [1:1:size(NeuralData, 1) * size(NeuralData, 2)];
%     end
%     
%     % Get the neural timestamps
%     %SdT : Reshapes the NeuralData matrix as a single vector so to match
%     %the indices provided by the previous "find" search
%     NeuralData_reshaped = reshape(NeuralData, size(NeuralData, 1) * size(NeuralData, 2), 1);
%     %Gives you the actual timestamps of each spikes for each unit,
%     %organized in a single, confusing, vector. I would run the previous
%     %find on each column seperatly, and register the timestamps, one column
%     %per unit.
%     timestamps(cChannel) = {NeuralData_reshaped(ind)};
%     
%     % Match the neural events with their unit ID
%     %SdT : This line assumes that each units have the same number of
%     %samples, which is not always the case
%     temp = ones(length(ind) / length(list), 1) * [NeuralInfo(:).SourceUnitID];
%     %SdT : This unit matrix is suppose to tell you to which units each row
%     %of the timestamps array refers. Here, this is not correct because we
%     %do not have the same number of samples in for each unit
%     units(cChannel) = {temp(:)};
%     % Remember how many neural events were found
%     ItemCount(cChannel) = length(timestamps{cChannel});
%     
%     NeuralData_reshaped(:) = [];
end
%clear NeuralData SegmentLabels NeuralLabels NeuralInfo;
%EntityInfo = rmfield(EntityInfo, 'ItemCount');

% Close data file. Should be done by the library but just in case.
ns_CloseFile(hfile);

% Unload DLL
clear mexprog;

%SdT : End of neural data extraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%Function to plot some waveforms, raster plot, ISI
for cChannel = 1 : 1 : length(channel),
    % Plot neural activity on that channel
    figure;
    subplot(4, 1, 1);
    % Plot the first 100 waveforms on that channel and apply different
    % colors for different unitIDs
    
    % Unclassified units
    %SdT : It seems here that it is assumed that Unclassified units have
    %unitIDs of 0
    ind = find(unitIDs(:, cChannel) == 0);
    if (~isempty(ind))
        plot(waveforms(:, ind, cChannel), 'Color', [.8 .8 .8]); % in grey
        hold on;
    end
    % Classified Units
    Colors = ['y' 'm' 'c' 'r' 'g'];
    %SdT : Classified units have UnitIDs of 1 to 5
    for cUnit = 1 : 1 : 5,
        ind = find(unitIDs(:, cChannel) == cUnit);
        if (~isempty(ind))
            plot(waveforms(:, ind, cChannel), Colors(cUnit)); %cChannel here beacuse it could be a multiple dimensions array if multiple channels were selected previously for the analyses
            hold on;
        end
    end
    
    axis tight;
    title(['Waveforms - ' EntityInfo(SegmentList(channel(cChannel))).EntityLabel]);
    ylabel('Voltage (uV)');
    
    subplot(4, 1, 2);
    % This creates lines (raster plot) for the unit
    % SdT : repmat replicates the same matrix on number of row and columns
    % defined by users. Here, just replicate it on another line.
    a = repmat(timestamps{cChannel}', 2, 1);          % X Values don't change so just two rows with the same values
    b = repmat([0.6; 1.4], 1, ItemCount(cChannel));   % Start and end points of each line
    
    % This plots the lines for each unit
    % Unclassified units
    ind = find(units{cChannel} == 0);   %SdT : Possible error here : should be unitIDs instead of units
    hold on;
    if (~isempty(ind))
        plot(a(1, ind), b(:, ind), 'Color', [.8 .8 .8]); % in grey
    end
    
    % Classified Units
    for cUnit = 1 : 1 : 5
        ind = find(unitIDs(:, cChannel) == cUnit);
        if (~isempty(ind))
            plot(a(:, ind), b(:, ind), Colors(cUnit));    %ERROR : Index exceeds ,matrix dimensions ::::: Some values of ind exceeds matrix dimension of a and b
        end
    end
    
%      % Classified Units SDT
%     for cUnit = 1 : 1 : size(NeuralData, 2)
%         if (~isempty(NeuralData))
%             plot(NeuralData(:, cUnit));    %ERROR : Index exceeds ,matrix dimensions ::::: Some values of ind exceeds matrix dimension of a and b
%         end
%     end
    
    title(['Raster plot - ' EntityInfo(SegmentList(channel(cChannel))).EntityLabel]);
    xlabel('Time [s]');
    set(gca, 'ytick', []);
    xlim([0 round(totaltime)]);
    %clear a b;
    
    % Use a sliding window 100 ms (delta t) long (equivalent to linear kernel)
    dt = 0.1; % ms
    for i = totaltime : -stepsize : 0,
        gaussdata(round(i / stepsize + 1)) = sum(exp(-(i - timestamps{cChannel}(find((timestamps{cChannel} > i - 4 * dt) & (timestamps{cChannel} < i + 4 * dt)))').^2 / (2 * dt^2))) / sqrt(2 * 3.1415) / dt;
    end
    
    % Throw away firing rates higher than 200 because they are not real and show up at the beginning
    gaussdata(find(gaussdata > 200)) = 0;
    
    subplot(4, 1, 3);
    plot(time, gaussdata(1 : length(time)));
    title(['Approximated firing rate (Gaussian) - ' EntityInfo(SegmentList(channel(cChannel))).EntityLabel]);
    xlabel('Time [s]');
    ylabel('Rate [Hz]');
    clear gaussdata;
    
    % Interspike Interval Distribution (1 ms resolution)
    if (length(timestamps{cChannel}) > 1)  % Have to make sure that we have at least two data points
        nbins = 500; % ms
        dtime = (timestamps{cChannel}(2 : end) - timestamps{cChannel}(1 : end - 1)) * 1000;
        
        subplot(4, 1, 4);
        bar(histc(dtime, 1 : 1 : nbins));
        axis tight;
        title(['Interspike Interval Histogram - ' EntityInfo(SegmentList(channel(cChannel))).EntityLabel]);
        xlabel('Time [ms]');
        ylabel('# Spikes');
        clear dtime;
    end
end
