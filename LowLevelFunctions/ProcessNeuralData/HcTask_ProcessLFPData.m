function HcTask_ProcessLFPData(XMazeStruct,LFPFile,OutputFolder)

% if cAnalog==0;
% No analog channels included in the file contianing the spike data
% (.nev or .plx) file. The Cerebus system, splits recordings in to
% two types (spike vs continuous). If there are no analog channels
% found on this recording, then the recording selected in previous
% sections will not contain your SyncPulse information. Load
% the appropriate file manually here.


%There should be a ns6 file with the same name as the nev
%     LFPFile = dir([InputFolder '*.ns6']);
%
%     if size(LFPFile,1)~=1
%
%         %Get nev file name
%         NevName = File(max(strfind(File, Slash))+1:end);
%
%         %Replace extension to ns6
%         NevName(end-2:end)='ns6';
%         FileIndex = strcmp({LFPFile.name}, NevName);
%
%         if ~sum(FileIndex) || sum(FileIndex)~=1
%
%             [s,v] = listdlg('PromptString','Select a file:',...
%                 'SelectionMode','single',...
%                 'ListString',{LFPFile.name});
%
%             if v
%                 LFPFile = LFPFile(s);
%             else
%                 warning('Error during file acquisition, self-destruct sequence initiated');
%                 LFPFile = [];
%             end
%
%
%         elseif sum(FileIndex)==1
%             LFPFile = LFPFile(FileIndex);
%
%         else
%             warning('Error during file acquisition, self-destruct sequence initiated');
%             LFPFile = [];
%         end
%     end
%

%     if ~isempty(LFPFile)


%To correct for the aliasing problem due to the fact that LFPs weren't
%filtered prior to digitizing, we need to load the ns6 raw file, using
%blackrock's NPMK library. We will then proceed to filter properly and
%downsample to a uniform 1kHz for all sessions.

%Balckrock's NPMK toolbox is used instead of neural share because there was
%a problem in accessing the 30 kHz signal with ns. 

%Get the basic channel properties
NS6 = openNSx(LFPFile, 't:1:2');

%Get recording info:
SR = NS6.MetaTags.SamplingFreq;
FinalSR = 1000;

NbrChans = size(NS6.ElectrodesInfo,2);

lfpwaitbar = waitbar(0, 'Extracting LFP data for all channels');

%Loop for all channels
for ch = 1: NbrChans
    
    waitbar(ch/NbrChans,lfpwaitbar);
    
    NS6 = openNSx(LFPFile, ['c:' num2str(ch)]);
    
    channel = deblank(NS6.ElectrodesInfo.Label);
    
    % added sept 2 2015 extra step because some old recordings have a dash
    % (-) in the names which is incompatible with matlab structure name
    channel = strrep(channel, '-', '_');
    
    %If there was a pause in the recording, the LFPs will be in two
    %separate cells.
    if iscell(NS6.Data)
        NS6.Data = cell2mat(NS6.Data);
    end
    
    %Removes zeros padding to have the first sample at time zero
    ZeroIndex = find(NS6.Data(1,:)~=0, 1, 'first')-1;
    
    NS6.Data(:,1:ZeroIndex) = [];
    
    NS6.Data = double(NS6.Data);
    LFPData = struct();
    
    % Before filtering, we need to convert the data from digital values
    % back to voltage values
    digiRange = double(NS6.ElectrodesInfo.MaxDigiValue);
    analogRange = double(NS6.ElectrodesInfo.MaxAnalogValue);
    
    NS6.Data = NS6.Data * analogRange / digiRange;
    
    %High-Pass filtering to correct for phase shifts
    %1 is butterworth, 0 is none
    if NS6.ElectrodesInfo.HighFilterType
        
        %High pass filtering to correct for phase shifts
        %Hardware filter data, for some reason need to divide by 1000 to
        %have frequency in Hz
        corner = double(NS6.ElectrodesInfo.HighFreqCorner)/1000;
        filtorder = NS6.ElectrodesInfo.HighFreqOrder;
        
        [b,a] = butter(filtorder, corner/(SR/2), 'high');
        NS6.Data = filter(b, a, fliplr(NS6.Data));
        
        %Low pass filtering to correct for phase shifts
        %Hardware filter data, for some reason need to divide by 1000 to
        %have frequency in Hz (no need for sos,g here since filter
        %order is usually low (i.e. 1).
        corner = double(NS6.ElectrodesInfo.LowFreqCorner)/1000;
        filtorder = NS6.ElectrodesInfo.LowFreqOrder;
        
        [b,a] = butter(filtorder, corner/(SR/2), 'low');
        
        %at this point the signal is still backwards
        NS6.Data = fliplr(filter(b, a, NS6.Data));
    end
    
    % low-pass from 250 before downsampling
    lpcorner = 250;
    filtorder = 4;
    
    %used a different syntax because of filtering problem (see matlab help
    %on IIR (butter):
    %Numerical Instability of Transfer Function Syntax
    %In general, use the [z,p,k] syntax to design IIR filters. To analyze
    %or implement your filter, you can then use the [z,p,k] output with zp2sos.
    %If you design the filter using the [b,a] syntax, you might encounter numerical problems.
    %These problems are due to round-off errors and can occur for n as low as 4.
    [z, p, k] = butter(filtorder,...
        lpcorner/(SR/2),...
        'low');
    [sos, g] = zp2sos(z, p, k);
    
    NS6.Data = filtfilt(sos, g, NS6.Data);
    
    %Power line removal (Zanos J neurophysiol 2012).
    %3rd order elliptical band-stop filter
    %60 HZ
    [z, p, k] = ellip(3, 1, 80, [58/(SR/2) 62/(SR/2)], 'stop');
    [sos, g] = zp2sos(z, p, k);
    NS6.Data = filtfilt(sos, g, NS6.Data);
    %120 HZ
    [z, p, k] = ellip(3, 1, 80, [118/(SR/2) 122/(SR/2)], 'stop');
    [sos, g] = zp2sos(z, p, k);
    NS6.Data = filtfilt(sos, g, NS6.Data);
    %180 HZ
    [z, p, k] = ellip(3, 1, 80, [178/(SR/2) 182/(SR/2)], 'stop');
    [sos, g] = zp2sos(z, p, k);
    NS6.Data = filtfilt(sos, g, NS6.Data);
    
    %Computer DownSampling factor
    if ~mod(SR, FinalSR)
        DSFactor = SR/FinalSR;
    else
        warning(['Sampling rate not compatible with ' num2str(FinalSR) ' Hz downsampling, no downsampling will occur']);
        DSFactor = 1;
    end
    
    LFPData.Channels.(channel).LFP.Voltage = NS6.Data(1:DSFactor:end);
    LFPData.Channels.(channel).LFP.Timestamps = 0 : DSFactor/SR : (numel(LFPData.Channels.(channel).LFP.Voltage)-1)*(DSFactor/SR);
    LFPData.Channels.(channel).LFPInfo = NS6.ElectrodesInfo;
    LFPData.Channels.(channel).LFPInfo.SampleRate = FinalSR;
    LFPData.Channels.(channel).LFPInfo.SampleCount = numel(LFPData.Channels.(channel).LFP.Voltage);
    
    LFPData.SessionInfo.MonkeyName = XMazeStruct.SessionInfo.MonkeyName;
    LFPData.SessionInfo.FileName = strrep(XMazeStruct.SessionInfo.FileName,'XMaze',[channel '_LFPData'] );
    
    savestr = [OutputFolder LFPData.SessionInfo.FileName];
    
    save(savestr,'LFPData', '-v7.3');
    
end

close(lfpwaitbar);
display(char(10));


%     end

% end