function [smoothedData,reformatted] = removeSections(smoothedData,BadData,reformatted)
% removes the x and y data from smoothedData and reformatted at the
% timestamps that match BadData, because these are not analyzable.
% It also carries out the PSOcalc because there are post saccadic
% oscillations that have to be removed. 
blkStartOffset = 40;
blkEndOffset   = 70;
velocityThreshold = 40;

%first, get the indices that match the timeStamps of the badData
BadDataInds = ismember(smoothedData(:,3),BadData);
%Make sure the first point is not a BadData Index, if so, shift it
%to 2, and both will get removed by the cleaning method
if BadDataInds(1) == 1
    BadDataInds(1) = 0;
    BadDataInds(2) = 1;
end
%now, there are sometimes BadData that is too close, so generate an
%array that matches up the end of one event with the start of the
%next like so
%[283 1638 1864  5153  5260  5267  5271  5339    0]
%[0   329  1670  1902  5156  5261  5268  5330  5342]

%However, first have to check that it actually gets back to good
%Data
%             if length(find(diff(BadDataInds)==1))==length(find(diff(BadDataInds)==-1))
%                 forRemoval = [find(diff(BadDataInds)==1)' 0;...
%                     -50 find(diff(BadDataInds)==-1)' ];
%             elseif length(find(diff(BadDataInds)==1))>length(find(diff(BadDataInds)==-1))
%                 forRemoval = [find(diff(BadDataInds)==1)' 0;...
%                     -50 find(diff(BadDataInds)==-1)' length(BadDataInds)];
%             else
%                 forRemoval = [1 find(diff(BadDataInds)==1)' 0;...
%                     -50 find(diff(BadDataInds)==-1)' ];
%             end


% Had issue with Wo 20150410 where the for removal matrix was:
%   [ 1 10 0; -50 2 12]
% It ends up being : [1;-50] since all columns are < 25
% samples. TO fix will replace the -50 and 0 by NaN.
if length(find(diff(BadDataInds)==1))==length(find(diff(BadDataInds)==-1))
    forRemoval = [find(diff(BadDataInds)==1)' NaN;...
        NaN find(diff(BadDataInds)==-1)' ];
elseif length(find(diff(BadDataInds)==1))>length(find(diff(BadDataInds)==-1))
    forRemoval = [find(diff(BadDataInds)==1)' NaN;...
        NaN find(diff(BadDataInds)==-1)' length(BadDataInds)];
else
    forRemoval = [1 find(diff(BadDataInds)==1)' NaN;...
        NaN find(diff(BadDataInds)==-1)' ];
end

%take all the columns that have less than 25 samples between them,
%and get rid of them
forRemoval(:,abs(diff(forRemoval,1))<25) = [];
%             forRemoval(2) = 0;


%now, StartInds is all the nonzero values in the first row, and
startInds = forRemoval(1,~isnan(forRemoval(1,:)));
endInds = forRemoval(2,~isnan(forRemoval(2,:)));

%Initialize the array for where to start NaNing Data
NaNStart = zeros(length(startInds),1);
NaNEnd = zeros(length(startInds),1);
Inds2NaN = zeros(length(smoothedData),1);

% does a cleaning for each section of BadData
for section = 1:length(startInds)
    %generate the speed for before the bad section (the start ind -
    %the offset
    
    if (startInds(section) - blkStartOffset) < 1 %if there isn't enough space for the
        %onset, start at the first index
        OnsetVelocity = sqrt(diff(smoothedData(1:startInds(section),1)).^2 +...
            diff(smoothedData(1:startInds(section),2)).^2 )./...
            diff(smoothedData(1:startInds(section),3));
    else %otherwise, start at the the section start - blink start offset
        OnsetVelocity = sqrt(diff(smoothedData(startInds(section)-...sqrt of x from the section before the
            blkStartOffset:startInds(section),1)).^2 +...
            diff(smoothedData(startInds(section)-...
            blkStartOffset:startInds(section),2)).^2 )./...
            diff(smoothedData(startInds(section)-...
            blkStartOffset:startInds(section),3));
    end
    
    if isempty(find(OnsetVelocity<velocityThreshold,1,'last'))
        NaNStart(section) = startInds(section)-length(OnsetVelocity);
    else
        NaNStart(section) =max(1,startInds(section)-blkStartOffset+find(OnsetVelocity<velocityThreshold,1,'last'));
    end
    
    %generate the velocity profile for after the period. Check that
    %the offset does not go beyond the end of the samples
    if blkEndOffset+endInds(section)<length(smoothedData(:,1))
        OffsetVelocity = [sqrt(diff(smoothedData(endInds(section):...
            blkEndOffset+endInds(section),1)).^2 +...
            diff(smoothedData(endInds(section):...
            blkEndOffset+endInds(section),2)).^2 )./...
            diff(smoothedData(endInds(section):...
            blkEndOffset+endInds(section),3)); 51];
    else % if it does, just go to the end of the  indices
        OffsetVelocity = [sqrt(diff(smoothedData(endInds(section):...
            end,1)).^2 +...
            diff(smoothedData(endInds(section):...
            end,2)).^2 )./...
            diff(smoothedData(endInds(section):...
            end,3)); 51];
    end
    %goes through and finds where it stays below velocityThreshold for at
    %least 6 samples, uses that as the beginning of the stable
    %period
    fast = [1;find(OffsetVelocity>velocityThreshold)];
    blocks = diff(fast);
    startOfGoodBlock = sum(blocks(1:find(blocks>=7,1)-1));
    start = min([startOfGoodBlock blkEndOffset]);
    %run the PSO offset test and then use that as the actualy
    %offset. This will remove PSO's from saccades from off
    %screen.
    if endInds(section)+start+31>= length(smoothedData(:,1))
        PSOEnd = length(smoothedData(:,1));
    else
        PSOX = PSOCalc(smoothedData(endInds(section)+start:...
            endInds(section)+start+31,1:2:3));
        PSOXind = find(smoothedData(:,3)==PSOX);
        
        PSOY = PSOCalc(smoothedData(endInds(section)+start:...
            endInds(section)+start+31,2:3));
        PSOYind = find(smoothedData(:,3)==PSOY);
        
        PSOEnd = max(PSOYind,PSOXind);
        
    end
    NaNEnd(section) = PSOEnd;
    Inds2NaN(NaNStart(section):NaNEnd(section)) = 1;
end
% NaN all the indices of BadData and the leadup and cool down as
% well
smoothedData(Inds2NaN==1,1:2)=NaN;
if ~isempty('reformatted')
    reformatted(Inds2NaN==1,1:2) = NaN;
    
end
end