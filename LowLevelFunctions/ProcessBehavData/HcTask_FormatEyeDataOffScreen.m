%% Function - Cleans Eye-Data and converts it to Degrees
% 1) Takes the eyesignal from each session, and identifies each portion where
% there is a blink (-7936) or beyond the range of the eye Tracker (0,0) and
% saves the time stampes in BadData for later identification.

% 2) Then the eyeSamples are converted to degrees using the calibration
% matrix.The conversion from units to degrees depends on eye calibration
% saved in the SessionInfo

% 3) Sample issues are addressed by removing double samples, and then filling
% in missing samples with a mean value of adjacent samples. This is also
% run on BadData. GoalsAppear and GoalRot indices are saved via their
% timestamp, and then updated to the new index, based on the saved
% timestamp.

% 4) Identify sample where the velocity is above 2000 degrees/s, and add it
% to MoreBadData list of indices (this is already corrected, and doesn't
% need to be interpolated.

% 5) Smoothing is done via a Savitzky Golay filter with order = 2 and window =
% 11. This assumes saccades are at least 10 ms.

% 6) Finally, the saved timestamps in BadData are identified in groups, for
% each section, the pre- and proceeding 20 ms are analyzed for the closest
% index where the velocity was below 30 degrees/s. This is assumed to be
% the last/first time where a valid, stable eyesample can be calculated.
% The blinks and bad samples are NaN'd. This replaces previous blink
% removal scripts

% 7) Create the PositionMatrix that combines the Unreal and EyeData. Also
% removes eye data that did not fall within the dva of the monitor. 

%%
 function [MLStruct] = HcTask_FormatEyeDataOffScreen(MLStruct, PosMatrix)

% Loops for every trial in structure
List_Trials = fieldnames(MLStruct.Trials);
blink = -7936;
blkEndOffset   = 70;
velocityThreshold = 40;
physiologicalCutOff = 4000;
AccelPhysiologicalCutoff = 200000;
MaxSpikes = cell(size(List_Trials, 1),2);
BadStart = [];


hWaitBar = waitbar(0, 'Formatting Eye Data.');
%%
for k=1:size(List_Trials, 1)
    if ~isempty(MLStruct.Trials.(List_Trials{k}).EyeSamples) && ...
            length(MLStruct.Trials.(List_Trials{k}).EyeSamples)>75;
        %% 1) Create BadData - the timestamps for all of the bad samples
        BadData = MLStruct.Trials.(List_Trials{k}).EyeSamples(... get eyesameple times where
            ((MLStruct.Trials.(List_Trials{k}).EyeSamples(:,1)==blink &... both x and y == blink
            MLStruct.Trials.(List_Trials{k}).EyeSamples(:,2)==blink)| ...
            (MLStruct.Trials.(List_Trials{k}).EyeSamples(:,1)==0 & ...or both x and y ==0
            MLStruct.Trials.(List_Trials{k}).EyeSamples(:,2)==0)),3);
        %% 2) Convert eye signal to degrees
        % Eye X = OffsetX + a * HREFX + b * HREFX^2 + c * HREFY + d * HREFY^2
        % Eye Y = OffsetY + f * HREFX + g * HREFX^2 + h * HREFY + i * HREFY^2
        
        tempEyeX = MLStruct.Trials.(List_Trials{k}).EyeSamples(:,1);
        tempEyeY = MLStruct.Trials.(List_Trials{k}).EyeSamples(:,2);
        
        %Convert to degrees
        tempEyeDegX = [ones(size(tempEyeX,1),1) tempEyeX tempEyeX.^2 tempEyeY tempEyeY.^2] * MLStruct.SessionInfo.EyeCalibration(1,:)';
        tempEyeDegY = [ones(size(tempEyeX,1),1) tempEyeX tempEyeX.^2 tempEyeY tempEyeY.^2] * MLStruct.SessionInfo.EyeCalibration(2,:)';
        
        for w=1:size(tempEyeDegX,1)
            
            %Quadrant specific correction factor
            if tempEyeDegX(w) < 0 && tempEyeDegY(w) > 0
                %Quadrant 1
                tempEyeDegX(w) = tempEyeDegX(w) + (MLStruct.SessionInfo.QuadrantCorrect(1,1) * tempEyeDegX(w) * tempEyeDegY(w));
                tempEyeDegY(w) = tempEyeDegY(w) + (MLStruct.SessionInfo.QuadrantCorrect(2,1) * tempEyeDegX(w) * tempEyeDegY(w));
                
            elseif tempEyeDegX(w) > 0 && tempEyeDegY(w) > 0
                %Quadrant 2
                tempEyeDegX(w) = tempEyeDegX(w) + (MLStruct.SessionInfo.QuadrantCorrect(1,2) * tempEyeDegX(w) * tempEyeDegY(w));
                tempEyeDegY(w) = tempEyeDegY(w) + (MLStruct.SessionInfo.QuadrantCorrect(2,2) * tempEyeDegX(w) * tempEyeDegY(w));
                
            elseif tempEyeDegX(w) < 0 && tempEyeDegY(w) < 0
                %Quadrant 3
                tempEyeDegX(w) = tempEyeDegX(w) + (MLStruct.SessionInfo.QuadrantCorrect(1,3) * tempEyeDegX(w) * tempEyeDegY(w));
                tempEyeDegY(w) = tempEyeDegY(w) + (MLStruct.SessionInfo.QuadrantCorrect(2,3) * tempEyeDegX(w) * tempEyeDegY(w));
                
            elseif tempEyeDegX(w) > 0 && tempEyeDegY(w) < 0
                %Quadrant 4
                tempEyeDegX(w) = tempEyeDegX(w) + (MLStruct.SessionInfo.QuadrantCorrect(1,4) * tempEyeDegX(w) * tempEyeDegY(w));
                tempEyeDegY(w) = tempEyeDegY(w) + (MLStruct.SessionInfo.QuadrantCorrect(2,4) * tempEyeDegX(w) * tempEyeDegY(w));
                
            else
                %Either terms = 0
                tempEyeDegX(w) = tempEyeDegX(w);
                tempEyeDegY(w) = tempEyeDegY(w);
                
            end
        end
        TempMatrix = [tempEyeDegX tempEyeDegY MLStruct.Trials.(List_Trials{k}).EyeSamples(:,3)];
        %% 3) Deal with Sampling problems - Double and missed samples
        % Sample issues are addressed by removing double samples, and then filling
        % in missing samples with a mean value of adjacent samples. This is also
        % run on BadData.
        
        
        TempMatrix(diff(TempMatrix(:,3))==0,:)=[];%remove double samples
        
        indices = find(diff(TempMatrix(:,3))>=.003); %gets the missing rows
        
        if ~isempty(indices)
            MissingRows = zeros(length(indices),3); %initializes a matrix to hold them
            
            for tLost = 1:length(indices)
                MissingRows(tLost,:) = mean([TempMatrix(indices(tLost) ,:); ...
                    TempMatrix(indices(tLost) +1,:)],1); %copies the previous sample's Location coordinates
                %and averages the eyesignal and timestamp
                
            end
        else
            MissingRows = [];
        end
        
        
        if ~isempty(MissingRows) %puts the missing row info into the output structure
            reformatted = insertrows(TempMatrix, MissingRows, indices);
        else
            reformatted = TempMatrix;
        end
        
        if ~isempty(BadData)
            %Also carried out on BadData
            BadData(diff(BadData)==0,:)=[];
            BDindices = find(diff(BadData)>=.003 & diff(BadData)<=.005); %gets the missing rows
            BDMissingRows = zeros(length(BDindices),1); %initializes a matrix to hold them
            for tLost = 1:length(BDindices)
                BDMissingRows(tLost,:) = mean([BadData(BDindices(tLost) ,:); ...
                    BadData(BDindices(tLost) +1,:)],1); %copies the previous sample's Location coordinates
                %and averages the eyesignal and timestamp
            end
            if ~isempty(BDMissingRows) %puts the missing row info into the output structure
                BadData = insertrows(BadData, BDMissingRows(:,1), BDindices);
            end
        end
        %% 4) Find periods where velocity is above physiological thresholds
        MoreBadData = reformatted(sqrt(diff(reformatted(:,1)).^2 + ...
            diff(reformatted(:,2)).^2)./diff(reformatted(:,3))>physiologicalCutOff,3);
        WorseBadData = reformatted(diff(sqrt(diff(reformatted(:,1)).^2 + ...
            diff(reformatted(:,2)).^2)./diff(reformatted(:,3)))./diff(reformatted(2:end,3))...
            >AccelPhysiologicalCutoff,3);
         % also find out if the trial starts on a saccade, and if so, NaN
        % that data as well
        StartVel = sqrt(diff(reformatted(1:40,1)).^2 + ...
            diff(reformatted(1:40,2)).^2)./diff(reformatted(1:40,3));
        if sum(StartVel(1:15)>40)>4
            fast = find([StartVel;100]>velocityThreshold);
                blocks = diff(fast);
                start = min([find(blocks>=4,1) blkEndOffset]);
                
                BadStart = reformatted(1:start,3);
        end
        %% 5) Find where the eye moves away and then back to the same spot
        % there are still some one or two sample spikes, despite the
        % algorithm executed by the EyeLink. Therefore, we use a couple of
        % measures to identify these regions. They are caused by loss of
        % corneal reflection, and sometimes it is found and lost repeatedly,
        % and so this shaky/spike region needs to be removed as well. 
        
        % Method is based on Larsson, L., Nystrom, M., & Stridh, M. (2013). 
        % Detection of saccades and postsaccadic oscillations in the presence 
        % of smooth pursuit. IEEE Transactions on Biomedical Engineering, 60(9), 
        % 2484?2493. http://doi.org/10.1109/TBME.2013.2258918
        [badSpikes] = HcTask_EyePositionSpikeChecker(reformatted);
        
   
%         figure
%         hold on
%         plot(reformatted(:,1))
%         plot(despiked(:,1))
%         pause
        %% 6) Smoothing
        % smoothing is done via a Savitzky Golay filter with order = 2 and window =
        % 11. This assumes saccades are at least 10 ms.
        % citation from: Nyström, M., & Holmqvist, K. (2010). An adaptive algorithm for fixation, 
        % saccade, and glissade detection in eyetracking data. Behavior Research Methods, 42(1), 
        % 188?204. http://doi.org/10.3758/BRM.42.1.188
        order = 2;
        window = 11;
        
        filtCleaned = sgolayfilt(reformatted(:,1:2),order,window);
        
        smoothedData = [filtCleaned,reformatted(:,3)];%outputs
        %% 7) Blink and Out Of Range removal
        % Finally, the saved timestamps in BadData are identified in groups, for
        % each section, the pre- and proceeding 20 ms are analyzed for the closest
        % index where the velocity was below 30 degrees/s. This is assumed to be
        % the last/first time where a valid, stable eyesample can be calculated.
        % The blinks and bad samples are NaN'd. This replaces previous blink
        % removal scripts. Because spikes cannot be treated this way, and
        % spikes are increased in size by smoothing, include the time
        % periods that are NaNed already in smoothed data
        AlreadyNaNd= smoothedData(isnan(smoothedData(:,1)),3);
        OffScreen = smoothedData(abs(smoothedData(:,1))>17|abs(smoothedData(:,2))>13,3);
        BadData = unique([BadData;OffScreen;BadStart;AlreadyNaNd]);
        if ~isempty(BadData) || ~isempty(MoreBadData)
            BadData = unique([BadData;MoreBadData]);
            [smoothedData,despiked] = removeSections(smoothedData,BadData,reformatted);
        else
            despiked = reformatted;
        end
        % Now we can safely despike the data
        if ~isempty(badSpikes)
            %just check that these spikes were not just part of blinks, and
            %removed by other processes, because we want to konw if this
            %was a noisy session
            fakeSpikes = ismember(badSpikes,smoothedData(isnan(smoothedData(:,1)),3));
            %remove the spikes already accounted for 
            badSpikes(fakeSpikes) = [];
            
            %remove the unaccounted for spikes
            despiked(ismember(despiked(:,3),badSpikes),1:2) = NaN;
            smoothedData(ismember(despiked(:,3),badSpikes),1:2) = NaN;
        end
        %% Plot figures for testing
        %         hFig = figure;
        %         set(hFig, 'Position', [20 300 1800 600])
        %         plot(reformatted(:,1),'r')
        %         % plotyy(reformatted(:,3),reformatted(:,1),...
        %         %     reformatted(1:end-1,3),diff(reformatted(:,1))./diff(reformatted(:,3)))
        %         hold on
        %         plot(smoothedData(:,1))
        %
        %
        %
        %         line([0 length(reformatted(:,1))], [17 17])
        %         line([0 length(reformatted(:,1))], [-17 -17])
        %         %                 line([MLStruct.Trials.(List_Trials{k}).GoalsAppearInd ...
        %         %                     MLStruct.Trials.(List_Trials{k}).GoalsAppearInd],...
        %         %                     [-30 30])
        %         ylim([-50 50])
        %         plot((diff(reformatted(:,1))./diff(reformatted(:,3)))/300)
        %         pause
        %
        %         close all
        %% 8)find the sections that are too short
        NanDiffs = diff(isnan([NaN;smoothedData(:,1)]));
        goodStart = find(NanDiffs==-1);
        goodEnd = find(NanDiffs==1);
        if~isempty(goodEnd)&& ~isempty(goodStart)
            if goodEnd(end)>goodStart(end)
                goodMat = [goodStart,goodEnd];
            else
                goodMat = [goodStart,[goodEnd;length(smoothedData)]];
            end
            goodMat(:,3) = diff(goodMat,1,2);
            goodMat(goodMat(:,3)>20,:) = [];
            
            if ~isempty(goodMat)
                for NanSection = 1:size(goodMat,1)
                    smoothedData(goodMat(NanSection,1):goodMat(NanSection,2),1:2)=NaN;
                    despiked(goodMat(NanSection,1):goodMat(NanSection,2),1:2)=NaN;
                end
            end
        end
        %
        %               subplot(2,2,1)
        %         plot(medFilteredResid)
        %         subplot(2,2,2)
        %         plot(smoothedData(:,1))
        %         subplot(2,2,3)
        %         plot(reformatted(:,1))
        %         subplot(2,2,4)
        %         plot(accResiduals(:,1))
        %         pause
        %% 9) generate position matrix
        if (PosMatrix) && ~isempty(smoothedData)
            % Loop through every eye sample to find corresponding Player position
            % and rotation in unreal space
            
            % Predefines PositionMatrix size
            PositionMatrix = NaN(size(smoothedData, 1), 6);
            
            % Creates new matrix with columns 1 - 6: PlayerPosX PlayerPosY Rotation EyePosX EyePosY  Time
            %     PositionMatrix(:, 4:6) = [Eye_PixX Eye_PixY MLStruct.Trials.(List_Trials{k}).EyeSamples(:,3)]; %Pixel Values
            PositionMatrix(:, 4:6) = smoothedData; %Degrees
            
            for m = 1 : size(PositionMatrix, 1)
                
                % Finds index of closer time point for Unreal data still < than
                % EyeSignal
                UnrealTimeIndex = find(MLStruct.Trials.(List_Trials{k}).Unreal_Times <= (PositionMatrix(m, 6)),1,'last');
                
                
                if ~isempty(UnrealTimeIndex)
                    % Adds the values at this position to the PositionMatrix
                    PositionMatrix(m, 1) = MLStruct.Trials.(List_Trials{k}).Unreal_PosXPosYRot(UnrealTimeIndex,1);
                    PositionMatrix(m, 2) = MLStruct.Trials.(List_Trials{k}).Unreal_PosXPosYRot(UnrealTimeIndex,2);
                    PositionMatrix(m, 3) = MLStruct.Trials.(List_Trials{k}).Unreal_PosXPosYRot(UnrealTimeIndex,3);
                end
            end
            
            
            MLStruct.Trials.(List_Trials{k}).PositionMatrix_deg = PositionMatrix;
            MLStruct.Trials.(List_Trials{k}).NonSmoothedEyes = despiked;
        else
            
            MLStruct.Trials.(List_Trials{k}).EyeDegrees = smoothedData;
            MLStruct.Trials.(List_Trials{k}).NonSmoothedEyes = despiked;

        end
        
        waitbar(k/size(List_Trials,1));
        % Compare the number of spikes found in this trial with the most
        % found so far in the session
        if length(badSpikes)>5
            MaxSpikes{k} = length(find(diff(badSpikes)>.003))+1;
            MaxSpikes{k,2} = length(find(diff(badSpikes)>.08))+1;
        end
    end
end
%if there are a bunch of bad spikes from corneal loss, note that in
%sessionInfo
    if max([MaxSpikes{:,2}])>10 
        MLStruct.SessionInfo.CornealLoss = 1;
    end
    close(hWaitBar);
%     figure
%     hold on
%     plot(reformatted(:,1))
%     plot(despiked(:,1))
%     plot(reformatted(:,2))
%     plot(despiked(:,2))
