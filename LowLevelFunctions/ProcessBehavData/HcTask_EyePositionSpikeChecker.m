function [badSpikes] = HcTask_EyePositionSpikeChecker(reformatted)
% there are still some one or two sample spikes, despite the
% algorithm executed by the EyeLink. Therefore, we use a couple of
% measures to identify these regions. They are caused by loss of
% corneal reflection, and sometimes it is found and lost repeatedly,
% and so this shaky/spike region needs to be removed as well. In
% our data, there are periods where the velocity is still below the
% physiological threshold, but the variable is not properly
% smoothed.
% Therefore, we find the spikes that are of at least .3 degrees,
% and then check whether they are over two or four samples, and
% then check the velocity before the spike, so that we are not
% removing PSOs

% Method is based on Larsson, L., Nystrom, M., & Stridh, M. (2013). 
        % Detection of saccades and postsaccadic oscillations in the presence 
        % of smooth pursuit. IEEE Transactions on Biomedical Engineering, 60(9), 
        % 2484?2493. http://doi.org/10.1109/TBME.2013.2258918
spikesExist=1;
badSpikes = [];
while spikesExist ==1
    wind = 20*[ -1 -1 -1 -1 0 1 1 1 1];
    badSpikes = [];
    
    Vel = [0,0;diff(reformatted(:,1:2))./...
        [diff(reformatted(:,3)),...
        diff(reformatted(:,3))]];
    Acc = abs([conv(Vel(:,1),wind,'same'),conv(Vel(:,2),wind,'same')]);
    Threshdiff = [3,3];
    PeakThresh = [20000,20000];
    for var = 1:2
        
        while Threshdiff(var)>1
            MeanVel(var) = mean(Acc(Acc(:,var)<PeakThresh(var),var));
            STD(var) = std(Acc(Acc(:,var)<PeakThresh(var),var));
            NewThresh = MeanVel(var)+(4*STD(var));
            Threshdiff(var) = PeakThresh(var)-NewThresh;
            PeakThresh(var) = NewThresh;
            %             count(var) = count(var)+1;
            
        end %calculate Threshold
        
    end
    
    %Take the periods above the threshold, and test the endpoints vs
    % the distance travelled
    PeakEyeVelThreshold = PeakThresh;
    OnsetThresh = MeanVel+(3*STD);
    %define the portions that are above the threshold for each component
    
    OnsetThresh = MeanVel+(3*STD);
    %define the portions that are above the threshold for each component
    isSpike = Acc(:,1)>PeakThresh(1)|Acc(:,2)>PeakThresh(2);
    if isempty(isSpike)
        spikesExist = 0;
        continue
    end
    if diff(isSpike(1:2))==-1||isSpike(1)==1
        %NonSpikeStart is all the points where isSpike changes from 1 to 0, but
        %we have to check that it doesn't do this at the beginning
        NonSpikeStart = find(diff(isSpike)==-1);
    else     NonSpikeStart = [1;find(diff(isSpike)==-1)];
    end
    %NonSpikeEnd is at the points where a sacade starts, so where isSpike
    %starts to go to 1
    NonSpikeEnd = [find(diff(isSpike)==1);length(isSpike)];
    %make a matrix and find the length of the subthreshold periods
    if isempty(NonSpikeStart)||isempty(NonSpikeEnd)
        spikesExist = 0;
        continue
    end
    d = [NonSpikeStart NonSpikeEnd(1:length(NonSpikeStart))];
    d(:,3) = diff(d,1,2);
    
    if length(NonSpikeStart)<2||length(NonSpikeEnd)<2
        spikesExist = 0;
        continue
    end
    e = [NonSpikeEnd(1:length(NonSpikeStart)-1) NonSpikeStart(2:end)];
    e(:,3) = diff(e,1,2);
    %get the short spike times
    SpikeTimes = find(e(:,3)<5);
    if isempty(SpikeTimes)
        spikesExist = 0;
        continue
    end
    %use these as the start and end of the bad data
    
    SpikeEnd = e(SpikeTimes,2);
    SpikeStart = e(SpikeTimes,1);
        
    % Calculate the sample to sample change, and the net displacement for each
    % movement
    
    for spike = 1:length(SpikeStart)
        %check if there are any displacememnts greater than .3degrees
        if any(sqrt(diff(reformatted(SpikeStart(spike):SpikeEnd(spike),1)).^2 + ...
                diff(reformatted(SpikeStart(spike):SpikeEnd(spike),2)).^2) >.3)
            %if yes, could be a spike, check the total displacement
            totalDisplacement = sqrt((reformatted(SpikeEnd(spike),1)-reformatted(SpikeStart(spike),1)).^2 +...
                (reformatted(SpikeEnd(spike),1)-reformatted(SpikeStart(spike),1)).^2);
            if totalDisplacement <.3
                % if the displacement is smaller, we have a likely spike. now just
                % check that the velocity is less than the preceeding velocities,
                % to ensure that this is not simply a post saccadic oscillation
                velSample = max(sqrt(diff(reformatted(SpikeStart(spike):SpikeEnd(spike),1)).^2 + ...
                    diff(reformatted(SpikeStart(spike):SpikeEnd(spike),2)).^2)...
                    ./diff(reformatted(SpikeStart(spike):SpikeEnd(spike),3)));
                firstIndex = max(SpikeStart(spike)-5,1);
                previousSamples = max(sqrt(diff(reformatted(firstIndex:SpikeStart(spike),1)).^2 + ...
                    diff(reformatted(firstIndex:SpikeStart(spike),2)).^2)...
                    ./diff(reformatted(firstIndex:SpikeStart(spike),3)));
                if velSample>previousSamples
                    badSpikes = [badSpikes;reformatted(SpikeStart(spike):SpikeEnd(spike),3)];
                end
            end
        end
    end
    spikesExist=0;
end
%% previous spike removal
% 
%         %define the portions that are above the threshold for each component
%         xdiff = [0;diff(reformatted(:,1))];
% %         xreversals = [0;(xdiff(1:end-1).*xdiff(2:end))<0];
%         xampP = xdiff>.3;
%         xampN = xdiff<-.3;
%         
%         maybeSpikes = xampP & [0;xampN(1:end-1)] | xampP & [0;0;xampN(1:end-2)]...
%             | xampP & [xampN(2:end);0] | xampP & [xampN(3:end);0;0];
%         
%        
%         spikes = maybeSpikes;%(Vel(maybeSpikes-1)<Vel(maybeSpikes))
%         
%         
%         HighAccelStart = [find(diff(spikes)==1)];
%         HighAccelEnd = [find(diff(spikes)==-1);length(reformatted)];
%         %make a matrix and find the length of the subthreshold periods
%         d = [HighAccelStart HighAccelEnd(1:length(HighAccelStart))];
%         
%         if~isempty(d)
%         for fastSec = 1:size(d,1)
%          
%             % for each high accel region, check whether a short period or a
%             % long period has a smaller total difference, and then check
%             % that it has a lower velocity than the samples before it(so
%             % that this doesn't catch Post saccadic oscillations)
% %             s2s(fastSec) = sum(degreesTravelled(HighAccelStart(fastSec):...
% %                 HighAccelEnd(fastSec)));
%             shortSample = sqrt((reformatted(HighAccelStart(fastSec),1)-...
%                 reformatted(HighAccelEnd(fastSec),1)).^2+...
%                 (reformatted(HighAccelStart(fastSec),2)-...
%                 reformatted(HighAccelEnd(fastSec),2)).^2);
%             longSample = sqrt((reformatted(max(HighAccelStart(fastSec)-1,1),1)-...
%                 reformatted(min(HighAccelEnd(fastSec)+1,length(reformatted)),1)).^2+...
%                 (reformatted(max(HighAccelStart(fastSec)-1,1),2)-...
%                 reformatted(min(HighAccelEnd(fastSec)+1,length(reformatted)),2)).^2);
% %             %calculate the smaller of the total displacements for the
% %             %potential spike samples
% %             totalDeg(fastSec) = min(shortSample,longSample);
% %             ratios(fastSec) = ((s2s(fastSec)*100)/(totalDeg(fastSec)*100))...
% %                 /(HighAccelEnd(fastSec)-HighAccelStart(fastSec));
% %             %if the data has a lot more inter sample variation per time
% %             %period than it does end to end difference, and the
% %             %variability is moer than one degree, its a spike and
% %             %should be removed
%             % if ratios(fastSec)>1 && s2s(fastSec)>1
%             % check which sample has less overall displacement, and then
%             % compare it's velocity to the velocity of the preceeding 5
%             % samples, if it is faster than the previous samples, it is a
%             % spike
%             if shortSample<longSample
%                 if HighAccelStart(fastSec)<5 || max(diff(reformatted(HighAccelStart(fastSec):...
%                     HighAccelEnd(fastSec))))>max(diff(reformatted(max(HighAccelStart(fastSec)-5,1):...
%                     HighAccelStart(fastSec)))) 
%                 badSpikes = [badSpikes;reformatted(HighAccelStart(fastSec):...
%                     HighAccelEnd(fastSec),3)];
%                 end
%             else
%                 if HighAccelStart(fastSec)<6 || max(diff(reformatted(max(HighAccelStart(fastSec)-1,1):...
%                     HighAccelEnd(fastSec)+1)))>max(diff(reformatted(max(HighAccelStart(fastSec)-6,1):...
%                     HighAccelStart(fastSec)-1))) 
% 
%                     badSpikes = [badSpikes;reformatted(max(HighAccelStart(fastSec)-1,1):...
%                     min(HighAccelEnd(fastSec)+1,length(reformatted)),3)];
%                 end
%             end
%             %end
%         end
%         %save the badSpikes timestamps, and NaN after the bad sections have
%         %been removed
%         end