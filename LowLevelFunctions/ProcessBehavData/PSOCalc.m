function [PSOendpoint,REnd] = PSOCalc(Input)
%takes a section of the post saccaddic period( X or Y) and finds if there
%is an oscillation in it, and returns the time of start and finish.
PPSO = Input(:,1);
Time = Input(:,2);
%initialize the variables
isPSO = 1; %changed if the section is not a PSO, or to end the algorithm
ThreshAngle = .08; %the minimum difference between the modeled signal and the decay
ender = min(20,length(PPSO)-7); %intially, set endpoint at 40 ms, or 20 samples, may be extended later
GNended = 0; %the variable for the while loop that checks if an end is found for GN
Vender = ender; %the variable ender that moves backwards until a linear signal is found in GN
PSOendpoint = Time(1);
REnd = 'good';
removed = 0; %logs any changes in the beginning of the signal (such as cutting it to make the model fit)
%check that the signal has evened out - so that the slope before and after
%have the same sign
while isPSO==1
    TLrefslope = polyfit([1:5]'*.002,PPSO(ender:ender+4),1);
    if ender-7>0 %check that ender is long enough to have a slope, if it is
        % too short, there must be not enough samples, so don't check the
        % slope, just keep ender there
        TLtestslope = polyfit([1:5]'*.002,PPSO(ender-4:ender),1);
        %check if there is any NaN in the sample
        if any(isnan(PPSO(1:ender)))
            isPSO = 0;
            continue
        end
        
        %check if we have to extend the sample( if there are opposite signs of
        %the slopes of the lines)
        if TLtestslope(1)*TLrefslope(1)<0
            ender = min(ender+10,length(PPSO));
        end
    end
    %check if there is any NaN in the sample
    if any(isnan(PPSO(1:ender)))
        isPSO = 0;
        continue
    end
    %Go through and find out where the linear part starts
    Vender = ender;
    while GNended == 0
        GNrefslope = polyfit([1:ender-Vender+3]'*.002, PPSO(Vender-2:ender),1)';
        GNtestslope = polyfit([1:2]'*.002, PPSO(Vender-3:Vender-2),1)';
        [Vender abs(GNrefslope(1)-GNtestslope(1))];%>30;
        if abs(GNrefslope(1)-GNtestslope(1))>17
            GNended = 1;
        else
            Vender = Vender-1;
        end
        if Vender<=5
            GNended = 1;
            isPSO = 0;
            continue
        end
        
    end
    GN = ones(ender,1)*PPSO(Vender-1);
    GN(1:Vender-1) = PPSO(1:Vender-1);
    GN = GN-PPSO(Vender-1);
    % check if there is anything left of the signal
    if isempty(find(GN,1))
        isPSO = 0;
        continue
    end
    %Check if there is enough of an amplitude to warrant a PSO
    if max(abs(GN))<.15
        isPSO = 0;
        continue
    end
    %% initialize variables to hold the models (G) and the RMSE
    
    
    adequateRMSE = 0; % variable to stop while loop
    
    while adequateRMSE ==0;
        GI = zeros(length(GN),4); % initialize the matrix to hold the impulse decay patterns
        RMSE = [1,1,1,1];
        for test = 1:4
            [num{test},denom{test}] = prony(GN,0,test);
            GI(:,test) = impz(num{test},denom{test},length(GN));
            RMSE(test) = (sqrt(sum((GN(:)-GI(:,test)).^2)/numel(GN)))/(max(abs(GN)));
            %     plot(impz(num{test},denom{test},length(GN)));
            
        end
%                 figure
%                 hold on
%                 plot(GN)
%                 plot(GI)
        %if there is at least one RMSE that is below .15, then continue,
        %otherwise, remove one sample from the beginning of the signal.
        if any(RMSE<.15) && [any(any(abs(GI)>.2)) || max(abs(GN))<.15]
            %we have a good fit, now check the poles, so that it isn't
            %something that is not decaying
            best = 1;
            for p = 2:4
                change = RMSE(best)-RMSE(p);
                if change>.05
                    best =  p;
                end
            end
            [~,Poles] = prony(GN,0,best);
            Poles = Poles(2:end);
            RMax = max(Poles);
            if RMax <.89
                adequateRMSE = 1;
            else
                GN = GN(2:end);
                removed = removed +1;
                % make sure there is still relevant signal
                if max(abs(GN))<.15
                    adequateRMSE = 1;
                end
            end
        else
            GN = GN(2:end);
            removed = removed +1;
            % make sure there is still relevant signal
            if max(abs(GN))<.15
                adequateRMSE = 1;
            end
                
        end
    end
    % make sure there is still relevant signal
    if max(abs(GN))<.15
        isPSO = 0;
        REnd = 'lowAmp';
        continue
    end
    %compare the RMSE
    best = 1;
    for p = 2:4
        change = RMSE(best)-RMSE(p);
        if change>.05
            best =  p;
        end
    end
    [~,Poles] = prony(GN,0,best);
    Poles = Poles(2:end);
    RMax = max(Poles);
    % Check that the amplitude is large enough to be a PSO
    Amp = max(abs(GI(:,best)));
    if Amp<.15
        isPSO = 0;
        REnd = 'lowAmp';
        continue
    end
%     %make sure that the pole is less than .89
%     if RMax>.89
%         isPSO = 0;
%         REnd = 'highPole'
%         continue
%     end
    %calculate the decay,f(n) = Amp(RMax^n) which is the amplitude times (RMax raised to the
    %power of the postion being calculated) starting at 0
    decay = ones(length(GI(:,best))-1,1);
    for nn = 0:length(GI(:,best))-1
        decay(nn+1) = max(abs(GI(:,best)))*(abs(RMax).^nn);
    end
    % to calculate the difference between the signal and the decay, the
    % decay is subtracted from the
    %get indices where the model is below zero
    negative = GI(:,best)<0;
    % subtract the difference between the decay and the model (
    signal = GI(:,best)+decay;
    negsignal = GI(:,best)-decay;
    signal(negative) = negsignal(negative);
    %gets all the points that are outside the decaying signal by at least
    %.08
    belowThresh = [0;abs(signal)<ThreshAngle;0];
    changes = find(belowThresh>0);
    %finds the point where the signal was below the threshold for more than
    %6 ms, that is the endpoint
    DipsBelowSignal =[];
    if~isempty(changes)
        DipsBelowSignal = max(find(diff(changes)>3),changes(end));
    end
    AlwaysBelowSignal = find(belowThresh==0);
    AlwaysBelowSignal = AlwaysBelowSignal(end-1);
    if isempty(DipsBelowSignal)
        endpoint = AlwaysBelowSignal;
    else
        
        endpoint = min(DipsBelowSignal,AlwaysBelowSignal);
    end
    endpoint = removed + endpoint(1) ;
    if isempty(endpoint);
        endpoint = length(GN);
    end
    endpoint = min(length(Time),endpoint);
    %     PSOendpoint = Time(endpoint);
    %      isPSO = 2;
    %     hold on
    %     plot(GI(:,best))
    %Now calculate the ratio of the amplitudes, and make sure that it is
    %not just a slow movement, by adding up the max positive and negative
    %amplitudes, and dividing by the time.
    S =  (max(PPSO) - min(PPSO)) ./ (endpoint*.002);
%     figure
%     hold on
%     plot(GN)
% plot(decay)
%     plot(GI(:,best))
    
    if S< 15
        isPSO = 2;
        continue
    end
    PSOendpoint = Time(endpoint);
    isPSO = 2;
    
    %     %if S is less than 400/2*Sampling Frequency, (in this case it is probably a slow
    %     %movement
    %     200/500
    %
end
a = 1;
end


