function [MLStruct] = Remove1000hzsampling(MLStruct,PosMatrix)
TrlNames = fieldnames(MLStruct.Trials);
for trl = 1:length(TrlNames)
    if PosMatrix ==1
        matrix = MLStruct.Trials.(TrlNames{trl}).PositionMatrix_deg;
        TimeIndex = 6;
    else
        matrix = MLStruct.Trials.(TrlNames{trl}).EyeDegrees;
        TimeIndex = 3;
    end
    % find if there are 1000hz samples
    sampleIntervals = diff(matrix(:,TimeIndex));
    FastSamples = find(sampleIntervals<.0015);
    % remove the samples, this will remove the even ones, may end up with a
    % difference of a couple ms between indices and trial time, but this is
    % because the sampling can switch from even time to odd time stamps
    % once it has the pupil.
    if ~isempty(FastSamples)
        length(FastSamples)
        evenSamples = FastSamples(find(mod(FastSamples,2)));
        matrix(evenSamples,:) = [];
    end
    
    if PosMatrix ==1
        MLStruct.Trials.(TrlNames{trl}).PositionMatrix_deg = matrix ;
        
    else
        MLStruct.Trials.(TrlNames{trl}).EyeDegrees = matrix;
        
    end
end
end