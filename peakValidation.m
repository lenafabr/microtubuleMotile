function [ validPeaks] = peakValidation(Th_rel_Intensity, Th_Intensity, Th_Angle, TH_AngleDiff, Candidate_Peaks, nPeaks, prev_angle)

if (Candidate_Peaks == 0) %never occurs
    validPeaks = [0;0];
else
    
    % threshold peaks (maybe normalize values)
    Intensity_validated_anglePeaks = Candidate_Peaks(:,(Candidate_Peaks(2,:)>Th_Intensity));
    %     IMPLEMENT NORMALIZED PEAKS
    if (isempty(Intensity_validated_anglePeaks))
%         disp('peaks too weak, done with path!');
        validPeaks = [0;0];
        return;
    end
    %% Th check on relative intensity
    norm_Peaks = Intensity_validated_anglePeaks(2,:)/max(Intensity_validated_anglePeaks(2,:));
    Intensity_validated_anglePeaks = Intensity_validated_anglePeaks(:,(abs(norm_Peaks)>Th_rel_Intensity));
    %% Th_angle check to remove large incremental angles 
    n=size(Intensity_validated_anglePeaks,2);
    good_angle_index=[];
    ii=1;
    for i=1:n
        delta_Angle = Intensity_validated_anglePeaks(1,i) - prev_angle;
        if (abs(delta_Angle)<Th_Angle) || (prev_angle==0)
            good_angle_index(ii)=i;
            ii=ii+1;
        end
    end
    if isempty(good_angle_index)
        validPeaks = [0;0];
        return;
    end
    if (good_angle_index(1)==0)
        validPeaks = [0;0];
        return;
    end
    angle_validated_Peaks = Intensity_validated_anglePeaks(:,good_angle_index);
%     angle_validated_Peaks = Intensity_validated_anglePeaks;
    
    
    [~,n] = size(angle_validated_Peaks);
    if isempty(angle_validated_Peaks)
        validPeaks=[0 0]';
        return
    end
    %% filter out similar angle peaks
    
    valToRemove = zeros;
    k=1;
    [m,n] = size(angle_validated_Peaks);
    if (n > 1)
        for i = 1:length(angle_validated_Peaks)-1
            diff = abs(angle_validated_Peaks(1,i) - angle_validated_Peaks(1,i+1));
            if diff < TH_AngleDiff
                minBrightness = min(angle_validated_Peaks(2,i), angle_validated_Peaks(2,i+1));
                valToRemove(1,k) = minBrightness;
                k = k+1;
            end
        end
        for i = 1:length(valToRemove)
            aIdx = find(angle_validated_Peaks(2,:) == valToRemove(i));
            angle_validated_Peaks(:,aIdx) = [];
        end
    end    
    
    %only keep nPeaks values 
    if (n > 0)   
        [m,n] = size(angle_validated_Peaks);
        
        if nPeaks > n
            nPeaks = n;
        end
        
        [d1,d2] = sort(angle_validated_Peaks(2,:));
        sortedAnglePeaks = angle_validated_Peaks(:,d2);
        
        %%sorted for brightness
        angle_validated_Peaks = sortedAnglePeaks(:,1:nPeaks);
        
        [Y,I] = sort(angle_validated_Peaks(1,:));
        validPeaks = angle_validated_Peaks(:,I);
    else
        validPeaks = [0;0];
    end
end


