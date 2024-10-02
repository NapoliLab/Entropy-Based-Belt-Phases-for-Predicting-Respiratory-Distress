function [zeroCrossingsCleaned, finalPeakValVec, finalPeakLocVec] = peak_detector(resp_belt, Fs, P)
%===================================
% Function Details 
%-----------------------------
% This function finds attributes of the respiratory belt transducer signal, specifically
% the zero crossing points. The algorithm also provides the peak or trough instance of the respiratory belt waveform. 
% 
% % Assumptions: 
% Assumption A) This algorithm assumes no breaths are faster than 30
% breaths per minute (30 bpm)
%=====================================
%     Inputs
%------------------------------------
% resp_belt: A vector 1 by N that contains the respiratory belt data
% Fs: Sampling frequency of the data
% P: A binary value that is 0 or 1. 0 indicates no plotting and 1 indicates plotting the figure.
%=====================================
% Outputs 
%-------------------------------------
% zeroCrossingsCleaned: A cleaned list of zero crossing points.
% finalPeakValVec: A vector of final peak values.
% finalPeakLocVec: A vector of final peak locations.
%
%======================================

%----------------------------------------------------
% Peak Instances for Expiratory and Inspiratory Times  
%----------------------------------------------------
% Statistical Assumed Constraints on the Breathing 
%-------------------------------------------------
% Min Peak Width 
BreathingRate = 30; 
MinTime = round(((BreathingRate / 60) / 2) * Fs); % Note divided by 2 b/c we are only looking at one phase of the breath

% Find Peaks 
[InspPeak, InspPeakLoc] = findpeaks(resp_belt, 'MinPeakWidth', MinTime);

InspPeakLoc = InspPeakLoc(InspPeak > 0);
InspPeak = InspPeak(InspPeak > 0);

[ExpPeak, ExpPeakLoc] = findpeaks(-1 * resp_belt, 'MinPeakWidth', MinTime);
ExpPeak = -1 * ExpPeak;

ExpPeakLoc = ExpPeakLoc(ExpPeak < 0);
ExpPeak = ExpPeak(ExpPeak < 0);

%-----------------
% Zero Crossings  
%-----------------
crossings = [];
for i = 2:length(resp_belt)
    if (resp_belt(i) >= 0 && resp_belt(i - 1) < 0) || (resp_belt(i) <= 0 && resp_belt(i - 1) > 0)
        crossings = [crossings; i];
    end
end

% Combine the two vectors and sort by index
combineLocs = [InspPeakLoc; ExpPeakLoc];
[sortedLocs, I_sortedLocs] = sort(combineLocs);

combinePeaks = [InspPeak; ExpPeak];
combinePeaksSorted = combinePeaks(I_sortedLocs);

% Iterate through each of the values in the combinePeaksSorted vector ensuring the values alternate between positive and negative
newPeakValVec = combinePeaksSorted(1);
newPeakLocVec = sortedLocs(1);

for i = 2:length(combinePeaksSorted)
    if (newPeakValVec(end) * combinePeaksSorted(i)) < 0
        newPeakValVec = [newPeakValVec; combinePeaksSorted(i)];
        newPeakLocVec = [newPeakLocVec; sortedLocs(i)];
    else
        if abs(newPeakValVec(end)) > abs(combinePeaksSorted(i))
            newPeakValVec(end) = newPeakValVec(end);
            newPeakLocVec(end) = newPeakLocVec(end);
        else
            newPeakValVec(end) = combinePeaksSorted(i);
            newPeakLocVec(end) = sortedLocs(i);
        end
    end
end

finalPeakValVec = newPeakValVec(newPeakValVec ~= 0);
finalPeakLocVec = newPeakLocVec(newPeakValVec ~= 0);

firstZeroCross = crossings(1);
lastZeroCross = crossings(end);
mask = (finalPeakLocVec > firstZeroCross) & (finalPeakLocVec < lastZeroCross);

finalPeakValVec = finalPeakValVec(mask);
finalPeakLocVec = finalPeakLocVec(mask);

%-----------------------------------------------------
% Linking each peak/trough to zero crossings
%-----------------------------------------------------
PeaksTroughs = finalPeakLocVec;

zeroCrossingsCleaned = [];
for thisPeakTrough = 1:length(PeaksTroughs)
    % For each peak/trough, find the nearest two zero crossings
    diffVals = crossings - PeaksTroughs(thisPeakTrough);

    % Find the last positive and first negative indices
    lastPositiveIndex = find(diffVals < 0, 1, 'last');
    firstNegativeIndex = find(diffVals > 0, 1, 'first');

    zeroCrossingsCleaned = [zeroCrossingsCleaned crossings(lastPositiveIndex) crossings(firstNegativeIndex)];
end

% Sort the remaining zero crossings
zeroCrossingsCleaned = sort(zeroCrossingsCleaned);

% Reshape into a 2 by N array where each pair is a zero crossing pair with a trough/peak in the breath
zeroCrossingsCleaned = reshape(zeroCrossingsCleaned, [2, length(zeroCrossingsCleaned) / 2]);

for i = 1:length(zeroCrossingsCleaned) - 1
    if zeroCrossingsCleaned(2, i) ~= zeroCrossingsCleaned(1, i + 1)
        zeroCrossingsCleaned(2, i) = NaN;
    end
end

zeroCrossingsCleaned = unique(zeroCrossingsCleaned(:)');
zeroCrossingsCleaned = zeroCrossingsCleaned(~isnan(zeroCrossingsCleaned));

%--------------------------------------------------------
% Plotting Option for Data Visualization
%---------------------------------------------------------
if P == 1
    figure;
    plot(1:length(resp_belt), resp_belt);
    hold on;
    plot(zeroCrossingsCleaned, resp_belt(zeroCrossingsCleaned), 'ro', 'MarkerSize', 8);
    hold on;
    plot(InspPeakLoc, InspPeak, 'o');
    hold on;
    plot(ExpPeakLoc, ExpPeak, 'o');
    plot(finalPeakLocVec, finalPeakValVec, 'o');
    title('Flow with Zero Crossings Marked');
    xlabel('Sample Index');
    ylabel('Flow Value');
    legend('Flow', 'Zero Crossings');
    hold off;
end 
end
