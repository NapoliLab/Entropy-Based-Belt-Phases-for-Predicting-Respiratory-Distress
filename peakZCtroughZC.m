% Closes all open figures, clears all variables, clears command window
clearvars;
close all;
clc;

% Loads the data as well as some coefficients for filters used later
load('PulseTrainData_HypoxiaVsSeaLevel.mat')

% Initialize structures for non-hypoxic and hypoxic data by participants
% TI - Times of Inspiration split by zero crossing. Values stored in
% chronological order [Ia, Ib, Ia, Ib, ...]
% TE - Times of Expiration split by zero crossing. Values stored in
% chronological order [Ea, Eb, Ea, Eb, ...]
% TT - vector of values stored in chronological order. Combination of TI
% and TE [Ia, Ib, Ea, Eb, Ia, Ib ...]
% TB - Total times of Inspiration and Expiration
PZTZresults = struct('TI', cell(2,49), 'TE', cell(2,49), 'TT', cell(2,49), 'TB', cell(2,49));
time_interval_names = {'fullLength', 'first5', 'last5', 'last2'};

for subnum = 1:49
    for condition = 1:2 % 1 for HYPOXIC, 2 for NONHYPOXIC
        if condition == 1
            tempData = CardRespDataH{subnum};
        else
            tempData = CardRespDataNH{subnum};
        end

        t = 0:1/256:length(tempData)/256 - 1/256;
        bre = tempData(:,3);
        smoothed_bre = smooth(bre, 256);
        
        %separate the signal into different time intervals
        %1. Full Signal Length
        %2. First 5 minutes
        %3. Last 5 minutes
        %4. Last 2 minutes
        intervals = calculate_intervals(length(smoothed_bre));
        [zeroCrossingsCleaned, PeakTroughVals, PeakTroughLocs] = peak_detector(smoothed_bre, 256, 0);
        [combinedIdxVec, I] = sort([PeakTroughLocs,zeroCrossingsCleaned]);
        
        % Analysis over different time intervals
        for time_interval = 1:size(intervals, 1)
            PZTZresults = analyze_intervals(t, zeroCrossingsCleaned, PeakTroughLocs, PeakTroughVals, combinedIdxVec, intervals(time_interval, :), time_interval, time_interval_names, PZTZresults, condition, subnum);
        end
    end
end

save("PZTZresults.mat","PZTZresults");


%% SUPPORTING FUNCTIONS

function intervals = calculate_intervals(length_sm)
    % Calculate the intervals for analysis
    mid = ceil(length_sm/2);
    last_2 = length_sm - ceil(length_sm/5);
    intervals = [1, length_sm; 1, mid; mid+1, length_sm; last_2+1, length_sm];
end

function results = analyze_intervals(t, zeroCrossings, PeakTroughLocs,PeakTroughVals, combinedIdxVec, interval, time_interval_iteration, time_interval_names, results, condition, subnum)
        indicator = combinedIdxVec > interval(1) & combinedIdxVec < interval(2);
        indices = find(indicator);

        %check whether it is ZC or PT val in order to modify the
        %starting_index
        if any(zeroCrossings == combinedIdxVec(indices(1)))
            indices = indices(2:end);
        end

        if any(zeroCrossings == combinedIdxVec(indices(end)))
            indices = indices(1:(end-1));
        end

        %find the index for the first peak
        if any(PeakTroughLocs == combinedIdxVec(indices(1)))
            idx_peak = find(PeakTroughLocs == combinedIdxVec(indices(1)));
            firstPeakVal = PeakTroughVals(idx_peak);
        end
    
    for i = 1:length(indices)-4
        if mod(i,4) == 1 
            Ia = t(combinedIdxVec(indices(i+1))) - t(combinedIdxVec(indices(i))); 
            Ib = t(combinedIdxVec(indices(i+2))) - t(combinedIdxVec(indices(i+1)));
            Ea = t(combinedIdxVec(indices(i+3))) - t(combinedIdxVec(indices(i+2)));
            Eb = t(combinedIdxVec(indices(i+4))) - t(combinedIdxVec(indices(i+3)));
            time_interval_field = time_interval_names{time_interval_iteration};

            % Check if the field exists, if not initialize it
            if ~isfield(results(condition, subnum).TT, time_interval_field)
                results(condition, subnum).TT.(time_interval_field) = [];
                results(condition, subnum).TB.(time_interval_field) = [];
                results(condition, subnum).TI.(time_interval_field) = [];
                results(condition, subnum).TE.(time_interval_field) = [];
            end

            results(condition, subnum).TB.(time_interval_field) = [results(condition, subnum).TB.(time_interval_field), Ia+Ib+Ea+Eb];
            results(condition, subnum).TT.(time_interval_field) = [results(condition, subnum).TT.(time_interval_field), Ia, Ib, Ea, Eb];
            if firstPeakVal < 0
                results(condition, subnum).TI.(time_interval_field) = [results(condition, subnum).TI.(time_interval_field), Ia, Ib];
                results(condition, subnum).TE.(time_interval_field) = [results(condition, subnum).TE.(time_interval_field), Ea, Eb];
            else
                results(condition, subnum).TI.(time_interval_field) = [results(condition, subnum).TI.(time_interval_field), Ea, Eb];
                results(condition, subnum).TE.(time_interval_field) = [results(condition, subnum).TE.(time_interval_field), Ia, Ib];
            end
        end
    end
end

