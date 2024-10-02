% Closes all open figures, clears all variables, clears command window
clearvars;
close all;
clc;

% Loads the data as well as some coefficients for filters used later
load('PulseTrainData_HypoxiaVsSeaLevel.mat')

% Initialize structures for non-hypoxic and hypoxic data by participants
% InterB_pos - The positive zero crossing to zero crossing interval
% InterB_neg - The negative zero crossing to zero crossing interval
% TT - vector of the InterB_pos and InterB_neg values stored in
% chronological order
% TB - vector of the addition of the InterB_pos and and InterB_neg values
% O2 - vector of the Oxygen Saturation during total breath (TB) windows
ZtZresults = struct('InterB_pos', cell(2,49), 'InterB_neg', cell(2,49), 'TT', cell(2,49), 'TB', cell(2,49), 'O2', cell(2,49));
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
        O2 = tempData(:,5);
        smoothed_bre = smooth(bre, 256);
        
        %separate the signal into different time intervals
        %1. Full Signal Length
        %2. First 5 minutes
        %3. Last 5 minutes
        %4. Last 2 minutes
        intervals = calculate_intervals(length(smoothed_bre));
        [zeroCrossingsCleaned, PeakTroughVals, PeakTroughLocs] = peak_detector(smoothed_bre, 256, 0);
        
        % Analysis over different time intervals
        for time_interval = 1:size(intervals, 1)
            ZtZresults = analyze_intervals(t, O2, zeroCrossingsCleaned, PeakTroughVals, PeakTroughLocs, intervals(time_interval, :), time_interval, time_interval_names, ZtZresults, condition, subnum);
        end
    end
end

save("ZtZresults.mat","ZtZresults");


%% SUPPORTING FUNCTIONS

function intervals = calculate_intervals(length_sm)
    % Calculate the intervals for analysis
    mid = ceil(length_sm/2);
    last_2 = length_sm - ceil(length_sm/5);
    intervals = [1, length_sm; 1, mid; mid+1, length_sm; last_2+1, length_sm];
end

function results = analyze_intervals(t, O2, zeroCrossings, PeakTroughVals, PeakTroughLocs, interval, time_interval_iteration, time_interval_names, results, condition, subnum)
    indicator = zeroCrossings > interval(1) & zeroCrossings < interval(2);
    indices = find(indicator);

    %find the first peak/trough value to understand orientation of breath
    indicator_peak = PeakTroughLocs > interval(1) & PeakTroughLocs < interval(2);
    idx_peak = find(indicator_peak);
    
    for i = 1:length(indices)-2
        if mod(i,2) == 1 
            O2_window = O2(zeroCrossings(indices(i)):zeroCrossings(indices(i+2)));
            pos_interval = t(zeroCrossings(indices(i+1))) - t(zeroCrossings(indices(i)));
            neg_interval = t(zeroCrossings(indices(i+2))) - t(zeroCrossings(indices(i+1)));
            time_interval_field = time_interval_names{time_interval_iteration};

            % Check if the field exists, if not initialize it
            if ~isfield(results(condition, subnum).TT, time_interval_field)
                results(condition, subnum).TT.(time_interval_field) = [];
                results(condition, subnum).TB.(time_interval_field) = [];
                results(condition, subnum).InterB_pos.(time_interval_field) = [];
                results(condition, subnum).InterB_neg.(time_interval_field) = [];
                results(condition, subnum).O2.(time_interval_field) = [];
            end

            results(condition, subnum).TT.(time_interval_field) = [results(condition, subnum).TT.(time_interval_field), pos_interval, neg_interval];
            results(condition, subnum).TB.(time_interval_field) = [results(condition, subnum).TB.(time_interval_field), pos_interval + neg_interval];
            results(condition, subnum).O2.(time_interval_field) = [results(condition, subnum).O2.(time_interval_field), mean(O2_window(O2_window~=0))];
            if PeakTroughVals(idx_peak(1)) > 0
                results(condition, subnum).InterB_pos.(time_interval_field) = [results(condition, subnum).InterB_pos.(time_interval_field), pos_interval];
                results(condition, subnum).InterB_neg.(time_interval_field) = [results(condition, subnum).InterB_neg.(time_interval_field), neg_interval];
            else
                results(condition, subnum).InterB_pos.(time_interval_field) = [results(condition, subnum).InterB_pos.(time_interval_field), neg_interval];
                results(condition, subnum).InterB_neg.(time_interval_field) = [results(condition, subnum).InterB_neg.(time_interval_field), pos_interval];
            end
        end
    end
end

