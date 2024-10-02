% Closes all open figures, clears all variables, clears command window
clearvars;
close all;
clc;

% Loads the data as well as some coefficients for filters used later
load('PulseTrainData_HypoxiaVsSeaLevel.mat')

% Initialize structures for non-hypoxic and hypoxic data by participants
% TI - Time of Inspirations (Time from Trough to Peak)
% TE - Time of Expirations (Time from Peak to Trough)
% TT - Time of Inspirations and Expirations saved in chronological order
% (TI, TE, TI, TE ...)
% TB - Time of one full breath Inspiration and Expiration (TI+TE, TI+TE ..)
% O2 - Oxygen Saturation for the duration of each TB window

PTresults = struct('TI', cell(2,49), 'TE', cell(2,49), 'TT', cell(2,49), 'TB', cell(2,49), 'O2', cell(2,49));
time_interval_names = {'fullLength', 'first5', 'last5', 'last2'};

for subnum = 1:49
    for condition = 1:2 % 1 for HYPOXIC, 2 for NORMOXIC
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
        [zeroCrossingsCleaned, PeakTroughVals, PeakTroughLocs] = peak_detector(smoothed_bre, 256, 1);
        
        % Analysis over different time intervals
        for time_interval = 1:size(intervals, 1)
            PTresults = analyze_intervals(t, O2, PeakTroughLocs, PeakTroughVals, intervals(time_interval, :), time_interval, time_interval_names, PTresults, condition, subnum);
        end
    end
end

save("PTresults.mat","PTresults");


%% SUPPORTING FUNCTIONS

function intervals = calculate_intervals(length_sm)
    % Calculate the intervals for analysis
    mid = ceil(length_sm/2);
    last_2 = length_sm - ceil(length_sm/5);
    intervals = [1, length_sm; 1, mid; mid+1, length_sm; last_2+1, length_sm];
end

function results = analyze_intervals(t,O2, PeakTroughLocs, PeakTroughVals, interval, time_interval_iteration, time_interval_names, results, condition, subnum)
    indicator = PeakTroughLocs > interval(1) & PeakTroughLocs < interval(2);
    indices = find(indicator);
    
    for i =  1:length(indices)-2
        if mod(i,2) == 1 
            O2_window = O2(PeakTroughLocs(indices(i)):PeakTroughLocs(indices(i+2)));
            te = t(PeakTroughLocs(indices(i+1))) - t(PeakTroughLocs(indices(i)));
            ti = t(PeakTroughLocs(indices(i+2))) - t(PeakTroughLocs(indices(i+1)));
            time_interval_field = time_interval_names{time_interval_iteration};

            % Check if the field exists, if not initialize it
            if ~isfield(results(condition, subnum).TT, time_interval_field)
                results(condition, subnum).TT.(time_interval_field) = [];
                results(condition, subnum).TB.(time_interval_field) = [];
                results(condition, subnum).TI.(time_interval_field) = [];
                results(condition, subnum).TE.(time_interval_field) = [];
                results(condition, subnum).O2.(time_interval_field) = [];
            end

            results(condition, subnum).TT.(time_interval_field) = [results(condition, subnum).TT.(time_interval_field), te, ti];
            results(condition, subnum).TB.(time_interval_field) = [results(condition, subnum).TB.(time_interval_field), te + ti];
            results(condition, subnum).O2.(time_interval_field) = [results(condition, subnum).O2.(time_interval_field),  mean(O2_window(O2_window~=0))];
            if PeakTroughVals(indices(1)) > 0
                results(condition, subnum).TE.(time_interval_field) = [results(condition, subnum).TE.(time_interval_field), te];
                results(condition, subnum).TI.(time_interval_field) = [results(condition, subnum).TI.(time_interval_field), ti];
            else
                results(condition, subnum).TE.(time_interval_field) = [results(condition, subnum).TE.(time_interval_field), ti];
                results(condition, subnum).TI.(time_interval_field) = [results(condition, subnum).TI.(time_interval_field), te];
            end
        end
    end

end
