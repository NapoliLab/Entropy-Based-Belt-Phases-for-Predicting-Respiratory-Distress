function [selectedFeaturesIdx, selectedFeaturesLabels] = stepwiseForwardSelection(DataPred, ResponseVariable, Labels, maxFeatureLim)

    % Initialization
    numTotalFeatures = length(Labels);
    itercount = 0;
    indicatorFeaturesMat = false(1, numTotalFeatures);  

    % Initialize matrix for AUC performance tracker
    AUCperformance = [];

    % Initial model performance
    basePerformance = -Inf;

    while itercount < maxFeatureLim
        bestPerformance = -Inf;
        featureIdx = find(~indicatorFeaturesMat);

        % Forward selection part
        for iter = 1:length(featureIdx)
            % Temporary model with current feature
            tempModel = [DataPred(:, indicatorFeaturesMat), DataPred(:, featureIdx(iter))];
            currentPerformance = LOOCV_function(tempModel, ResponseVariable);

            % Check if this feature improves the model
            if currentPerformance > bestPerformance
                bestPerformance = currentPerformance;
                bestAttribute = featureIdx(iter);
            end
        end

        itercount = itercount+1;
        AUCperformance(itercount) = bestPerformance;

        % Update the model with the best new feature
        if bestPerformance > basePerformance
            basePerformance = bestPerformance;
        end

        indicatorFeaturesMat(bestAttribute) = true;
    end

    
    selectedFeaturesIdx = find(indicatorFeaturesMat);
    selectedFeaturesLabels = Labels(1, find(indicatorFeaturesMat));


    %% SUPPORTING FUNCTIONS    
    function AUC = LOOCV_function(DataPred, responseVec)
        n = size(DataPred, 1);
        cv = cvpartition(n, 'LeaveOut');

        predicted_scores = zeros(n, 1);
        true_labels = zeros(n, 1);

        for i = 1:n
            % Training/testing indices
            trainIdx = cv.training(i);
            testIdx = cv.test(i);

            % Training/testing sets
            X_train = DataPred(trainIdx, :);
            y_train = responseVec(trainIdx, :);
            X_test = DataPred(testIdx, :);

            % Fit logistic regression model
            model = fitglm(X_train, y_train, 'Distribution', 'binomial');

            % Test model and store the predicted score
            prob = predict(model, X_test);
            predicted_scores(i) = prob;

            % Store the actual label
            true_labels(i) = responseVec(testIdx);
        end

        % Calculate AUC
        [~, ~, ~, AUC] = perfcurve(true_labels, predicted_scores, 1);
    end
end


function out = ifelse(condition, trueVal, falseVal)
    if condition
        out = trueVal;
    else
        out = falseVal;
    end
end
