classdef plsrLearner
    %plsrLearner Uses partial least square regression to learn associations
    % between the spectrum of a metabolite and a response variable

    properties
        xData; % the spectra for each meatbolite
        yData; % the reponse variable
    end

    methods
        function obj = plsrLearner(x, y)
            %plsrLearner Construct an instance of plsrLearner
            obj.xData = x;
            obj.yData = y;
        end

        function [BETA, Ypred, loss, sse, XL, YL, XS, YS, PCTVAR] = learn(obj, x, y, n)
            %learn does PLSR with n components
            [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(x, y, n);
            Ypred = BETA(1,:) + x * BETA(2:end,:);
            loss = (Ypred - y).^2;
            sse = sum(loss(:));
        end


        function [BETA, Ypred, loss, sse, XL, YL, XS, YS, PCTVAR] = learnWithAllData(obj, n)
            %learn does PCR with n components
            [BETA, Ypred, loss, sse XL, YL, XS, YS, PCTVAR] =...
                learn(obj, obj.xData, obj.yData, n);
        end



        function [Ypred, trainSse, leaveOneOutSse, betaDistribution] = leaveOneOutEvaluation(obj, n)
            %leaveOneOutEvaluation Evaluates a model with n components
            ndata = size(obj.xData, 1);
            for j = 1:ndata
                idxTrainig = ones(ndata, 1);
                idxTrainig(j) = 0;
                idxTrainig = logical(idxTrainig);
                xtrain = obj.xData(idxTrainig, :);
                ytrain = obj.yData(idxTrainig, :);
                xTest  = obj.xData(~idxTrainig, :);
                [BETA, ~, ~, sse] = obj.learn(xtrain, ytrain, n);
                sseModel(j) = sse;
                Ypred(j, :) = [ones(size(xTest, 1), 1), xTest] * BETA;
                betaDistribution(:, :, j) = BETA;
            end
            leaveOneOutSse = (Ypred(:) - obj.yData(:)).^2;
            trainSse = mean(sseModel);
        end


        function [Ypred, trainSse, kFoldSse, betaDistribution] = kFoldEvaluation(obj, n, k)
            % Evaluates a PLSR model with 'n' components using k-fold cross-validation.
            % Input:
            %   n: Number of components to use in the PLSR model.
            %   k: Number of folds for cross-validation.
            % Output:
            %   Ypred: Predicted values of the response variable.
            %   trainSse: Mean sum of squared errors for training data.
            %   kFoldSse: Sum of squared errors for k-fold cross-validation.

            % Set up k-fold cross-validation indices
            cvIndices = crossvalind('Kfold', size(obj.xData, 1), k);

            Ypred = zeros(size(obj.yData));
            trainSse = zeros(k, 1);

            for i = 1:k
                % Split the data into training and test sets
                trainIdx = (cvIndices ~= i);
                testIdx = (cvIndices == i);

                xTrain = obj.xData(trainIdx, :);
                yTrain = obj.yData(trainIdx, :);
                xTest = obj.xData(testIdx, :);

                % Train the PLSR model on the training data
                [BETA, ~, ~, sse] = obj.learn(xTrain, yTrain, n);

                % Predict the response variable for the test set
                YpredTest = [ones(size(xTest, 1), 1), xTest] * BETA;
                Ypred(testIdx, :) = YpredTest;
                % Keep the sum of squared errors for the training set
                trainSse(i) = sse;
                kFoldSse(i, 1) = sum((YpredTest - obj.yData(testIdx, :)).^2, "all");
                % build a distribution of coefficients
                idx = find(testIdx);
                for j = 1:length(idx)
                    betaDistribution(:, :, idx(j)) = BETA;
                end
            end
            % Calculate the mean sum of squared errors for the training data
            trainSse = mean(trainSse);
        end


        function [nopt, Ypred, trainSse, testSse] = optimizeComponentsAndLearn(obj, maxn, kfold)
            %optimizeComponentsAndLearn determine the number of components
            %that minimizes the loss of leave-on-out evaluations
            ncomps = 1:maxn;
            h = waitbar(0,'Optimizing PLSR components: Please wait...');
            for ncomp = ncomps
                waitbar(ncomp/length(ncomps),h)
                if nargin == 2 || kfold == 0
                    [Yp, trainSseO, testSseO] = obj.leaveOneOutEvaluation(ncomp);
                else
                    [Yp, trainSseO, testSseO] = obj.kFoldEvaluation(ncomp, kfold);
                end
                Ypred(:,:, ncomp) = Yp;
                trainSse(ncomp) = trainSseO;
                testSse(:, ncomp) = testSseO;
            end
            close(h)
            [~, nopt] = min(mean(testSse));
        end
    end
end