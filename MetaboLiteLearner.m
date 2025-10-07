classdef MetaboLiteLearner
    %plsrStandAloneLearner Uses partial least square regression to learn associations
    % between the spectrum of a metabolite and a response variable

    properties
        % input properties
        xFullData; % the spectra for each meatbolite
        yFullData; % the reponse variable
        kfold; % the kfold used to estimate loss and optimize components
        maxn; % the maximum number of components tested
        nrandomized; % the number of shuffiling iterations used to test
        % state variables
        cvIndices; % cross valiation indices
        nopt; % the optimum number of components
        testSse; % sum of squared errors of held-out examples (evaluation 
        BETA; % coefficients of linear model trained with all data
        PCTVAR; % percent variance explained
        Ypred; % fit results
        XL; % loadings of the X
        YL; % laodings of the Y
        XS; % scores of input data mapped onto nopt-dimensional latent space
        YS; % scores of output data mapped onto nopt-dimensional latent space
        stats; % the output of the plsregress
    end

    methods
        function obj = MetaboLiteLearner(x, y, kfold, maxn, nrandomized)
            %plsrLearner Construct an instance of plsrLearner
            obj.xFullData = x;
            obj.yFullData = y;
            obj.maxn = maxn;
            obj.nrandomized = nrandomized;
            % set up cross-validation indices
            if kfold == 0
                % Set up leave-one-out cross-validation indices
                obj.cvIndices = 1:size(obj.xFullData, 1);
            else
                % Set up k-fold cross-validation indices
                obj.cvIndices = crossvalind('Kfold', size(obj.xFullData, 1), kfold);
            end
            obj.kfold = max(obj.cvIndices);
            % optimize the components
            [obj.nopt, testSse, ~] = obj.optimizeComponentsAndLearn;
            % save the evaluation loss for the best model
            obj.testSse = mean(testSse(:, obj.nopt));
            % learn model with optimal number of components
            [obj.BETA, obj.Ypred, ~ , ~, obj.XL, obj.YL, obj.XS, obj.YS, obj.PCTVAR, obj.stats] =...
                obj.learn(obj.xFullData, obj.yFullData, obj.nopt);
            % plot the result
            % calculate correlation between predictions and data, as well as p-value
            [r, p] = corr(obj.Ypred(:), obj.yFullData(:));

            subplot(2, 2, 2)
            scatter(obj.yFullData, obj.Ypred,  10, 'filled', MarkerEdgeColor='k')
            xline(0,'k-');
            yline(0,'k-');
            legend({'Brain-homing cells' 'Lung-homing cells'}, 'Location', 'northwest')
            title({'True metabolite abundances vs. model prediction';...
                sprintf('(Learner with %d components: \\rho=%0.2f, P=%0.0e)', obj.nopt, r, p)})
            xlabel('log_2(FC) from actual data')
            ylabel({'log_2(FC) fit by model'})
            axis equal square
            grid on;
        end

        function [ls, predictedY] = mapToLatentSpace(obj, x)
            % mapp a set of spectra to the latent space
            X0 = x - mean(obj.xFullData,1);
            ls = X0 * obj.stats.W;
            predictedY = obj.BETA(1,:) + x * obj.BETA(2:end,:);
        end

        function [] = shufflingTest(obj)            
            % Shuffling experiment:
            % shuffles the rows and traines the model with optimal number of components
            h = waitbar(0,'Shuffle test. Please wait...');
            i = 1;
            while i<obj.nrandomized
                waitbar(i/obj.nrandomized,h)
                parfor j = i:(i+3)
                    % shuffle Y columns
                    Yrand = obj.yFullData(randperm(size(obj.yFullData, 1)), :);
                    % evaluate model
                    [~, ~, shuffleSse, ~] =...
                        obj.crossValidationEvaluation(obj.xFullData, Yrand, obj.nopt);
                    testSseRand(j) = mean(shuffleSse);
                end
                i = i+4;
            end
            close(h)


            % plot the distribution histogram compared with the real loss
            histogram(testSseRand, 100)
            xline(obj.testSse, '--',...
                sprintf('MSE of model (N=%d) trained in real data', obj.nopt));
            xlabel(sprintf('MSE of models (N=%d) trained in shuffled data', obj.nopt))
            ylabel('Distribution density')
            xlimits = get(gca, 'XLim');
            set(gca, 'XLim', [obj.testSse*0.5 xlimits(2)])
        end

        function [BETA, Ypred, loss, sse, XL, YL, XS, YS, PCTVAR, stats] = learn(obj, x, y, n)
            %learn does PLSR with n components
            [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(x, y, n);
            Ypred = BETA(1,:) + x * BETA(2:end,:);
            loss = (Ypred - y).^2;
            sse = sum(loss(:));
        end

        function [Ypred, trainSse, kFoldSse, betaDistribution] = crossValidationEvaluation(obj, x, y, n)
            % Evaluates a PLSR model with 'n' components using k-fold cross-validation.
            % Input:
            %   n: Number of components to use in the PLSR model.
            % Output:
            %   Ypred: Predicted values of the response variable.
            %   trainSse: Mean sum of squared errors for training data.
            %   kFoldSse: Sum of squared errors for k-fold cross-validation.
            %   betaDistribution: the distribution of betavalues accross
            %   cross validations

            % pre-allocate outputs
            Ypred = zeros(size(y));
            trainSse = zeros(obj.kfold, 1);
            kFoldSse = zeros(obj.kfold, 1);
            betaDistribution = zeros(size(x, 2) + 1, size(y, 2), obj.kfold);
            % TODO runs faster with parfor
            for i = 1:obj.kfold
                % Split the data into training and test sets
                trainIdx = (obj.cvIndices ~= i);
                testIdx = (obj.cvIndices == i);

                xTrain = x(trainIdx, :);
                yTrain = y(trainIdx, :);
                xTest = x(testIdx, :);

                % Train the PLSR model on the training data
                [BETA, ~, ~, sse] = obj.learn(xTrain, yTrain, n);

                % Predict the response variable for the test set
                YpredTest = [ones(size(xTest, 1), 1), xTest] * BETA;
                Ypred(testIdx, :) = YpredTest;
                % Keep the sum of squared errors for the training set
                trainSse(i) = sse / size(xTrain, 1);
                %kFoldSse(i, 1) = sum((YpredTest - y(testIdx, :)).^2, "all") * obj.kfold;
                kFoldSse(i) = sum((YpredTest - y(testIdx, :)).^2, "all") / size(xTest, 1);
                % build a distribution of coefficients
                idx = find(testIdx);
                for j = 1:length(idx)
                    betaDistribution(:, :, idx(j)) = BETA;
                end
            end
            % Calculate the mean sum of squared errors for the training data
            %trainSse = mean(trainSse);
        end


        function [nopt, testSse, trainSse] = optimizeComponentsAndLearn(obj)
            %optimizeComponentsAndLearn determine the number of components
            %that minimizes the evaluationi loss
            % TODO: can be optimized, by solving once for all components

            ncomps = 1:obj.maxn;
            h = waitbar(0,'Determining number of latent components...');

            % pre-allocate outputs
            Ypred = zeros(size(obj.yFullData, 1), size(obj.yFullData, 2), obj.kfold);
            trainSse = zeros(obj.kfold, obj.maxn);
            testSse = zeros(obj.kfold, obj.maxn);
            % iterate through the components to determine best
            for ncomp = ncomps
                waitbar(ncomp/length(ncomps),h)
                [Yp, trainSseO, testSseO] =...
                    obj.crossValidationEvaluation(obj.xFullData,...
                    obj.yFullData, ncomp);
                Ypred(:,:, ncomp) = Yp;
                trainSse(:, ncomp) = trainSseO;
                testSse(:, ncomp) = testSseO;
            end
            close(h)
            % the ooptimal number of components is the one with minimum
            % average loss + 1 standard error
            meanTestSse = mean(testSse);
            [~, nMin] = min(meanTestSse);
            % calculate standard error
            stderror = std(testSse) / sqrt( size(testSse, 1) );
            %
            minError = meanTestSse(nMin) +  stderror(nMin);
            candidateIdx = find(meanTestSse > minError & (1:numel(stderror)) < nMin);
            if isempty(candidateIdx)
                nopt = nMin;
            else
                nopt = max(candidateIdx);
            end
  
            % darw figure with evaluation and training
            figure
            % plot training loss
            subplot(1, 2, 1)
            hold on
            plot(ncomps, mean(trainSse), 'o-', 'LineWidth',2)
            %stderror2 = std(trainSse) / sqrt( size(trainSse, 1) );
            %errorbar(ncomps, mean(trainSse), stderror2,'k-')
            % plot mean evaluation loss
            errorbar(ncomps, mean(testSse), stderror, 'k-', 'Marker', 'none')
            plot(ncomps, mean(testSse), 'o-', 'LineWidth',2)
            hold off
            xlabel('Number of latent components')
            ylabel('MSE')
            % Draw a vertical line to show number of latent components that
            % minimizes the leave-on-out loss
            xline(ncomps(nMin),'k:', sprintf('Minimum (N=%d)', nMin));
            xline(ncomps(nopt),'k--', sprintf('Optimal (N=%d)', nopt));

            % plot the leave one out data 
            ypred = Ypred(:, :, nopt);
            % calculate correlation between predictions and data, as well as p-value
            [r, p] = corr(ypred(:), obj.yFullData(:));

            subplot(2, 2, 4)
            scatter(obj.yFullData, ypred,  10, 'filled', MarkerEdgeColor='k')
            xline(0,'k-');
            yline(0,'k-');
            legend({'Brain-homing cells' 'Lung-homing cells'}, 'Location', 'northwest')
            title({'True metabolite abundances vs. model prediction';...
                sprintf('(Learner with %d components: \\rho=%0.2f, P=%0.0e)', nopt, r, p)})
            xlabel('log_2(FC) from actual data')
            ylabel({'log_2(FC) predicted by model';'in leave-one-out cross-validation'})
            axis equal square
            %set(gca, 'XLim', [-2.5 2.5], 'YLim', [-2.5 2.5])
            grid on;

        end
    end
end