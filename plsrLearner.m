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

        function [BETA, Ypred, loss, sse] = learn(obj, x, y, n)
            %learn does PLSR with n components
            [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(x, y, n);
            Ypred = BETA(1,:) + x * BETA(2:end,:);
            loss = (Ypred - y).^2;
            sse = sum(loss(:));
        end


        function [Ypred, trainSse, leaveOneOutSse] = leaveOneOutEvaluation(obj, n)
            %leaveOneOutEvaluation Evaluates a model with n components
            ndata = size(obj.xData, 1);
            parfor j = 1:ndata
                idxTrainig = ones(ndata, 1);
                idxTrainig(j) = 0;
                idxTrainig = logical(idxTrainig);
                xtrain = obj.xData(idxTrainig, :);
                ytrain = obj.yData(idxTrainig, :);                
                [BETA, ~, ~, sse] = obj.learn(xtrain, ytrain, n);
                sseModel(j) = sse;
                Ypred(j, :) = BETA(1,:) + obj.xData(~idxTrainig, :) * BETA(2:end,:);
            end
            leaveOneOutSse = sum((Ypred(:) - obj.yData(:)).^2);
            trainSse = mean(sseModel);
        end

        function [nopt, Ypred, trainSse, leaveOneOutSse] = optimizeComponentsAndLearn(obj, maxn)
            %optimizeComponentsAndLearn determine the number of components
            %that minimizes the loss of leave-on-out evaluations
            ncomps = 1:maxn;
            h = waitbar(0,'PLSR: Please wait...');
            for ncomp = ncomps
                waitbar(ncomp/length(ncomps),h)
                [YpredOut, trainSseOut, leaveOneOutSseOut] = obj.leaveOneOutEvaluation(ncomp);
                Ypred(:,:, ncomp) = YpredOut;
                trainSse(ncomp) = trainSseOut;
                leaveOneOutSse(ncomp) = leaveOneOutSseOut;
            end
            close(h)
            [~, nopt] = min(leaveOneOutSse);
        end
    end
end