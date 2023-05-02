%% load spectra for all peaks
tblSpectra = readtable('extractedPeaks/tblSpectra.csv');
KFOLDLEARN = 0; % k-fold for learning; use 0 for leave-one-out.
MAXN = 30; % maximum number of components to test
NRANDOMIZED = 1000;

%% load fold changes for all peaks
tblFoldChanges = readtable('folds/peakFoldChanges.csv');
% remove spectra that are not in foldchange table (compounds unchanged)
tblSpectra(~ismember(tblSpectra.peakId, tblFoldChanges.peakID), :) = [];

%% sort the two tables to make sure they have the same order
tblSpectra = sortrows(tblSpectra, "peakId");
tblFoldChanges = sortrows(tblFoldChanges, "peakID");

%% determine the optimum number of components (latent variables in model)
mSpectra = tblSpectra{:, 2:end};
X = mSpectra./ sqrt(sum(mSpectra.^2, 2));
Y = [tblFoldChanges.B, tblFoldChanges.L];

obj = plsrLearner(X, Y);
[nopt, Ypred, trainSse, testSse] = obj.optimizeComponentsAndLearn(MAXN, KFOLDLEARN);

%% Plot the training and evaluation losses
xComponents = 1:MAXN;

figure(1)
% plot training loss
subplot(2, 1, 1)
plot(xComponents, trainSse, 'o-', 'LineWidth',2)
xlabel('Number of latent components')
ylabel('Training loss')
xline(xComponents(nopt),'-', sprintf('Optimal (n=%d)', nopt));

% plot mean evaluation loss
subplot(2, 1, 2)
plot(xComponents, mean(testSse), 'o-', 'LineWidth',2)
hold off
xlabel('Number of latent components')
ylabel('Evaluation loss')
% Draw a vertical line to show number of latent components that
% minimizes the leave-on-out loss
xline(xComponents(nopt),'-', sprintf('Optimal (n=%d)', nopt));


%% compare predictions of model with optimal number of components and the real data

ypred = Ypred(:, :, nopt);
% calculate correlation between predictions and data, as well as p-value
[r, p] = corr(ypred(:), Y(:));

figure(2)
scatter(Y, ypred,  10, 'filled', MarkerEdgeColor='k')
xline(0,'k-');
yline(0,'k-');
grid on;
legend({'Brain-homing cells' 'Lung-homing cells'}, 'Location', 'northwest')
title({'True metabolite abundances vs. best prediction';...
    sprintf('(Learner with %d components: \\rho=%0.2f, P=%0.0e)', nopt, r, p)})
xlabel('log_2(abundance)')
ylabel('Prediction Score')
axis square


%% compare the loss of real data with loss from randomized data
% Shuffling experiment:
% shuffles the rows and traines the model with optimal number of components
h = waitbar(0,'Building loss distribution from shuffled Y: Please wait...');
i = 1;
while i<NRANDOMIZED
    waitbar(i/NRANDOMIZED,h)
    parfor j = i:(i+3)
        Yrand = Y(randperm(size(Y, 1)), :);
        objRand = plsrLearner(X, Yrand);
        [~, ~, leaveOneOutSse, ~] = objRand.leaveOneOutEvaluation(nopt);
        testSseRand(j) = sum(leaveOneOutSse);
    end
    i = i+4;
end

close(h)

%% plot the distribution histogram compared with the real loss
figure(3)
histogram(testSseRand, 100)
xline(sum(testSse(:, nopt)),'-', sprintf('Evaluation loss of model trained in real data', nopt));
xlabel('Evaluation loss of model trained in shuffled data')
ylabel('Distribution density')

%% train the optimal model with all the data
[BETA, Ypred, loss, sse, XL, YL, XS, YS, PCTVAR] = obj.learnWithAllData(nopt);

% calculate correlation between fit and data, as well as p-value
[r, p] = corr(Ypred(:), Y(:));

figure(4)
subplot(2, 1, 1)
scatter(Y, Ypred,  10, 'filled', MarkerEdgeColor='k')
xline(0,'k-');
yline(0,'k-');
grid on;
legend({'Brain-homing cells' 'Lung-homing cells'}, 'Location', 'northwest')
title({'True metabolite abundances vs. model fit with all data';...
    sprintf('(Learner with %d components: \\rho=%0.2f, P=%0.0e)', nopt, r, p)})
xlabel('log_2(abundance)')
ylabel('Prediction Score')
axis square

subplot(2, 1, 2)
gscatter(YL(1, :), YL(2, :), 1:size(YL, 2), jet(size(YL, 2)))
legend('AutoUpdate','off')
xline(0)
yline(0)
refline(1, 0)
grid on
axis square equal tight
xlabel('c weights in brain-homing' )
ylabel('c weights in lung-homing' )
title('The contribution of latent components to brain- and/or lung-homing')

%% interpret the learned model by looking at identified metabolites
% load FiehLib
importMsl;

% keep only the derivatized compounds withih FiehnLib
idx = (massSpectralLibrary.abundance73 == 0);
fprintf('FiehLib has %d non-derivatized compounds (abundance73 ==0)\n',...
    sum(idx));
disp('removed those.')
massSpectralLibrary(idx, :) = [];

% get the spectra for FiehnLib
xFiehnLib = cell2mat(massSpectralLibrary.abundance);
% normalize the spectra as done for X in our data
xFiehnLib = xFiehnLib./ sqrt(sum(xFiehnLib.^2, 2));

%% address the questions in paper dicussion
%what is C1?
mzArray = 50:599;
fprintf('C1 = (%f,%f)\n', YL(1,1), YL(2,1))

%what is C2?
fprintf('C2 = (%f,%f)\n', YL(1,2), YL(2,2))

figure(7)
stem(mzArray, XL(:, 2), 'Marker','none')
xlabel('m/z')
ylabel('W')
title('The weight vector for latent component 2');

%% get the score for component 2 of all identified metabolites
tblIndentified = readtable('identifiedFiehnLib/tblIdentity.csv');
tblW2Score = table(tblFoldChanges.peakID, XS(:, 2));
tblW2Score.Properties.VariableNames = {'peakID' 'w2Score'};
tblW2Score = innerjoin(tblW2Score, tblIndentified, "LeftKeys","peakID", "RightKeys","peakId");
tblW2Score(~tblW2Score.bestMatch95Rt, :) = [];
tblW2Score = sortrows(tblW2Score, "w2Score","descend")

%% plot their spectra
meanX = mean(X,1);
X0 = X - meanX;

n = 10;
for i = 1:n
    idx = find(tblFoldChanges.peakID == tblW2Score.peakID(i));
    subplot(ceil(n/2), 2, i);
    %stem(mzArray, X(idx, :), 'Marker','none');
    stem(mzArray, X0(idx, :), 'Marker','none');
    hold on
    plot(218, X0(idx, mzArray==218), 'r*')
    hold off
    title(sprintf('%s w2score=%f',tblW2Score.bestMatchFiehnLib{i},...
        tblW2Score.w2Score(i)))
end


%% This code illustrates how to reconstruct XS
X = obj.xData;
meanX = mean(X,1);
X0 = X - meanX;
Xtemp = X0*XL;

XSreconstructed = pinv(X0') * XL;

figure(5)
scatter(XS(:), XSreconstructed(:))
refline(1, 0)
xlabel('XS')
xlabel('XSreconstructed')

