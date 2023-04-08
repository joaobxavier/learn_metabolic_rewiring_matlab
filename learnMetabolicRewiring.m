%% load spectra for all peaks
tblSpectra = readtable('extractedPeaks/tblSpectra.csv');

%% load fold changes for all peaks
tblFoldChanges = readtable('folds/peakFoldChanges.csv');

figure(1)
errorbar(tblFoldChanges.B, tblFoldChanges.L,...
    tblFoldChanges.L - tblFoldChanges.lowerL, tblFoldChanges.upperL - tblFoldChanges.L,...
    tblFoldChanges.B - tblFoldChanges.lowerB, tblFoldChanges.upperB - tblFoldChanges.B, '.', 'Color', [1 1 1]*0.8)
hold on;
scatter(tblFoldChanges.B, tblFoldChanges.L, [], -log10(tblFoldChanges.pValB .* tblFoldChanges.pValL), 'filled', 'MarkerEdgeColor','k')
hold off
xlabel('Brain-homing vs. parental [log_2(FC)]')
ylabel('Lung-homing vs. parental [log_2(FC)]')
grid on
axis equal square
h = refline(1, 0)
set(h, 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--')
set(gca, 'CLim', [1.3 2])
% correlation between brain-homing and lung-homing 
[r, p] = corr(tblFoldChanges.B, tblFoldChanges.L);
title(sprintf('%d compounds; \\rho = %0.2f; P-value = %0.0e',...
    height(tblFoldChanges), r, p));

%% remove spectra that are not in foldchange table (compounds unchanged)
tblSpectra(~ismember(tblSpectra.peakId, tblFoldChanges.peakID), :) = [];

%% sort the two tables to make sure they have the same order
tblSpectra = sortrows(tblSpectra, "peakId");
tblFoldChanges = sortrows(tblFoldChanges, "peakID");

%% determine the optimum number of components (latent variables in model)
mSpectra = tblSpectra{:, 2:end};
X = mSpectra./ sqrt(sum(mSpectra.^2, 2));
Y = [tblFoldChanges.B, tblFoldChanges.L];
%obj = pcaRegressionLearner(X, Y);
obj = plsrLearner(X, Y);
maxn = 30;
[nopt, Ypred, trainSse, testSse] = obj.optimizeComponentsAndLearn(maxn, 5);

%% plot the results
figure(3)
% plot training and evaluation loss
subplot(2, 2, 1)
plot(1:maxn, trainSse, 'o-')
xlabel('Number of latent components')
ylabel('Training loss')

subplot(2, 2, 3)
plot(1:maxn, testSse, 'o-')
xlabel('Number of latent components')
ylabel('Evaluation loss')
% Draw a vertical line to show number of latent components that
% minimizes the leave-on-out loss
xline(nopt,'-', sprintf('Best (n=%d)', nopt));

% compare predictions from optimal model with real data
subplot(1, 2, 2)
ypred = Ypred(:, :, nopt);

[r, p] = corr(ypred(:), Y(:));
scatter(Y, ypred,  30, 'filled', MarkerEdgeColor='k')
hold on
plot(Y(:, :)', ypred(:, :)', 'k-')
hold off
xline(0,'k-');
yline(0,'k-');
% refline(1, 0)
grid on;
legend({'Brain-homing cells' 'Lung-homing cells'}, 'Location', 'northwest')
title({'True metabolite abundances vs. best prediction';...
    sprintf('(Learner with %d components: \\rho=%0.2f, P=%0.0e)', nopt, r, p)})
xlabel('log_2(abundance)')
ylabel('Prediction Score')
axis square

%% plot the credible intervals for the effects
%[~, ~, ~, betaDistribution] = obj.leaveOneOutEvaluation(nopt);
[~, ~, ~, betaDistribution] = obj.kFoldEvaluation(nopt, 5);

modeBeta = mode(betaDistribution, 3);
quantileBeta = quantile(betaDistribution, [0.025 0.975],3);

figure(4)
subplot(1, 3, 1)
errorbar(modeBeta(:, 1), modeBeta(:, 2),...
    modeBeta(:, 2) - quantileBeta(:, 2, 1),...
    quantileBeta(1, 2, 2) - modeBeta(:, 2),...
    modeBeta(:, 1) - quantileBeta(:, 1, 1),...
    quantileBeta(:, 1, 2) - modeBeta(:, 2),...
    'o', 'Color', [1 1 1]*0.5,...
    'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k')
hold on
scatter(modeBeta(1, 1), modeBeta(1, 2), 'filled')
scatter(modeBeta(2:end, 1), modeBeta(2:end, 2), 'filled',...
    'MarkerEdgeColor', 'k')
hold off
xlabel('\beta Brain-homing')
ylabel('\beta Lung-homing')
grid on
refline(1, 0)
axis equal tight
legend('\beta_0', '\beta_{mz}')
title('Effect sizes of ion fragments on metabolic rewiring')

% plot the effect spectra
mzArray = 50:599;

subplot(2, 3, 2)
stem(mzArray, modeBeta(2:end, 1), 'Marker','none')
xlabel('mz of ion fragment')
ylabel('\beta Brain-homing')

subplot(2, 3, 5)
stem(mzArray, modeBeta(2:end, 2), 'Marker','none')
xlabel('mz of ion fragment')
ylabel('\beta Lung-homing')

corrArray = corr(X, Y);
subplot(2, 3, 3)
stem(mzArray, corrArray(:, 1), 'Marker','none')
xlabel('mz of ion fragment')
ylabel('\rho Brain-homing')

subplot(2, 3, 6)
stem(mzArray, corrArray(:, 2), 'Marker','none')
xlabel('mz of ion fragment')
ylabel('\rho Lung-homing')

%% to interpret what the model is doing we will now focus on the 3 metabolites
% most depleted in both B and L (they are the same in both
[~, idxFocus] = sortrows(tblFoldChanges.B);
idxFocus = idxFocus(1:3);

% plot the distrubution of the model predictions for each of these 3 compounds
for j = 1:3
    i = idxFocus(j);
    xTest = X(i, :);
    yTest = Y(i, :)
    BETA = betaDistribution(:, :, i);
    Ypred = [ones(size(xTest, 1), 1), xTest] * BETA

    figure(5)
    subplot(3, 3, j)
    stem(mzArray, BETA(2:end, 1) .* xTest', 'Marker','none')
    xlabel('mz')
    ylabel('\beta_BI')
    subplot(3, 3, j+3)
    stem(mzArray, BETA(2:end, 2) .* xTest', 'Marker','none')
    xlabel('mz')
    ylabel('\beta_LI')
    subplot(3, 3, j+6)
    scatter(BETA(2:end, 1) .* xTest', BETA(2:end, 2) .* xTest', 'filled')
    xlabel('\beta_BI')
    ylabel('\beta_LI')
    grid on
    refline(1, 0)
    axis equal tight
    title(sprintf('PeakId=%d', i))
end

%% important mzs: 305 and 214
m = find(mzArray == 305)
d = 2;
b = betaDistribution(m+1, d, :);

figure(6)
subplot(2, 1, 1)
histogram(b(:), 100)
hold on
plot(squeeze(betaDistribution(m+1, d, idxFocus)), [0 0 0], 'r*')
hold off
title(sprintf('mz = %d', mzArray(m)))
xlabel('\beta_BI')
ylabel('Density')

subplot(2, 1, 2)
scatter(X(:, m), Y(:, d))
h = lsline;
set(h, 'Color', 'r')
xlabel(sprintf('Intensity mz = %d', mzArray(m)))
ylabel('log2FC')
[r, p] = corr(X(:, m), Y(:, d));
title(sprintf('\\rho = %0.2f; P-value = %0.0e', r, p));
