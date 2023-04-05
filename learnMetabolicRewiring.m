%% load spectra for all peaks
tblSpectra = readtable('extractedPeaks/tblSpectra.csv');

%% load fold changes for all peaks
tblFoldChanges = readtable('folds/peakFoldChanges.csv');

%% remove spectra that are not in foldchange table (compounds unchanged)
tblSpectra(~ismember(tblSpectra.peakId, tblFoldChanges.peakID), :) = [];

%% sort the two tables to make sure they have the same order
tblSpectra = sortrows(tblSpectra, "peakId");
tblFoldChanges = sortrows(tblFoldChanges, "peakID");

%% determine the optimum number of components (latent variables in model)
X = mSpectra./ sqrt(sum(mSpectra.^2, 2));
Y = [tblFoldChanges.B, tblFoldChanges.L];
obj = plsrLearner(X, Y);
maxn = 25;
[nopt, Ypred, trainSse, leaveOneOutSse] = obj.optimizeComponentsAndLearn(maxn);

%% plot the results
figure(2)
% plot training and evaluation loss
subplot(2, 1, 1)
plot(1:maxn, trainSse, 'o-', 1:maxn, leaveOneOutSse, 'o-')
xlabel('Number of latent components')
ylabel('Loss')
% Draw a vertical line to show number of latent components that
% minimizes the leave-on-out loss
xline(nOptimal,'-', sprintf('Best (n=%d)', nopt));
legend('Training (mean SSE)', 'Leave-one-out (SSE)')
title('Number of components vs. learning loss')

% compare predictions from optimal model with real data
subplot(2, 1, 2)
ypred = Ypred(:, :, nopt);
[r, p] = corr(ypred(:), Y(:));
scatter(Y, ypred,  10, 'filled', MarkerEdgeColor='k')
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

