%% load spectra for all peaks
tblSpectra = readtable('extractedPeaks/tblSpectra.csv');
KFOLDLEARN = 0; % k-fold for learning
MAXN = 30; % maximum number of components to test
NRANDOMIZED = 1000;

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

%% explore the structure in the explanatory variable Y
D = pdist(X, 'cosine');
[cmdsScores, eigvals] = cmdscale(D, 2);

figure(2);
subplot(1, 2, 1)
scatter(cmdsScores(:,1),cmdsScores(:,2),'o', 'filled', 'MarkerEdgeColor','k')
xlabel('PCoA 1')
ylabel('PCoA 2')
grid on
axis equal square
title(sprintf('%d peaks from GC/MS data', height(tblFoldChanges)))

subplot(2, 2, 2)
scatter(cmdsScores(:,1),cmdsScores(:,2), [], tblFoldChanges.B, 'o', 'filled', 'MarkerEdgeColor','k')
colorbar
xlabel('PCoA 1')
ylabel('PCoA 2')
grid on
axis equal
title('log_2(FC) in brain-homing')
colormap redbluecmap
set(gca, 'CLim', [-2 2])
axis equal square

subplot(2, 2, 4)
scatter(cmdsScores(:,1),cmdsScores(:,2), [], tblFoldChanges.L, 'o', 'filled', 'MarkerEdgeColor','k')
colorbar
xlabel('PCoA 1')
ylabel('PCoA 2')
grid on
axis equal
colormap redbluecmap
title('log_2(FC) in lung-homing')
set(gca, 'CLim', [-2 2])
axis equal square