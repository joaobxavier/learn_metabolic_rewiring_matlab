%% load peak data
tlbSampleTicPeaksIntegrated = readtable('extractedPeaks/tblPeaksIntegrated.csv');

%% extract sample type
cellType = tlbSampleTicPeaksIntegrated.Properties.VariableNames(2:end)';
for i = 1:length(cellType)
    cellType{i} = cellType{i}(1);
end

%% Remove peaks that are not different from media
pValues = ones(height(tlbSampleTicPeaksIntegrated), 1);
for i = 1:height(tlbSampleTicPeaksIntegrated)
    y = log(tlbSampleTicPeaksIntegrated{i, 2:end});
    x = cellType;
    m = fitlm(categorical(x), y);
    aov = anova(m);
    pValues(i) = aov.pValue(1);
end

%% PCA before rescalling 
figure(1)
m = tlbSampleTicPeaksIntegrated{pValues<0.05, 2:end}';
[coeff, score, latent, tsquared, explained] = pca(m);
gscatter(score(:, 1), score(:, 2), cellType)
xlabel(sprintf('PC1 (%0.1f%%EV)', explained(1)))
ylabel(sprintf('PC2 (%0.1f%%EV)', explained(2)))
grid on;
axis equal;

%% keep only the metabolites singificant in cells
cols = tlbSampleTicPeaksIntegrated.Properties.VariableNames;
tblBeforeScalling = tlbSampleTicPeaksIntegrated(pValues<0.05, ~contains(cols, 'M'));

%% rescale the samples with intracellular metabolome 
% log-transform to make multiplicative model of scalling
tblScalled = stack(tblBeforeScalling, 2:width(tblBeforeScalling),...
    'NewDataVariableName','area', 'IndexVariableName','sample');
% add 'batch'
tblScalled.sample = string(tblScalled.sample);
tblScalled.batch = contains(tblScalled.sample, 'plate1') +...
    2*contains(tblScalled.sample, 'plate2') +...
    3*contains(tblScalled.sample, 'plate3');
tblScalled.batch = categorical(tblScalled.batch);
% log transfor to make multiplicative model
tblScalled.area = log(tblScalled.area);

% create a categorical for the peak ID and add empty refference category
tblScalled.peakId = categorical(tblScalled.peakId);
c = categories(tblScalled.peakId);
tblScalled.peakId = addcats(tblScalled.peakId,...
    'ref', 'Before', c{1});

% get the first letter of the sample name, which is the type of cell (P, B
% or L)
tblScalled.cell =...
    cellfun(@(x) x(1), tblScalled.sample, 'UniformOutput',false);

% fit model and calculate corrected peak areas (fold changes relative to corrected mean)
mdl = fitlme(tblScalled,...
       'area ~ peakId - 1 + (1|sample) + (1|peakId:batch)');

%% tlbCellPeaksScalled.areaModel = residuals(mdl);
tblScalled.areaModel = mdl.Residuals.Raw;
tblScalledUnstack = tblScalled(:, {'peakId' 'sample' 'areaModel'});
tblScalledUnstack = unstack(tblScalledUnstack, 'areaModel', 'sample');

%% do the PCA to compare before and after correction

cellType2 = tblScalledUnstack.Properties.VariableNames(2:end)';
for i = 1:length(cellType2)
    cellType2{i} = cellType2{i}(1);
end

batch = tblScalledUnstack.Properties.VariableNames(2:end)';
batch = contains(batch, 'plate1') +...
    2*contains(batch, 'plate2') +...
    3*contains(batch, 'plate3');

figure(1)
% PCA before correction
m = tblBeforeScalling{:, 2:end}';
m = zscore(log(m), [], 2);
[coeff, score, latent, tsquared, explained] =...
    pca(m);
subplot(2, 2, 1)
gscatter(score(:, 1), score(:, 2), cellType2)
xlabel(sprintf('PC1 (%0.1f%%EV)', explained(1)))
ylabel(sprintf('PC2 (%0.1f%%EV)', explained(2)))
grid on;
axis equal;
title('Cell type: after correction')

subplot(2, 2, 2)
gscatter(score(:, 1), score(:, 2), batch')
xlabel(sprintf('PC1 (%0.1f%%EV)', explained(1)))
ylabel(sprintf('PC2 (%0.1f%%EV)', explained(2)))
grid on;
axis equal;
title('Batch: after correction')

% PCA after correction
m = tblScalledUnstack{:, 2:end}';
[coeff, score, latent, tsquared, explained] =...
    pca(m);
subplot(2, 2, 3)
gscatter(score(:, 1), score(:, 2), cellType2)
xlabel(sprintf('PC1 (%0.1f%%EV)', explained(1)))
ylabel(sprintf('PC2 (%0.1f%%EV)', explained(2)))
grid on;
axis equal;
title('Cell type: after correction')

subplot(2, 2, 4)
gscatter(score(:, 1), score(:, 2), batch')
xlabel(sprintf('PC1 (%0.1f%%EV)', explained(1)))
ylabel(sprintf('PC2 (%0.1f%%EV)', explained(2)))
grid on;
axis equal;
title('Batch: after correction')

%% compute the fold change for each metabolite in B and L
tblScalled.cell = categorical(tblScalled.cell, {'P' 'B' 'L'});
mdl2 = fitlme(tblScalled,...
       'area ~ peakId + cell:peakId - 1 + (1|sample) + (1|peakId:batch)');
coefficientTable = mdl2.Coefficients;
coefficientTableP = coefficientTable(~contains(coefficientTable.Name, ':'), :);
coefficientTableBL = dataset2table(coefficientTable(contains(coefficientTable.Name, ':'), :));
for i = 1:height(coefficientTableBL)
    s = strsplit(coefficientTableBL.Name{i}, ':');
    coefficientTableBL.peakID{i} = strrep(s{1}, 'peakId_', '');
    coefficientTableBL.cell{i} = strrep(s{2}, 'cell_', '');
end

folds = unstack(coefficientTableBL(:, {'peakID' 'cell' 'Estimate'}), 'Estimate', 'cell');
% covert to base 2
folds.B = folds.B ./ log(2); 
folds.L = folds.L ./ log(2); 
% add the p-vales
pVals = unstack(coefficientTableBL(:, {'peakID' 'cell' 'pValue'}), 'pValue', 'cell');
folds.pValB = pVals.B;
folds.pValL = pVals.L;
% add the lower confidence intervals
lowerCI = unstack(coefficientTableBL(:, {'peakID' 'cell' 'Lower'}), 'Lower', 'cell');
folds.lowerB = lowerCI.B./ log(2);
folds.lowerL = lowerCI.L./ log(2);
% add the upper confidence intervals
upperCI = unstack(coefficientTableBL(:, {'peakID' 'cell' 'Upper'}), 'Upper', 'cell');
folds.upperB = upperCI.B./ log(2);
folds.upperL = upperCI.L./ log(2);


% plot the fold-changes in a scatter plot with each metabolite colored by
% p-value
figure(2)
errorbar(folds.B, folds.L,...
    folds.L - folds.lowerL, folds.upperL - folds.L,...
    folds.B - folds.lowerB, folds.upperB - folds.B, '.', 'Color', [1 1 1]*0.8)
hold on;
scatter(folds.B, folds.L, [], -log10(folds.pValB .* folds.pValL), 'filled', 'MarkerEdgeColor','k')
hold off
xlabel('Brain-homing vs. parental [log_2(FC)]')
ylabel('Lung-homing vs. parental [log_2(FC)]')
grid on
axis equal square
h = refline(1, 0)
set(h, 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--')
set(gca, 'CLim', [1.3 2])
% correlation between brain-homing and lung-homing 
[r, p] = corr(folds.B, folds.L);
title(sprintf('%d compounds; \\rho = %0.2f; P-value = %0.0e',...
    height(folds), r, p));
% hcb=colorbar;
% title(hcb,'-log10(PValue)');

%% save the folds table
if exist('folds', 'dir') == 0
    disp('creating folds directory');
    mkdir('folds');
end
writetable(folds, 'folds/peakFoldChanges.csv',...
    "FileType","text", 'Delimiter', ',');


