%% Set variables to false to skip workflow steps
DOWNLOAD = false;
CONVERT_AGILENT_TO_CSVS = false;
EXTRACT_SPECTRA_AND_INTEGRATE = false;
SHUFFLE_TEST = false;

%% Download the raw data from Zenodo
if (DOWNLOAD) 
    dataUrl = ['https://zenodo.org/record/8193580/files/'...
        'metastasis_lineages_3replicates_3runs_10302020.zip?download=1'];
    if ~isdir('rawAgilentData')
        mkdir('rawAgilentData');
    end
    cd('rawAgilentData')
    disp('Dowloading raw data from Zenodo...')
    imageFileFullPath =...
        websave('metastasis_lineages_3replicates_3runs_10302020.zip',...
        dataUrl);
    disp('Done. Unzipping...')
    unzip('metastasis_lineages_3replicates_3runs_10302020.zip')
    cd('..')
    disp('Done.')
end

%% import Agilent, bin and export as cvs (convertAgilentToCvs.m)
% requires chromatography-master in the path
if (CONVERT_AGILENT_TO_CSVS)
    disp('Converting Agilent files to CSVs.')
    convertAgilentToCvs(...
        'rawAgilentData/metastasis_lineages_3replicates_3runs_10302020',...
        'gcmsCSVs');
    disp('Done.')
    % remove the ethyl acetate files
    delete('gcmsCSVs/ethyl_acetate*.csv');
end

%% extract spectra and integrate (extractSpectraAndIntegrate.m)
if (EXTRACT_SPECTRA_AND_INTEGRATE)
    [tblPeaksIntegrated, tblSpectra] =...
        extractSpectraAndIntegrate('gcmsCSVs', 'extractedPeaks');
end

%% load peak data
tlbSampleTicPeaksIntegrated =...
    readtable('extractedPeaks/tblPeaksIntegrated.csv');
tblSpectra = readtable('extractedPeaks/tblSpectra.csv');

%% extract sample type from sample name
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
% keep only the metabolites singificant in cells
cols = tlbSampleTicPeaksIntegrated.Properties.VariableNames;
tblBeforeScalling = tlbSampleTicPeaksIntegrated(pValues<0.05,...
    ~contains(cols, 'M'));

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

%% compute the fold change for each metabolite in B and L
tblScalled.areaModel = mdl.Residuals.Raw;
tblScalledUnstack = tblScalled(:, {'peakId' 'sample' 'areaModel'});
tblScalledUnstack = unstack(tblScalledUnstack, 'areaModel', 'sample');

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
tblFoldChanges = unstack(coefficientTableBL(:, {'peakID' 'cell' 'Estimate'}), 'Estimate', 'cell');
% covert to base 2
tblFoldChanges.B = tblFoldChanges.B ./ log(2); 
tblFoldChanges.L = tblFoldChanges.L ./ log(2); 
% add the p-vales
pVals = unstack(coefficientTableBL(:, {'peakID' 'cell' 'pValue'}), 'pValue', 'cell');
tblFoldChanges.pValB = pVals.B;
tblFoldChanges.pValL = pVals.L;
% add the lower confidence intervals
lowerCI = unstack(coefficientTableBL(:, {'peakID' 'cell' 'Lower'}), 'Lower', 'cell');
tblFoldChanges.lowerB = lowerCI.B./ log(2);
tblFoldChanges.lowerL = lowerCI.L./ log(2);
% add the upper confidence intervals
upperCI = unstack(coefficientTableBL(:, {'peakID' 'cell' 'Upper'}), 'Upper', 'cell');
tblFoldChanges.upperB = upperCI.B./ log(2);
tblFoldChanges.upperL = upperCI.L./ log(2);

% save the folds table
if exist('folds', 'dir') == 0
    disp('creating folds directory');
    mkdir('folds');
end
writetable(tblFoldChanges, 'folds/peakFoldChanges.csv',...
    "FileType","text", 'Delimiter', ',');

%% remove spectra that are not in foldchange table (compounds unchanged)
tblSpectra = readtable('extractedPeaks/tblSpectra.csv');
tblSpectra.peakId = arrayfun(@(x) num2str(x), tblSpectra.peakId, 'UniformOutput', false);
tblSpectra(~ismember(tblSpectra.peakId, tblFoldChanges.peakID), :) = [];
% sort the two tables to make sure they have the same order
tblSpectra = sortrows(tblSpectra, "peakId");
tblFoldChanges = sortrows(tblFoldChanges, "peakID");

%% determine number of latent components (metaboLiteLearner.m)
KFOLDLEARN = 0; % k-fold for learning; use 0 for leave-one-out.
MAXN = 30; % maximum number of components to test
NRANDOMIZED = 1000;

mSpectra = tblSpectra{:, 2:end};
X = mSpectra./ sqrt(sum(mSpectra.^2, 2));
Y = [tblFoldChanges.B, tblFoldChanges.L];

obj = MetaboLiteLearner(X, Y, KFOLDLEARN, MAXN, NRANDOMIZED);

% Shuffle the rows of the response variable NRANDOMIZED times
if SHUFFLE_TEST
    subplot(1, 2, 2)
    obj.shufflingTest
end

%% Plot the percentage of variance explained by each component for X and Y
cumulativeVarX = cumsum(obj.PCTVAR(1,:));
cumulativeVarY = cumsum(obj.PCTVAR(2,:));

figure
set(gcf, 'Position', [1000        1062         278         275])
% Plotting for X (predictor matrix)
subplot(2, 1, 1)
bar(1:obj.nopt, obj.PCTVAR(1,:)*100);
xlabel('Latent component');
ylabel({'Explained variance' 'in X (%)'} );
%title('Explained variance for predictors (m/z of metabolite spectra)');
grid on;

% Plotting for Y (response matrix) as well, you can uncomment the lines below
subplot(2, 1, 2)
bar(1:obj.nopt, obj.PCTVAR(2,:)*100);
xlabel('Latent component');
ylabel({'Explained variance' 'in Y (%)'} );
%title('Explained variance for responses (log_2(FC) of metabolite)');
grid on;

%% plot the loadings of the latent compoenet
figure
set(gcf, 'Position', [955   555   839   782]);

% import kegg compounds with spectrum
load kegg/keggCompoundsWithFiehlibSpectrum.mat;

% remove antibiotics, bufanoil and cofactors so that we are left with 7
% categories only
tblKegg3(strcmp(tblKegg3.metClassLevel1, 'Antibiotics'), :) = [];
tblKegg3(strcmp(tblKegg3.metClassLevel1, 'Bufanolide derivatives [Fig]'), :) = [];
tblKegg3(strcmp(tblKegg3.metClassLevel1, 'Vitamins and cofactors'), :) = [];

% normalize the spectra
xKegg = cell2mat(tblKegg3.abundance);
xKegg = xKegg./ sqrt(sum(xKegg.^2, 2));

% mapp Kegg to latent space
[lsKegg, predictedY] = obj.mapToLatentSpace(xKegg);
tblKegg3.latent1 = lsKegg(:,1);
tblKegg3.latent2 = lsKegg(:,2);
tblKegg3.latent3 = lsKegg(:,3);
tblKegg3.latent4 = lsKegg(:,4);
tblKegg3.latent5 = lsKegg(:,5);
tblKegg3.predictedY1 = predictedY(:,1);
tblKegg3.predictedY2 = predictedY(:,2);



% plot the destibution of latent compnent inputes per kegg category
tblKegg3.metClassLevel1 = categorical(tblKegg3.metClassLevel1);

subplot(2, 2, 1)
% plot the beta values for all the mzs
plot([obj.BETA(:, 1)*0 obj.BETA(:, 1)]', [obj.BETA(:, 2)*0 obj.BETA(:, 2)]', 'k.-')
hold on
nBetaLAbels = 2;
mzArray =  50:599;
[~, iSort] = sort((obj.BETA(:, 1).^2 + obj.BETA(:, 2).^2), 'descend');
for i = 1:nBetaLAbels
    text(obj.BETA(iSort(i), 1), obj.BETA(iSort(i), 2),...
        num2str(mzArray(iSort(i))));
end
% plot the fold change predicted by model
%scatter(obj.Ypred(:, 1), obj.Ypred(:, 2), 'ko', 'MarkerFaceColor', [1 0 0], 'MarkerFaceAlpha', 0.2)
xline(0,'k--');
yline(0,'k--');
% plot the fold change predicted by model for KEGG compounds with biolgical
% activity
 h = gscatter(tblKegg3.predictedY1, tblKegg3.predictedY2,...
     tblKegg3.metClassLevel1);
 l  = legend('Location', 'SouthEast')
hold off
xlabel('log_2(FC) in brain-homing cells')
ylabel('log_2(FC) in lung-homing cells')
axis  equal square
grid on
xlim([-2 2]), ylim([-3 3]);


subplot(2, 2, 3)
hMarkers = gscatter(obj.YL(1, :), obj.YL(2, :), 1:obj.nopt, [], [], 40);
legend('AutoUpdate','off')
hold on
hLines = plot([obj.YL(1, :)*0 ;obj.YL(1, :)],...
    [obj.YL(2, :)*0; obj.YL(2, :)], 'k-', 'LineWidth', 2);
for i = 1:length(hLines)
    hLines(i).Color = hMarkers(i).Color;
end
hold off
xlabel('log_2(FC) in brain-homing cells')
ylabel('log_2(FC) in lung-homing cells')
xline(0,'k--');
yline(0,'k--');
axis  equal square
grid on
l  = legend('Location', 'SouthEast')
title(l, 'Latent component','FontSize',12);
xlim([-2 2]), ylim([-3 3]);

% Plot the distribution of latent components in Kegg compounds of
% biological activity
cmap = lines(5);
for i = 1:5
    subplot(5, 2, (i*2))
    latentComponent = i;
    field = sprintf('latent%d', latentComponent);
    swarmchart(tblKegg3, 'metClassLevel1', field,...
        'filled', 'MarkerEdgeColor','none',...
        'MarkerFaceColor', cmap(i, :), 'MarkerFaceAlpha', 0.2);
    hold on
    tblMeans = grpstats(tblKegg3, 'metClassLevel1', 'mean', 'DataVars', field);
    %plot(1:7, tblMeans{:, 3}, '+')
    h = plot(tblMeans, 'metClassLevel1', 3);
    set(h, 'Marker', '_', 'MarkerSize', 20, 'LineWidth', 4,...
        'Color', cmap(i, :), 'LineStyle', 'none')
    hold off
    title('')
    xlabel('');
    ylabel(sprintf('Input into\nlatent component %d', i));
    grid on
    yline(0,'k-');
    ylim([-0.4 0.6])
    labels = get(gca, 'XTickLabels');
    set(gca, 'XTickLabels', []);
end
set(gca, 'XTickLabels', labels);