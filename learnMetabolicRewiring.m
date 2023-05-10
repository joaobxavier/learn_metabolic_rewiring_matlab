%% load spectra for all peaks
tblSpectra = readtable('extractedPeaks/tblSpectra.csv');
KFOLDLEARN = 0; % k-fold for learning; use 0 for leave-one-out.
MAXN = 20; % maximum number of components to test
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

obj = plsrStandAloneLearner05102023(X, Y, KFOLDLEARN, MAXN, NRANDOMIZED);
obj.shufflingTest

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
fprintf('C1 = (%f,%f)\n', obj.YL(1,1), obj.YL(2,1))

%what is C2?
fprintf('C2 = (%f,%f)\n', obj.YL(1,2), obj.YL(2,2))

%% get the score for component 2 of all identified metabolites
tblIndentified = readtable('identifiedFiehnLib/tblIdentity.csv');
tblW2Score = table(tblFoldChanges.peakID, obj.XS(:, 2));
tblW2Score.Properties.VariableNames = {'peakID' 'w2Score'};
tblW2Score = innerjoin(tblW2Score, tblIndentified, "LeftKeys","peakID", "RightKeys","peakId");
tblW2Score(~tblW2Score.bestMatch95Rt, :) = [];
tblW2Score = sortrows(tblW2Score, "w2Score","descend")


