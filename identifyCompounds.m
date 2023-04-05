% load spectra
tblSpectra = readtable("extractedPeaks/tblSpectra.csv");

% load fold changes
tblFoldChanges = readtable("folds/peakFoldChanges.csv");

% remove peaks that are absent in fold change table (peaks unchanges in
% media vs cells)
tblSpectra(~ismember(tblSpectra.peakId, tblFoldChanges.peakID), :) = [];

% make sure both tables are ordered the same way
tblFoldChanges = sortrows(tblFoldChanges, "peakID");
tblSpectra = sortrows(tblSpectra, "peakId");

%% load FiehnLib
importMsl;
mz = 50:599;
mFiehLib = cell2mat(massSpectralLibrary.abundance);
rtFiehnLib = massSpectralLibrary.RT;

%% find best match by cosine distance
mGcms = tblSpectra{:, 2:end};
cosineDist = pdist2(mGcms, mFiehLib, 'cosine');
cosineSimilarity = 1 - cosineDist;
[cosineSimilarityOfBestMatch, bestMatch] = max(cosineSimilarity');

%%
tblFoldChanges.cosineSimilarityOfBestMatch = cosineSimilarityOfBestMatch';
tblFoldChanges.bestMatchFiehnLib = massSpectralLibrary.NAME(bestMatch);
tblFoldChanges.bestMatchCasNo = massSpectralLibrary.CASNO(bestMatch);

%%
tblFoldChangesRanked = sortrows(tblFoldChanges, "cosineSimilarityOfBestMatch", "descend");
n = 31;

figure(1)
subplot(2, 1, 1)
plot(tblFoldChangesRanked.cosineSimilarityOfBestMatch, 'o-')
hold on
plot(tblFoldChangesRanked.cosineSimilarityOfBestMatch(1:n), 'ro-')
hold off
title('Matching GC/MS peaks with FiehnLib');
xlabel('Peak')
ylabel('Cosine similarity')
grid on

tblFoldChangesRanked = tblFoldChangesRanked(1:n, :);
subplot(2, 1, 2)
barh([tblFoldChangesRanked.B tblFoldChangesRanked.L]);
set(gca, 'YTick', 1:height(tblFoldChangesRanked),...
    'YTickLabel', tblFoldChangesRanked.bestMatchFiehnLib)

%% save the identified table
if exist('identifiedFiehnLib', 'dir') == 0
    disp('creating folds directory');
    mkdir('identifiedFiehnLib');
end
writetable(tblFoldChanges, 'identifiedFiehnLib/tblFoldChanges.csv',...
    "FileType","text", 'Delimiter', ',');
