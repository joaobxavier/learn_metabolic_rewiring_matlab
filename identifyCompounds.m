% load spectra
tblSpectra = readtable("extractedPeaks/tblSpectra.csv");

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

%% create a new table with those best matches
tblIdentity = table();
tblIdentity.peakId = tblSpectra.peakId;
tblIdentity.cosineSimilarityOfBestMatch = cosineSimilarityOfBestMatch';
tblIdentity.bestMatchFiehnLib = massSpectralLibrary.NAME(bestMatch);
tblIdentity.bestMatchCasNo = massSpectralLibrary.CASNO(bestMatch);
tblIdentity.bestMatchRt = massSpectralLibrary.RT(bestMatch);

%% compare peak retention time with best match retention time


% get indices of peaks with excellent matches (>95% similarity)
idx95 = tblIdentity.cosineSimilarityOfBestMatch>=0.95;

% calibrate retention time using robust lienar fit with excellect mathces
mdlRt = fitlm(tblIdentity, 'bestMatchRt ~ peakId', 'Exclude',~idx95, 'RobustOpts','on');
[ypred, ypredci] = predict(mdlRt, tblIdentity, 'Alpha', 0.01);

% add a flag stating if best hit is >=95% similar *and* is inside 
% the 95% confidence retention time window
tblIdentity.bestMatch95Rt = idx95 &...
    tblIdentity.bestMatchRt >= ypredci(:, 1) &...
    tblIdentity.bestMatchRt <= ypredci(:, 2);

% plot the model and matches
figure(1)
subplot(2, 1, 1)


scatter(tblIdentity.peakId(~idx95), tblIdentity.bestMatchRt(~idx95),...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor','w');
hold on
scatter(tblIdentity.peakId(idx95&~tblIdentity.bestMatch95Rt),...
    tblIdentity.bestMatchRt(idx95&~tblIdentity.bestMatch95Rt),...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor','y');
scatter(tblIdentity.peakId(tblIdentity.bestMatch95Rt),...
    tblIdentity.bestMatchRt(tblIdentity.bestMatch95Rt),...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor','r');
lgd = legend({'Weak match', 'Outside RT window' 'Good match'},...
    'Location','northwest');
set(lgd, 'AutoUpdate', 'off');
plot(tblIdentity.peakId, predict(mdlRt, tblIdentity), 'r:')
plot(tblIdentity.peakId, ypredci(:, 1), 'r-')
plot(tblIdentity.peakId, ypredci(:, 2), 'r-')
xlabel('Peak retention time [min]');
ylabel('Retention time of best FiehnLib match [min]')
grid on;

% draw the linear fit
hold off;
title('Calibrating the retention time of GC/MS with FiehnLib');

% plot the number of peks with good match
tblIdentityRanked = sortrows(tblIdentity, "cosineSimilarityOfBestMatch", "descend");
tblIdentityRanked.rank = (1:height(tblIdentityRanked))';
idx95 = tblIdentityRanked.cosineSimilarityOfBestMatch>=0.95 % indices of good matches

subplot(2, 1, 2)
scatter(tblIdentityRanked.rank(~idx95),...
    tblIdentityRanked.cosineSimilarityOfBestMatch(~idx95),...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor','w');
hold on
scatter(tblIdentityRanked.rank(idx95&~tblIdentity.bestMatch95Rt),...
    tblIdentityRanked.cosineSimilarityOfBestMatch(idx95&~tblIdentity.bestMatch95Rt),...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor','y');
scatter(tblIdentityRanked.rank(tblIdentityRanked.bestMatch95Rt),...
    tblIdentityRanked.cosineSimilarityOfBestMatch(tblIdentityRanked.bestMatch95Rt),...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor','r');
hold off
%
title(sprintf('Only %d out of %d peaks match well with FiehnLib',...
    sum(tblIdentityRanked.bestMatch95Rt), height(tblIdentityRanked)));
xlabel('Peak')
ylabel('Cosine similarity')
grid on

%% save the identified table
if exist('identifiedFiehnLib', 'dir') == 0
    disp('creating folds directory');
    mkdir('identifiedFiehnLib');
end
writetable(tblIdentity, 'identifiedFiehnLib/tblIdentity.csv',...
    "FileType","text", 'Delimiter', ',');