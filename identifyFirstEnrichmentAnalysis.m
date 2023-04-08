% load the fold changes
tblFolds = readtable('folds/peakFoldChanges.csv');

% load the peak identities
tblPeakIdentities = readtable('identifiedFiehnLib/tblIdentity.csv');

% keep only the best matches
tblPeakIdentities(tblPeakIdentities.bestMatch95Rt==0, :) = [];

% do an inner join to keep only the peaks identified
tblFoldsIdentified = innerjoin(tblFolds, tblPeakIdentities,...
    "LeftKeys","peakID", "RightKeys","peakId");

%% cleanup metabolite names
tblFoldsIdentified.name = tblFoldsIdentified.bestMatchFiehnLib;
for i = 1:height(tblFoldsIdentified)
    str = tblFoldsIdentified.name{i};
    idxSpaces = strfind(str, ' ');
    str = str((idxSpaces(1)+1):(idxSpaces(end)-1));
    tblFoldsIdentified.name{i} = str;
end

%% make a bar plot to show fold changes in identified compounds in
% brain-hominc and lung-homing cells
figure(1)
barh([tblFoldsIdentified.B tblFoldsIdentified.L])
set(gca, 'YTick', 1:height(tblFoldsIdentified), 'YTickLabel', tblFoldsIdentified.name)
legend({'Bone-homing' 'Lung-homing'})
xlabel('log_2(FC)')

%% plot the fold-changes in a scatter plot with each metabolite colored by
% p-value
figure(2)
errorbar(tblFoldsIdentified.B, tblFoldsIdentified.L,...
    tblFoldsIdentified.L - tblFoldsIdentified.lowerL, tblFoldsIdentified.upperL - tblFoldsIdentified.L,...
    tblFoldsIdentified.B - tblFoldsIdentified.lowerB, tblFoldsIdentified.upperB - tblFoldsIdentified.B, '.', 'Color', [1 1 1]*0.8)
hold on;
plot(tblFolds.B, tblFolds.L, 'ko')
scatter(tblFoldsIdentified.B, tblFoldsIdentified.L, [], -log10(tblFoldsIdentified.pValB .* tblFoldsIdentified.pValL), 'filled', 'MarkerEdgeColor','k')
hold off
xlabel('Brain-homing vs. parental [log_2(FC)]')
ylabel('Lung-homing vs. parental [log_2(FC)]')
grid on
axis equal square
h = refline(1, 0)
set(h, 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--')
set(gca, 'CLim', [1.3 2])
% correlation between brain-homing and lung-homing 
[r, p] = corr(tblFoldsIdentified.B, tblFoldsIdentified.L);
title(sprintf('%d compounds; \\rho = %0.2f; P-value = %0.0e',...
    height(tblFoldsIdentified), r, p));

%% add text labels to metabolites significantly changed
idx = tblFoldsIdentified.pValB .* tblFoldsIdentified.pValL < 0.05;
text(tblFoldsIdentified.B(idx), tblFoldsIdentified.L(idx),...
    tblFoldsIdentified.name(idx))


