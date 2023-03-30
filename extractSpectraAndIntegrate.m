function [tblPeaksIntegrated, tblSpectra] = extractSpectraAndIntegrate(csvDataDir, outputSpectraDir)

%extractSpectraAndIntegrate(csvDataDir, outputSpectraDir)
%   Analyzes the GCMS data in the directory csvDataDir. The data should be
% converted into CSV files. Indetifies peaks in the total ion chromatogram,
% extract their spectra and integrates the peaks to make a peak area table.

% initialize tblSpectra
time = (6:0.01:30)';
mz   = (50:599);
[mzMesh, timeMesh] = meshgrid(mz, time);
tblSpectra = table();
tblSpectra.mz = mzMesh(:);
tblSpectra.time = timeMesh(:);

%% load the CSV files
fileList = dir([csvDataDir '/*.csv']);
if isempty(fileList)
    error(['No CSV files in ' csvDataDir])
end

h = waitbar(0,'Loading CSV...');
sampleNames = {};
for i = 1:length(fileList)
    xic = readmatrix( [csvDataDir '/' fileList(i).name]);
    % add a new columns to the tblSpectra
    sName = fileList(i).name(1:end-4);
    tblSpectra{:, {sName}} = xic(:);
    sampleNames{end+1} = sName;
    %
    waitbar(i/length(fileList),h)
end
close(h);

%% make individual samples TIC
tblSampleTic = grpstats(tblSpectra,'time','sum','DataVars',sampleNames);
tblSampleTic.GroupCount = [];
tblSampleTic.Properties.VariableNames(2:end) = sampleNames;

%% make the bulk TIC table
tblTic = tblSampleTic(:,1);
tblTic.tic = sum(tblSampleTic{:,sampleNames}, 2);

%% find peaks in bulk TIC
[peakTicList, peakTicWHH, ~] = mspeaks(tblTic.time, tblTic.tic);
tblTicPeaks = table(peakTicList(:, 1), peakTicList(:, 2), peakTicWHH(:, 1), peakTicWHH(:, 2));
tblTicPeaks.Properties.VariableNames = {'rtPeak' 'intensity' 'rtStart' 'rtEnd' };

% create a peak Id
tblTicPeaks.peakId = (1:height(tblTicPeaks))';

%% assign peakIds to all the retention time values in the TICtable
a = bsxfun(@gt,tblTic.time,tblTicPeaks.rtStart');
b = bsxfun(@lt,tblTic.time,tblTicPeaks.rtEnd');
tblRtPeakID = tblTic(:, 1);
tblRtPeakID.peakId = sum((a&b) .* tblTicPeaks.peakId', 2);
% remove RTs outside any peak
tblRtPeakID(tblRtPeakID.peakId == 0, :) = [];

%% integrate the TIC peaks of each sample
tblPeaksIntegrated = grpstats(innerjoin(tblSampleTic, tblRtPeakID),...
    'peakId','sum','DataVars',sampleNames);
tblPeaksIntegrated.GroupCount = [];
tblPeaksIntegrated.Properties.VariableNames(2:end) = sampleNames;

%% extract the spectra for each peak
% make the bulk chromatogram
tblSummedSpectra = tblSpectra(:, 1:2);
tblSummedSpectra.allSamples = sum(tblSpectra{:,sampleNames}, 2);
tblSummedSpectra = unstack(tblSummedSpectra, 'allSamples', "mz");

tblSpectra = grpstats(innerjoin(tblSummedSpectra, tblRtPeakID),...
    'peakId','sum','DataVars',2:width(tblSummedSpectra));
tblSpectra.GroupCount = [];
columnNames = tblSummedSpectra.Properties.VariableNames(2:end);
columnNames = strrep(columnNames, 'x', 'mz');
tblSpectra.Properties.VariableNames(2:end) = columnNames;

%% check if outputDir exists. If not, create it.
if exist(outputSpectraDir, 'dir') == 0
    disp(['creating ' outputSpectraDir]);
    mkdir(outputSpectraDir);
end

%% save output
writetable(tblPeaksIntegrated, [outputSpectraDir '/tblPeaksIntegrated.csv'],...
    "FileType","text", 'Delimiter', ',');
writetable(tblSpectra, [outputSpectraDir '/tblSpectra.csv'],...
    "FileType","text", 'Delimiter', ',');

end