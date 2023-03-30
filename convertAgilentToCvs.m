function [] = convertAgilentToCvs(dataDir, outputDir)

% [] = convertAgilentToCvs(dataDir, outputDir)
% Convert Agilent .D files into CVS text files. The GCMS data is 
% represented as a matrix of 2,401 x 550, where the 2,401 rows represent 
% the time steps from 6 to 30 minutes (in 0.01 minute intervals) and the 
% 500 columns represent the mz steps from 50 to 599 (in unit intervals).
% 
% dataDir - path to the directory containing .D files
% outputDir - path to the directory to save the CVS files
%
% Requires chromatography-master in the path
% https://www.mathworks.com/matlabcentral/fileexchange/47696-chromatography-toolbox
%
% Joao Xavier (xavierj@mskcc.org) March 2022


% check if chromatography-master exists
if  isempty(which('ImportAgilent'))
    error('Requires chromatography-master in the path');
end

% check if outputDir exists. If not, create it.
if exist(outputDir, 'dir') == 0
    disp(['creating ' outputDir]);
    mkdir(outputDir);
end

% list the .D files in a directory
fileList = dir([dataDir '/*.D']);

% load each file
h = waitbar(0,'Converting to CSV...');
for i = 1:length(fileList)
    fn     = [dataDir '/' fileList(i).name];
    sample = importAgilentAndResample(fn);
    sample.name = fileList(i).name(1:end-2);
    writematrix(sample.xic, [outputDir '/' sample.name '.csv']);
    disp(['saved ' outputDir '/' sample.name '.csv'])
    waitbar(i/length(fileList),h)
end
close(h);
