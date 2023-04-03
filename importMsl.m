% import a MSL file (mass spectral library)

fn = 'Library/Fiehn-2013.MSL';

% get the number of compounds in MSL file
[status, cmdout] =  system(['cat "' fn '" | grep "NAME:" | wc -l']);
nCompounds = str2double(cmdout);

%% create table and prealocate number of rows
massSpectralLibrary = table('Size', [nCompounds 12], 'VariableTypes',...
    repmat({'cell'}, [1 12]));
massSpectralLibrary.Properties.VariableNames{1} = 'NAME';
massSpectralLibrary.Properties.VariableNames{2} = 'FORM';
massSpectralLibrary.Properties.VariableNames{3} = 'MW';
massSpectralLibrary.Properties.VariableNames{4} = 'CASNO';
massSpectralLibrary.Properties.VariableNames{5} = 'RI';
massSpectralLibrary.Properties.VariableNames{6} = 'RT';
massSpectralLibrary.Properties.VariableNames{7} = 'RF';
massSpectralLibrary.Properties.VariableNames{8} = 'RSN';
massSpectralLibrary.Properties.VariableNames{9} = 'COMMENT';
massSpectralLibrary.Properties.VariableNames{10} = 'numberOfPeaks';
massSpectralLibrary.Properties.VariableNames{11} = 'mz';
massSpectralLibrary.Properties.VariableNames{12} = 'abundance';

%% Iterate through every compound in the text file and import spectrum
% open the text file
fid = fopen(fn);
% iterate through each compound
tline = fgetl(fid);
h = waitbar(0,'Importing the spectra from Fiehn library...');
for i = 1:nCompounds
    waitbar(i/nCompounds, h);
    % skip until next compound
    while ~contains(tline, 'NAME:')
        tline = fgetl(fid);
    end
    % cycle through the first 9 lines of the compound
    for j = 1:9
        variableName = massSpectralLibrary.Properties.VariableNames{j};
        massSpectralLibrary{i, variableName} = {tline(length(variableName)+3:end)};
        tline = fgetl(fid);
    end
    % Get number of peaks
    massSpectralLibrary.numberOfPeaks{i} = tline(13:end);
    npeaks = str2double(tline(13:end));
    % Get rows containing ions and intensities
    v = [];
    while size(v, 1) < npeaks
        tline = fgetl(fid);
        tline = strrep(tline, ') (', ';');
        tline = strrep(tline, '(', '[');
        tline = strrep(tline, ')', ']');
        v = [v; eval(tline)];
    end
    % some compounds have repat peaks. Add dum here:
    v2 = [];
    v2(:, 1) = unique(v(:, 1));
    v2(:, 2) = grpstats(v(:, 2), v(:, 1), 'sum');
    % save the spectra as arrays of abundances between mz 50 and 599
    massSpectralLibrary.mz{i} = 50:599;
    massSpectralLibrary.abundance{i} =  zeros(size(50:599));
    massSpectralLibrary.abundance{i}(ismember(50:599, v2(:, 1))) = v2(:, 2); 
end
fclose(fid);
close(h)

%% convert some of the variables to number
massSpectralLibrary.MW = str2double(massSpectralLibrary.MW);
massSpectralLibrary.RI = str2double(massSpectralLibrary.RI);
massSpectralLibrary.RT = str2double(massSpectralLibrary.RT);
massSpectralLibrary.RF = [];
massSpectralLibrary.numberOfPeaks = str2double(massSpectralLibrary.numberOfPeaks);

%% add another variable to designate if this is a derivatized compound
m = cell2mat([massSpectralLibrary.abundance]);
idx = find(massSpectralLibrary.mz{1} == 73);
massSpectralLibrary.abundance73 = m(:, idx);




