# learn_metabolic_rewiring_matlab
MATLAB code to learn metabolic rewiring from GCMS data.

## convertAgilentToCvs(dataDir, outputDir)
Converts Agilent .D files into CVS text files. Requires [chromatography-toolbox](https://www.mathworks.com/matlabcentral/fileexchange/47696-chromatography-toolbox) in the path.

## extractSpectraAndIntegrate(csvDataDir, outputSpectraDir)
Detects peaks in the total ion chromatogram, extract their spectra and integrates the peaks to make a peak area table.

## scaleGcmsPeaks
This is a script, not a function. Therefore, it doesn't have input arguments. Loads the peak areas table, scales each sample using a mixed effects model, does a PCA to compare the data before and after scaling, then uses another mixed effects model to calculate the log-2 fold change of each peak detected. Saves the result as a table in a new directory 'folds'.

## indentifyCompounds
This script loads two data tables: spectra data from a CSV file (tblSpectra.csv) and fold changes data from another CSV file (peakFoldChanges.csv). It imports a mass spectral library (FiehnLib) and finds the best match for each spectrum by maximizing cosine similarity and stores the results in the tblFoldChanges table. Plots the results in a two-part figure: The first plot displays the cosine similarity for each peak, with the top 31 matches highlighted (cosine similarity 0.95). The second plot shows the fold changes of brain-homing and lung-homing for the top 31 matches as horizontal bars. Finally, it saves the identified table as a CSV file in the 'identifiedFiehnLib' directory.

## importMsl
This script imports a mass spectral library (MSL) file (Fiehn-2013.MSL) and creates a table, *massSpectralLibrary*, with information about each compound in the library. *massSpectralLibrary* has 12 columns, which store various properties of the compounds like name, molecular weight, CAS number, retention index, retention time, number of peaks, m/z values, and abundance values. The spectra are represented as arrays of ion abundances between m/z 50 and 599. Column *abundance73* represents the abundance of each compound at m/z 73.
