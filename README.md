# learn_metabolic_rewiring_matlab
MATLAB code to learn metabolic rewiring from GCMS data.

<b>convertAgilentToCvs(dataDir, outputDir)</b>
Converts Agilent .D files into CVS text files. Requires <a href="https://www.mathworks.com/matlabcentral/fileexchange/47696-chromatography-toolbox">chromatography-toolbox</a> in the path.


<b>extractSpectraAndIntegrate(csvDataDir, outputSpectraDir)</b>
detifies peaks in the total ion chromatogram, extract their spectra and integrates the peaks to make a peak area table.
