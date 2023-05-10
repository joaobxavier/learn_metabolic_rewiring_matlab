# MetaboLiteLearner: a lightweight machine learning algorithm to learn metabolic rewiring from GC/MS data

## Functions

### convertAgilentToCvs(dataDir, outputDir)

Converts Agilent .D files into CSV text files. The GCMS data is represented as a 2,401 x 550 matrix, with 2,401 rows representing the time steps from 6 to 30 minutes (in 0.01-minute intervals) and 550 columns representing the m/z steps from 50 to 599 (in unit intervals). The function requires the chromatography-toolbox in the path.
* `dataDir`: Path to the directory containing .D files.
* `outputDir`: Path to the directory to save the CSV files.

### extractSpectraAndIntegrate(csvDataDir, outputSpectraDir)

Analyzes GCMS data in the provided directory containing CSV files. It identifies peaks in the total ion chromatogram (TIC), extracts their spectra, and integrates the peaks to create a peak area table. The output tables, tblPeaksIntegrated and tblSpectra, are saved as CSV files in the specified outputSpectraDir directory.
- `csvDataDir`: The directory containing the GCMS data in CSV format. Each CSV file should represent a sample.
- `outputSpectraDir`: The directory where the output files, `tblPeaksIntegrated.csv` and `tblSpectra.csv`, will be saved.

## Scripts

### scaleGcmsPeaks

This is a script, not a function. Therefore, it doesn't have input arguments. Loads the peak areas table, scales each sample using a mixed effects model, does a PCA to compare the data before and after scaling, then uses another mixed effects model to calculate the log-2 fold change of each peak detected. Saves the result as a table in a new directory 'folds'.

### unsupervisedAnalysis.m

This script performs unsupervised analysis on the metabolite data, exploring the structure within the explanatory variable X and examining correlations within the response variable Y. It generates two figures:

1. Correlation between brain-homing and lung-homing fold changes.
2. PCoA plots for GC/MS peaks, log2(FC) in brain-homing, and log2(FC) in lung-homing.



### identifyCompounds 
Identify Compounds with GC/MS Data and FiehnLib. Matches GC/MS spectra from the `extractedPeaks` folder against the FiehnLib mass spectral library. It assesses best matches based on cosine similarity and the retention time window.

#### Script Workflow

1. Load the spectra from the `tblSpectra.csv` file in the `extractedPeaks` folder.
2. Import the FiehnLib mass spectral library.
3. Compute the cosine similarity between the input spectra and the FiehnLib spectra.
4. Identify the best matches based on cosine similarity.
5. Create a new table with the best matches and their properties (name, CAS number, retention time).
6. Calibrate the retention time using a robust linear fit model, considering only matches with excellent similarity (>= 95%).
7. Add a flag to the table indicating if a match has >= 95% similarity and is within the 95% confidence retention time window.
8. Generate two plots:
   - A scatter plot showing the calibrated retention times of the input spectra and the best FiehnLib matches.
   - A scatter plot ranking the cosine similarity of the identified peaks, highlighting good matches.
9. Save the identified compound table in the `identifiedFiehnLib` folder as `tblIdentity.csv`.

### importMsl

This script imports a mass spectral library (MSL) file (Fiehn-2013.MSL) and creates a table, *massSpectralLibrary*, with information about each compound in the library. *massSpectralLibrary* has 12 columns, which store various properties of the compounds like name, molecular weight, CAS number, retention index, retention time, number of peaks, m/z values, and abundance values. The spectra are represented as arrays of ion abundances between m/z 50 and 599. Column *abundance73* represents the abundance of each compound at m/z 73.

### learnMetabolicRewiring.m

The `learnMetabolicRewiring.m` script analyzes metabolic rewiring in cancer cells using partial least squares regression (PLSR). It processes metabolomics data by loading and organizing spectra and fold changes, determining the optimal number of latent variables for the model, and comparing predictions with data. The script also explores the significance of the model by comparing it to shuffled data, visualizes the learned model by examining identified metabolites, and investigates latent components and their associated spectra. The analysis is performed using the `plsrLearner` class for PLSR.

## Classes
### PLSR Standalone Learner (plsrStandAloneLearner05102023)

This MATLAB class uses Partial Least Squares Regression (PLSR) to establish associations between the spectrum of a metabolite and a response variable.

#### Features

- **xFullData**: The spectra for each metabolite.
- **yFullData**: The response variable.
- **kfold**: The k-fold used to estimate loss and optimize components.
- **maxn**: The maximum number of components tested.
- **nrandomized**: The number of shuffling iterations used to test.
- **cvIndices**: Cross-validation indices.
- **nopt**: The optimal number of components.
- **testSse**: Sum of squared errors of held-out examples (evaluation).
- **BETA**: Coefficients of linear model trained with all data.
- **Ypred**: Fit results.
- **XL, YL**: Loadings of the X and Y.
- **XS, YS**: Scores of input and output data mapped onto nopt-dimensional latent space.

#### Usage

```MATLAB
plsrObject = plsrStandAloneLearner05102023(x, y, kfold, maxn, nrandomized)
```
The class constructor plsrStandAloneLearner05102023 takes as input parameters:

x: The input spectra for each metabolite.
y: The response variable.
kfold: The k-fold used to estimate loss and optimize components.
maxn: The maximum number of components tested.
nrandomized: The number of shuffling iterations used to test.
The class has methods to perform shuffling tests, learn the PLSR model, evaluate the model using k-fold cross-validation and optimize the model's components.

Please refer to the code documentation for more details on each method.


