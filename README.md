# Metabolic rewiring discovery with MetaboLiteLearner, a lightweight and efficient PLSR-based statistical learning algorithm.

MATLAB code to learn metabolic rewiring from GCMS data.

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

### learnMetabolicRewiring

Uses the `plsrLearner` class to analyze metabolic rewiring in cells based on their metabolite spectra and fold changes. It determines the optimum number of components (latent variables) for the PLSR model and plots the learning loss and the predictions of the model with the optimal number of components.

#### Data used

- `tblSpectra.csv`: A table containing the spectra for all peaks.
- `peakFoldChanges.csv`: A table containing fold changes for all peaks.

#### Workflow

1. Load data.
2. Determine the optimum number of components (latent variables) for the PLSR model.
3. Plot the training and evaluation loss, and compare predictions from the optimal model with real data.

#### Output

- A plot showing the training and evaluation loss versus the number of latent components in the PLSR model.
- A scatter plot comparing the true metabolite abundances with the best model predictions.

### identifyFirstEnrichmentAnalysis
This script generates bar and scatter plots to visualize the fold changes of identified metabolites in bone-homing and lung-homing cells. 


## Classes
### PLSR Learner

The `plsrLearner` class is a MATLAB class for learning associations between the spectrum of a metabolite and a response variable using Partial Least Square Regression (PLSR). The class is designed to handle the training, evaluation, and optimization of PLSR models.

#### Properties

- `xData`: The spectra for each metabolite.
- `yData`: The response variable.

#### Methods

##### plsrLearner(x, y)

Constructor for the `plsrLearner` class. Initializes an instance with the given spectra `x` and response variable `y`.

##### learn(x, y, n)

Performs PLSR with `n` components using input `x` and response variable `y`. Returns the following:

- `BETA`: Regression coefficients.
- `Ypred`: Predicted values of the response variable.
- `loss`: Squared difference between predicted and true values.
- `sse`: Sum of squared errors.

##### leaveOneOutEvaluation(n)

Evaluates a PLSR model with `n` components using leave-one-out cross-validation. Returns the following:

- `Ypred`: Predicted values of the response variable.
- `trainSse`: Mean sum of squared errors for training data.
- `leaveOneOutSse`: Sum of squared errors for leave-one-out evaluation.

##### optimizeComponentsAndLearn(maxn)

Determines the optimal number of components (up to `maxn`) that minimizes the loss of leave-one-out evaluations. Returns the following:

- `nopt`: Optimal number of components.
- `Ypred`: Predicted values of the response variable for the optimal model.
- `trainSse`: Mean sum of squared errors for training data for each tested number of components.
- `leaveOneOutSse`: Sum of squared errors for leave-one-out evaluation for each tested number of components.

