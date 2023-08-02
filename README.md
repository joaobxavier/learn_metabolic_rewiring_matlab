# MetaboLiteLearner: a lightweight machine learning algorithm to learn metabolic rewiring from GC/MS data

This repository contains a MATLAB script `scMeboLiteLearnerWorkFlow.m` which implements the entire MetaboLiteLearner workflow. Below is the detailed description of the script.

## How to use the script
This script can be executed in MATLAB environment. The script contains the steps to preprocess the data, compute the fold change for each metabolite, and plot the percentage of variance explained by each component. The default settings are set to skip workflow steps, but you can change this by setting variables to true. For example, if you set `DOWNLOAD = true;`, the script will download the raw data from Zenodo.

The script is divided into several steps:

1. **Download the raw data from Zenodo:** This step downloads the raw Agilent data from the Zenodo repository and saves it in a directory named 'rawAgilentData'. 
   
2. **Convert Agilent files to CSV:** This step imports the Agilent data, bins it, and exports it as CSV files. It requires the chromatography-master in the path.

3. **Extract Spectra and Integrate:** This step extracts the spectra and integrates it, generating two tables: `tblPeaksIntegrated` and `tblSpectra`.

4. **Load peak data:** This step loads the previously generated tables and extracts the sample type from the sample name.

5. **Remove peaks that are not different from media:** This step performs ANOVA tests to remove peaks that are not significantly different from media.

6. **Rescale samples:** This step rescales the samples using the intracellular metabolome. 

7. **Fit model and calculate corrected peak areas:** This step uses linear mixed effects models to calculate corrected peak areas.

8. **Compute the fold change for each metabolite in B and L:** This step computes the fold change for each metabolite in B and L and saves the folds table.

9. **Remove spectra that are not in foldchange table (compounds unchanged):** This step removes spectra that are not present in the fold change table.

10. **Determine number of latent components (metaboLiteLearner.m):** This step determines the number of latent components using the MetaboLiteLearner.

11. **Plot the percentage of variance explained by each component for X and Y:** This step generates a bar plot to visualize the percentage of variance explained by each component for X and Y.

12. **Plot the loadings of the latent component:** This step plots the loadings of the latent component. 

## Required Dependencies
To use this script, the following dependencies are required:
- MATLAB: A high-level language and interactive environment developed by MathWorks.
- Chromatography-Master (for the conversion of Agilent data to CSV).

## Data
The data used in this script is from Zenodo repository: https://zenodo.org/record/8193580/

## Output
The output of this script includes a variety of graphs and CSV files related to the MetaboLiteLearner workflow, including the percentage of variance explained by each component and the fold changes of each metabolite in B and L.

## Functions

### convertAgilentToCvs(dataDir, outputDir)

Converts Agilent .D files into CSV text files. The GCMS data is represented as a 2,401 x 550 matrix, with 2,401 rows representing the time steps from 6 to 30 minutes (in 0.01-minute intervals) and 550 columns representing the m/z steps from 50 to 599 (in unit intervals). The function requires the chromatography-toolbox in the path.
* `dataDir`: Path to the directory containing .D files.
* `outputDir`: Path to the directory to save the CSV files.

### extractSpectraAndIntegrate(csvDataDir, outputSpectraDir)

Analyzes GCMS data in the provided directory containing CSV files. It identifies peaks in the total ion chromatogram (TIC), extracts their spectra, and integrates the peaks to create a peak area table. The output tables, tblPeaksIntegrated and tblSpectra, are saved as CSV files in the specified outputSpectraDir directory.
- `csvDataDir`: The directory containing the GCMS data in CSV format. Each CSV file should represent a sample.
- `outputSpectraDir`: The directory where the output files, `tblPeaksIntegrated.csv` and `tblSpectra.csv`, will be saved.

## Classes
### MetaboLiteLearner
The `MetaboLiteLearner` class employs partial least square regression (PLSR) to study associations between the spectrum of a metabolite and a corresponding response variable.

##### Properties

###### Input Properties
- **xFullData**: Spectra for each metabolite.
- **yFullData**: Response variable.
- **kfold**: k-fold used for estimating loss and optimizing components.
- **maxn**: Maximum number of components tested.
- **nrandomized**: Number of shuffling iterations for testing.

###### State Variables
- **cvIndices**: Cross-validation indices.
- **nopt**: Optimum number of components.
- **testSse**: Sum of squared errors of held-out examples for evaluation.
- **BETA**: Coefficients of the linear model trained with all data.
- **PCTVAR**: Percent variance explained.
- **Ypred**: Fit results.
- **XL**: Loadings of the X.
- **YL**: Loadings of the Y.
- **XS**: Scores of input data mapped onto `nopt`-dimensional latent space.
- **YS**: Scores of output data mapped onto `nopt`-dimensional latent space.
- **stats**: Output of the `plsregress`.

##### Methods

###### Constructor
- `MetaboLiteLearner(x, y, kfold, maxn, nrandomized)`: Constructs an instance of the class.

###### Public Methods
- `mapToLatentSpace(x)`: Maps a set of spectra to the latent space.
- `shufflingTest()`: Performs a shuffling test.
- `learn(x, y, n)`: Executes PLSR with `n` components.
- `crossValidationEvaluation(x, y, n)`: Evaluates a PLSR model with 'n' components using k-fold cross-validation.
- `optimizeComponentsAndLearn()`: Determines the number of components that minimizes the evaluation loss.
