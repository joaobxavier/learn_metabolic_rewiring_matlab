# MetaboLiteLearner: Python port of the metabolic rewiring GC/MS workflow

This repository now includes a Python implementation of the original MATLAB-based MetaboLiteLearner workflow. The Python code mirrors the same top-level stages:

1. convert Agilent sample directories into the repository's matrix CSV format,
2. extract TIC peaks and bulk spectra,
3. remove peaks that are not significantly different from media,
4. fit mixed-effects models to estimate corrected fold changes,
5. learn latent components with partial least squares regression, and
6. generate the same workflow summary artifacts in the `folds/` directory.

The original MATLAB files are still present as reference implementations, but the supported workflow entry point is now the Python package in `metabolite_learner/`.

## Python package layout

- `metabolite_learner/agilent.py`: Agilent conversion step with an injectable raw-data loader.
- `metabolite_learner/extract.py`: TIC peak detection, peak integration, and spectrum export.
- `metabolite_learner/workflow.py`: end-to-end analysis pipeline, including ANOVA, mixed models, and plotting.
- `metabolite_learner/pls.py`: `MetaboLiteLearner` partial least squares implementation.
- `metabolite_learner/cli.py`: command-line interface.

## Installation

The Python port targets Python 3.10+ and depends on:

- `numpy`
- `pandas`
- `scipy`
- `scikit-learn`
- `statsmodels`
- `matplotlib`

Install the project in editable mode:

```bash
python3 -m pip install -e .
```

## Usage

### 1. Extract peaks and spectra from the existing GC/MS matrix CSV files

```bash
python3 -m metabolite_learner extract-peaks gcmsCSVs extractedPeaks
```

This writes:

- `extractedPeaks/tblTicPeaks.csv`
- `extractedPeaks/tblPeaksIntegrated.csv`
- `extractedPeaks/tblSpectra.csv`

### 2. Run the full workflow

```bash
python3 -m metabolite_learner run-workflow \
  --gcms-csv-dir gcmsCSVs \
  --extracted-peaks-dir extractedPeaks \
  --folds-dir folds
```

This produces:

- `folds/peakFoldChanges.csv`
- `folds/variance_explained.png`
- `folds/loadings.png`

### 3. Agilent raw-data conversion

```bash
python3 -m metabolite_learner convert-agilent rawAgilentData/... gcmsCSVs
```

Unlike the MATLAB version, the Python implementation does not hard-code a dependency on `chromatography-master`. Instead, `convert_agilent_to_csv` accepts a pluggable loader because Agilent `.D` directories are proprietary vendor containers. Out of the box, the function looks for a pre-exported matrix CSV inside each `.D` directory.

## Notes on parity with MATLAB

- The statistical workflow follows the same sequence as `scMeboLiteLearnerWorkFlow.m`.
- `statsmodels` is used for the OLS ANOVA and mixed-effects models.
- `scikit-learn` is used for the PLS learner and cross-validation.
- The Python peak extractor uses `scipy.signal.find_peaks` and half-height peak widths as a replacement for MATLAB's `mspeaks`.
- KEGG plotting remains best-effort because the `.mat` file contains a MATLAB table structure whose representation can vary by SciPy version.

## Data

The repository still contains the processed GC/MS matrix CSVs in `gcmsCSVs/`, the extracted peak tables in `extractedPeaks/`, and the reference fold-change table in `folds/`.

## License

This project is licensed under the MIT License. See `LICENSE` for details.
