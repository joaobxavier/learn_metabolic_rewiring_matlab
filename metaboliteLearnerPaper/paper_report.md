# MetaboLiteLearner Paper Report

Paper PDF: [metaboliteLearnerPaper.pdf](/Users/jxavier/dev/metaboliteLearner/metaboliteLearnerPaper/metaboliteLearnerPaper.pdf)

## Executive Summary

The paper presents **MetaboLiteLearner** as a supervised-learning workflow for studying **metabolic rewiring** from GC/MS scan-mode data without requiring explicit metabolite identification. The central idea is to treat each detected GC/MS peak as a metabolite-like object represented by a **550-dimensional electron-impact fragmentation spectrum** and to learn a mapping from that spectrum to **two response variables**: the log2 fold changes of that peak in **brain-homing** and **lung-homing** metastatic breast cancer cells relative to the parental line.

The method is intentionally lightweight. It avoids spectral deconvolution and avoids dependence on curated pathway structure during training. Instead, it builds a peak table from a virtual bulk total-ion chromatogram, removes peaks not distinguishable from blank media, estimates corrected fold changes with mixed-effects models, and then trains a **partial least squares regression (PLSR)** model to connect fragment patterns to metabolic adaptation.

The paper’s main claim is not that the model replaces metabolite identification, but that useful biological signal about metabolic rewiring can be recovered from **fragmentation structure alone**. In that sense, the method treats EI spectra as compressed structural fingerprints that carry enough information to generalize from observed metabolites to withheld ones.

## What the Paper Says

From the paper, the workflow is:

1. culture parental, lung-homing, brain-homing, and media samples
2. extract and derivatize metabolites for GC/MS
3. convert raw Agilent `.D` files into matrix-form CSVs
4. create a “virtual bulk sample” and detect TIC peaks directly, without deconvolution
5. extract a summed spectrum for each peak and compute per-sample peak areas
6. remove peaks that are not significantly different from blank media
7. compute log2 fold changes for brain- and lung-homing cells
8. train a PLSR model from spectra to fold changes
9. choose the latent dimensionality by leave-one-out cross-validation
10. interpret the learned loadings and latent projections, including KEGG/Fiehn reference spectra

The paper explicitly frames this as an answer to two common metabolomics bottlenecks:

- untargeted workflows often produce many peaks that are hard to identify confidently
- supervised metabolomics often suffers from small sample sizes and weak interpretability

MetaboLiteLearner addresses that by using a low-dimensional latent PLS model that remains interpretable through:

- fragment coefficients (`BETA`)
- latent component loadings
- predicted fold-change projections for biologically annotated compounds

## Interpretation of the Method

The method is best understood as a **structure-to-rewiring learner**. It does not infer biochemical mechanisms directly. Instead, it learns statistical associations between:

- molecular fragmentation patterns, which indirectly encode structural chemistry
- phenotypic abundance shifts associated with metastatic adaptation

That is a narrow but useful target. If structurally similar metabolites tend to be up- or down-regulated together under a given adaptive program, then EI spectra should contain enough regularity for a supervised model to exploit.

This is also why the method can work without compound identification. Identification is useful for interpretation, but the learning problem itself only requires:

- a consistent spectrum representation
- reliable peak-level abundance changes
- enough coherence between chemical structure and biological response

The paper’s strongest conceptual contribution is showing that **fragmentation-space learning** can act as a surrogate for metabolite annotation when the practical question is prediction of rewiring rather than formal metabolite naming.

## Strengths of the Paper

1. The method is operationally simple.
   It uses standard GC/MS scan-mode data and a linear supervised learner rather than a heavy end-to-end black box.

2. The model is interpretable.
   PLS loadings, latent components, and projected KEGG/Fiehn spectra make the results biologically discussable.

3. It is pragmatic about metabolite identification.
   The paper does not require perfect annotation before extracting signal from the data.

4. It uses a relevant biological benchmark.
   The MDA-MB-231 parental/brain/lung system is a meaningful metastasis model with prior metabolomics characterization.

## Limits and Cautions

1. The training labels still come from a fairly specific experimental design.
   The model is predictive within this biological system, but the paper does not establish broad out-of-domain generalization.

2. The method is peak-centric, not compound-centric.
   One chromatographic peak is treated as a learning unit. That is practical, but it can blur the distinction between true single metabolites and unresolved signal mixtures.

3. No-deconvolution is both a strength and a tradeoff.
   It reduces pipeline complexity and dependence on deconvolution quality, but it can also leave mixed spectral content in the learned representation.

4. PLSR is interpretable but still sensitive to preprocessing.
   Peak detection windows, normalization, and fold-change estimation materially affect the learned latent space.

## How the MATLAB Implementation Encodes the Paper

The MATLAB workflow is the closest direct implementation of the paper.

Main orchestrator:
- [scMeboLiteLearnerWorkFlow.m](/Users/jxavier/dev/metaboliteLearner/scMeboLiteLearnerWorkFlow.m)

### Step 1: Raw Agilent conversion

MATLAB converts `.D` folders to fixed-grid CSV matrices with:
- [convertAgilentToCvs.m](/Users/jxavier/dev/metaboliteLearner/convertAgilentToCvs.m)

Important details:
- uses `ImportAgilent` from `chromatography-master`
- bins to a `2401 x 550` matrix
- interpolates onto the fixed retention-time grid `6:0.01:30`

### Step 2: Peak detection and spectra extraction

Implemented in:
- [extractSpectraAndIntegrate.m](/Users/jxavier/dev/metaboliteLearner/extractSpectraAndIntegrate.m)

Important details:
- builds a bulk TIC from all sample matrices
- detects peaks with MATLAB `mspeaks`
- integrates peak areas per sample
- sums spectra per peak window

Persistent outputs:
- [`extractedPeaks/tblPeaksIntegrated.csv`](/Users/jxavier/dev/metaboliteLearner/extractedPeaks/tblPeaksIntegrated.csv)
- [`extractedPeaks/tblSpectra.csv`](/Users/jxavier/dev/metaboliteLearner/extractedPeaks/tblSpectra.csv)

### Step 3: Statistical filtering and fold-change estimation

Implemented directly in:
- [scMeboLiteLearnerWorkFlow.m](/Users/jxavier/dev/metaboliteLearner/scMeboLiteLearnerWorkFlow.m#L54)

Important details:
- ANOVA removes peaks not distinguishable from media
- `fitlme` is used twice:
  - first for scaling / residual correction
  - second for estimating B and L fold changes plus p-values and confidence intervals

Persistent output:
- [`folds/peakFoldChanges.csv`](/Users/jxavier/dev/metaboliteLearner/folds/peakFoldChanges.csv)

### Step 4: PLS learning and diagnostics

Implemented in:
- [MetaboLiteLearner.m](/Users/jxavier/dev/metaboliteLearner/MetaboLiteLearner.m)

Important details:
- uses MATLAB `plsregress`
- uses leave-one-out when `KFOLDLEARN = 0`
- selects `nopt` with the “minimum + 1 standard error” rule
- exposes loadings, scores, coefficients, and shuffle testing

### Step 5: Interpretation figures

The MATLAB workflow generates three figures interactively:

- model diagnostics / model-fit figure
- explained-variance figure
- loadings + KEGG projection figure

Reference exported MATLAB figures from this repo:
- [figure_1.png](/Users/jxavier/dev/metaboliteLearner/testoutputs/matlab_reference/figure_1.png)
- [figure_2.png](/Users/jxavier/dev/metaboliteLearner/testoutputs/matlab_reference/figure_2.png)
- [figure_3.png](/Users/jxavier/dev/metaboliteLearner/testoutputs/matlab_reference/figure_3.png)

## How the Python Implementation Maps to the Paper

The Python implementation now reproduces most of the MATLAB workflow structure and artifacts.

Python entry points:
- [README.md](/Users/jxavier/dev/metaboliteLearner/README.md)
- [metabolite_learner/cli.py](/Users/jxavier/dev/metaboliteLearner/metabolite_learner/cli.py)

### Step 1: Raw Agilent conversion

Implemented in:
- [metabolite_learner/agilent.py](/Users/jxavier/dev/metaboliteLearner/metabolite_learner/agilent.py)

Interpretation:
- the Python version is more self-contained than MATLAB because it can decode `data.ms` directly
- this removes the hard dependency on `chromatography-master`
- this is a real engineering improvement over the paper’s original MATLAB dependency chain

### Step 2: Peak extraction

Implemented in:
- [metabolite_learner/extract.py](/Users/jxavier/dev/metaboliteLearner/metabolite_learner/extract.py)

Difference from MATLAB:
- Python uses `scipy.signal.find_peaks` and `peak_widths`
- MATLAB uses `mspeaks`

Interpretation:
- the representation is conceptually the same
- exact peak boundaries can differ slightly
- so parity should be judged at the level of workflow outputs, not bitwise identity of intermediate windows

### Step 3: Filtering and fold-change estimation

Implemented in:
- [metabolite_learner/workflow.py](/Users/jxavier/dev/metaboliteLearner/metabolite_learner/workflow.py#L363)

Interpretation:
- Python follows the same logic as the paper and MATLAB:
  - ANOVA filtering
  - mixed-effects scaling
  - mixed-effects fold-change estimation
- the fitted **B** and **L** estimates match MATLAB essentially exactly in the parity runs
- p-values and confidence intervals remain slightly different because `statsmodels` mixed-model inference is not numerically identical to MATLAB `fitlme`

### Step 4: PLS learner

Implemented in:
- [metabolite_learner/pls.py](/Users/jxavier/dev/metaboliteLearner/metabolite_learner/pls.py)

Interpretation:
- the current Python version uses a MATLAB-aligned SIMPLS-style implementation rather than the earlier `scikit-learn` PLS shortcut
- this was necessary to match the MATLAB latent-component selection behavior
- after parity fixes, both MATLAB and Python select `nopt = 5` on the repo dataset

### Step 5: KEGG/Fiehn interpretation layer

Implemented in:
- [metabolite_learner/workflow.py](/Users/jxavier/dev/metaboliteLearner/metabolite_learner/workflow.py#L158)

Important implementation note:
- the KEGG/Fiehn reference data in MATLAB lives in a `.mat` table that SciPy cannot reliably decode in this repo
- the Python implementation now uses a fallback CSV export:
  - [keggCompoundsWithFiehlibSpectrum.csv](/Users/jxavier/dev/metaboliteLearner/kegg/keggCompoundsWithFiehlibSpectrum.csv)

That fix is important because otherwise the right-hand KEGG analysis panels in the Python `figure_3` equivalent would render empty.

## What Outputs the Python Implementation Produces

Current Python parity outputs:
- [peakFoldChanges.csv](/Users/jxavier/dev/metaboliteLearner/testoutputs/python_parity/peakFoldChanges.csv)
- [learner_diagnostics.png](/Users/jxavier/dev/metaboliteLearner/testoutputs/python_parity/learner_diagnostics.png)
- [variance_explained.png](/Users/jxavier/dev/metaboliteLearner/testoutputs/python_parity/variance_explained.png)
- [loadings.png](/Users/jxavier/dev/metaboliteLearner/testoutputs/python_parity/loadings.png)
- [figure_1.png](/Users/jxavier/dev/metaboliteLearner/testoutputs/python_parity/figure_1.png)
- [figure_2.png](/Users/jxavier/dev/metaboliteLearner/testoutputs/python_parity/figure_2.png)
- [figure_3.png](/Users/jxavier/dev/metaboliteLearner/testoutputs/python_parity/figure_3.png)

These map approximately to the MATLAB exported references:
- [figure_1.png](/Users/jxavier/dev/metaboliteLearner/testoutputs/matlab_reference/figure_1.png)
- [figure_2.png](/Users/jxavier/dev/metaboliteLearner/testoutputs/matlab_reference/figure_2.png)
- [figure_3.png](/Users/jxavier/dev/metaboliteLearner/testoutputs/matlab_reference/figure_3.png)

## Assessment of MATLAB vs Python in This Repo

### MATLAB implementation

Strengths:
- closest direct realization of the paper
- uses the exact MATLAB functions cited by the article, especially `plsregress` and `mspeaks`
- naturally handles the MATLAB `.mat` table used for KEGG/Fiehn context

Weaknesses:
- depends on external toolbox code (`ImportAgilent`)
- less convenient for reproducible automation outside MATLAB
- figures are interactive by default rather than exported by the paper script itself

### Python implementation

Strengths:
- easier to automate and test
- no hard dependency on MATLAB toolboxes for routine execution
- direct `.ms` decoding for Agilent files
- explicit on-disk output figures and parity artifacts

Weaknesses:
- mixed-model inferential statistics still differ slightly from MATLAB
- exact peak detection behavior differs because `find_peaks` is not `mspeaks`
- `.mat` table interoperability required a fallback export to keep the KEGG analysis working

## Overall Interpretation

The paper is methodologically modest and that is part of its value. It does not introduce a deep model or a broad metabolomics foundation model. Instead, it shows that one can get meaningful, biologically interpretable rewiring predictions from a disciplined linear workflow built on:

- consistent peak extraction
- careful fold-change estimation
- a low-dimensional supervised latent model

That makes the paper useful in two ways:

1. as a practical workflow for GC/MS scan-mode data when identification is incomplete
2. as a conceptual argument that **fragmentation patterns themselves contain enough structural signal to support prediction of adaptive metabolic behavior**

In this repository, the MATLAB implementation remains the canonical expression of the paper. The Python implementation has matured into a practical parity port that reproduces the same major outputs and figures, with the main remaining discrepancy confined to mixed-model p-values and confidence intervals rather than the core fold-change estimates or the latent-space interpretation.

## Recommended Reading of the Repo

If you want to understand the paper through the code, read in this order:

1. [scMeboLiteLearnerWorkFlow.m](/Users/jxavier/dev/metaboliteLearner/scMeboLiteLearnerWorkFlow.m)
2. [extractSpectraAndIntegrate.m](/Users/jxavier/dev/metaboliteLearner/extractSpectraAndIntegrate.m)
3. [MetaboLiteLearner.m](/Users/jxavier/dev/metaboliteLearner/MetaboLiteLearner.m)
4. [metabolite_learner/workflow.py](/Users/jxavier/dev/metaboliteLearner/metabolite_learner/workflow.py)
5. [metabolite_learner/pls.py](/Users/jxavier/dev/metaboliteLearner/metabolite_learner/pls.py)
6. [testoutputs/python_parity/comparison_report.md](/Users/jxavier/dev/metaboliteLearner/testoutputs/python_parity/comparison_report.md)

## Bottom Line

The paper’s core claim is credible and well matched to the code in this repository: spectral structure alone can be used to learn predictable aspects of metabolic rewiring. The MATLAB implementation expresses the original method directly. The Python implementation now reproduces the main workflow, the key tables, and the figure suite closely enough to serve as a practical, testable port for further work.
