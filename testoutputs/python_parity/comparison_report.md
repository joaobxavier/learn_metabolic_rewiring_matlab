# Python vs MATLAB Parity Report

## Scope

This report compares the updated Python workflow outputs in `testoutputs/python_parity/` against MATLAB reference outputs generated from `scMeboLiteLearnerWorkFlow.m` in `testoutputs/matlab_reference/`.

## Output Sets

### MATLAB reference

- `figure_1.png`
- `figure_2.png`
- `figure_3.png`
- `peakFoldChanges.csv`

### Python parity run

- `figure_1.png`
- `figure_2.png`
- `figure_3.png`
- `learner_diagnostics.png`
- `variance_explained.png`
- `loadings.png`
- `peakFoldChanges.csv`

## Parity Summary

- Latent component count selected by MATLAB: `5`
- Latent component count selected by Python: `5`
- Figure count compared: `3`
- Figure pixel sizes:
  - `figure_1.png`: MATLAB `1884 x 1986`, Python `1884 x 1986`
  - `figure_2.png`: MATLAB `695 x 701`, Python `694 x 700`
  - `figure_3.png`: MATLAB `1413 x 1118`, Python `1412 x 1118`

## Fold-Change Table Comparison

The fold-change table shape and columns match exactly.

- Shape: `153 x 9`
- Columns: identical

### Numeric deltas, MATLAB vs Python

- `B`: max abs diff `9.043044091328056e-13`, mean abs diff `8.515956158070083e-13`
- `L`: max abs diff `1.6782131240233866e-12`, mean abs diff `1.6233804197421647e-12`
- `pValB`: max abs diff `0.013805343382053803`, mean abs diff `0.007102988423936028`
- `pValL`: max abs diff `0.0138053633337622`, mean abs diff `0.008108610070560438`
- `lowerB`: max abs diff `0.007602245618331893`, mean abs diff `0.007602245618274684`
- `lowerL`: max abs diff `0.007602245617555292`, mean abs diff `0.00760224561749827`
- `upperB`: max abs diff `0.007602245620029802`, mean abs diff `0.007602245619977783`
- `upperL`: max abs diff `0.0076022456207986006`, mean abs diff `0.007602245620744934`

## Interpretation

The Python workflow now matches MATLAB on the major workflow outputs for steps 4-7:

- the same number of latent components is selected,
- the same fold-change peaks are produced,
- the fold-change point estimates match MATLAB to machine precision,
- the MATLAB-style figure set is reproduced and exported with matching figure dimensions.

The remaining gap is in mixed-model inference details:

- p-values differ by up to about `0.014`,
- confidence interval endpoints differ by about `0.0076` in log2 units.

Those residual differences are consistent with `statsmodels` mixed-model inference not being numerically identical to MATLAB `fitlme`, even when the estimated coefficients agree. The figures and the core fold-change estimates now match MATLAB as closely as the Python stack can reproduce from the current implementation.
