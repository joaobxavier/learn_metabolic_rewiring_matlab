# Agilent Conversion Comparison Report

## Test Case

- Sample: `B_plate1_run1_10302020`
- Raw input: `rawAgilentData/metastasis_lineages_3replicates_3runs_10302020/B_plate1_run1_10302020.D`
- Generated CSV: `testoutputs/B_plate1_run1_10302020.csv`
- Reference CSV copy: `testoutputs/B_plate1_run1_10302020.reference.csv`

## Conversion Command

The conversion was run from the Python package in the local virtual environment at `~/dev/venvs/metabolite-learner` using `metabolite_learner.agilent.convert_agilent_to_csv(...)` with `sample_dirs=[sample]` so only this single `.D` directory was processed.

## Matrix Comparison Summary

- Generated shape: `2401 x 550`
- Reference shape: `2401 x 550`
- Same shape: `yes`
- NaN count, generated: `0`
- NaN count, reference: `0`
- Exact equality including NaNs: `no`
- `np.allclose(..., atol=1e-9, rtol=0)`: `no`
- `np.allclose(..., atol=1e-6, rtol=0)`: `no`
- `np.allclose(..., atol=1e-3, rtol=0)`: `yes`
- Finite entries compared: `1,320,550`
- Max absolute difference: `1.7294660210609436e-06`
- Mean absolute difference: `6.200057905807928e-11`
- Median absolute difference: `0.0`
- RMSE: `3.618852419222505e-09`
- Entries with abs diff `> 1e-9`: `5,610`
- Entries with abs diff `> 1e-6`: `3`
- Entries with abs diff `> 1e-3`: `0`

## File-Level Summary

- Generated file size: `13,942,061` bytes
- Reference file size: `11,516,705` bytes
- Generated SHA-256: `38b9f27b59f85939ac49138c1ca129eac4b8767db61727d012d93b4957434720`
- Reference SHA-256: `b73cf842b2f63449479d08a71cae7075fedf79d9f57c8e8150aeaa9d4fb91522`
- Generated matrix sum: `3252204990.5386596`
- Reference matrix sum: `3252204990.5386496`
- Generated min/max: `0.0 / 8388096.0`
- Reference min/max: `0.0 / 8388096.0`

## Largest Numeric Differences

| Row | Column | Generated | Reference | Abs diff |
| --- | --- | ---: | ---: | ---: |
| 1120 | 155 | 5483825.834667459 | 5483825.83466573 | 1.7294660210609436e-06 |
| 1113 | 23 | 4889894.127659852 | 4889894.12765822 | 1.6316771507263184e-06 |
| 1126 | 269 | 3159866.02666755 | 3159866.02666654 | 1.0100193321704865e-06 |
| 1120 | 23 | 2939913.5573337837 | 2939913.5573328 | 9.83942300081253e-07 |
| 1113 | 269 | 4677599.31914877 | 4677599.31914975 | 9.806826710700989e-07 |

## Interpretation

The generated matrix is not byte-identical to the checked-in CSV, but the numerical agreement is extremely close. All entries match within an absolute tolerance of `1e-3`, and the worst observed difference is about `1.7e-6`. The larger generated file size and differing file hash are consistent with different float serialization rather than a material change in the signal matrix.
