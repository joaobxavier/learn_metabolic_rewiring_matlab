from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.io import loadmat
from statsmodels.stats.anova import anova_lm

from .extract import extract_spectra_and_integrate
from .pls import MetaboLiteLearner


@dataclass(slots=True)
class WorkflowResult:
    fold_changes: pd.DataFrame
    spectra: pd.DataFrame
    learner: MetaboLiteLearner


def _infer_cell_type(sample_name: str) -> str:
    return sample_name[0]


def _ols_anova_p_value(values: np.ndarray, cell_types: list[str]) -> float:
    frame = pd.DataFrame({"value": np.log(values), "cell": pd.Categorical(cell_types)})
    model = smf.ols("value ~ C(cell)", data=frame).fit()
    return float(anova_lm(model).iloc[0]["PR(>F)"])


def _prepare_scaled_table(peaks_integrated: pd.DataFrame) -> pd.DataFrame:
    melted = peaks_integrated.melt(id_vars="peakId", var_name="sample", value_name="area")
    melted["sample"] = melted["sample"].astype(str)
    melted["batch"] = pd.Categorical(
        np.select(
            [melted["sample"].str.contains("plate1"), melted["sample"].str.contains("plate2")],
            ["1", "2"],
            default="3",
        )
    )
    melted["cell"] = pd.Categorical(melted["sample"].map(_infer_cell_type), categories=["P", "B", "L"])
    melted["peakId"] = melted["peakId"].astype(str)
    melted["peak_batch"] = melted["peakId"] + ":" + melted["batch"].astype(str)
    melted["area"] = np.log(melted["area"].astype(float))
    return melted


def _fit_scaling_model(frame: pd.DataFrame):
    return smf.mixedlm(
        "area ~ 0 + C(peakId)",
        frame,
        groups=frame["sample"],
        vc_formula={"peak_batch": "0 + C(peak_batch)"},
    ).fit(reml=False, method="lbfgs", disp=False)


def _fit_fold_change_model(frame: pd.DataFrame):
    return smf.mixedlm(
        "area ~ 0 + C(peakId) + C(peakId):C(cell)",
        frame,
        groups=frame["sample"],
        vc_formula={"peak_batch": "0 + C(peak_batch)"},
    ).fit(reml=False, method="lbfgs", disp=False)


def _coefficient_frame(model) -> pd.DataFrame:
    conf = model.conf_int()
    params = model.params
    pvalues = model.pvalues
    records = []
    for name, estimate in params.items():
        if not name.startswith("C(peakId)") or ":C(cell)" not in name:
            continue
        peak_part, cell_part = name.split(":")
        peak_id = peak_part.split("[")[-1].rstrip("]")
        cell = cell_part.split("[")[-1].rstrip("]").replace("T.", "")
        records.append(
            {
                "peakID": peak_id,
                "cell": cell,
                "Estimate": float(estimate),
                "pValue": float(pvalues[name]),
                "Lower": float(conf.loc[name, 0]),
                "Upper": float(conf.loc[name, 1]),
            }
        )
    return pd.DataFrame.from_records(records)


def _plot_variance_explained(learner: MetaboLiteLearner, output_dir: Path) -> Path:
    figure, axes = plt.subplots(2, 1, figsize=(4, 5), constrained_layout=True)
    axes[0].bar(range(1, learner.nopt + 1), learner.pctvar[0] * 100)
    axes[0].set_xlabel("Latent component")
    axes[0].set_ylabel("Explained variance\nin X (%)")
    axes[0].grid(True)
    axes[1].bar(range(1, learner.nopt + 1), learner.pctvar[1] * 100)
    axes[1].set_xlabel("Latent component")
    axes[1].set_ylabel("Explained variance\nin Y (%)")
    axes[1].grid(True)
    path = output_dir / "variance_explained.png"
    figure.savefig(path, dpi=200)
    plt.close(figure)
    return path


def _safe_extract_kegg_table(mat_path: Path) -> pd.DataFrame | None:
    if not mat_path.exists():
        return None
    payload = loadmat(mat_path, squeeze_me=True, struct_as_record=False)
    table = payload.get("tblKegg3")
    if table is None:
        return None
    try:
        records = []
        iterable = table if np.ndim(table) else [table]
        for row in iterable:
            met_class = getattr(row, "metClassLevel1", None)
            abundance = getattr(row, "abundance", None)
            if met_class is None or abundance is None:
                continue
            records.append({"metClassLevel1": str(met_class), "abundance": np.asarray(abundance, dtype=float)})
        return pd.DataFrame.from_records(records)
    except Exception:
        return None


def _plot_loadings(learner: MetaboLiteLearner, output_dir: Path, kegg_mat: Path) -> Path:
    figure, axes = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)
    axes[0].plot(learner.beta[:, 0], learner.beta[:, 1] if learner.beta.shape[1] > 1 else learner.beta[:, 0], "k.-")
    axes[0].axhline(0, color="black", linestyle="--", linewidth=1)
    axes[0].axvline(0, color="black", linestyle="--", linewidth=1)
    axes[0].set_xlabel("Brain-homing coefficient")
    axes[0].set_ylabel("Lung-homing coefficient")
    axes[0].grid(True)

    if learner.y_loadings.shape[1] >= 2:
        axes[1].scatter(learner.y_loadings[0], learner.y_loadings[1], c=np.arange(learner.nopt), cmap="viridis")
    else:
        axes[1].scatter(np.arange(learner.nopt), learner.y_loadings[0], c=np.arange(learner.nopt), cmap="viridis")
    axes[1].axhline(0, color="black", linestyle="--", linewidth=1)
    axes[1].set_xlabel("Component")
    axes[1].set_ylabel("Y loading")
    axes[1].grid(True)

    kegg_table = _safe_extract_kegg_table(kegg_mat)
    if kegg_table is not None and not kegg_table.empty:
        kegg_table = kegg_table.loc[
            ~kegg_table["metClassLevel1"].isin(
                ["Antibiotics", "Bufanolide derivatives [Fig]", "Vitamins and cofactors"]
            )
        ].copy()
        spectra = np.vstack(kegg_table["abundance"].to_numpy())
        spectra = spectra / np.linalg.norm(spectra, axis=1, keepdims=True)
        _, predicted = learner.map_to_latent_space(spectra)
        axes[0].scatter(predicted[:, 0], predicted[:, 1], alpha=0.35, s=10)

    path = output_dir / "loadings.png"
    figure.savefig(path, dpi=200)
    plt.close(figure)
    return path


def run_workflow(
    gcms_csv_dir: str | Path = "gcmsCSVs",
    extracted_peaks_dir: str | Path = "extractedPeaks",
    folds_dir: str | Path = "folds",
    kegg_mat_path: str | Path = "kegg/keggCompoundsWithFiehlibSpectrum.mat",
    *,
    regenerate_peaks: bool = False,
    kfold_learn: int = 0,
    max_components: int = 30,
    nrandomized: int = 1000,
) -> WorkflowResult:
    gcms_csv_dir = Path(gcms_csv_dir)
    extracted_peaks_dir = Path(extracted_peaks_dir)
    folds_dir = Path(folds_dir)
    folds_dir.mkdir(parents=True, exist_ok=True)

    peaks_path = extracted_peaks_dir / "tblPeaksIntegrated.csv"
    spectra_path = extracted_peaks_dir / "tblSpectra.csv"
    if regenerate_peaks or not peaks_path.exists() or not spectra_path.exists():
        peaks_integrated, spectra = extract_spectra_and_integrate(gcms_csv_dir, extracted_peaks_dir)
    else:
        peaks_integrated = pd.read_csv(peaks_path)
        spectra = pd.read_csv(spectra_path)

    sample_columns = list(peaks_integrated.columns[1:])
    cell_types = [_infer_cell_type(name) for name in sample_columns]
    p_values = [
        _ols_anova_p_value(peaks_integrated.loc[row_index, sample_columns].to_numpy(dtype=float), cell_types)
        for row_index in range(len(peaks_integrated))
    ]
    filtered = peaks_integrated.loc[np.array(p_values) < 0.05, ["peakId", *[name for name in sample_columns if not name.startswith("M")]]].copy()

    scaled = _prepare_scaled_table(filtered)
    scaling_model = _fit_scaling_model(scaled)
    scaled["areaModel"] = scaling_model.resid

    fold_change_model = _fit_fold_change_model(scaled)
    coefficients = _coefficient_frame(fold_change_model)
    if coefficients.empty:
        raise RuntimeError("The mixed-effects fold-change model did not produce any B/L coefficients.")

    fold_changes = coefficients.pivot(index="peakID", columns="cell", values="Estimate").reset_index().rename_axis(columns=None)
    pvals = coefficients.pivot(index="peakID", columns="cell", values="pValue").reset_index().rename_axis(columns=None)
    lower = coefficients.pivot(index="peakID", columns="cell", values="Lower").reset_index().rename_axis(columns=None)
    upper = coefficients.pivot(index="peakID", columns="cell", values="Upper").reset_index().rename_axis(columns=None)

    fold_changes = fold_changes.merge(pvals, on="peakID", suffixes=("", "_p"))
    fold_changes = fold_changes.merge(lower, on="peakID", suffixes=("", "_lower"))
    fold_changes = fold_changes.merge(upper, on="peakID", suffixes=("", "_upper"))

    log2_factor = np.log(2.0)
    for column in ["B", "L", "B_lower", "L_lower", "B_upper", "L_upper"]:
        if column in fold_changes:
            fold_changes[column] = fold_changes[column] / log2_factor

    fold_changes = fold_changes.rename(
        columns={
            "B_p": "pValB",
            "L_p": "pValL",
            "B_lower": "lowerB",
            "L_lower": "lowerL",
            "B_upper": "upperB",
            "L_upper": "upperL",
        }
    )
    fold_changes = fold_changes[["peakID", "B", "L", "pValB", "pValL", "lowerB", "lowerL", "upperB", "upperL"]]
    fold_changes = fold_changes.sort_values("peakID").reset_index(drop=True)
    fold_changes.to_csv(folds_dir / "peakFoldChanges.csv", index=False)

    spectra = spectra.copy()
    spectra["peakId"] = spectra["peakId"].astype(str)
    fold_changes["peakID"] = fold_changes["peakID"].astype(str)
    spectra = spectra.loc[spectra["peakId"].isin(fold_changes["peakID"]), :].sort_values("peakId").reset_index(drop=True)
    fold_changes = fold_changes.sort_values("peakID").reset_index(drop=True)

    matrix_spectra = spectra.iloc[:, 1:].to_numpy(dtype=float)
    x = matrix_spectra / np.linalg.norm(matrix_spectra, axis=1, keepdims=True)
    y = fold_changes[["B", "L"]].to_numpy(dtype=float)
    learner = MetaboLiteLearner(x, y, kfold=kfold_learn, max_components=max_components, nrandomized=nrandomized)

    _plot_variance_explained(learner, folds_dir)
    _plot_loadings(learner, folds_dir, Path(kegg_mat_path))

    return WorkflowResult(fold_changes=fold_changes, spectra=spectra, learner=learner)
