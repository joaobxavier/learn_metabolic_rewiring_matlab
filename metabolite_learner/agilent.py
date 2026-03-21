from __future__ import annotations

from pathlib import Path
from typing import Callable, Iterable

import numpy as np
import pandas as pd

TIME_GRID = np.round(np.arange(6.0, 30.0 + 0.01, 0.01), 2)
MZ_GRID = np.arange(50, 600, 1)
EXPECTED_SHAPE = (TIME_GRID.size, MZ_GRID.size)

Loader = Callable[[Path], np.ndarray]


def _identity_matrix_loader(sample_path: Path) -> np.ndarray:
    candidates = [
        sample_path / "matrix.csv",
        sample_path / "xic.csv",
        sample_path / f"{sample_path.stem}.csv",
    ]
    for candidate in candidates:
        if candidate.exists():
            data = pd.read_csv(candidate, header=None).to_numpy(dtype=float)
            return data
    raise FileNotFoundError(
        "No matrix CSV was found inside the Agilent sample directory. "
        "The Python port expects a pluggable loader that can decode the proprietary "
        ".D directory into the 2401x550 GC/MS matrix used elsewhere in the repo."
    )


def convert_agilent_to_csv(
    data_dir: str | Path,
    output_dir: str | Path,
    loader: Loader | None = None,
    sample_dirs: Iterable[Path] | None = None,
) -> list[Path]:
    """Convert Agilent sample directories into the repo's CSV matrix format.

    The original MATLAB implementation relied on `ImportAgilent` from the external
    chromatography-master toolbox. The Python port exposes the same workflow step but
    makes the raw-data decoder injectable because Agilent `.D` directories are a
    proprietary vendor format.

    Parameters
    ----------
    data_dir:
        Directory that contains Agilent `.D` sample directories.
    output_dir:
        Destination for the generated CSV matrix files.
    loader:
        Callable that receives a sample directory and returns a 2401x550 matrix.
        When omitted, the function looks for a pre-exported matrix CSV inside each
        sample directory.
    sample_dirs:
        Optional explicit iterable of sample directories to process.
    """

    input_root = Path(data_dir)
    output_root = Path(output_dir)
    output_root.mkdir(parents=True, exist_ok=True)

    if sample_dirs is None:
        sample_dirs = sorted(path for path in input_root.iterdir() if path.suffix == ".D")
    else:
        sample_dirs = sorted(sample_dirs)

    if not sample_dirs:
        raise FileNotFoundError(f"No Agilent sample directories (*.D) found in {input_root}.")

    matrix_loader = loader or _identity_matrix_loader
    written_files: list[Path] = []

    for sample_dir in sample_dirs:
        matrix = np.asarray(matrix_loader(sample_dir), dtype=float)
        if matrix.shape != EXPECTED_SHAPE:
            raise ValueError(
                f"Expected matrix shape {EXPECTED_SHAPE} for {sample_dir.name}, "
                f"but received {matrix.shape}."
            )
        sample_name = sample_dir.stem
        destination = output_root / f"{sample_name}.csv"
        pd.DataFrame(matrix).to_csv(destination, header=False, index=False)
        written_files.append(destination)

    return written_files
