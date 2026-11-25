"""Handle writing of outputs"""

import logging
import os
from pathlib import Path

import pandas as pd

from prymer.models import Result, Settings

logger = logging.getLogger(__name__)

def write_output(result: Result, output_dir: Path) -> Path | None:
    """Write results of one Job to disk. Returns the path or None if no primers."""
    logger.info("%s: %s", result.job.name, result.explain)
    if not result.pairs:
        logger.warning("%s: Primer3 found no valid primers", result.job.name)
        return None
    path = output_dir / f"{result.job.name}.csv"
    pd.DataFrame.from_dict(result.pairs).to_csv(path)
    return path

def write_outputs(results: list[Result], settings: Settings) -> list[Path]:
    paths = (write_output(r, settings.output_dir) for r in results)
    return [p for p in paths if p is not None]

def validate_output_path(path_str: str) -> Path:
    """Validate the output directory exists and is writable.
    Parameters:
        path_str: String containing user-specified output directory path
    Returns:
          Output directory path.
    Raises:
        FileNotFoundError: If output directory does not exist.
        PermissionError: If output directory is not writable.
        """
    path = Path(path_str)
    if not path.exists():
        raise FileNotFoundError(f"Output path does not exist: {path}")
    if not os.access(path, os.W_OK):
        raise PermissionError(f"Output path is not writable: {path}")
    return path