"""Central project locations and read-only resolution of historical paths."""
from __future__ import annotations

from pathlib import Path
import sys

PROJECT = Path(__file__).resolve().parents[1]
PYTHON = PROJECT / ".venv" / "bin" / "python"
HUMANIZE = PROJECT / "core" / "humanization" / "humanize.py"
BACK = PROJECT / "core" / "backmutation" / "back.py"
INPUTS = PROJECT / "data" / "inputs"
RUNS = PROJECT / "results" / "analysis"
HUMAN_DATABASE = PROJECT / "data" / "databases" / "human_imgt" / "human_germline_frameworks"
MOUSE_DATABASE = PROJECT / "data" / "databases" / "mouse_ncbi" / "mouse_V"
CATALOGS = PROJECT / "data" / "databases" / "catalogs.json"

# This maps locations only; saved metadata, XML and original results are not rewritten.
_RELOCATIONS = (
    ("lucianommartins/imgt_germline_database", "data/databases/human_imgt"),
    ("lucianommartins/imgt_download", "data/databases/human_imgt/sources"),
    ("200-back/ncbi-igblast/internal_data/mouse_imgt_preserved_20261008", "data/databases/mouse_imgt"),
    ("200-back/ncbi-igblast/internal_data/mouse", "data/databases/mouse_ncbi"),
    ("200-back/mouse_imgt_download", "data/databases/mouse_imgt/sources"),
    ("app_runs_integrado", "results/analysis"),
    ("000-input", "data/inputs"),
    ("runs", "results/manual"),
    ("lucianommartins", "core/humanization"),
    ("200-back", "core/backmutation"),
)


def resolve_path(value: str | Path) -> Path:
    """Remap only known project locations, never arbitrary external paths."""
    path = Path(value).expanduser()
    absolute = path if path.is_absolute() else PROJECT / path
    try:
        relative = absolute.relative_to(PROJECT).as_posix()
    except ValueError:
        return absolute.resolve()
    for previous, current in _RELOCATIONS:
        if relative == previous or relative.startswith(previous + "/"):
            suffix = relative[len(previous):].lstrip("/")
            return (PROJECT / current / suffix).resolve()
    return absolute.resolve()


def normalize_metadata(metadata: dict | None) -> dict | None:
    if metadata is None:
        return None
    result = dict(metadata)
    for field in ("input_file", "output_dir", "log_file", "database_path"):
        if result.get(field):
            result[field] = str(resolve_path(result[field]))
    return result


def require_virtual_environment() -> None:
    if Path(sys.prefix).resolve() != PYTHON.parent.parent.resolve():
        raise SystemExit("Activate the project virtual environment first: source .venv/bin/activate")
