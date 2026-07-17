#  path resolution for all data-acquisition scripts

from pathlib import Path


_THIS = Path(__file__).resolve()
PROJECT_ROOT = _THIS.parent.parent

INPUT_DIR = PROJECT_ROOT / "data_input"
OUTPUT_DIR = PROJECT_ROOT / "data_output" / "from_data_acquisition"

# subfolders of the output tree
SATELLITE_RAW_DIR = OUTPUT_DIR / "satellite_data"       # downloaded granules
SATELLITE_MATCHUP_DIR = OUTPUT_DIR / "satellite_matchup"  # extracted weeklies
AUDIT_DIR = SATELLITE_MATCHUP_DIR / "audit"


def ensure_dirs():
    """Create the output tree if missing. Inputs are never auto-created."""
    for d in (OUTPUT_DIR, SATELLITE_RAW_DIR, SATELLITE_MATCHUP_DIR, AUDIT_DIR):
        d.mkdir(parents=True, exist_ok=True)


# convenience helpers for the named artifacts
def input_file(name):
    return INPUT_DIR / name


def output_file(name):
    return OUTPUT_DIR / name
