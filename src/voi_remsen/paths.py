"""Path resolution for the voi-remsen world-model pipeline.

Single shared paths module (import everywhere). Assumes this file lives at
`src/voi_remsen/paths.py`, so `parents[2]` is the project root:

    <project_root>/src/voi_remsen/paths.py
    parents[0] = voi_remsen   parents[1] = src   parents[2] = <project_root>

Data-stage convention (cookiecutter-data-science):
    input   raw, immutable   (calHABMAP csv, CDPH xlsx)
    interim cleaned/derived  (joint_panel.csv, toxicity_panel.csv)
    output  final artifacts  (fitted Q, model-selection tables, figures)
"""

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]

INPUT_DIR = PROJECT_ROOT / "data" / "input"
INTERIM_DIR = PROJECT_ROOT / "data" / "interim"
OUTPUT_DIR = PROJECT_ROOT / "data" / "output"

for _d in (INPUT_DIR, INTERIM_DIR, OUTPUT_DIR):
    _d.mkdir(parents=True, exist_ok=True)


def input_file(name):
    return INPUT_DIR / name


def interim_file(name):
    return INTERIM_DIR / name


def output_file(name):
    return OUTPUT_DIR / name
