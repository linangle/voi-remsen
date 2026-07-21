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
