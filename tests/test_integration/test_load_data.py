"""Test loading datasets."""

from pathlib import Path

from cmfa.data_preparation import (
    import_fluxomics_dataset_from_json,
    load_dataset_from_csv,
)

HERE = Path(__file__).parent
DATA_DIR = HERE / ".." / ".." / "data" / "test_data"
FLUX_MEASUREMENT_FILE = DATA_DIR / "flux.csv"
MID_MEASUREMENT_FILE = DATA_DIR / "ms_measurements.csv"
REACTION_FILE = DATA_DIR / "reactions.csv"
TRACER_FILE = DATA_DIR / "tracers.csv"
MODEL_FILE = DATA_DIR / "model.json"


def test_load_dataset_from_csv():
    """Test loading a dataset from csv."""
    ds_from_csv = load_dataset_from_csv(
        tracer_file=TRACER_FILE,
        flux_measurement_file=FLUX_MEASUREMENT_FILE,
        mid_measurement_file=MID_MEASUREMENT_FILE,
        reaction_file=REACTION_FILE,
    )
    ds_from_json = import_fluxomics_dataset_from_json(MODEL_FILE)
    assert ds_from_csv == ds_from_json


test_load_dataset_from_csv()
