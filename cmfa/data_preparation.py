"""Provides functions prepare_data_x.

These functions should take in a dataframe of measurements and return a
PreparedData object.
"""
import json
import os

import pandas as pd
import pandera as pa
from pandera.typing import DataFrame
from pydantic import BaseModel

from cmfa import util

NAME_FILE = "name.txt"
COORDS_FILE = "coords.json"
MEASUREMENTS_FILE = "meausurements.csv"
DIMS: dict = {}
HERE = os.path.dirname(__file__)
DATA_DIR = os.path.join(HERE, "..", "data")
RAW_DIR = os.path.join(DATA_DIR, "raw")
PREPARED_DIR = os.path.join(DATA_DIR, "prepared")
RAW_DATA_FILES: dict = {}


class MeasurementsDF(pa.SchemaModel):
    """A PreparedData should have a measurements dataframe like this."""


class PreparedData(BaseModel, arbitrary_types_allowed=True):
    """What prepared data looks like in this analysis."""

    name: str
    coords: util.CoordDict
    measurements: DataFrame[MeasurementsDF]


def load_prepared_data(directory: str) -> PreparedData:
    """Load prepared data from files in directory."""
    with open(os.path.join(directory, COORDS_FILE), "r") as f:
        coords = json.load(f)
    with open(os.path.join(directory, NAME_FILE), "r") as f:
        name = f.read()
    measurements = pd.read_csv(os.path.join(directory, MEASUREMENTS_FILE))
    return PreparedData(
        name=name,
        coords=coords,
        measurements=DataFrame[MeasurementsDF](measurements),
    )


def write_prepared_data(prepped: PreparedData, directory):
    """Write prepared data files to a directory."""
    if not os.path.exists(directory):
        os.mkdir(directory)
        prepped.measurements.to_csv(os.path.join(directory, MEASUREMENTS_FILE))
    with open(os.path.join(directory, COORDS_FILE), "w") as f:
        json.dump(prepped.coords, f)
    with open(os.path.join(directory, NAME_FILE), "w") as f:
        f.write(prepped.name)


def prepare_data_main(
    measuremements: DataFrame[MeasurementsDF],
) -> PreparedData:
    return PreparedData(
        name="",
        coords=util.CoordDict({"": [""]}),
        measurements=DataFrame[MeasurementsDF](),
    )


def prepare_data():
    """Run main function."""
    print("Reading raw data...")
    raw_data = {
        k: pd.read_csv(v, index_col=None) for k, v in RAW_DATA_FILES.items()
    }
    data_preparation_functions_to_run = [prepare_data_main]
    print("Preparing data...")
    for dpf in data_preparation_functions_to_run:
        print(f"Running data preparation function {dpf.__name__}...")
        prepared_data = dpf(raw_data["raw_measurements"])
        output_dir = os.path.join(PREPARED_DIR, prepared_data.name)
        print(f"\twriting files to {output_dir}")
        if not os.path.exists(PREPARED_DIR):
            os.mkdir(PREPARED_DIR)
        write_prepared_data(prepared_data, output_dir)
