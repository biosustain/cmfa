"""Integration tests for functions in src/data_preparation.py."""

from typing import Callable

import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from cmfa.data_preparation import prepare_data_main
from cmfa.util import CoordDict

EXAMPLE_RAW_MEASUREMENTS = pd.DataFrame(
    {
        "X1": [1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2],
        "X2": [1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2],
    }
)


@pytest.mark.parametrize(
    "prepare_data_function,name,raw_measurements,expected_measurements,"
    "expected_coords",
    [
        (
            prepare_data_main,
            "main",
            EXAMPLE_RAW_MEASUREMENTS,
            pd.DataFrame(
                {
                    "X1": [1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2],
                    "X2": [1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2],
                }
            ),
            {"": [""]},
        ),
    ],
)
def test_prepare_data_function(
    prepare_data_function: Callable,
    name: str,
    raw_measurements: pd.DataFrame,
    expected_measurements: pd.DataFrame,
    expected_coords: CoordDict,
):
    """Check that a prepare data function behaves as expected."""
    prepped = prepare_data_function(raw_measurements)
    assert prepped.name == name
    assert prepped.coords == expected_coords
    assert_frame_equal(prepped.measurements, expected_measurements)
