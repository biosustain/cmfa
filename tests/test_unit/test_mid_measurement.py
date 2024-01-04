"""Unit tests for models of fluxomics measurements."""

from cmfa.fluxomics_data.mid_measurement import MIDMeasurement

EXAMPLE_MID_MEASUREMENT_INPUT = {
    "experiment_id": "e1",
    "compound_id": "F",
    "fragment_id": "F1",
    "mass_isotopomer_id": 0,
    "measured_intensity": 0.2,
}


def test_mid_measurement():
    """Test good case of loading a mid measurement."""
    MIDMeasurement.model_validate(EXAMPLE_MID_MEASUREMENT_INPUT)
