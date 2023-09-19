"""Unit tests for models of fluxomics measurements."""

from cmfa.fluxomics_data.measurement import FluxMeasurement, MIDMeasurement

EXAMPLE_FLUX_MEASUREMENT_INPUT = {
    "reaction_id": "v6",
    "experiment_id": "e1",
    "measured_flux": 1.2,
    "measurement_error": 0.1,
}
EXAMPLE_MID_MEASUREMENT_INPUT = {
    "experiment_id": "e1",
    "compound_id": "F",
    "fragment_id": "F1",
    "mass_isotopomer_id": 0,
    "measured_intensity": 0.2,
}


def test_flux_measurement():
    FluxMeasurement.model_validate(EXAMPLE_FLUX_MEASUREMENT_INPUT)


def test_mid_measurement():
    MIDMeasurement.model_validate(EXAMPLE_MID_MEASUREMENT_INPUT)
