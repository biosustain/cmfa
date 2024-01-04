"""Unit tests for model of flux measurements."""

from cmfa.fluxomics_data.flux_measurement import FluxMeasurement

EXAMPLE_FLUX_MEASUREMENT_INPUT = {
    "reaction_id": "v6",
    "experiment_id": "e1",
    "replicate": 1,
    "measured_flux": 1.2,
    "measurement_error": 0.1,
}


def test_flux_measurement():
    """Test good case of loading a flux measurement."""
    FluxMeasurement.model_validate(EXAMPLE_FLUX_MEASUREMENT_INPUT)
