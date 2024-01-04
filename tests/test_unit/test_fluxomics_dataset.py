"""Unit tests for the FluxomicsDataset data model."""

from cmfa.fluxomics_data.fluxomics_dataset import FluxomicsDataset

from .test_flux_measurement import EXAMPLE_FLUX_MEASUREMENT_INPUT
from .test_mid_measurement import EXAMPLE_MID_MEASUREMENT_INPUT
from .test_reaction_network import EXAMPLE_NETWORK_INPUT
from .test_tracer import EXAMPLE_TRACER_EXPERIMENT_INPUT, EXAMPLE_TRACER_INPUT

EXAMPLE_FLUXOMICS_DATASET_INPUT = {
    "id": "My fluxomics dataset",
    "method": "labelling",
    "reaction_network": EXAMPLE_NETWORK_INPUT,
    "tracers": [EXAMPLE_TRACER_INPUT],
    "tracer_experiments": [EXAMPLE_TRACER_EXPERIMENT_INPUT],
    "flux_measurements": [EXAMPLE_FLUX_MEASUREMENT_INPUT],
    "mid_measurements": [EXAMPLE_MID_MEASUREMENT_INPUT],
}


def test_fluxomics_dataset():
    """Test good case of loading a fluxomics dataset."""
    FluxomicsDataset.model_validate(EXAMPLE_FLUXOMICS_DATASET_INPUT)
