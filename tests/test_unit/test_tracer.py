"""Unit tests for models of tracers."""

from cmfa.fluxomics_data.tracer import Tracer, TracerExperiment

EXAMPLE_TRACER_INPUT = {
    "id": "[1,2-13C]glucose",
    "compound_id": "A",
    "labelled_atom_positions": [2],
    "purity": 0.95,
}
EXAMPLE_TRACER_EXPERIMENT_INPUT = {
    "experiment_id": "e1",
    "tracer_id": "[1,2-13C]glucose",
    "enrichment": 1.0,
}


def test_tracer():
    Tracer.model_validate(EXAMPLE_TRACER_INPUT)


def test_tracer_experiment():
    TracerExperiment.model_validate(EXAMPLE_TRACER_EXPERIMENT_INPUT)
