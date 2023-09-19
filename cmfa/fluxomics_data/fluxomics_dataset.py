from typing import List

from pydantic import BaseModel

from cmfa.fluxomics_data.measurement import FluxMeasurement, MIDMeasurement
from cmfa.fluxomics_data.reaction_network import ReactionNetwork
from cmfa.fluxomics_data.tracer import Tracer, TracerExperiment


class FluxomicsDataset(BaseModel):
    id: str
    reaction_network: ReactionNetwork
    tracers: List[Tracer]
    tracer_experiments: List[TracerExperiment]
    flux_measurements: List[FluxMeasurement]
    mid_measurements: List[MIDMeasurement]
