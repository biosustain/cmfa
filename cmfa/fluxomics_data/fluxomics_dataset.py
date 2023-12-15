"""fluxomics_dataset.py includes the classes of a fluxomics dataset."""
from typing import List

from pydantic import BaseModel, field_validator

from cmfa.fluxomics_data.flux_measurement import FluxMeasurement
from cmfa.fluxomics_data.mid_measurement import MIDMeasurement
from cmfa.fluxomics_data.reaction_network import ReactionNetwork
from cmfa.fluxomics_data.tracer import Tracer, TracerExperiment


class FluxomicsDataset(BaseModel):
    """
    A class that encapsulates all data relevant to a fluxomics dataset.

    This class holds information about a specific fluxomics study, including
    the reaction network, tracers used, tracer experiments, and measurements
    of fluxes and metabolite isotopomer distributions (MIDs).

    Parameters
    ----------
    id : str
        A unique identifier for the fluxomics dataset.
    reaction_network : ReactionNetwork
        The reaction network model associated with the dataset. This includes
        details about reactions and compounds involved in the metabolic study.
    tracers : List[Tracer]
        A list of tracers used in the fluxomics experiments. Each tracer is
        represented with details like tracer ID, compound ID, labeled atom
        positions, and purity.
    tracer_experiments : List[TracerExperiment]
        A list of tracer experiments performed. Each experiment includes
        details about the experiment ID and the tracers used with their
        corresponding enrichments.
    flux_measurements : List[FluxMeasurement]
        A list of flux measurements. Each entry in the list represents a
        measured flux value in the metabolic network, typically obtained
        from experimental data.
    mid_measurements : List[MIDMeasurement]
        A list of Metabolite Isotopomer Distribution (MID) measurements.
        Each MID measurement includes details about the metabolite, sample
        source, and isotopomer distribution.

    Methods
    -------
    __repr__
        Returns a string representation of the fluxomics data.

    load_data()
        Given a set of data in csv or xlsx, read it into fluxomics data class.
    """

    reaction_network: ReactionNetwork
    tracers: List[Tracer]
    tracer_experiments: List[TracerExperiment]
    flux_measurements: List[FluxMeasurement]
    mid_measurements: List[MIDMeasurement]

    def __repr__(self):
        """Return a string representation of the fluxomics data."""
        return (
            f"<FluxomicsDataset, "
            f"num_reactions={len(self.reaction_network.reactions)}, "
            f"num_compounds={len(self.reaction_network.compounds)}, "
            f"num_tracers={len(self.tracers)}, "
            f"num_tracer_experiments={len(self.tracer_experiments)}, "
            f"num_flux_measurements={len(self.flux_measurements)}, "
            f"num_mid_measurements={len(self.mid_measurements)}>"
        )

    def __eq__(self, other):
        """Check equality with another FluxomicsDataset instance."""
        if not isinstance(other, FluxomicsDataset):
            return NotImplemented

        return (
            self.reaction_network == other.reaction_network
            and self.tracers == other.tracers
            and self.tracer_experiments == other.tracer_experiments
            and self.flux_measurements == other.flux_measurements
            and self.mid_measurements == other.mid_measurements
        )

    @field_validator("flux_measurements")
    def check_unique_replicates(cls, v) -> List[FluxMeasurement]:
        """Check flux measurement ids are unique for each experiment."""
        experiment_replicates = {}
        for measurement in v:
            if measurement.experiment_id in experiment_replicates:
                if (
                    measurement.replicate
                    in experiment_replicates[measurement.experiment_id]
                ):
                    raise ValueError(
                        f"Duplicate replicate {measurement.replicate} found for experiment_id {measurement.experiment_id}"
                    )
                experiment_replicates[measurement.experiment_id].add(
                    measurement.replicate
                )
            else:
                experiment_replicates[measurement.experiment_id] = {
                    measurement.replicate
                }

        return v
