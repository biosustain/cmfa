"""tracer.py includes the description of tracers and experiments."""

from typing import Dict

from pydantic import BaseModel, Field, PositiveInt, field_validator


class Tracer(BaseModel):
    """
    Tracer Class records the trace compound information.

    Attributes
    ----------
    isotope: str
        Identifier of the tracer used. For instance, [2-13C]A indicates the metabolite A is labeled with 13C at the 2nd carbon position.

    compound: str
        tracer compound identifier, e.g. "Glu__L"

    labelled_atom_positions: set[PositiveInt]
        Positions of the carbon atoms in the tracer.

    purity: float
        The purity of the tracer, a float between 0 and 1, can be called as auto_mdv.
    """

    isotope: str
    compound: str
    labelled_atom_positions: set[PositiveInt]  # If empty is allowed
    purity: float = Field(default=1, gt=0, le=1)

    def __repr__(self):
        """Return a string representation of the Tracer instance."""
        atom_positions = ", ".join(map(str, self.labelled_atom_positions))
        return (
            f"<Tracer isotope={self.isotope}, compound={self.compound}, "
            f"labelled_atom_positions=[{atom_positions}], purity={self.purity}>"
        )

    def __hash__(self) -> int:
        """Return a unique hash of the compound."""
        return hash((self.isotope, self.compound, self.purity))


class TracerExperiment(BaseModel):
    """
    The tracer experiment that takes are dictionary of all tracers and each of the enrichment.

    Attributes
    ----------
    experiment_id:
        Identifier for the experiment.

    tracer_enrichments: Dict[str, float]
        A dictionary of tracers id and enrichment value, e.g. 90% of [1,2]A and 10% of [2]B is in the medium.

    """

    experiment_id: str
    tracer_enrichments: Dict[str, float] = dict()

    def __repr__(self):
        """Return a string representation of the TracerExperiment instance."""
        enrichments_repr = ", ".join(
            [
                f"{isotope}: {enrichment}"
                for isotope, enrichment in self.tracer_enrichments.items()
            ]
        )
        return (
            f"<TracerExperiment experiment_id={self.experiment_id}, "
            f"tracer_enrichments={{ {enrichments_repr} }}>"
        )

    @field_validator("tracer_enrichments")
    def validate_tracers(cls, v):
        """Validate the tracers in the experiment."""
        if not all(
            isinstance(isotope, str) and isinstance(enrichment, float)
            for isotope, enrichment in v.items()
        ):
            raise ValueError(
                "Tracers must be a dictionary with tracer ID strings as keys and float as values"
            )
        # TODO: Check if tracer is found.
        return v

    def add_tracer(self, isotope: str, enrichment: float):
        """Add a new tracer to the experiment."""
        self.tracers[isotope] = enrichment

    def update_tracer_enrichment(self, isotope: str, enrichment: float):
        """Update tracer information."""
        if isotope in self.tracers:
            self.tracers[isotope] = enrichment
        else:
            raise ValueError(f"Tracer {isotope} not found in experiment")
