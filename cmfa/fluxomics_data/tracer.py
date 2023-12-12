"""tracer.py includes the description of tracers and experiments."""

from typing import Dict

from pydantic import BaseModel, Field, PositiveInt, field_validator


class Tracer(BaseModel):
    """
    Tracer Class records the trace compound information.

    Attributes
    ----------
    tracer_id: str
        Identifier of the tracer used. For instance, [2-13C]A indicates the metabolite A is labeled with 13C at the 2nd carbon position.

    compound_id: str
        tracer compound identifier, e.g. "Glu__L"

    labelled_atom_positions: set[PositiveInt]
        Positions of the carbon atoms in the tracer.

    purity: float
        The purity of the tracer, a float between 0 and 1, can be called as auto_mdv.
    """

    tracer_id: str
    compound_id: str
    labelled_atom_positions: set[PositiveInt]  # If empty is allowed
    purity: float = Field(default=1, gt=0, le=1)

    def __repr__(self):
        """Return a string representation of the Tracer instance."""
        atom_positions = ", ".join(map(str, self.labelled_atom_positions))
        return (
            f"<Tracer tracer_id={self.tracer_id}, compound_id={self.compound_id}, "
            f"labelled_atom_positions=[{atom_positions}], purity={self.purity}>"
        )

    def __hash__(self) -> int:
        """Return a unique hash of the compound."""
        return hash((self.tracer_id, self.compound_id, self.purity))


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
                f"{tracer_id}: {enrichment}"
                for tracer_id, enrichment in self.tracer_enrichments.items()
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
            isinstance(tracer_id, str) and isinstance(enrichment, float)
            for tracer_id, enrichment in v.items()
        ):
            raise ValueError(
                "Tracers must be a dictionary with tracer ID strings as keys and float as values"
            )
        # TODO: Check if tracer is found.
        return v

    def add_tracer(self, tracer_id: str, enrichment: float):
        """Add a new tracer to the experiment."""
        self.tracers[tracer_id] = enrichment

    def update_tracer_enrichment(self, tracer_id: str, enrichment: float):
        """Update tracer information."""
        if tracer_id in self.tracers:
            self.tracers[tracer_id] = enrichment
        else:
            raise ValueError(f"Tracer {tracer_id} not found in experiment")


# t = Tracer.model_validate(
#     {
#         "tracer_id": "[1,2-13C]glucose",
#         "compound_id": "A",
#         "labelled_atom_positions": [1, 2],
#         "purity": 0.95,
#     }
# )
# te = TracerExperiment.model_validate(
#     {
#         "experiment_id": "e1",
#         "tracer_enrichments": {"[1,2-13C]glucose": 1.0},
#     }
# )
# print(t)
# print(te)
