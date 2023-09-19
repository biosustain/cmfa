from pydantic import BaseModel, Field, PositiveInt


class Tracer(BaseModel):
    id: str
    compound_id: str
    labelled_atom_positions: list[PositiveInt]
    purity: float = Field(gt=0, le=1)


class TracerExperiment(BaseModel):
    experiment_id: str
    tracer_id: str
    enrichment: float = Field(gt=0, le=1)


t = Tracer.model_validate(
    {
        "id": "[1,2-13C]glucose",
        "compound_id": "A",
        "labelled_atom_positions": [1, 2],
        "purity": 0.95,
    }
)
te = TracerExperiment.model_validate(
    {"experiment_id": "e1", "tracer_id": "[1,2-13C]glucose", "enrichment": 1.0}
)
