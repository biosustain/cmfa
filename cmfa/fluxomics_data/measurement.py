from pydantic import BaseModel, NonNegativeInt, PositiveFloat


class FluxMeasurement(BaseModel):
    reaction_id: str
    experiment_id: str
    measured_flux: float
    measurement_error: PositiveFloat


class MIDMeasurement(BaseModel):
    experiment_id: str
    compound_id: str
    fragment_id: str
    mass_isotopomer_id: NonNegativeInt
    measured_intensity: PositiveFloat
