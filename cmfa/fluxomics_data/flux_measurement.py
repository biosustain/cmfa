"""flux_measurement.py includes the tracer flux data coming out of experiments."""

from typing import Optional

from pydantic import BaseModel, PositiveFloat


class FluxMeasurement(BaseModel):
    """
    A class to represent a single flux measurement in a metabolic network.

    This class holds the details of a flux measurement including the reaction
    ID, experiment ID, replicate information, measured flux value, and the
    associated measurement error.

    Parameters
    ----------
    experiment_id : str
        Identifier of the experiment under which the flux measurement was taken.
    reaction_id : str
        Identifier of the reaction for which the flux is measured.
    replicate : int
        An integer representing the replicate number of the experiment.
        Unique for each measurement.
    measured_flux : float
        The measured flux value for the reaction in the specified experimental
        condition and replicate.
    measurement_error : PositiveFloat
        The error or uncertainty associated with the measured flux.

    """

    experiment_id: str
    reaction_id: str
    replicate: int  # Unique int pydantic
    measured_flux: Optional[float] = None
    measurement_error: Optional[PositiveFloat] = None

    def __repr__(self):
        """Return a string representation of the FluxMeasurement instance."""
        return (
            f"<FluxMeasurement experiment_id={self.experiment_id}, "
            f"reaction_id={self.reaction_id}, replicate={self.replicate}, "
            f"measured_flux={self.measured_flux}, "
            f"measurement_error={self.measurement_error}>"
        )
