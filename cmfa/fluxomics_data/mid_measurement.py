"""mid_measurement.py includes MID data from mass spectrometry for isotopomer fragments."""

import hashlib
from typing import List, Optional, Set

from pydantic import (
    BaseModel,
    Field,
    NonNegativeFloat,
    PositiveFloat,
    field_validator,
)


class MIDMeasurementComponent(BaseModel):
    """
    A class to represent a component of a Mass Isotopomer Distribution (MID) measurement.

    This class holds the details of a single mass isotopomer within an MID measurement,
    including its identifier, measured intensity, standard deviation, and normalized intensity.

    Parameters
    ----------
    mass_isotopomer_id : str
        Identifier of the mass isotopomer. It represents the mass shift due to isotopic labeling.
    measured_intensity : PositiveFloat
        The measured intensity (abundance) of the mass isotopomer.
    measured_std_dev : PositiveFloat
        The standard deviation of the measured intensity.
    normalized_intensity : PositiveFloat
        The normalized intensity of the mass isotopomer, calculated as a proportion of the total intensity of all isotopomers.

    Methods
    -------
    __repr__()
        Returns a string representation of each isotopomer measurement.

    """

    mass_isotopomer_id: str
    measured_intensity: PositiveFloat
    measured_std_dev: PositiveFloat
    normalized_intensity: NonNegativeFloat = 0

    def __repr__(self):
        """Return a string representation of each isotopomer measurement."""
        return (
            f"<MIDMeasurementComponent mass_isotopomer_id={self.mass_isotopomer_id}, "
            f"measured_intensity={self.measured_intensity}, "
            f"measured_std_dev={self.measured_std_dev}, "
            f"normalized_intensity={self.normalized_intensity}>"
        )

    # def __hash__(self):
    #     return hashlib.md5()

    # def __eq__(self, other):
    #     if not isinstance(other, MIDMeasurementComponent):
    #         return False
    #     return (
    #         self.mass_isotopomer_id == other.mass_isotopomer_id and
    #         self.measured_intensity == other.measured_intensity and
    #         self.measured_std_dev == other.measured_std_dev and
    #         self.normalized_intensity == other.normalized_intensity
    #     )


class MIDMeasurement(BaseModel):
    """
    A class to represent a Mass Isotopomer Distribution (MID) measurement.

    This class is used to document the details of an MID measurement,
    including the experiment ID, compound ID, fragment ID, and the
    individual mass isotopomer components.

    Parameters
    ----------
    experiment_id : str
        Identifier of the experiment under which the MID measurement was taken.
    compound_id : str
        Identifier of the compound for which the MID is measured.
    fragment_id : str
        Identifier of the fragment of the compound measured.
    measured_components : List[MIDMeasurementComponent]
        A list of MIDMeasurementComponent instances representing individual
        mass isotopomers and their measured properties.

    Methods
    -------
    __repr__()
        Returns a string representation of the MID measurement.

    normalize_components()
        Normalizes the intensities of the MIDMeasurementComponent instances
        to sum up to 1.
    """

    experiment_id: str
    compound_id: str
    fragment_id: str
    measured_components: List[MIDMeasurementComponent] = Field(
        default_factory=list
    )

    def __repr__(self):
        """Return a string representation of the MID measurement."""
        components_repr = ", ".join(
            [repr(comp) for comp in self.measured_components]
        )
        return (
            f"<MIDMeasurement experiment_id={self.experiment_id}, "
            f"compound_id={self.compound_id}, fragment_id={self.fragment_id}, "
            f"measured_components=[{components_repr}]>"
        )

    def normalize_components(self):
        """Ensure all the components are normalized."""
        total_intensity = sum(
            comp.measured_intensity for comp in self.measured_components
        )
        for comp in self.measured_components:
            comp.normalized_intensity = (
                comp.measured_intensity / total_intensity
                if total_intensity
                else 0
            )
