"""reaction_network.py includes the description of compounds."""

from operator import gt, lt
from typing import List, Optional

from pydantic import BaseModel, PositiveFloat


class Compound(BaseModel):
    """
    A class to represent a compound in a reaction network.

    Attributes
    ----------
    id : str
        A unique identifier for the compound.
    name : str
        The name of the compound.
    formula : Optional[str]
        The chemical formula of the compound.
    carbon_label: Optional[str]
        The labelled carbon in the compound.

    Methods
    -------
    __repr__()
        Returns a string representation of the compound.
    """

    id: str
    name: Optional[str] = None
    formula: Optional[str] = None
    carbon_label: Optional[str] = None
    # For calculating the atom transition, and largest EMU carbon
    # n_labellable_atoms: int = Field(ge=0)
    fragment_id: Optional[str] = None
    m_z: Optional[PositiveFloat] = None

    def __repr__(self):
        """Return a self-description of the compound."""
        return f"Compound id: {self.id}, name: {self.name}, carbon_label: {self.carbon_label}"

    def __hash__(self):
        """Return a unique hash of the compound."""
        return hash((self.id, self.carbon_label))
