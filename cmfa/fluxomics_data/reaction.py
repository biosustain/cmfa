"""reaction_network.py includes the description of reactions."""

from operator import gt, lt
from typing import Dict, List

from pydantic import BaseModel, Field, computed_field, model_validator

from cmfa.fluxomics_data.compound import Compound


class Reaction(BaseModel):
    """
    A class to represent a reaction in a reaction network.

    Attributes
    ----------
    id : str
        A unique identifier for the reaction.
    name : str
        The name of the reaction.
    compounds : Dict[Compounds, float]
        A dictionary with compound objects as keys and their stoichiometric coefficients as values.
    direction : str
        <-> means reversible reaction; <-- means backward reaction; --> means forward reaction
    atom_transition:


    Methods
    -------
    __repr__()
        Returns a string representation of the reaction.
    """

    id: str
    name: str
    compounds: Dict[Compound, float]
    direction: str = "<->"
    atom_transition: Dict[Compound, str]

    def __repr__(self):
        """<Return a string representation of the reaction."""
        return (
            f"Reaction id: {self.id}, name: {self.name},",
            f"number of compounds: {len(self.compounds.keys)},direction: {self.direction} >",
        )

    @model_validator(mode="after")
    def check_atom_balance(self) -> "Reaction":
        """Check if the atoms are balanced after the class is created."""
        atoms_in_str, atoms_out_str = (
            "".join(
                sorted(t.atom_pattern for t in self.transitions if f(t.coef, 0))
            )
            for f in (lt, gt)
        )
        msg = (
            f"Reaction {self.id} has unbalanced atom transitions."
            f"\n\tAtoms in: {atoms_in_str}"
            f"\n\tAtoms out: {atoms_out_str}\n"
        )
        assert atoms_in_str == atoms_out_str, msg
        return self
