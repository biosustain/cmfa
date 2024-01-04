"""reaction_network.py includes the description of reactions."""

import hashlib
from enum import Enum
from operator import gt, lt
from typing import Dict, List, Optional

from pydantic import BaseModel, Field, computed_field, model_validator

from cmfa.fluxomics_data.compound import Compound


class ReactionDirection(Enum):
    """A class that defines reversible and forward reactions."""

    REVERSIBLE = 0  # "<->"
    FORWARD = 1  # "-->"


class Reaction(BaseModel):
    """
    A class to represent a reaction in a reaction network.

    Attributes
    ----------
    id : str
        A unique identifier for the reaction.
    name : str
        The name of the reaction.
    stoichiometry : Dict[str, float]
        A dictionary with compound objects as keys and their stoichiometric coefficients as values.
    reversible : bool
        Whether or not the reaction is reversible
    atom_transition:
        Dict{str (compound id), str (atom pattern)}

    Methods
    -------
    __repr__()
        Returns a string representation of the reaction.

    check_atom_balance()
        Check if an reaction is mass balanced.

    """

    id: str
    name: Optional[str] = None
    stoichiometry: Dict[str, float]
    reversible: bool = True
    atom_transition: Dict[str, list] = dict()

    def __repr__(self):
        """Return a string representation of the reaction."""
        return (
            f"Reaction id={self.id}, name={self.name},"
            f"stoichiometry={self.stoichiometry},"
            f"direction={self.reversible}, "
            f"atom_transtions={self.atom_transition}>"
        )

    def __hash__(self) -> int:
        """Return a unique hash of the reaction."""
        return hash(self.id)

    @model_validator(mode="after")
    def check_atom_balance(self):
        """
        Check if the atoms are balanced in the reaction based on atom transitions.

        Parameters
        ----------
        atom_transition : Dict[str, str]
            A dictionary mapping each compound to its atom transition pattern.

        Returns
        -------
        Dict[str, str]
            The validated atom transition dictionary.

        Raises
        ------
        ValueError
            If the atoms are not balanced.
        """
        if self.stoichiometry == dict():
            return self.atom_transition  # Cannot validate without compounds

        lhs_atoms, rhs_atoms = "", ""
        for compound, coeff in self.stoichiometry.items():
            transition = self.atom_transition.get(compound, "")
            if coeff < 0:  # Reactant
                lhs_atoms += "".join(transition)
            else:  # Product
                rhs_atoms += "".join(transition)

        if sorted(lhs_atoms) != sorted(rhs_atoms):
            raise ValueError(
                f"Unbalanced atoms in reaction {self.id}: {lhs_atoms} != {rhs_atoms}"
            )

        return self
