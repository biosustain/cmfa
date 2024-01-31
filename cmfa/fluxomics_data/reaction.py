"""reaction_network.py includes the description of reactions."""

import hashlib
import warnings
from typing import Dict, Optional

from pydantic import BaseModel, computed_field, field_validator, model_validator

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
    stoichiometry : Dict[str, float]
        A dictionary with compound objects as keys and their stoichiometric coefficients as values.
        For example, {'a': -1, 'b': 1}
    reversible : bool
        Whether or not the reaction is reversible

    Methods
    -------
    __repr__()
        Returns a string representation of the reaction.

    __hash__()
        Returns an unique id for a reaction.

    check_atom_balance()
        Check if an reaction is mass balanced.

    """

    id: str
    name: Optional[str] = None
    reversible: bool = True
    stoichiometry: Dict[str, Dict[tuple[int, ...], float]]

    def __repr__(self):
        """Return a string representation of the reaction."""
        return (
            f"Reaction id={self.id}, name={self.name},"
            f"stoichiometry={self.stoichiometry},"
            f"direction={self.reversible}, "
        )

    def __hash__(self) -> int:
        """Return a unique hash of the reaction."""
        return hash(self.id)

    @field_validator("stoichiometry", mode="before")
    def convert_atom_transition_in_stoichiometry(cls, v):
        """Convert string form of atom transition into integer representation. e.g. "abdc" will become (1,2,4,3)."""
        converted_stoichiometry = {}

        for compound, transitions in v.items():
            converted_transitions = {}
            for transition_str, value in transitions.items():
                if transition_str is None or transition_str == "":
                    warnings.warn(
                        f"No atom transition provided for compound {compound}."
                    )
                    continue

                if any(
                    char.isdigit() or char.isupper() for char in transition_str
                ):
                    raise ValueError(
                        f"Invalid atom transition '{transition_str}' for compound {compound}. Only lowercase alphabets are allowed."
                    )

                # Convert transition string to list of integers
                transition_tuple = tuple(
                    [ord(char) - ord("a") + 1 for char in transition_str]
                )

                converted_transitions[transition_tuple] = value
            converted_stoichiometry[compound] = converted_transitions

        return converted_stoichiometry

    @model_validator(mode="after")
    def check_atom_balance(self):
        """Check if the atoms are balanced in the reaction based on atom transitions."""
        if self.stoichiometry == dict():
            raise ValueError("Stoichiometry is invalid")

        lhs_atoms, rhs_atoms = 0, 0
        for compound, transitions in self.stoichiometry.items():
            for transition, coeff in transitions.items():
                if coeff < 0:  # Reactant
                    print(sum(transition) * abs(coeff))
                    lhs_atoms += sum(transition) * abs(coeff)
                else:  # Product
                    rhs_atoms += sum(transition) * abs(coeff)

        if lhs_atoms != rhs_atoms:
            raise ValueError(
                f"Unbalanced atoms in reaction {self.id}: {lhs_atoms} != {rhs_atoms}"
            )

        return self
