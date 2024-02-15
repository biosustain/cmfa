"""reaction_network.py includes the description of reactions."""

import hashlib
import warnings
from typing import Optional

from pydantic import BaseModel, field_validator, model_validator

from cmfa.fluxomics_data.compound import AtomPattern, Compound

type ReactionStoichiometry = dict[str, dict[AtomPattern, float]]


class Reaction(BaseModel):
    """
    A class to represent a reaction in a reaction network.

    Attributes
    ----------
    id : str
        A unique identifier for the reaction.
    name : str
        The name of the reaction.
    stoichiometry_input : dict[str, float]
        A dictionary with compound names as keys, whose values are maps of atom
        patterns to stoichiometric coefficients.

        For example, {'A': {"abc": -0.5, "bca": -0.5}, "B": {"acb": 1}}
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
    stoichiometry_input: dict[str, dict[str, float]]

    @property
    def stoichiometry(self) -> ReactionStoichiometry:
        """Get the stoichiometry in the right form."""
        return {
            compound: {
                AtomPattern(pattern=pattern): coef
                for pattern, coef in compound_stoich.items()
            }
            for compound, compound_stoich in self.stoichiometry_input.items()
        }

    def __repr__(self):
        """Return a string representation of the reaction."""
        return (
            f"Reaction id={self.id}, name={self.name},"
            f"stoichiometry={self.stoichiometry},"
            f"reversible={self.reversible}, "
        )

    def __hash__(self) -> int:
        """Return a unique hash of the reaction."""
        return hash(self.id)

    @field_validator("stoichiometry_input")
    def check_stoichiometry(
        cls, v: dict[str, dict[str, float]]
    ) -> dict[str, dict[str, float]]:
        """Check that the stoichiometry input is non-empty."""
        assert len(v.keys()) > 1, "Stoichiometry must be non-empty."
        return v

    @model_validator(mode="after")
    def check_atom_balance(self):
        """Check if the atoms are balanced in the reaction."""
        lhs_atoms, rhs_atoms = 0, 0
        for compound, transitions in self.stoichiometry.items():
            for transition, coeff in transitions.items():
                if coeff < 0:  # Reactant
                    lhs_atoms += sum(transition.pattern_tuple) * abs(coeff)
                else:  # Product
                    rhs_atoms += sum(transition.pattern_tuple) * abs(coeff)
        if lhs_atoms != rhs_atoms:
            raise ValueError(
                f"Unbalanced atoms in reaction {self.id}: {lhs_atoms} != {rhs_atoms}"
            )
        return self
