"""reaction_network.py includes the description of reactions."""

import hashlib
import warnings
from typing import Dict, Optional

from pydantic import (
    BaseModel,
    Field,
    computed_field,
    field_validator,
    model_validator,
)

from cmfa.fluxomics_data.compound import Compound

type ReactionStoichiometry = Dict[str, Dict[AtomPattern, float]]


class AtomPattern(BaseModel):
    """A string representing the order of labellable atoms in a compound.

    For example "abcd".

    Atom patterns are understood in the context of a reaction. For example if
    the reaction turns compound A with atom pattern "ab" into compound B with
    atom pattern "ba", that means that the reaction swaps the order of the two
    labellable atoms in these compounds.

    Note that atom patterns rely on alphabetical ordering, so only alphabetic
    characters are allowed, and that each atom is represented by a single unique
    character. For simplicity upper case letters are also not allowed.

    Atom patterns can also be represented as tuples of integers.

    """

    pattern_string: str = Field(alias="pattern")

    @field_validator("pattern_string")
    def check_characters(cls, v: str) -> str:
        """Check that the atom pattern is valid."""
        for char in v:
            assert char.isalpha(), f"Found non-alphabetic character {char}."
            assert not char.isupper(), f"Found upper case character {char}."
            assert v.count(char) == 1, f"Found duplicate character {char}."
        return v

    @computed_field
    def pattern_tuple(self) -> tuple[int, ...]:
        """Convert atom pattern to integer representation.

        e.g. "abdc" will become (1,2,4,3).
        """
        return tuple(ord(l) - ord("a") + 1 for l in self.pattern_string)

    def __hash__(self) -> int:
        """Hash an atom pattern."""
        return hash(self.pattern_string)

    def __repr__(self) -> str:
        """Hash an atom pattern."""
        return self.pattern_string


class Reaction(BaseModel):
    """
    A class to represent a reaction in a reaction network.

    Attributes
    ----------
    id : str
        A unique identifier for the reaction.
    name : str
        The name of the reaction.
    stoichiometry_input : Dict[str, float]
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
    stoichiometry_input: Dict[str, Dict[str, float]]

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
            f"direction={self.reversible}, "
        )

    def __hash__(self) -> int:
        """Return a unique hash of the reaction."""
        return hash(self.id)

    @field_validator("stoichiometry_input")
    def check_stoichiometry(
        cls, v: Dict[str, Dict[str, float]]
    ) -> Dict[str, Dict[str, float]]:
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
                    print(sum(transition.pattern_tuple) * abs(coeff))
                    lhs_atoms += sum(transition.pattern_tuple) * abs(coeff)
                else:  # Product
                    rhs_atoms += sum(transition.pattern_tuple) * abs(coeff)
        if lhs_atoms != rhs_atoms:
            raise ValueError(
                f"Unbalanced atoms in reaction {self.id}: {lhs_atoms} != {rhs_atoms}"
            )
        return self
