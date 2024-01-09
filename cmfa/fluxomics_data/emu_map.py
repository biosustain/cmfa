"""emu_map.py includes a class for maps that generated from EMU algorithm."""

from typing import Dict, List, Optional

from pydantic import BaseModel, PositiveInt, field_validator, model_validator

from cmfa.fluxomics_data.compound import Compound
from cmfa.fluxomics_data.reaction import Reaction


class EMUReaction(BaseModel):
    """
    A class to represent a single EMU reaction.

    Attributes
    ----------
    reaction : Reaction
        Identifier of the reaction.
    emu_size: PositiveInt
        The size of EMU map.
    stoichiometry : Dict[Compound, PositiveInt]
        A dictionary that representing the reaction stoichiometry.
    atom_transition :  Dict[Compound, List]
        A dictionary representing the atom transition in the reaction.

    Methods
    -------
    __repr__()
        Return a string representation of the EMU reaction.

    check_emu_size_balance()
        Check if the emu reaction is balanced.

    """

    emu_size: PositiveInt
    reaction: Reaction
    stoichiometry: Dict[Compound, PositiveInt]
    atom_transition: Dict[Compound, List]

    def __repr__(self):
        """Return a string representation of the EMU reaction."""
        stoich_repr = ", ".join(
            f"{k}: {v}" for k, v in self.stoichiometry.items()
        )
        atom_trans_repr = ", ".join(
            f"{k}: {v}" for k, v in self.atom_transition.items()
        )
        return (
            f'EMUReaction(emu_size={self.emu_size}, reaction="{self.reaction}", '
            f"stoichiometry={{{stoich_repr}}}, atom_transition={{{atom_trans_repr}}})"
        )

    @model_validator(mode="after")
    def check_emu_size_balance(self):
        """Check if the emu reaction is balanced."""
        if self.stoichiometry.keys() != self.atom_transition.keys():
            raise ValueError(
                "EMU stoichiometry and atom_transition must have the same keys"
            )

        total_sum = sum(
            len(self.stoichiometry[key]) * self.atom_transition[key]
            for key in self.stoichiometry
        )
        if total_sum != 0:
            raise ValueError(
                "The EMU size balance is not 0, check your EMU map again."
            )

        return self


class EMUMap(BaseModel):
    """
    A class representing a map of EMU reactions.

    Attributes
    ----------
    emu_map: List[EMUReaction]
        A list of EMUReaction objects that make up the EMU map.

    Methods
    -------
    __repr__()
        Return a string representation of the EMU map.

    calculate_mid()
        calculate the mid for each EMU.

    """

    emu_reactions: List[EMUReaction]

    def __repr__(self):
        """Return a string representation of the EMU map."""
        reactions_repr = ", ".join(
            repr(reaction) for reaction in self.emu_reactions
        )
        return f"EMUMap(emu_reactions=[{reactions_repr}])"

    def calculate_mid(self):
        """Calculate the mid for each EMU."""
        return self
