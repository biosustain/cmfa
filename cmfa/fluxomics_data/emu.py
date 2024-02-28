"""emu_map.py includes a class for maps that generated from EMU algorithm."""

from pydantic import (
    BaseModel,
    Field,
    PositiveInt,
    computed_field,
    model_validator,
)


class EMU(BaseModel):
    """
    A class to represent a EMU. Define EMU (i.e., elementary metabolite unit) object and its operations.

    EMUs in the same metabolite and with the same atom NOs are considered as identical,
    while metabolites which they derived from could be different.

    EMUs can be compared based self.compound_id and self.atom_nos.
    EMU and iterable object of EMUs can also be compared. In this case EMU will be put into
    the same iterable object with single item, and comparison between two iterables are performed.

    Attributes
    ----------
    id: str
        A unique EMU ID.
    compound_id: str
        The associated compound id.
    atom_number: str
        The atom numbers of a compound involved in a reaction, e.g. if first three atoms in A('acbd') is used in EMU map, it will be recorded as '123'.

    Methods
    -------
    __repr__()
        Return a string representation of the EMU reaction.

    """

    id: str
    compound: str
    atom_number_input: str

    def __repr__(self):
        """Return a self-description of the EMU."""
        return f"EMU id: {self.id}, compound id: {self.compound}, atom pattern: {self.atom_number}"

    def __hash__(self):
        """Return a unique hash of the EMU."""
        return hash((self.compound, self.atom_number))

    @computed_field
    @property
    def atom_number(self) -> tuple[int, ...]:
        """Convert atom pattern to integer representation.

        e.g. "234" will become (2,3,4).
        """
        if self.atom_number_input.isdigit():
            return tuple(map(int, self.atom_number_input))
        else:
            raise ValueError("The EMU string contains non-numeric characters.")

    @computed_field
    @property
    def size(self) -> PositiveInt:
        """Return the size of EMU."""
        return len(self.atom_number_input)


class EMUReaction(BaseModel):
    """
    EMU Reaction defined as a subset of atoms that participated in a given reaction of the model.

    Attributes
    ----------
    reaction_id: str
        The reaction id of the corresponding EMU reaction
    emu_stoichiometry: dict[str, dict[str, float]]
        The stoichiometry that describes the atoms that participated in the reaction.
    """

    reaction_id: str
    emu_stoichiometry: dict[str, dict[str, float]]

    def __repr__(self):
        """Return a self-description of the EMU reaction."""
        return f"reaction id: {self.reaction_id}, EMU stoichiometry: {self.emu_stoichiometry}"

    def __hash__(self):
        """Return a unique hash of the EMU."""
        return hash((self.reaction_id, self.emu_stoichiometry))

    @model_validator(mode="after")
    def check_emu_balance(self):
        """Check if the atoms are balanced in the EMU reaction."""
        lhs_atoms, rhs_atoms = 0, 0
        for compound, transitions in self.emu_stoichiometry.items():
            for atom, coeff in transitions.items():
                if coeff < 0:  # Reactant
                    lhs_atoms += len(atom) * abs(coeff)
                else:  # Product
                    rhs_atoms += len(atom) * abs(coeff)
        if lhs_atoms != rhs_atoms:
            raise ValueError(
                f"Unbalanced atoms in the EMU reaction {self.id}: {lhs_atoms} != {rhs_atoms}"
            )
        return self
