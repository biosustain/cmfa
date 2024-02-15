"""reaction_network.py includes the description of compounds."""

from operator import gt, lt
from typing import List, Optional

from pydantic import (
    BaseModel,
    Field,
    PositiveFloat,
    computed_field,
    field_validator,
    model_validator,
)


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
        """Representation of the atom pattern."""
        return self.pattern_string

    def __len__(self) -> int:
        """Legnth of the atom pattern."""
        return len(self.pattern_string)
