"""emu_map.py includes a class for maps that generated from EMU algorithm."""

from pydantic import BaseModel, Field, PositiveInt, computed_field

from cmfa.fluxomics_data.compound import AtomPattern


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
    atom_pattern: AtomPattern
        The atom numbers, e.g. 'acb' is recorded as '(1, 3, 2)'

    Methods
    -------
    __repr__()
        Return a string representation of the EMU reaction.

    """

    id: str
    compound: str
    atom_pattern: AtomPattern = Field(alias="pattern")

    def __repr__(self):
        """Return a self-description of the EMU."""
        return f"EMU id: {self.id}, compound id: {self.compound_id}, atom pattern: {self.atom_pattern}"

    def __hash__(self):
        """Return a unique hash of the EMU."""
        return hash((self.compound_id, self.atom_pattern))

    @computed_field
    @property
    def size(self) -> PositiveInt:
        """Return the size of EMU."""
        return len(self.atom_pattern)
