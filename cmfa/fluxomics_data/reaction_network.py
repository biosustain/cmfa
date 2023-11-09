"""reaction_network.py includes the description of compound, reaction and reaction_network."""

from operator import gt, lt
from typing import List

from pydantic import BaseModel, Field, computed_field, model_validator


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
    carbon_label: str
        The labelled carbon in the compound.

    Methods
    -------
    __repr__()
        Returns a string representation of the compound.
    """

    id: str
    name: str
    formula: Optional[str] = None
    carbon_label: str
    n_labellable_atoms: int = Field(ge=1)

    def __repr__(self):
        """Return a self-description of the compound."""
        return f"Compound id: {self.id}, name: {self.name}, carbon_label: {self.atom_numbers}"


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
        <-> means reversible reaction; <- means backward reaction; -> means forward reaction
    atom_transition:


    Methods
    -------
    __repr__()
        Returns a string representation of the reaction.
    """

    id: str
    name: str
    compounds: Dict[Compound, float]
    direction: str
    atom_transition: Dict[Compound, str]

    def __repr__(self):
        """Return a string representation of the reaction."""
        return f"Reaction id: {self.id}, name: {self.name}: "

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


class ReactionNetwork(BaseModel):
    """
    The reaction network incorporates the compounds and reactions into a small model.

    -----------

    Parameters
    ----------
    id: str
        A unique identifier for the model.
    name: str
        The name of the model.
    reactions: List[Reaction]
        A list of reactions in the model.
    compounds: List[Compound]
        A list of compounds in the model.

    Methods
    -------
    add_reaction(reaction: Reaction)
        Adds a reaction to the model.
    add_compound(compound: Compound)
        Adds a compound to the model.
    """

    id: str
    name: str
    reactions: List[Reaction]
    compounds: List[Compound]

    @computed_field
    def compounds(self) -> List[Compound]:
        """Create all the compounds for the model."""
        out = []
        for r in self.reactions:
            for t in r.transitions:
                new = Compound(
                    id=t.compound_id, n_labellable_atoms=len(t.atom_pattern)
                )
                for cpd in out:
                    if cpd.id == new.id:
                        assert cpd.n_labellable_atoms == new.n_labellable_atoms
                if not any(t.compound_id == c.id for c in out):
                    out += [new]
        return out

    @model_validator(mode="after")
    def check_compounds(self) -> "ReactionNetwork":
        """To check if the reaction network has all the compounds."""
        return self
