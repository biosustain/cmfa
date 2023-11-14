"""reaction_network.py includes the description of compound, reaction and reaction_network."""

from copy import deepcopy
from operator import gt, lt
from typing import List

from pydantic import (
    BaseModel,
    Field,
    computed_field,
    model_validator,
    validator,
)

from cmfa.fluxomics_data.compound import Compound
from cmfa.fluxomics_data.reaction import Reaction


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
    reactions: List[Reaction] = []
    compounds: List[Compound] = []

    def __repr__(self):
        """Return a string representation of the reaction network."""
        return (
            f"<ReactionNetwork id={self.id}, name={self.name}, "
            f"num_reactions={len(self.reactions)}, num_compounds={len(self.compounds)}>"
        )

    def generate_compounds(self) -> List[Compound]:
        """Create all the compounds for the model."""
        # Original:
        # out = []
        # for r in self.reactions:
        #     for t in r.transitions:
        #         new = Compound(
        #             id=t.compound_id, n_labellable_atoms=len(t.atom_pattern)
        #         )
        #         for cpd in out:
        #             if cpd.id == new.id:
        #                 assert cpd.n_labellable_atoms == new.n_labellable_atoms
        #         if not any(t.compound_id == c.id for c in out):
        #             out += [new]
        # return out

        compound_set = set()
        for r in self.reactions:
            for c in r.compounds:
                c_new = c.deepcopy()
                compound_set.add(c_new)

        return list(compound_set)

    @validator("compounds", each_item=True, pre=True, always=True)
    def check_all_compounds(self) -> "ReactionNetwork":
        """To check if the reaction network has all the compounds."""
        reaction_compounds = set(
            c for r in self.reactions for c in r.compounds()
        )
        model_compounds = set(self.compounds)
        if not reaction_compounds.issubset(model_compounds):
            missing = reaction_compounds - model_compounds
            raise ValueError(f"Missing compounds in the model: {missing}")

        return self

    def add_reaction(self, reaction: Reaction):
        """To add a reaction to the model. Duplicate reactions (based on id) are not added."""
        if reaction.id not in [r.id for r in self.reactions]:
            self.reactions.append(reaction)

    def add_compound(self, compound: Compound):
        """To add a compound to the model. Duplicate compounds (based on id) are not added."""
        if compound.id not in [c.id for c in self.compounds]:
            self.compounds.append(compound)

    @validator("reactions", each_item=True)
    def check_reaction_uniqueness(cls, reaction, values):
        """Ensure each reaction is unique within the network."""
        if (
            len([r for r in values.get("reactions", []) if r.id == reaction.id])
            > 1
        ):
            raise ValueError(f"Duplicate reaction ID found: {reaction.id}")
        return reaction

    @validator("compounds", each_item=True)
    def check_compound_uniqueness(cls, compound, values):
        """Ensure each compound is unique within the network."""
        if (
            len([c for c in values.get("compounds", []) if c.id == compound.id])
            > 1
        ):
            raise ValueError(f"Duplicate compound ID found: {compound.id}")
        return compound
