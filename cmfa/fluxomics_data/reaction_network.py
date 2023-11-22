"""reaction_network.py includes the description of compound, reaction and reaction_network."""

# TODO: Set or List
from copy import deepcopy
from operator import gt, lt
from typing import List, Optional, Set

from pydantic import (
    BaseModel,
    Field,
    computed_field,
    field_validator,
    model_validator,
)

from cmfa.fluxomics_data.compound import Compound
from cmfa.fluxomics_data.reaction import Reaction


class ReactionNetwork(BaseModel):
    """
    A class representing a reaction network in a metabolic model.

    This class encapsulates all the necessary details of a metabolic reaction network,
    including its compounds and reactions. It allows for the addition of reactions
    and compounds to the network.

    Parameters
    ----------
    id : str
        A unique identifier for the model.
    name : Optional[str], default: None
        The name of the model, optional.
    reactions : Set[Reaction], default: empty set
        A set of reactions in the model. Ensures that each reaction is unique.
    compounds : Set[Compound], default: empty set
        A set of compounds in the model. Ensures that each compound is unique.

    Attributes
    ----------
    id : str
        Unique identifier of the reaction network.
    name : Optional[str]
        Name of the reaction network.
    reactions : Set[Reaction]
        Set of reactions in the network.
    compounds : Set[Compound]
        Set of compounds in the network.

    Methods
    -------
    add_reaction(reaction: Reaction)
        Adds a unique reaction to the network.
        If the reaction is already present, it is ignored.
    add_compound(compound: Compound)
        Adds a unique compound to the network.
        If the compound is already present, it is ignored.
    """

    id: str
    name: Optional[str] = None
    reactions: Set[Reaction] = set()  # List or Set?
    compounds: Set[Compound] = set()  # List or Set?

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

        return compound_set

    @field_validator("compounds", each_item=True, pre=True, always=True)
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
