"""reaction_network.py includes the description of compound, reaction and reaction_network."""

import warnings
from copy import deepcopy
from operator import gt, lt
from typing import List, Optional, Set

import pandas as pd
from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    computed_field,
    field_serializer,
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
    """

    id: str
    name: str = Field("")
    reactions: Set[Reaction] = Field(default_factory=set)
    user_compounds: Set[Compound] = Field(
        default_factory=set, alias="compounds"
    )
    model_config = ConfigDict(arbitrary_types_allowed=True)

    def __repr__(self):
        """Return a string representation of the reaction network."""
        return (
            f"<ReactionNetwork id={self.id}, name={self.name}, "
            f"num_reactions={len(self.reactions)}, num_compounds={len(self.compounds)}>"
        )

    def __eq__(self, other):
        """Check equality with another ReactionNetwork instance."""
        if not isinstance(other, ReactionNetwork):
            return NotImplemented

        return (
            self.reactions == other.reactions
            and self.compounds == other.compounds
        )

    @computed_field
    @property
    def compounds(self: "ReactionNetwork") -> Set[Compound]:
        """Add the compounds field."""
        compounds = {Compound.model_validate(c) for c in self.user_compounds}
        for reaction in self.reactions:
            for compound_id in reaction.stoichiometry.keys():
                new_compound = Compound(id=compound_id)
                if not any(c.id == compound_id for c in compounds):
                    warnings.warn(
                        f"adding auto-generated compound {new_compound}"
                    )
                    compounds.add(new_compound)
        return compounds

    @model_validator(mode="after")
    def check_all_compounds(self) -> "ReactionNetwork":
        """Check if the reaction network has all the compounds."""
        reaction_compounds = set()
        for reaction in self.reactions:
            # Access compound IDs
            reaction_compounds.update(reaction.stoichiometry.keys())
        model_compounds = {compound.id for compound in self.compounds}
        missing = reaction_compounds - model_compounds
        if missing != set():
            raise ValueError(f"Missing compounds in the model: {missing}")
        return self

    @property
    def reaction_adjacency_matrix(self: "ReactionNetwork") -> pd.DataFrame:
        """
        Convert ReactionNetwork into an adjacency matrix.

        Parameters
        ----------
        reaction_network : ReactionNetwork
            The reaction network to convert.

        Returns
        -------
        pd.DataFrame
            The adjacency matrix representing the reaction network. The row are representing reactants, and columns are products. The value is the reaction id.
        """
        # Extract all unique compounds
        all_compounds = set()
        for reaction in self.reactions:
            for cpd in reaction.stoichiometry.keys():
                all_compounds.add(cpd)
        ix = pd.MultiIndex.from_tuples(
            set(
                [
                    (cpd, pat.pattern_string)
                    for r in self.reactions
                    for cpd, tr in r.stoichiometry.items()
                    for pat in tr.keys()
                ]
            )
        ).sort_values()
        adj = pd.DataFrame("", index=ix, columns=ix).apply(
            lambda c: c.str.split()
        )
        for reaction in self.reactions:
            rid = reaction.id
            rid_rev = rid + "_rev"
            stoich = reaction.stoichiometry
            substrates = {
                c for c, t in stoich.items() if any(s < 0 for s in t.values())
            }
            products = {
                c for c, t in stoich.items() if any(s > 0 for s in t.values())
            }
            for sub in substrates:
                for spat in stoich[sub].keys():
                    for prod in products:
                        for ppat in stoich[prod].keys():
                            intersection = set(spat.pattern_string) & set(
                                ppat.pattern_string
                            )
                            if len(intersection) > 0:
                                adj.at[
                                    (sub, spat.pattern_string),
                                    (prod, ppat.pattern_string),
                                ].append(rid)
                                if reaction.reversible:
                                    adj.at[
                                        (prod, ppat.pattern_string),
                                        (sub, spat.pattern_string),
                                    ].append(rid_rev)
        return pd.DataFrame(adj)
