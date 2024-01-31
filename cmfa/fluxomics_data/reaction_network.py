"""reaction_network.py includes the description of compound, reaction and reaction_network."""

import warnings
from copy import deepcopy
from operator import gt, lt
from typing import List, Optional, Set

import pandas as pd
from pydantic import BaseModel, Field, computed_field, model_validator

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
            for compound, transitions in reaction.stoichiometry.items():
                all_compounds.add(compound)

        all_compounds_list = sorted(list(all_compounds))
        adjacency_matrix = pd.DataFrame(
            index=all_compounds_list, columns=all_compounds_list
        )

        for reaction in self.reactions:
            reactants = {
                compound
                for compound, transitions in reaction.stoichiometry.items()
                if any(coeff < 0 for coeff in transitions.values())
            }
            products = {
                compound
                for compound, transitions in reaction.stoichiometry.items()
                if any(coeff > 0 for coeff in transitions.values())
            }

            print(reactants, products)
            for reactant in reactants:
                for product in products:
                    adjacency_matrix.at[reactant, product] = reaction.id

            # Handle reversible reactions
            if reaction.reversible:
                for reactant in reactants:
                    for product in products:
                        adjacency_matrix.at[product, reactant] = (
                            reaction.id + "_rev"
                        )

        return adjacency_matrix


reactions = {
    Reaction(
        id="R1",
        reversible=True,
        stoichiometry={"A": {"ab": -1}, "B": {"ab": 1}},
    ),
    Reaction(
        id="R2",
        reversible=False,
        stoichiometry={"B": {"ab": -1}, "C": {"a": 1, "b": 1}},
    ),
    Reaction(
        id="R3",
        reversible=False,
        stoichiometry={
            "D": {"ab": -1},
            "A": {"ab": 0.5},
            "B": {"ab": 0.5},
            "Z": {"": 1},
        },
    ),
}

reaction_network = ReactionNetwork(id="a", reactions=reactions)
matrix = reaction_network.reaction_adjacency_matrix()
print(matrix)
