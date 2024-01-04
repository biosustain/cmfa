"""Unit tests for the reaction network data model."""

from copy import deepcopy

import pytest

from cmfa.fluxomics_data.compound import Compound
from cmfa.fluxomics_data.reaction import Reaction
from cmfa.fluxomics_data.reaction_network import ReactionNetwork

EXAMPLE_NETWORK_INPUT = {
    "id": "Example metabolic network",
    "compounds": {
        Compound(id="B", name="compound B", formula="b2o"),
        Compound(id="C", name="compound C", formula="c2o"),
        Compound(id="D", name="compound D", formula="d2o"),
        Compound(id="E", name="compound E", formula="e2o"),
        Compound(id="F", name="compound F", formula="f2o"),
    },
    "reactions": {
        Reaction(
            id="v1",
            stoichiometry={"A": -1, "B": 1},
            atom_transition={"A": ["abc"], "B": ["abc"]},
        ),
        Reaction(
            id="v2",
            stoichiometry={"B": -1, "D": 1},
            atom_transition={"B": ["abc"], "D": ["abc"]},
        ),
        Reaction(
            id="v3",
            stoichiometry={"D": -1, "B": 1},
            atom_transition={"D": ["abc"], "B": ["abc"]},
        ),
        Reaction(
            id="v4",
            stoichiometry={"B": -1, "C": 1, "E": 1},
            atom_transition={"B": ["abc"], "C": ["ab"], "E": ["c"]},
        ),
        Reaction(
            id="v5",
            stoichiometry={"B": -1, "C": -1, "D": 1, "E": 1, "F": 1},
            atom_transition={
                "B": ["abc"],
                "C": ["de"],
                "D": ["bcd"],
                "E": ["a"],
                "F": ["e"],
            },
        ),
        Reaction(
            id="v6",
            stoichiometry={"D": -1, "F": 1},
            atom_transition={"D": ["abc"], "F": ["abc"]},
        ),
    },
}


def test_example_network():
    """Test that loading a reaction networks works in a simple case."""
    ReactionNetwork.model_validate(EXAMPLE_NETWORK_INPUT)
