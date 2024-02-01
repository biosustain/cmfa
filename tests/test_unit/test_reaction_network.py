"""Unit tests for the reaction network data model."""

from copy import deepcopy

import pandas as pd
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
            stoichiometry={"A": {"abc": -1}, "B": {"abc": 1}},
        ),
        Reaction(
            id="v2",
            stoichiometry={"B": {"abc": -1}, "D": {"abc": 1}},
        ),
        Reaction(
            id="v3",
            stoichiometry={"D": {"abc": -1}, "B": {"abc": 1}},
        ),
        Reaction(
            id="v4",
            stoichiometry={"B": {"abc": -1}, "C": {"ab": 1}, "E": {"c": 1}},
        ),
        Reaction(
            id="v5",
            stoichiometry={
                "B": {"abc": -1},
                "C": {"de": -1},
                "D": {"bcd": 1},
                "E": {"a": 1},
                "F": {"e": 1},
            },
        ),
        Reaction(
            id="v6",
            stoichiometry={"D": {"abc": -1}, "F": {"abc": 1}},
        ),
    },
}


def test_example_network():
    """Test that loading a reaction networks works in a simple case."""
    rn = ReactionNetwork.model_validate(EXAMPLE_NETWORK_INPUT)
    assert rn.reaction_adjacency_matrix.loc[("D", "abc"), ("B", "abc")] == [
        "v2_rev",
        "v3",
    ]
