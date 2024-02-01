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
EXPECTED_REACTION_ADJACENCY_MATRIX = pd.DataFrame(
    [
        [pd.NA, "v1", pd.NA, pd.NA, pd.NA, pd.NA],
        ["v1_rev", pd.NA, "v4", "v2", "v4", "v5"],
        [pd.NA, "v4_rev", pd.NA, "v5", "v5", "v5"],
        [pd.NA, "v2_rev", "v5_rev", pd.NA, pd.NA, "v6"],
        [pd.NA, "v4_rev", "v5_rev", pd.NA, pd.NA, pd.NA],
        [pd.NA, "v5_rev", "v5_rev", "v6_rev", pd.NA, pd.NA],
    ],
    index=["A", "B", "C", "D", "E", "F"],
    columns=["A", "B", "C", "D", "E", "F"],
)


def test_example_network():
    """Test that loading a reaction networks works in a simple case."""
    rn = ReactionNetwork.model_validate(EXAMPLE_NETWORK_INPUT)
    pd.testing.assert_frame_equal(
        rn.reaction_adjacency_matrix, EXPECTED_REACTION_ADJACENCY_MATRIX
    )
