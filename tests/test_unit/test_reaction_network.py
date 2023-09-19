"""Unit tests for the reaction network data model."""

from copy import deepcopy

import pytest

from cmfa.fluxomics_data.reaction_network import ReactionNetwork

EXAMPLE_NETWORK_INPUT = {
    "id": "Example metabolic network",
    "reactions": [
        {
            "id": "v1",
            "transitions": [
                {"compound_id": "A", "atom_pattern": "abc", "coef": -1},
                {"compound_id": "B", "atom_pattern": "abc", "coef": 1},
            ],
        },
        {
            "id": "v2",
            "transitions": [
                {"compound_id": "B", "atom_pattern": "abc", "coef": -1},
                {"compound_id": "D", "atom_pattern": "abc", "coef": 1},
            ],
        },
        {
            "id": "v3",
            "transitions": [
                {"compound_id": "D", "atom_pattern": "abc", "coef": -1},
                {"compound_id": "B", "atom_pattern": "abc", "coef": 1},
            ],
        },
        {
            "id": "v4",
            "transitions": [
                {"compound_id": "B", "atom_pattern": "abc", "coef": -1},
                {"compound_id": "C", "atom_pattern": "ab", "coef": 1},
                {"compound_id": "E", "atom_pattern": "c", "coef": 1},
            ],
        },
        {
            "id": "v5",
            "transitions": [
                {"compound_id": "B", "atom_pattern": "abc", "coef": -1},
                {"compound_id": "C", "atom_pattern": "de", "coef": -1},
                {"compound_id": "D", "atom_pattern": "bcd", "coef": 1},
                {"compound_id": "E", "atom_pattern": "a", "coef": 1},
                {"compound_id": "E", "atom_pattern": "e", "coef": 1},
            ],
        },
        {
            "id": "v6",
            "transitions": [
                {"compound_id": "D", "atom_pattern": "abc", "coef": -1},
                {"compound_id": "F", "atom_pattern": "abc", "coef": 1},
            ],
        },
    ],
}
EXAMPLE_NETWORK_INPUT_BAD = deepcopy(EXAMPLE_NETWORK_INPUT)
EXAMPLE_NETWORK_INPUT_BAD["reactions"][1]["transitions"][0][
    "atom_pattern"
] = "abcd"


def test_example_network():
    ReactionNetwork.model_validate(EXAMPLE_NETWORK_INPUT)


@pytest.mark.xfail
def test_example_network_bad():
    ReactionNetwork.model_validate(EXAMPLE_NETWORK_INPUT_BAD)
