"""Functions performing Elementary Metabolite Unit calculations.

For the original paper see here:

    Antoniewicz, M. R., Kelleher, J. K., et al. (2007),
    Elementary Metabolite Units (EMU): A Novel Framework for Modeling Isotopic
    Distributions. doi:10.1016/j.ymben.2006.09.001


"""

from collections import deque

import pandas as pd

from cmfa.fluxomics_data.compound import AtomPattern
from cmfa.fluxomics_data.fluxomics_dataset import FluxomicsDataset
from cmfa.fluxomics_data.reaction import Reaction
from cmfa.fluxomics_data.reaction_network import ReactionNetwork


def decompose_network(
    target_compound: str, reaction_network: ReactionNetwork
) -> pd.DataFrame:
    """
    Decompose the Reaction network based on an initial EMU.

    Parameters
    ----------
    initial_compound:
        The starting compound id` of decomposing network for EMUs.

    reaction_network : ReactionNetwork
        The reaction network from which EMUs are generated.

    Returns
    -------
    EMUMap: pd.DataFrame

    """
    MAM = reaction_network.reaction_adjacency_matrix()

    queue = deque([target_compound])
    visited = set()
    emu_maps = {}

    while queue:
        # Keep poping new EMU added from previous round.
        cur_emu = queue.popleft()
        print(f"Current EMU: {cur_emu}, have visited {visited}")
        if cur_emu in visited:
            continue
        visited.add(cur_emu)

        # Pattern matching
        _cur_emu = EMU(
            id="_".join(cur_emu),
            compound=cur_emu[0],
            pattern=AtomPattern(pattern=cur_emu[1]),
        )

        emu_size = _cur_emu.size
        if emu_size not in emu_maps:
            emu_maps[emu_size] = {}

        # Find the reaction that is participated in the current emu
        col = MAM[cur_emu]
        reactants = col[col.apply(lambda x: len(x) > 0)]

        # Finding multiple reactants
        _reactants_str = reactants.apply(lambda x: str(x))
        multi_reactants = _reactants_str.duplicated(keep=False)

        if multi_reactants.any():
            # Use the duplicates mask to filter the original Series and get the multi-index keys
            duplicate_keys = reactants[multi_reactants].index.tolist()

            overlap_pat = _find_matchable_parts(duplicate_keys, cur_emu[1])

            print("Duplicate multi-index keys:", overlap_pat)
        else:
            print("No duplicates found.")

        # Add reactants to the queue and emu map
        print(
            reactants,
            "\n",
        )
        for i in reactants.index:
            print(i)
            match_cpd_pat = [key for key in MAM.index if key[0] == i[0]]
            for m in match_cpd_pat:
                queue.append(tuple(m))
            # Add a condition to stop when the number of atoms are not matching -> Stop when this happens(B23 + C1)
        for next_emu in reactants:
            emu_maps[emu_size].setdefault(cur_emu, []).append(next_emu)

    return emu_maps


def _find_matchable_parts(keys, match):
    """Return new EMU reactants that maps partial atoms if available."""
    new_keys = []
    for key in keys:
        # Unpack the tuple
        letter, sequence = key

        # Check for matchable parts and keep characters that are found in the matchable part
        matched_sequence = "".join([char for char in sequence if char in match])

        # Construct a new tuple with the matched parts
        new_key = (letter, matched_sequence)
        new_keys.append(new_key)

    return new_keys


def determine_emus(reaction: Reaction) -> list:
    """
    Determine the EMUs in a given reaction.

    Parameters
    ----------
    reaction : Reaction
        The reaction to analyze.

    Returns
    -------
    List
        A list of EMUs involved in the reaction.
    """
    return []


def create_emu_reaction(reaction: Reaction, emu) -> EMUReaction:
    """
    Create an EMUReaction based on a reaction and an EMU.

    Parameters
    ----------
    reaction : Reaction
        The reaction from which the EMUReaction is derived.
    emu
        The EMU for which the EMUReaction is created.

    Returns
    -------
    EMUReaction
        The created EMUReaction.
    """
    return EMUReaction(...)


def emu_simulate(df: FluxomicsDataset):
    """Perform the second part of the EMU algorithm.

    Specifically, given a set of fluxes, a list of EMU reactions, and a set of
    known mass isotopomer distributions (i.e. tracers), find the steady state
    mass isotopomer distribution vector for any EMU.

    Questions:

     - what should the return value be?
    """
