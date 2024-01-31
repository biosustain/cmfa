"""Functions performing Elementary Metabolite Unit calculations.

For the original paper see here:

    Antoniewicz, M. R., Kelleher, J. K., et al. (2007),
    Elementary Metabolite Units (EMU): A Novel Framework for Modeling Isotopic
    Distributions. doi:10.1016/j.ymben.2006.09.001


"""

from typing import Dict, List

import pandas as pd

from cmfa.fluxomics_data.emu_map import EMUMap, EMUReaction
from cmfa.fluxomics_data.fluxomics_dataset import FluxomicsDataset
from cmfa.fluxomics_data.reaction import Reaction
from cmfa.fluxomics_data.reaction_network import ReactionNetwork


def decompose_network(
    initial_emu: Dict[str, list[int]], reaction_network: ReactionNetwork
):
    """
    Decompose the Reaction network based on an initial EMU.

    Parameters
    ----------
    initial_emu : Dict[str]
        The starting point of decomposing network for EMUs.

    reaction_network : ReactionNetwork
        The reaction network from which EMUs are generated.

    Returns
    -------
    EMUMap
    """
    return []


def determine_emus(reaction: Reaction) -> List:
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
