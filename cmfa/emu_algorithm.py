"""Functions performing Elementary Metabolite Unit Network Decomposition.

For the original paper see here:

    Antoniewicz, M. R., Kelleher, J. K., et al. (2007),
    Elementary Metabolite Units (EMU): A Novel Framework for Modeling Isotopic
    Distributions. doi:10.1016/j.ymben.2006.09.001


"""

from collections import deque

import pandas as pd

from cmfa.fluxomics_data.emu import EMU, EMUReaction
from cmfa.fluxomics_data.fluxomics_dataset import FluxomicsDataset
from cmfa.fluxomics_data.reaction import AtomPattern, Reaction
from cmfa.fluxomics_data.reaction_network import ReactionNetwork


def decompose_network(
    target_EMU: EMU, reaction_network: ReactionNetwork
) -> pd.DataFrame:
    """Return a EMU map with a list of EMUs based on the target compound from the reaction network."""
    MAM = reaction_network.reaction_adjacency_matrix()
    queue = deque([target_EMU])
    visited = set()
    emu_maps = {}
    emu_cpd_list = {}

    while queue:
        # Keep poping new EMU added from previous round.
        cur_emu = queue.popleft()
        # print(f"Current EMU: {cur_emu.compound}, have visited {visited}")
        if cur_emu in visited:
            continue
        visited.add(cur_emu)

        # Add new EMU map if the EMU size is not added
        emu_size = cur_emu.size
        if emu_size not in emu_maps:
            emu_maps[emu_size] = []
            emu_cpd_list[emu_size] = []

        # Find the reaction that is participated in the current emu
        participated_reactions = _find_participated_reactions(cur_emu, MAM)
        new_map = []
        emus = []
        for r in participated_reactions:
            new_emu_reactions, new_emus = _atom_mapping_to_reaction(
                cur_emu,
                reaction_network.find_reaction_by_id(r),
                reverse=False,
            )
            # print(new_emu_reactions, new_emus)
            new_map.append(new_emu_reactions)
            emus.append(new_emus)
            # TODO: Find equivalent EMUs
            for e in new_emus:
                if e not in visited and e not in queue:
                    queue.append(e)
                if e not in emu_cpd_list[emu_size]:
                    emu_cpd_list[emu_size] += [e]

        if len(new_map) > 0:
            for i in new_map:
                emu_maps[emu_size] += i

    return emu_maps, emu_cpd_list


def _find_participated_reactions(cur_emu, MAM):
    """
    Return the reactions in which a given EMU participates within the Metabolite Adjacency Matrix (MAM).

    Parameters
    ----------
    MAM: DataFrame representing the Metabolite Adjacency Matrix with multi-index columns.
    cur_emu: The current EMU object, expected to have at least a 'compound' attribute.

    Returns
    -------
    A DataFrame filtered to show only the reactions (columns) in which the current EMU's compound participates.
      The cells are boolean, indicating whether the reaction is non-empty (True) or not (False).
    """
    # Extract reactions for the current EMU's compound
    participated_reactions = MAM.xs(cur_emu.compound, level=0, axis=1).map(
        lambda x: x if isinstance(x, list) and len(x) > 0 else []
    )

    all_participated_reactions = set(
        item
        for sublist in participated_reactions.values.flatten()
        for item in sublist
    )
    return all_participated_reactions


def _atom_mapping_to_reaction(
    emu: EMU, reaction: Reaction, reverse: bool = False
) -> list:
    """Given an EMU and reaction information, returns an EMUReaction object."""
    try:
        atom_transitions = reaction.stoichiometry_input[emu.compound].keys()
    except:
        raise ValueError(
            f"Cannot find the EMU compound in the corresponding reaction"
        )

    emu_reactions = []
    # Use emu to find matching atoms in given compounds.
    # TODO: What if the emu matches two same compounds in reaction
    for transition in atom_transitions:
        adjusted_indices = [i - 1 for i in emu.atom_number]
        try:
            matched_atoms = "".join(transition[i] for i in adjusted_indices)
        except:
            raise ValueError(
                f"Cannot find atoms in compound {emu.compound} with indices {emu.atom_number_input}"
            )

        # Extracting the stoichiometry for the specific EMU
        emu_stoichiometry = {}
        emu_reactants = {emu}
        for cpd, stoich in reaction.stoichiometry_input.items():
            # new_stoich = {}
            for pattern, coeff in stoich.items():
                if cpd == emu.compound:  # Add the original EMU
                    # new_stoich[emu.atom_number_input] = -coeff if reverse else coeff
                    # emu_stoichiometry[cpd] = new_stoich
                    emu_stoichiometry[emu.id] = -coeff if reverse else coeff

                elif (not reverse and coeff < 0) or (
                    reverse and reaction.reversible and coeff > 0
                ):
                    # For matched reactants or products in reverse reaction
                    reactant_matched_atom_number = _find_matchable_atoms(
                        pattern, matched_atoms
                    )
                    if len(reactant_matched_atom_number) != 0:
                        emu_id = f"{cpd}_{reactant_matched_atom_number}"
                        emu_reactants.add(
                            EMU(
                                id=emu_id,
                                compound=cpd,
                                atom_number_input=reactant_matched_atom_number,
                            )
                        )
                        # emu_stoichiometry.setdefault(cpd, {})[reactant_matched_atom_number] = -coeff if reverse else coeff
                        emu_stoichiometry[emu_id] = -coeff if reverse else coeff

        emu_reactions.append(
            EMUReaction(
                reaction_id=reaction.id, emu_stoichiometry=emu_stoichiometry
            )
        )

    return emu_reactions, set(emu_reactants)


def _find_matchable_atoms(pattern: str, matched_atoms: str):
    """Return matchable atoms between a pattern and matched_atoms, returning the positions as a string."""
    matched_pattern = [
        str(pattern.index(char) + 1)
        for char in matched_atoms
        if char in pattern
    ]
    return "".join(sorted(matched_pattern))
