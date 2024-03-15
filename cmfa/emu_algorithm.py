"""Functions performing Elementary Metabolite Unit calculations.

For the original paper see here:

    Antoniewicz, M. R., Kelleher, J. K., et al. (2007),
    Elementary Metabolite Units (EMU): A Novel Framework for Modeling Isotopic
    Distributions. doi:10.1016/j.ymben.2006.09.001


"""

from collections import deque

import pandas as pd
from sympy import Mul, Symbol

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
            # Handle forward and reverse case
            #
            if r.endswith("_rev"):
                new_emu_reactions, new_emus = _atom_mapping_to_reaction(
                    cur_emu,
                    reaction_network.find_reaction_by_id(r[:-4]),
                    reverse=True,
                )
            else:
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
            # emu_maps[emu_size].setdefault(cur_emu, []).append(new_map)
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


def create_emu_stoichiometry_matrix(emu_network, emu_size):
    """
    Return a stoichiometry matrix for a specified EMU size from a given EMU network.

    Parameters
    ----------
    emu_network : dict
        A dictionary representing the EMU network where keys are EMU sizes and values are lists of reaction dictionaries.
    emu_size : int
        The specific size of EMU for which the stoichiometry matrix is to be created.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame representing the stoichiometry matrix, where rows and columns are indexed by reactants and products.
    """
    # Create matrix indices
    reactants = set()
    products = set()

    for reaction in emu_network[emu_size]:
        # Directly extend the sets without checking the length
        reactants.update([" * ".join(reaction.reactants)])
        products.update([" * ".join(reaction.products)])

    # Combine reactants and products into a single tuple for output
    indices = tuple(reactants.union(products))

    # Initialize the DataFrame with indices for both rows and columns
    matrix = pd.DataFrame(index=indices, columns=indices, data=[])
    # Populate the matrix
    for reaction in emu_network[emu_size]:
        row_key, col_key = None, None  # Initialize keys

        for emu_id, stoich in reaction.emu_stoichiometry.items():
            if stoich < 0:  # Reactant
                for idx in indices:
                    # Find the matching reactant row
                    if emu_id in idx.split(" * "):
                        row_key = idx

            elif stoich > 0:  # Product
                for idx in indices:
                    # Find the matching product column
                    if emu_id in idx.split(" * "):
                        col_key = idx

            # Update the matrix cell with reaction ID and stoichiometry
            if row_key and col_key:
                cell_value = Mul(stoich, reaction.flux_symbol)
                matrix.at[row_key, col_key] = cell_value

    return matrix


def process_and_split_matrix(matrix):
    """
    Return A and B matrices according to the EMU algorithm.

    Parameters
    ----------
    matrix : pd.DataFrame
        The matrix to be processed and split, typically representing stoichiometries or interactions within an EMU network.

    Returns
    -------
    tuple of pd.DataFrame
        A tuple containing two pandas DataFrames after splitting and processing the input matrix. The first matrix's diagonal is adjusted to include non-NaN values from corresponding rows.
    """
    # Remove columns with all NaN values
    cleaned_matrix = matrix.dropna(axis=1, how="all")

    # Identify the row indices to split the matrix
    split_indices = [
        idx for idx in matrix.index if idx not in cleaned_matrix.columns
    ]

    # Split the matrix into two based on the identified indices
    matrix_A = (cleaned_matrix.loc[~cleaned_matrix.index.isin(split_indices)]).T
    matrix_B = (cleaned_matrix.loc[cleaned_matrix.index.isin(split_indices)]).T

    # Adjust the diagonal of the first matrix
    for i in range(len(matrix_A)):
        non_nan_values = (
            matrix_A.iloc[i].dropna().tolist()
            + matrix_B.iloc[i].dropna().tolist()
        )
        updated_values = 0
        for item in non_nan_values:
            # Copy the list to avoid modifying the original matrices
            updated_values += -1 * item.copy()

        # Assign the updated list to the diagonal element
        matrix_A.iat[i, i] = updated_values

    return matrix_A, matrix_B


def emu_simulate(df: FluxomicsDataset):
    """Perform the second part of the EMU algorithm.

    Specifically, given a set of fluxes, a list of EMU reactions, and a set of
    known mass isotopomer distributions (i.e. tracers), find the steady state
    mass isotopomer distribution vector for any EMU.

    Questions:

     - what should the return value be?
    """
