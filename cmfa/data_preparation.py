"""Provides functions prepare_data_x.

These functions should take in a dataframe of measurements and return a
PreparedData object.
"""

import json
import logging
import re
import warnings
from pathlib import Path
from typing import Dict, List, Set, Tuple

import pandas as pd

from cmfa.fluxomics_data.compound import Compound
from cmfa.fluxomics_data.flux_measurement import FluxMeasurement
from cmfa.fluxomics_data.fluxomics_dataset import FluxomicsDataset
from cmfa.fluxomics_data.mid_measurement import (
    MIDMeasurement,
    MIDMeasurementComponent,
)
from cmfa.fluxomics_data.reaction import Reaction
from cmfa.fluxomics_data.reaction_network import ReactionNetwork
from cmfa.fluxomics_data.tracer import Tracer, TracerExperiment


def parse_tracer_table(
    tracer_table: pd.DataFrame,
) -> Tuple[List[Tracer], List[TracerExperiment]]:
    """
    Read tracer data from a CSV file.

    Parameters
    ----------
    tracer_table : pd.DataFrame
        Result of running pd.read_csv against a suitable file, e.g.
        https://github.com/biosustain/cmfa/blob/main/data/test_data/tracers.csv

    Returns
    -------
    List[Tracer]
        A list of Tracer objects from the file.
    List[TracerExperiment]
        A list of TracerExperiment objects loaded from the file.
    """
    te_dict = dict()
    tracers = []
    for _, row in tracer_table.iterrows():
        if not any(t.isotope == row["tracer_id"] for t in tracers):
            new_tracer = Tracer(
                isotope=row["tracer_id"],
                compound=row["met_id"],
                labelled_atom_positions=set(json.loads(row["atom_ids"])),
                purity=row["ratio"],
            )
            tracers.append(new_tracer)
        if row["experiment_id"] not in te_dict.keys():
            te_dict[row["experiment_id"]] = dict()
        te_dict[row["experiment_id"]][row["tracer_id"]] = row["enrichment"]
    tracer_experiments = [
        TracerExperiment(experiment_id=k, tracer_enrichments=v)
        for k, v in te_dict.items()
    ]
    return tracers, tracer_experiments


def parse_flux_measurements(
    measurement_table: pd.DataFrame,
) -> List[FluxMeasurement]:
    """
    Parse flux measurements from a CSV file.

    Parameters
    ----------
    measurement_table : pd.DataFrame
        Result of running pd.read_csv against a suitable table, e.g.
        https://github.com/biosustain/cmfa/blob/main/data/test_data/flux.csv

    Returns
    -------
    List[FluxMeasurement]
        A list of FluxMeasurement objects loaded from the file.
    """
    flux_measurements = []
    for _, row in measurement_table.iterrows():
        flux_measurement = FluxMeasurement(
            experiment_id=row["experiment_id"],
            reaction_id=row["rxn_id"],
            replicate=row["replicate"],
            measured_flux=float(row["flux"]),
            measurement_error=float(row["flux_std_error"]),
        )
        flux_measurements.append(flux_measurement)
    return flux_measurements


def parse_mid_measurements(
    measurements_table: pd.DataFrame,
) -> List[MIDMeasurement]:
    """
    Load MID measurements from a CSV file.

    Parameters
    ----------
    file_path : pd.DataFrame
        Result of running pd.read_csv against a suitable table, e.g.
        https://github.com/biosustain/cmfa/blob/main/data/test_data/ms_measurements.csv

    Returns
    -------
    List[MIDMeasurement]
        A list of MIDMeasurement objects loaded from the file.
    """
    out = list()
    measurements_dict = dict()
    for _, row in measurements_table.iterrows():
        experiment_id = row["experiment_id"]
        compound_id = row["met_id"]
        fragment_id = row["ms_id"]
        key = (experiment_id, compound_id, fragment_id)
        if key not in measurements_dict.keys():
            measurements_dict[key] = {
                "experiment_id": experiment_id,
                "compound_id": compound_id,
                "fragment_id": fragment_id,
                "measured_components": [],
            }
        measurements_dict[key]["measured_components"].append(
            MIDMeasurementComponent(
                mass_isotopomer_id=str(row["mass_isotope"]),
                measured_intensity=float(row["intensity"]),
                measured_std_dev=float(row["intensity_std_error"]),
            )
        )
    for m in measurements_dict.values():
        out.append(MIDMeasurement.model_validate(m))
    return out


def parse_reaction_equation(
    equation: str,
) -> Tuple[Dict[str, Dict[str, float]], bool]:
    """
    Parse a reaction equation into a dictionary combining stoichiometry and atom transitions.

    Parameters
    ----------
    equation : str
        The reaction equation string.

    Returns
    -------
    Tuple[Dict[str, Dict[Tuple[int, ...], float]], bool]
        A tuple containing:
        - A dictionary mapping compounds to another dictionary with tuples of atom transitions
          and stoichiometric coefficients.
        - A boolean indicating if the reaction is reversible.
    """
    # Remove spaces and identify the reaction direction
    equation = equation.replace(" ", "")
    if "<->" in equation:
        lhs, rhs = equation.split("<->")
        reversible = True
    elif "->" in equation:
        lhs, rhs = equation.split("->")
        reversible = False
    else:
        raise ValueError(
            "Invalid reaction equation format, direction is missing"
        )

    stoichiometry = {}
    # Define the regex pattern for parsing
    pattern = r"(\d*\.*\d*)\**([A-Za-z0-9]+)(\(([^)]*)\))*"
    lhs, rhs = re.split("<->|->", equation)

    # Function to parse each side of the equation
    def parse_side(side: str, sign: int):
        for match in re.finditer(pattern, side):
            coeff_str, compound_id, _, atom_transition = match.groups()
            coeff = float(coeff_str) if coeff_str else 1.0

            if compound_id not in stoichiometry:
                stoichiometry[compound_id] = {}

            if atom_transition not in stoichiometry[compound_id]:
                stoichiometry[compound_id][atom_transition] = 0
            elif atom_transition == ():
                if atom_transition == None:
                    atom_transition = ""
                stoichiometry[compound_id][atom_transition] = 0

            stoichiometry[compound_id][atom_transition] += coeff * sign

    # Parse left-hand side and right-hand side
    parse_side(lhs, -1)  # Negative coefficients for LHS
    parse_side(rhs, 1)  # Positive coefficients for RHS

    return stoichiometry, reversible


def parse_reaction_table(
    reaction_table: pd.DataFrame,
    network_id: str,
    network_name: str = "",
) -> ReactionNetwork:
    """
    Load the reaction network from a CSV file.

    Parameters
    ----------
    reaction_table : pd.DataFrame
        Result of running pd.read_csv against a suitable table, e.g.
        https://github.com/biosustain/cmfa/blob/main/data/test_data/reactions.csv

    Returns
    -------
    ReactionNetwork
        A reaction network that consists of reactions and compounds.
    """
    reactions_set: Set[Reaction] = set()

    for _, row in reaction_table.iterrows():
        stoichiometry, reversible = parse_reaction_equation(str(row["rxn_eqn"]))
        reaction = Reaction(
            id=str(row["rxn_id"]),
            name=str(row["rxn_id"]),
            stoichiometry_input=stoichiometry,
            reversible=reversible,
        )
        reactions_set.add(reaction)
    return ReactionNetwork(
        id=network_id,
        name=network_name,
        reactions=reactions_set,
    )


def load_dataset_from_csv(
    tracer_file: Path,
    flux_measurement_file: Path,
    mid_measurement_file: Path,
    reaction_file: Path,
) -> FluxomicsDataset:
    """
    Load all existing data in a single model.

    Parameters
    ----------
    data_id : str
        The id for the dataset.

    Returns
    -------
    FluxomicsDataset
        A fluxomics dataset model that consists of fluxes, tracers, mid measurements, and a reaction network.

    """
    logging.info("Reading raw data...")
    tracer_table = pd.read_csv(tracer_file)
    flux_measurements_table = pd.read_csv(flux_measurement_file)
    mid_measurements_table = pd.read_csv(mid_measurement_file)
    reactions_table = pd.read_csv(reaction_file)
    logging.info("Parsing tables...")
    tracers, tracer_experiments = parse_tracer_table(tracer_table)
    flux_measurements = parse_flux_measurements(flux_measurements_table)
    mid_measurements = parse_mid_measurements(mid_measurements_table)
    reaction_network = parse_reaction_table(reactions_table, "a")
    logging.info("Aggregating...")
    FD = FluxomicsDataset(
        reaction_network=reaction_network,
        tracers=tracers,
        tracer_experiments=tracer_experiments,
        flux_measurements=flux_measurements,
        mid_measurements=mid_measurements,
    )
    logging.info("Created fluxomics dataset:\n" + repr(FD))
    return FD


def export_fluxomics_dataset_to_json(dataset: FluxomicsDataset, filename: str):
    """
    Export a FluxomicsDataset instance to a JSON file.

    Parameters
    ----------
    dataset : FluxomicsDataset
        The FluxomicsDataset instance to be exported.
    filename : str
        The path of the file where the JSON data will be saved.
    """
    with open(filename, "w", encoding="utf-8") as file:
        json_data = dataset.model_dump_json(indent=4)
        file.write(json_data)


def import_fluxomics_dataset_from_json(filename: Path) -> FluxomicsDataset:
    """
    Import a FluxomicsDataset instance from a JSON file.

    Parameters
    ----------
    filename : str
        The path of the JSON file to be imported.

    Returns
    -------
    FluxomicsDataset
        The imported FluxomicsDataset instance.
    """
    with open(filename, "r", encoding="utf-8") as file:
        json_data = file.read()
        dataset = FluxomicsDataset.model_validate_json(json_data)
        return dataset


def extract_atom_pattern(compound_pattern: str):
    """Extract the atom pattern from a compound pattern string.

    Parameters
    ----------
    compound_pattern : str
        The compound pattern string, e.g., "Glu(abcde)".

    Returns
    -------
    str
        The extracted atom pattern, e.g., "abcde".
    """
    match = re.search(r"(.+)\((.+)\)", compound_pattern)
    if match:
        compound, pattern = match.groups()
        return compound, pattern
    else:
        raise ValueError(f"No atom pattern found in the input: {compound}")
