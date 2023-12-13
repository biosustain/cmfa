"""Provides functions prepare_data_x.

These functions should take in a dataframe of measurements and return a
PreparedData object.
"""
import csv
import json
import os
import re
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd
import pandera as pa
from pandera.typing import DataFrame
from pydantic import BaseModel

from cmfa.fluxomics_data.compound import Compound
from cmfa.fluxomics_data.flux_measurement import FluxMeasurement
from cmfa.fluxomics_data.fluxomics_dataset import FluxomicsDataset
from cmfa.fluxomics_data.mid_measurement import (
    MIDMeasurement,
    MIDMeasurementComponent,
)
from cmfa.fluxomics_data.reaction import Reaction, ReactionDirection
from cmfa.fluxomics_data.reaction_network import ReactionNetwork
from cmfa.fluxomics_data.tracer import Tracer, TracerExperiment

FLUX_MEASUREMENTS_FILE = "flux.csv"
MS_MEASUREMENTS_FILE = "ms_measurements.csv"
REACTION_FILE = "reactions.csv"
TRACER_FILE = "tracers.csv"
MODEL_FILE = "model.json"
CUR_DIR = os.getcwd()
DATA_DIR = os.path.join(CUR_DIR, "data", "test_data")


def load_tracer_data(file_path: str) -> List[TracerExperiment]:
    """
    Load tracer data from a CSV file into Tracer and TracerExperiment objects.

    Parameters
    ----------
    file_path : str
        The path to the CSV file containing the tracer data.

    Returns
    -------
    List[Tracer]
        A list of Tracer objects from the file.
    List[TracerExperiment]
        A list of TracerExperiment objects loaded from the file.
    """
    tracer_experiments = {}
    tracers = {}
    with open(file_path, mode="r", encoding="utf-8") as file:
        reader = csv.DictReader(file)
        for row in reader:
            # Create Tracer object
            atom_positions = set(
                map(int, row["atom_ids"].strip("[]").split(","))
            )
            isotope = row["tracer_id"]
            if isotope not in tracers:
                tracers[isotope] = Tracer(
                    isotope=isotope,
                    compound=row["met_id"],
                    labelled_atom_positions=atom_positions,
                    purity=float(row["ratio"]),
                )
            print(tracers)
            # Add to TracerExperiment
            exp_id = row["experiment_id"]
            enrichment = float(row["enrichment"])
            if exp_id not in tracer_experiments:
                tracer_experiments[exp_id] = TracerExperiment(
                    experiment_id=exp_id,
                    tracer_enrichments={tracers[isotope].isotope: enrichment},
                )
            else:
                tracer_experiments[exp_id].tracer_enrichments[
                    tracers[isotope].isotope
                ] = enrichment
            print(tracer_experiments)
    return list(tracers.values()), list(tracer_experiments.values())


def load_flux_measurements(file_path: str) -> List[FluxMeasurement]:
    """
    Load flux measurements from a CSV file.

    Parameters
    ----------
    file_path : str
        The path to the CSV file containing flux measurements.

    Returns
    -------
    List[FluxMeasurement]
        A list of FluxMeasurement objects loaded from the file.
    """
    flux_measurements = []
    with open(file_path, mode="r", encoding="utf-8") as file:
        reader = csv.DictReader(file)
        for row in reader:
            flux_measurement = FluxMeasurement(
                experiment_id=row["experiment_id"],
                reaction_id=row["rxn_id"],
                replicate=row["replicate"],
                measured_flux=float(row["flux"]),
                measurement_error=float(row["flux_std_error"]),
            )
            flux_measurements.append(flux_measurement)
    print(flux_measurements)

    return flux_measurements


def load_mid_measurements(file_path: str) -> List[MIDMeasurement]:
    """
    Load MID measurements from a CSV file.

    Parameters
    ----------
    file_path : str
        The path to the CSV file containing MID measurements.

    Returns
    -------
    List[MIDMeasurement]
        A list of MIDMeasurement objects loaded from the file.
    """
    measurements = {}
    with open(file_path, mode="r", encoding="utf-8-sig") as file:
        reader = csv.DictReader(file)
        for row in reader:
            experiment_id = row["experiment_id"]
            compound_id = row["met_id"]
            fragment_id = row["ms_id"]

            # Group components into MIDMeasurements
            key = (experiment_id, compound_id, fragment_id)
            if key not in measurements:
                measurements[key] = MIDMeasurement(
                    experiment_id=experiment_id,
                    compound_id=compound_id,
                    fragment_id=fragment_id,
                )

            # Create MIDMeasurementComponent
            component = MIDMeasurementComponent(
                mass_isotopomer_id=row["mass_isotope"],
                measured_intensity=float(row["intensity"]),
                measured_std_dev=float(row["intensity_std_error"]),
            )
            measurements[key].measured_components.append(component)

    # Normalize intensities for each MIDMeasurement
    for mid_measurement in measurements.values():
        print(mid_measurement)
        mid_measurement.normalize_components()

    return list(measurements.values())


def parse_reaction_equation(
    equation: str, compounds_set: Set[Compound]
) -> Tuple[Dict[str, float], Dict[str, str], ReactionDirection]:
    """
    Parse a reaction equation into compounds, atom transitions, and direction.

    Parameters
    ----------
    equation : str
        The reaction equation string.

    Returns
    -------
    Tuple[Dict[str, float], Dict[str, str], ReactionDirection]
        A tuple containing:
        - A dictionary of compounds and their stoichiometric coefficients.
        - A dictionary of atom transitions for each compound.
        - The direction of the reaction (ReactionDirection).
    """
    compounds = {}
    atom_transitions = {}
    direction = ReactionDirection.REVERSIBLE

    # Remove spaces and identify the reaction direction
    equation = equation.replace(" ", "")
    if "<->" in equation:
        lhs, rhs = equation.split("<->")
        direction = ReactionDirection.REVERSIBLE
    elif "->" in equation:
        lhs, rhs = equation.split("->")
        direction = ReactionDirection.FORWARD
    else:
        raise ValueError(
            "Invalid reaction equation format, direction is missing"
        )

    # Define the regex pattern for parsing
    pattern = r"(\d*)([A-Za-z]+)\(([^)]+)\)"

    # Function to parse each side of the equation
    def parse_side(side: str, sign: int):
        for match in re.finditer(pattern, side):
            coeff_str, compound, atom_transition = match.groups()
            coeff = float(coeff_str) if coeff_str else 1.0
            coeff *= sign

            if compound not in compounds:
                compounds[compound] = 0
                atom_transitions[compound] = []

            compounds[compound] += coeff
            atom_transitions[compound].append(atom_transition)

            # Create and add new Compound objects
            new_compound = Compound(id=compound, carbon_label=atom_transition)
            compounds_set.add(new_compound)

    # Parse left-hand side and right-hand side
    parse_side(lhs, -1)  # Negative coefficients for LHS
    parse_side(rhs, 1)  # Positive coefficients for RHS

    return compounds, atom_transitions, direction


def load_reactions(
    file_path: str, network_id: str, network_name: Optional[str] = None
) -> ReactionNetwork:
    """
    Load the reaction network from a CSV file.

    Parameters
    ----------
    file_path : str
        The path to the CSV file containing the reaction network.

    Returns
    -------
    ReactionNetwork
        A reaction network that consists of reactions and compounds.
    """
    reactions_set: Set[Reaction] = set()
    compounds_set: Set[Compound] = set()

    with open(file_path, mode="r", encoding="utf-8") as file:
        reader = csv.DictReader(file)
        for row in reader:
            compounds, atom_transitions, direction = parse_reaction_equation(
                row["rxn_eqn"], compounds_set
            )
            reaction = Reaction(
                id=row["rxn_id"],
                name=row.get("rxn_name"),
                compounds=compounds,
                direction=direction,
                atom_transition=atom_transitions,
            )
            reactions_set.add(reaction)
    RN = ReactionNetwork(
        id=network_id,
        name=network_name,
        reactions=reactions_set,
        compounds=compounds_set,
    )
    print(RN)
    return RN


def prepare_data() -> FluxomicsDataset:
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
    print("Reading raw data...")

    print("\n----------\nloading tracer & experiments\n----------\n")
    tracers, tracer_experiments = load_tracer_data(f"{DATA_DIR}/{TRACER_FILE}")
    print("\n----------\nloading flux\n----------\n")
    flux_measurements = load_flux_measurements(
        f"{DATA_DIR}/{FLUX_MEASUREMENTS_FILE}"
    )
    print("\n----------\nloading mid measurements\n----------\n")
    mid_measurements = load_mid_measurements(
        f"{DATA_DIR}/{MS_MEASUREMENTS_FILE}"
    )
    print("\n----------\nloading reaction\n----------\n")
    reaction_network = load_reactions(f"{DATA_DIR}/{REACTION_FILE}", "test")
    print("\n----------\nmaking fluxomics dataset\n----------\n")

    FD = FluxomicsDataset(
        reaction_network=reaction_network,
        tracers=tracers,
        tracer_experiments=tracer_experiments,
        flux_measurements=flux_measurements,
        mid_measurements=mid_measurements,
    )
    print(repr(FD))
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


def import_fluxomics_dataset_from_json(filename: str) -> FluxomicsDataset:
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


FD = prepare_data()

export_fluxomics_dataset_to_json(FD, f"{DATA_DIR}/{MODEL_FILE}")

FD2 = import_fluxomics_dataset_from_json(f"{DATA_DIR}/{MODEL_FILE}")

print(FD == FD2)
