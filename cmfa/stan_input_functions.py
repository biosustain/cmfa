"""Functions for generating input to Stan from prepared data."""


from typing import Dict

from cmfa.data_preparation import PreparedData


# todo: complete this function
def get_stan_input(prepared_data: PreparedData) -> Dict:
    """General function for creating a Stan input."""
    measurements = prepared_data.measurements
    return {
        "N": len(measurements),
        "N_train": len(measurements),
        "N_test": len(measurements),
        "ix_train": [i + 1 for i in range(len(measurements))],
        "ix_test": [i + 1 for i in range(len(measurements))],
    }
