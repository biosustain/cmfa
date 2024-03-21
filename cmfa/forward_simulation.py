"""Functions performing Elementary Metabolite Unit Matrix Calculations.

For the original paper see here:

    Antoniewicz, M. R., Kelleher, J. K., et al. (2007),
    Elementary Metabolite Units (EMU): A Novel Framework for Modeling Isotopic
    Distributions. doi:10.1016/j.ymben.2006.09.001


"""

import numpy as np
import pandas as pd
import sympy as sp


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

    matrix_B = matrix_B.map(lambda x: -x if not pd.isnull(x) else x)

    return matrix_A, matrix_B


def _process_matrix(matrix_df, fluxes):
    """
    Return the matrix based on reaction fluxes.

    Parameters
    ----------
    matrix_df : pd.DataFrame
        The matrix to be processed, typically representing stoichiometries or interactions within an EMU network.
    fluxes : dict
        A dictionary of reaction sympy symbols and flux values.

    Returns
    -------
    numpy.ndarray
        A NumPy array of calculated reaction fluxes.
    """
    # Define the lambda function for substitution and conversion
    process_element = lambda x: (
        float(x.subs(fluxes)) if isinstance(x, sp.Expr) else 0
    )

    # Apply the lambda function to each element in the DataFrame
    processed_matrix = matrix_df.map(process_element)

    # Convert the processed DataFrame to a NumPy array
    return processed_matrix.to_numpy(dtype=float)


def _solve_matrix(matrix_A, matrix_B, Y):
    """Return the MID of each emu given the A and B matrices."""
    try:
        inverse = np.linalg.inv(matrix_A)
    except:
        raise ValueError("Cannot find inverse of A matrix")
    return inverse @ matrix_B @ Y


def _construct_mid_from_emu(emu):
    """Return the MID matrix of an EMU."""
    return None
