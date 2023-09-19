"""Functions performing Elementary Metabolite Unit calculations.

For the original paper see here:

    Antoniewicz, M. R., Kelleher, J. K., et al. (2007),
    Elementary Metabolite Units (EMU): A Novel Framework for Modeling Isotopic
    Distributions. doi:10.1016/j.ymben.2006.09.001


"""

from typing import List

from pydantic import BaseModel

from cmfa.fluxomics_data.fluxomics_dataset import FluxomicsDataset


class EMUTransition(BaseModel):
    """Stub for EMUTransition data model."""


class EMUReaction(BaseModel):
    """Stub for EMUReaction data model."""


def emu_decompose(ds: FluxomicsDataset) -> List[EMUReaction]:
    """Perform EMU decomposition."""
    return []


def emu_simulate(df: FluxomicsDataset):
    """Perform the second part of the EMU algorithm.

    Specifically, given a set of fluxes, a list of EMU reactions, and a set of
    known mass isotopomer distributions (i.e. tracers), find the steady state
    mass isotopomer distribution vector for any EMU.

    Questions:

     - what should the return value be?
    """
