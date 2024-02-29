"""Proof of concept script for compositional regression with Stan.

The script uses the package composition_stats to make a compositional dataset
from some hardcoded numbers, then fits this dataset using the Stan.

"""

from pathlib import Path

import numpy as np
from cmdstanpy import CmdStanModel
from composition_stats import clr, clr_inv

HERE = Path(__file__).parent
N_MEASUREMENT = 5
SIGMA = 0.2
PROBS = [
    (0.3, 0.7),
    (0.1, 0.3, 0.6),
    (0.1, 0.1, 0.2, 0.4, 0.2),
    (0.8, 0.2),
    (0.3, 0.5, 0.2),
]
y = []

for p_i in PROBS:
    trans_i = clr(p_i)
    y_trans_i = np.random.normal(trans_i, scale=SIGMA)
    y_i = list(clr_inv(y_trans_i))
    y += y_i

data = {
    "N": len(y),
    "N_measurement": N_MEASUREMENT,
    "y_sizes": [len(p_i) for p_i in PROBS],
    "stacked_y": y,
}
print(HERE)
model = CmdStanModel(stan_file=HERE / "stan" / "ragged_comp_demo.stan")
mcmc = model.sample(data=data)
