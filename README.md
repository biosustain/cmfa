CMFA
==============================

Compositional Metabolic Flux Analysis

# How to run the analysis

To run the analysis, run the command `make analysis` from the project root. This
will install a fresh virtual environment if one doesn't exist already, activate
it and install python dependencies and cmdstan, then run the analysis with the
following commands:

- `python cmfa/prepare_data.py`
- `python cmfa/sample.py`
- `jupyter execute cmfa/investigate.ipynb`

# How to create a pdf report

First make sure you have installed [quarto](https://https://quarto.org/).

Now run this command from the project root:

```
make docs
```




# How to run tests

Run this command from the project root:

```
python -m pytest
```
