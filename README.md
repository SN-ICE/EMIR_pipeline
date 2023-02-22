# EMIR_pipeline
A pipeline for reducing spectra from raw EMIR data.

The purpose of this repository is to create a file managing python interface for the PyEMIR analysis code. This pipeline is built on top of the [PyEMIR package](https://pyemir.readthedocs.io/en/v0.16.0/) maintained by Sergio Pascual and Nicolás Cardiel. 

## Installation

- Install PyEMIR as specified [here](https://pyemir.readthedocs.io/en/v0.16.0/installation/index.html#install-in-conda)
- add `path/to/EMIR_pipeline/src` to your $PYTHONPATH
- store `path/to/EMIR_pipeline/src` to a system variable called $EMIR_PIPE

When running any of these functions, be sure to activate your EMIR conda environment, probably with something like 
```
conda activate emir
```

## How to Use the Pipeline

Take a look at `EMIR_pipeline/examples/OB0011_example.ipynb` for an example of how to run this on one observation.
