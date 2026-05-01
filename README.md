# EMIR_pipeline
An interactive pipeline for reducing spectra from raw EMIR data.

The purpose of this repository is to create a file managing python interface for the PyEMIR analysis code. This pipeline is built on top of the [PyEMIR package](https://pyemir.readthedocs.io/en/v0.16.0/) maintained by Sergio Pascual and Nicolás Cardiel. 

## Installation

- Install PyEMIR __with python 3.11__ as specified [here](https://pyemir.readthedocs.io/en/stable/installation/index.html#install-in-conda). If you don't install with 3.11 or higher there will be issues with the interactive options.
- add `path/to/EMIR_pipeline/src` to your $PYTHONPATH
- store `path/to/EMIR_pipeline/src` to a system variable called $EMIR_PIPE
- `pip install PyQt5`
- `pip install mpl-interactions`

To run the notebooks in the examples in your new conda environment:
- `pip install jupyter`

When running any of these functions, be sure to activate your EMIR conda environment, probably with something like 
```
conda activate emir
```

## Known Issues / Required Patches

After installing PyEMIR (numina 0.35.2), two bugs in the numina package must be patched manually before the pipeline will run correctly.

### 1. Typo in `numina/array/interpolation.py`

numina 0.35.2 contains a typo (`'rigth'` instead of `'right'`) that causes a `ValueError` when running the wavelength calibration step.

Find the file inside your conda environment:
```
$CONDA_PREFIX/lib/python3.XX/site-packages/numina/array/interpolation.py
```

On line ~153, change:
```python
x_indices = np.searchsorted(self._x, v, side='rigth')
```
to:
```python
x_indices = np.searchsorted(self._x, v, side='right')
```

### 2. NumPy 2.0 incompatibility in `numina/array/wavecalib/fix_pix_borders.py`

`np.alltrue` was removed in NumPy 2.0. If you are using a recent NumPy, numina will crash with `AttributeError: module 'numpy' has no attribute 'alltrue'`.

Find the file inside your conda environment:
```
$CONDA_PREFIX/lib/python3.XX/site-packages/numina/array/wavecalib/fix_pix_borders.py
```

On line ~44, change:
```python
if not np.alltrue(sp == sought_value):
```
to:
```python
if not np.all(sp == sought_value):
```

## How to Use the Pipeline

There are a series of tutorials in the `examples/` file. It would be best to go through them in the following order:
1. [`OB0011_example.ipynb`](https://github.com/HOSTFLOWS/EMIR_pipeline/blob/main/examples/OB0011_example.ipynb)
2. [`OB0012_response_curve.ipynb`](https://github.com/HOSTFLOWS/EMIR_pipeline/blob/main/examples/OB0012_response_curve.ipynb)
3. [`OB0012_interactive_response_curve.ipynb`](https://github.com/HOSTFLOWS/EMIR_pipeline/blob/main/examples/OB0012_interactive_response_curve.ipynb)
4. [`OB0011_interactive_example.ipynb`](https://github.com/HOSTFLOWS/EMIR_pipeline/blob/main/examples/OB0011_interactive_example.ipynb)
