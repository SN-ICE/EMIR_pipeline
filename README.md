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

### Command-line scripts

The standard spectroscopy pipeline is [`reduce_spec.py`](/Users/lluisgalbany/EMIR_pipeline/reduce_spec.py). It is interactive by design: for each grism it lets you confirm the standard-star sensitivity-function anchor points and the SN trace position before flux calibration continues.

```bash
conda activate emir
export EMIR_PIPE=/path/to/EMIR_pipeline/src   # or add to your shell profile
export EMIR_DATA_DIR=/path/to/your/data       # optional: avoids --data-dir every time
```

---

#### `reduce_spec.py` — interactive spectroscopy reduction

Reduces a supernova OB and a standard-star OB end-to-end. Grisms are auto-detected from the QC files, standard-star magnitudes are fetched automatically from Simbad, and already-completed steps are skipped on re-runs. For each grism, the script opens two matplotlib windows:

1. **Sensitivity function fit** — drag the sliders to place the spline anchor points on the standard-star continuum.
2. **SN position finder** — drag the sliders to mark the positive and negative trace rows used for the final extraction.

Click **Finalize** in each window to save the chosen parameters and continue. On re-runs the sliders are pre-loaded with the previously saved values from `<OB>/results/interactive_parameters.pkl`.

```bash
python reduce_spec.py <SN_OB> <STD_OB> [options]
```

| Option | Description |
|---|---|
| `--data-dir PATH` | Root directory containing the OB folders. Defaults to `$EMIR_DATA_DIR` or cwd. |
| `--grisms G [G ...]` | Grism(s) to reduce, e.g. `YJ HK`. Defaults to all grisms common to both QC files. |
| `--clean` | Delete existing results and start from scratch. |
| `--no-plot` | Save the final figure but skip the final interactive spectrum window. |
| `--linear` | Use a linear y-axis scale instead of log. |
| `--out-dir PATH` | Output directory for the spectrum and figure. Defaults to `<SN_OB>/results/`. |
| `--skip-sens-interactive` | Skip the sensitivity-function window and reuse saved anchor points, or fall back to defaults if none exist. |

**Examples:**

```bash
# full interactive run
python reduce_spec.py OB0001 OB0002 --data-dir ~/Desktop/EMIR_REDUX

# skip the sensitivity-function window, but still confirm the SN trace interactively
python reduce_spec.py OB0001 OB0002 --skip-sens-interactive

# reduce a single grism on a linear scale
python reduce_spec.py OB0001 OB0002 --grisms YJ --linear

# with EMIR_DATA_DIR set, the flag can be omitted
python reduce_spec.py OB0001 OB0002
```

**Outputs** (written to `<SN_OB>/results/` by default):

- `SNNAME_DATEOBS_YJ.txt` — flux-calibrated YJ spectrum.
- `SNNAME_DATEOBS_HK.txt` — flux-calibrated HK spectrum.
- `SNNAME_DATEOBS.txt` — combined spectrum.
- `SNNAME_DATEOBS.png` — single combined plot with telluric bands marked and the propagated error band shown.

---

#### `reduce_imaging.py` — EMIR imaging reduction

For EMIR imaging OBs, use:

```bash
python reduce_imaging.py <OB_FOLDER> [--skip-check]
```

This expects an imaging-style OB containing at least `object/` and `flat/`, and writes the reduced imaging products into the OB’s reduction area.

---

### Jupyter notebook tutorials

For a step-by-step walkthrough of what each pipeline stage does, see the notebooks in `examples/`:

1. [`OB0011_interactive_example.ipynb`](https://github.com/SN-ICE/EMIR_pipeline/blob/main/examples/OB0011_interactive_example.ipynb)
2. [`OB0012_interactive_sensitivity_function.ipynb`](https://github.com/SN-ICE/EMIR_pipeline/blob/main/examples/OB0012_interactive_sensitivity_function.ipynb)
