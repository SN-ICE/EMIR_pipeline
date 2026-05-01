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

The quickest way to reduce a pair of observations is with the two driver scripts in the repository root. Both require the `emir` conda environment to be active and `$EMIR_PIPE` to be set.

```bash
conda activate emir
export EMIR_PIPE=/path/to/EMIR_pipeline/src   # or add to your shell profile
export EMIR_DATA_DIR=/path/to/your/data       # optional: avoids --data-dir every time
```

---

#### `reduce.py` — fully automatic reduction

Reduces a supernova OB and a standard-star OB end-to-end without any user interaction. Grisms are auto-detected from the QC files, standard-star magnitudes are fetched automatically from Simbad, and already-completed steps are skipped on re-runs.

```bash
python reduce.py <SN_OB> <STD_OB> [options]
```

| Option | Description |
|---|---|
| `--data-dir PATH` | Root directory containing the OB folders. Defaults to `$EMIR_DATA_DIR` or cwd. |
| `--grisms G [G ...]` | Grism(s) to reduce, e.g. `YJ HK`. Defaults to all grisms common to both QC files. |
| `--clean` | Delete existing results and start from scratch. |
| `--no-plot` | Save the figure but skip the interactive plot window. |
| `--linear` | Use a linear y-axis scale instead of log. |
| `--out-dir PATH` | Output directory for the spectrum and figure. Defaults to `<SN_OB>/results/`. |

**Examples:**

```bash
# reduce all common grisms, auto-fetch magnitudes, show plot
python reduce.py OB0001 OB0002 --data-dir ~/Desktop/EMIR_REDUX

# reduce YJ only, linear scale, no interactive window
python reduce.py OB0001 OB0002 --grisms YJ --linear --no-plot

# with EMIR_DATA_DIR set, the flag can be omitted
python reduce.py OB0001 OB0002
```

**Outputs** (written to `<SN_OB>/results/` by default):

- `<SN_OB>_spectrum.txt` — flux-calibrated spectrum, wavelength-sorted (YJ then HK), with a header including `OBJECT`, `GRISM`, `DATE OBS`, exposure time, and airmass.
- `<SN_OB>_spectrum.png` — plot with telluric bands marked.

---

#### `reduce_interactive.py` — reduction with interactive parameter tuning

Same as `reduce.py` but pauses twice per grism to open interactive matplotlib windows:

1. **Sensitivity function fit** — drag sliders to set the spline knot positions used to fit the standard-star continuum.
2. **SN position finder** — drag sliders to mark the A (positive) and B (negative) spectral trace rows in the ABBA image.

Click **Finalize** in each window to save the chosen parameters and continue. On re-runs the sliders are pre-loaded with the previously saved values.

```bash
python reduce_interactive.py <SN_OB> <STD_OB> [options]
```

All options from `reduce.py` are available, plus:

| Option | Description |
|---|---|
| `--skip-sens-interactive` | Skip the sensitivity function window; reuse saved parameters or fall back to defaults. Useful once you are happy with the sensitivity function from a prior run. |

**Examples:**

```bash
# full interactive run
python reduce_interactive.py OB0001 OB0002 --data-dir ~/Desktop/EMIR_REDUX

# skip the sensitivity function window, still pick the SN position interactively
python reduce_interactive.py OB0001 OB0002 --skip-sens-interactive

# reduce a single grism interactively, linear scale
python reduce_interactive.py OB0001 OB0002 --grisms YJ --linear
```

---

### Jupyter notebook tutorials

For a step-by-step walkthrough of what each pipeline stage does, see the notebooks in `examples/`:

1. [`OB0011_example.ipynb`](https://github.com/SN-ICE/EMIR_pipeline/blob/main/examples/OB0011_example.ipynb)
2. [`OB0012_response_curve.ipynb`](https://github.com/SN-ICE/EMIR_pipeline/blob/main/examples/OB0012_response_curve.ipynb)
3. [`OB0012_interactive_response_curve.ipynb`](https://github.com/SN-ICE/EMIR_pipeline/blob/main/examples/OB0012_interactive_response_curve.ipynb)
4. [`OB0011_interactive_example.ipynb`](https://github.com/SN-ICE/EMIR_pipeline/blob/main/examples/OB0011_interactive_example.ipynb)
