#!/usr/bin/env python3
"""
EMIR interactive pipeline driver script.

Mirrors the spectroscopy reduction workflow but adds two interactive
matplotlib windows per grism:
  1. Standard-star sensitivity function fit — drag sliders to position the
     spline knots used to fit the std-star continuum.
  2. SN position finder — drag sliders to mark the A (positive) and B
     (negative) spectral trace rows in the ABBA image.

Both windows show a "Finalize" button; clicking it saves the chosen
parameters to <OB>/results/interactive_parameters.pkl and closes the
window. Re-running the script will pre-load those saved values as
initial slider positions.

Usage
-----
    python reduce_spec.py <SN_OB> <STD_OB> [options]

Positional arguments
    SN_OB     Name or full path of the supernova observation directory.
    STD_OB    Name or full path of the standard-star observation directory.

Optional arguments
    --data-dir PATH   Root directory that contains the OB folders.
                      Defaults to $EMIR_DATA_DIR if set, otherwise cwd.
    --grisms  G ...   Grism(s) to reduce (e.g. YJ HK).
                      Defaults to all grisms common to both QC files.
    --clean           Delete existing results and start fresh.
    --no-plot         Do not display the final spectrum window.
    --linear          Use linear y-axis scale instead of log.
    --out-dir PATH    Where to write the final spectrum and figure.
                      Defaults to <SN_OB>/results/.
    --skip-sens-interactive
                      Skip the interactive sensitivity function fit and
                      use any previously saved parameters (or defaults).

Example
-------
    python reduce_spec.py OB0001 OB0002
    python reduce_spec.py OB0001 OB0002 --data-dir ~/Desktop/EMIR_REDUX --grisms YJ
"""

import argparse
import os
import shutil
import sys

# ------------------------------------------------------------------
# src/ on the path regardless of working directory
# ------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.join(SCRIPT_DIR, 'src')
sys.path.insert(0, SRC_DIR)
os.environ.setdefault('EMIR_PIPE', SRC_DIR)

# Interactive matplotlib backend must be set before importing pyplot
import matplotlib
matplotlib.use('TkAgg')

import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt
import matplotlib as mpl

from Observation import Observation, DEFAULT_FLAT_FIELD_FILE, FLAT_FIELD_FILENAME
from reduction_tools import get_std_magnitudes
from interactive_reduction_tools import (
    interactive_SN_position_finder,
    interactive_sens_function_fit,
    get_parameters,
)
from templates import OBS_RES_TEMPLATE


# ------------------------------------------------------------------
# Observation fallback patches
# ------------------------------------------------------------------

def _iter_sibling_ob_dirs(self):
    parent_dir = os.path.dirname(self.obs_dir.rstrip('/'))
    if not os.path.isdir(parent_dir):
        return

    for entry in sorted(os.listdir(parent_dir)):
        sibling_dir = os.path.join(parent_dir, entry)
        if sibling_dir == self.obs_dir:
            continue
        if os.path.isdir(sibling_dir) and entry.startswith('OB'):
            yield sibling_dir


def _collect_calibration_files_from_directory(path, discovered, expected_filters):
    if not os.path.isdir(path):
        return

    for fname in sorted(os.listdir(path)):
        if not fname.lower().endswith((".fits", ".fit", ".fts")):
            continue

        try:
            header = fits.getheader(os.path.join(path, fname))
        except OSError:
            continue

        filter_name = header.get('FILTER')
        if filter_name is None:
            continue
        if expected_filters is not None and filter_name not in expected_filters:
            continue

        discovered.setdefault(filter_name, [])
        if fname not in discovered[filter_name]:
            discovered[filter_name].append(fname)


def _discover_files_from_headers(self, path):
    """
    Groups FITS files in a directory by FILTER and supplements missing
    calibration grisms from sibling OB folders under the same parent.
    """
    discovered = {}
    expected_filters = set(
        row['FILTER'] for row in self.qc.object_table if row['GRISM'] != 'OPEN'
    )

    _collect_calibration_files_from_directory(path, discovered, expected_filters)

    missing_filters = expected_filters - set(discovered.keys())
    section_name = os.path.basename(path.rstrip('/'))
    for sibling_dir in self._iter_sibling_ob_dirs():
        if not missing_filters:
            break
        sibling_section = os.path.join(sibling_dir, section_name)
        _collect_calibration_files_from_directory(
            sibling_section,
            discovered,
            missing_filters,
        )
        missing_filters = expected_filters - set(discovered.keys())

    return discovered


def _resolve_data_file(self, analysis_type, fname):
    local_path = os.path.join(self.obs_dir, analysis_type, fname)
    if os.path.exists(local_path):
        return local_path

    for sibling_dir in self._iter_sibling_ob_dirs():
        candidate = os.path.join(sibling_dir, analysis_type, fname)
        if os.path.exists(candidate):
            return candidate

    return None


def _copy_files_with_sibling_fallback(self, analysis_type, grism_type):
    try:
        file_dict = getattr(self, "%s_files" % analysis_type)
        data_file_list = file_dict[grism_type]
        if analysis_type not in ['object', 'arc', 'flat']:
            raise ValueError(
                "'%s' is not a valid analysis type (options: 'object', 'arc', 'flat')"
                % analysis_type
            )

        for fname in data_file_list:
            source_path = self._resolve_data_file(analysis_type, fname)
            if source_path is None:
                print('WARNING: Qualty control file lists a file that does not exist (%s)' % fname)
                continue
            target_path = os.path.join(self.emir_path, "data", fname)
            shutil.copyfile(source_path, target_path)

    except KeyError:
        raise KeyError(
            "%s is not a valid grism type for %s analysis (options: %s)"
            % (grism_type, analysis_type, str(list(file_dict.keys())))
        )


def _generate_flat_field_with_sibling_fallback(self):
    """
    Takes flat field files and generates the flat field file for EMIR analysis,
    allowing missing flats to be sourced from sibling OB folders.
    """
    for fil in self.flat_files:
        ff_coadd = None
        for f in self.flat_files[fil]:
            flat_path = self._resolve_data_file('flat', f)
            if flat_path is None:
                raise FileNotFoundError(
                    "Missing flat-field file %s for observation %s" % (f, self.name)
                )
            fi = fits.open(flat_path)
            im = fi[-1].data
            ff_coadd = im if ff_coadd is None else ff_coadd + im
            fi.close()

        emir_defaults_dir = os.path.join(os.environ['EMIR_PIPE'], 'default_EMIR_files')
        default_fits = fits.open(os.path.join(emir_defaults_dir, DEFAULT_FLAT_FIELD_FILE))
        default_fits[0].data = ff_coadd
        default_fits.writeto(
            os.path.join(self.emir_path, 'data', FLAT_FIELD_FILENAME % fil),
            overwrite=True,
        )
        default_fits.close()


def _generate_obs_res_with_sibling_fallback(self, analysis_type, grism_type):
    """
    Generates EMIR obs-result YAML files and keeps borrowed arc/flat files in
    the run even when they live in another sibling OB folder.
    """
    file_list = ''
    with open(os.path.join(self.emir_path, '%s_obs_res_%s.yaml' % (analysis_type, grism_type)), 'w') as f:
        try:
            all_files = getattr(self, "%s_files" % analysis_type)[grism_type]
            for i, filename in enumerate(all_files):
                fpath = self._resolve_data_file(analysis_type, filename)
                if fpath is None:
                    print("WARNING: skipping file %s" % filename)
                    continue

                if analysis_type == 'arc':
                    file_list += " - %s\n" % filename
                else:
                    file_list = " - %s\n" % filename
                    text = OBS_RES_TEMPLATE % (self.name + "_" + filename.split('-')[0], file_list)
                    f.write(text)
                    if i < len(all_files) - 1:
                        f.write('---\n')

            if analysis_type == 'arc' and file_list:
                text = OBS_RES_TEMPLATE % (self.name + '_arc', file_list)
                f.write(text)

        except KeyError:
            pass


Observation._iter_sibling_ob_dirs = _iter_sibling_ob_dirs
Observation._discover_files_from_headers = _discover_files_from_headers
Observation._resolve_data_file = _resolve_data_file
Observation._copy_files = _copy_files_with_sibling_fallback
Observation._generate_flat_field = _generate_flat_field_with_sibling_fallback
Observation._generate_obs_res = _generate_obs_res_with_sibling_fallback


# ------------------------------------------------------------------
# Argument parsing
# ------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description='Interactively reduce a GTC/EMIR long-slit spectrum.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument('sn_ob',  metavar='SN_OB',  help='Supernova OB directory (name or full path)')
    p.add_argument('std_ob', metavar='STD_OB', help='Standard-star OB directory (name or full path)')
    p.add_argument('--data-dir', default=None,
                   help='Root directory containing the OB folders '
                        '(default: $EMIR_DATA_DIR or cwd)')
    p.add_argument('--grisms', nargs='+', default=None,
                   help='Grism(s) to reduce, e.g. YJ HK (default: all common grisms)')
    p.add_argument('--clean', action='store_true',
                   help='Delete existing results and start fresh')
    p.add_argument('--no-plot', action='store_true',
                   help='Skip the final interactive spectrum window')
    p.add_argument('--linear', action='store_true',
                   help='Use linear y-axis scale instead of log')
    p.add_argument('--out-dir', default=None,
                   help='Output directory for spectrum and figure '
                        '(default: <SN_OB>/results/)')
    p.add_argument('--skip-sens-interactive', action='store_true',
                   help='Skip interactive sensitivity function fit; '
                        'use saved parameters or defaults')
    return p.parse_args()


# ------------------------------------------------------------------
# Helpers (shared with reduce.py)
# ------------------------------------------------------------------

def resolve_ob_path(ob_arg, data_dir):
    if os.path.isabs(ob_arg):
        return ob_arg
    candidate = os.path.join(data_dir, ob_arg)
    if os.path.isdir(candidate):
        return candidate
    if os.path.isdir(ob_arg):
        return os.path.abspath(ob_arg)
    raise FileNotFoundError(
        "OB directory '%s' not found (tried '%s' and cwd)" % (ob_arg, candidate)
    )


def abba_done(ob_path, grism):
    return os.path.exists(
        os.path.join(ob_path, 'results', 'ABBA', 'ABBA_subtracted_%s.fits' % grism)
    )


def sens_done(ob_path, grism):
    return os.path.exists(
        os.path.join(ob_path, 'results', 'sens_function', '%s_sens_function.fits' % grism)
    )


def section(msg):
    print('\n' + '=' * 60)
    print('  ' + msg)
    print('=' * 60)


# ------------------------------------------------------------------
# Main reduction logic
# ------------------------------------------------------------------

def reduce(sn_path, std_path, grisms, clean, no_plot, out_dir,
           linear=False, skip_sens_interactive=False):

    section('Initialising observations')
    print('  SN  :', sn_path)
    print('  STD :', std_path)

    osn = Observation(sn_path)
    ost = Observation(std_path)

    if clean:
        print('  --clean: removing existing results')
        osn._clean_files()
        ost._clean_files()
        osn = Observation(sn_path)
        ost = Observation(std_path)

    def spectro_grisms_for(obs):
        return set(row['FILTER'] for row in obs.qc.object_table if row['GRISM'] != 'OPEN')

    sn_grisms  = spectro_grisms_for(osn)
    std_grisms = spectro_grisms_for(ost)
    common_grisms = sorted(sn_grisms & std_grisms)

    if not common_grisms:
        raise ValueError(
            'No common spectroscopic grisms found between %s and %s'
            % (sn_path, std_path)
        )

    if sn_grisms - std_grisms:
        print('  Warning: grism(s) %s present in SN but not in std star — skipping'
              % sorted(sn_grisms - std_grisms))

    if grisms is None:
        grisms = common_grisms
    else:
        unknown = [g for g in grisms if g not in common_grisms]
        if unknown:
            raise ValueError(
                "Grism(s) %s not available in both observations. Common grisms: %s"
                % (unknown, common_grisms)
            )

    print('  Grisms to reduce:', grisms)

    # ------ per-grism reduction + interactive steps ------
    sn_kwargs = {}   # keyed by grism

    for grism in grisms:

        # --- SN ABBA ---
        section('Grism %s — SN reduction' % grism)
        if abba_done(sn_path, grism):
            print('  ABBA already exists, skipping rectification')
        else:
            osn.initialize('object', grism)
            osn.rectify_and_analyze('object', grism)
            osn.ABBA_subtract(grism)

        # --- Std ABBA ---
        section('Grism %s — standard-star reduction' % grism)
        if abba_done(std_path, grism):
            print('  ABBA already exists, skipping rectification')
        else:
            ost.initialize('object', grism)
            ost.rectify_and_analyze('object', grism)
            ost.ABBA_subtract(grism)

        # --- Interactive sensitivity function fit ---
        if not skip_sens_interactive:
            section('Grism %s — interactive sensitivity function fit' % grism)
            print('  Adjust the spline knots to fit the std-star continuum.')
            print('  Click "Finalize" when done.\n')
            std_kwargs = interactive_sens_function_fit(ost, grism)
            if std_kwargs is None or 'sample_points' not in std_kwargs:
                raise RuntimeError(
                    'Interactive sensitivity fit for grism %s was closed before '
                    'clicking Finalize' % grism
                )
        else:
            # Use saved parameters (or defaults) to (re)build the sens function
            if sens_done(std_path, grism):
                print('  Sensitivity function already exists, skipping')
            else:
                section('Grism %s — sensitivity function (non-interactive)' % grism)
                std_kwargs = get_parameters(ost, grism)
                ost.make_sens_function(grism, **std_kwargs)

        # --- Interactive SN position finder ---
        section('Grism %s — interactive SN position finder' % grism)
        print('  Drag the sliders to align the extraction window with the SN trace.')
        print('  Click "Finalize" when done.\n')
        sn_result = interactive_SN_position_finder(osn, grism)
        if sn_result is None or 'SN_position' not in sn_result:
            raise RuntimeError(
                'Interactive SN position finder for grism %s was closed before '
                'clicking Finalize' % grism
            )
        sn_kwargs[grism] = sn_result
        print('  SN position saved:', sn_kwargs[grism].get('SN_position'))

    # ------ magnitude lookup ------
    section('Fetching standard-star 2MASS magnitudes from Simbad')
    magnitudes = get_std_magnitudes(ost)

    # ------ flux calibration ------
    section('Flux calibration')
    results = {}
    for grism in grisms:
        print('  Grism %s ...' % grism)
        wave, flux, flux_err, flux_err_stat = osn.get_reduced_spectrum(
            magnitudes, ost, grism, **sn_kwargs.get(grism, {})
        )
        results[grism] = (wave, flux, flux_err, flux_err_stat)

    # ------ save spectrum ------
    if out_dir is None:
        out_dir = os.path.join(sn_path, 'results')
    os.makedirs(out_dir, exist_ok=True)

    output_paths = osn.save_output_spectra(out_dir)
    combined_spec_path = output_paths['combined']
    print('\n  Spectra saved to')
    for grism in sorted(k for k in output_paths if k != 'combined'):
        print('   %s : %s' % (grism, output_paths[grism]))
    print('   combined : %s' % combined_spec_path)

    # ------ plot ------
    section('Plotting')
    colors = {'YJ': 'royalblue', 'HK': 'firebrick'}
    default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    mpl.rcParams.update({'font.size': 13, 'font.family': 'serif'})
    fig, ax = plt.subplots(figsize=(14, 5))

    for i, (grism, (wave, flux, flux_err, flux_err_stat)) in enumerate(results.items()):
        color = colors.get(grism, default_colors[i % len(default_colors)])
        ax.plot(wave, flux, color=color, lw=0.8, label=grism)
        upper = flux + flux_err
        lower = flux - flux_err
        if not linear:
            lower = np.where(lower > 0, lower, np.nan)
        ax.fill_between(wave, lower, upper, color=color, alpha=0.18, linewidth=0)

    # telluric bands
    telluric_bands = [(13400, 14500), (18000, 20000)]
    for lo, hi in telluric_bands:
        ax.axvspan(lo, hi, color='gray', alpha=0.2, zorder=0)
        ax.text((lo + hi) / 2, 0.97, r'$\oplus$', ha='center', va='top',
                fontsize=13, color='dimgray', transform=ax.get_xaxis_transform())

    ax.set_xlabel(r'Wavelength ($\AA$)')
    ax.set_ylabel(r'Flux (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)')
    plot_basename = osn.get_output_basename()
    ax.set_title('%s — GTC/EMIR NIR spectrum' % plot_basename.split('_')[0])
    ax.legend()
    if not linear:
        ax.set_yscale('log')

    fig_path = os.path.join(out_dir, plot_basename + '.png')

    if not no_plot:
        plt.show()

    fig.savefig(fig_path, dpi=150, bbox_inches='tight')
    print('  Figure saved to', fig_path)

    section('Done')
    print('  YJ/HK/combined spectra written in %s' % out_dir)
    print('  Combined : %s' % combined_spec_path)
    print('  Figure   : %s' % fig_path)


# ------------------------------------------------------------------
# Entry point
# ------------------------------------------------------------------

def main():
    args = parse_args()

    data_dir = (
        args.data_dir
        or os.environ.get('EMIR_DATA_DIR')
        or os.getcwd()
    )

    sn_path  = resolve_ob_path(args.sn_ob,  data_dir)
    std_path = resolve_ob_path(args.std_ob, data_dir)

    reduce(
        sn_path=sn_path,
        std_path=std_path,
        grisms=args.grisms,
        clean=args.clean,
        no_plot=args.no_plot,
        out_dir=args.out_dir,
        linear=args.linear,
        skip_sens_interactive=args.skip_sens_interactive,
    )


if __name__ == '__main__':
    main()
