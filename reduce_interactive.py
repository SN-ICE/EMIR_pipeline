#!/usr/bin/env python3
"""
EMIR interactive pipeline driver script.

Mirrors reduce.py but adds two interactive matplotlib windows per grism:
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
    python reduce_interactive.py <SN_OB> <STD_OB> [options]

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
    python reduce_interactive.py OB0001 OB0002
    python reduce_interactive.py OB0001 OB0002 --data-dir ~/Desktop/EMIR_REDUX --grisms YJ
"""

import argparse
import os
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
from matplotlib import pyplot as plt
import matplotlib as mpl

from Observation import Observation
from reduction_tools import get_std_magnitudes
from interactive_reduction_tools import (
    interactive_SN_position_finder,
    interactive_sens_function_fit,
    get_parameters,
)


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
            interactive_sens_function_fit(ost, grism)
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
        interactive_SN_position_finder(osn, grism)
        sn_kwargs[grism] = get_parameters(osn, grism)
        print('  SN position saved:', sn_kwargs[grism].get('SN_position'))

    # ------ magnitude lookup ------
    section('Fetching standard-star 2MASS magnitudes from Simbad')
    magnitudes = get_std_magnitudes(ost)

    # ------ flux calibration ------
    section('Flux calibration')
    results = {}
    for grism in grisms:
        print('  Grism %s ...' % grism)
        wave, flux = osn.get_reduced_spectrum(magnitudes, ost, grism,
                                              **sn_kwargs.get(grism, {}))
        results[grism] = (wave, flux)

    # ------ save spectrum ------
    if out_dir is None:
        out_dir = os.path.join(sn_path, 'results')
    os.makedirs(out_dir, exist_ok=True)

    sn_name = os.path.basename(sn_path.rstrip('/'))
    spec_path = os.path.join(out_dir, '%s_spectrum.txt' % sn_name)
    osn.save_spectrum(spec_path)
    print('\n  Spectrum saved to', spec_path)

    # ------ plot ------
    section('Plotting')
    colors = {'YJ': 'royalblue', 'HK': 'firebrick'}
    default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    mpl.rcParams.update({'font.size': 13, 'font.family': 'serif'})
    fig, ax = plt.subplots(figsize=(14, 5))

    for i, (grism, (wave, flux)) in enumerate(results.items()):
        color = colors.get(grism, default_colors[i % len(default_colors)])
        ax.plot(wave, flux, color=color, lw=0.8, label=grism)

    # telluric bands
    telluric_bands = [(13400, 14500), (18000, 20000)]
    for lo, hi in telluric_bands:
        ax.axvspan(lo, hi, color='gray', alpha=0.2, zorder=0)
        ax.text((lo + hi) / 2, 0.97, r'$\oplus$', ha='center', va='top',
                fontsize=13, color='dimgray', transform=ax.get_xaxis_transform())

    ax.set_xlabel(r'Wavelength ($\AA$)')
    ax.set_ylabel(r'Flux (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)')
    ax.set_title('%s — GTC/EMIR NIR spectrum' % sn_name)
    ax.legend()
    if not linear:
        ax.set_yscale('log')

    fig_path = os.path.join(out_dir, '%s_spectrum.png' % sn_name)
    fig.savefig(fig_path, dpi=150, bbox_inches='tight')
    print('  Figure saved to', fig_path)

    if not no_plot:
        plt.show()

    section('Done')
    print('  Spectrum : %s' % spec_path)
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
