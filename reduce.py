#!/usr/bin/env python3
"""
EMIR pipeline driver script.

Usage
-----
    python reduce.py <SN_OB> <STD_OB> [options]

Positional arguments
    SN_OB     Name or full path of the supernova observation directory.
    STD_OB    Name or full path of the standard-star observation directory.

Optional arguments
    --data-dir PATH   Root directory that contains the OB folders.
                      Defaults to $EMIR_DATA_DIR if set, otherwise the
                      current working directory.
    --grisms  G ...   One or more grisms to reduce (e.g. YJ HK).
                      Defaults to all grisms found in the QC file.
    --clean           Delete existing results and start fresh.
    --no-plot         Do not display an interactive plot window.
    --out-dir PATH    Where to write the final spectrum and figure.
                      Defaults to <SN_OB>/results/.

Example
-------
    python reduce.py OB0001 OB0002
    python reduce.py OB0001 OB0002 --data-dir ~/Desktop/EMIR_REDUX --grisms YJ
"""

import argparse
import os
import sys

# ------------------------------------------------------------------
# Make sure the src/ directory is on the path regardless of where
# the script is called from.
# ------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.join(SCRIPT_DIR, 'src')
sys.path.insert(0, SRC_DIR)
os.environ.setdefault('EMIR_PIPE', SRC_DIR)

import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import matplotlib as mpl

from Observation import Observation
from reduction_tools import get_std_magnitudes


# ------------------------------------------------------------------
# Argument parsing
# ------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description='Reduce a GTC/EMIR long-slit spectrum end-to-end.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument('sn_ob',  metavar='SN_OB',  help='Supernova OB directory (name or full path)')
    p.add_argument('std_ob', metavar='STD_OB', help='Standard-star OB directory (name or full path)')
    p.add_argument('--data-dir', default=None,
                   help='Root directory containing the OB folders '
                        '(default: $EMIR_DATA_DIR or cwd)')
    p.add_argument('--grisms', nargs='+', default=None,
                   help='Grism(s) to reduce, e.g. YJ HK (default: all in QC file)')
    p.add_argument('--clean', action='store_true',
                   help='Delete existing results and start fresh')
    p.add_argument('--no-plot', action='store_true',
                   help='Skip the interactive plot window')
    p.add_argument('--linear', action='store_true',
                   help='Use linear y-axis scale instead of log')
    p.add_argument('--out-dir', default=None,
                   help='Output directory for spectrum and figure '
                        '(default: <SN_OB>/results/)')
    return p.parse_args()


# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------

def resolve_ob_path(ob_arg, data_dir):
    """Return an absolute path for an OB argument."""
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

def reduce(sn_path, std_path, grisms, clean, no_plot, out_dir, linear=False):

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

    # Exclude acquisition images: those have GRISM='OPEN' in the QC table.
    # _get_available_grisms() returns all unique FILTER values, including J/H/K
    # from imaging frames, so we filter them out here.
    # Use the intersection of SN and std grisms — both must have the data.
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

    # ------ per-grism reduction ------
    for grism in grisms:

        section('Grism %s — SN reduction' % grism)
        if abba_done(sn_path, grism):
            print('  ABBA already exists, skipping rectification')
        else:
            osn.initialize('object', grism)
            osn.rectify_and_analyze('object', grism)
            osn.ABBA_subtract(grism)

        section('Grism %s — standard-star reduction' % grism)
        if abba_done(std_path, grism):
            print('  ABBA already exists, skipping rectification')
        else:
            ost.initialize('object', grism)
            ost.rectify_and_analyze('object', grism)
            ost.ABBA_subtract(grism)

        if sens_done(std_path, grism):
            print('  Sensitivity function already exists, skipping')
        else:
            section('Grism %s — sensitivity function' % grism)
            ost.make_sens_function(grism)

    # ------ magnitude lookup ------
    section('Fetching standard-star 2MASS magnitudes from Simbad')
    magnitudes = get_std_magnitudes(ost)

    # ------ flux calibration ------
    section('Flux calibration')
    results = {}
    for grism in grisms:
        print('  Grism %s ...' % grism)
        wave, flux = osn.get_reduced_spectrum(magnitudes, ost, grism)
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
    )


if __name__ == '__main__':
    main()
