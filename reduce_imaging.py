#!/usr/bin/env python3
"""
EMIR imaging pipeline entrypoint.

This is a thin wrapper around the imaging reducer implemented in
`/Users/lluisgalbany/Desktop/O2000_rex`, which was adapted to support
EMIR NIR imaging OB folders.

Usage
-----
    python reduce_imaging.py <OB_FOLDER> [options]

Example
-------
    python reduce_imaging.py /Users/lluisgalbany/Desktop/EMIR_redux/OB0002
"""

import argparse
import os
import sys


SCRIPT_DIR = "/Users/lluisgalbany/Desktop/O2000_rex"
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

from organize import process_day_folder
from python_reduction import run_python_reduction


def parse_args():
    parser = argparse.ArgumentParser(
        description="Reduce one EMIR imaging OB folder with the Python imaging pipeline."
    )
    parser.add_argument("ob_folder", help="EMIR imaging OB folder containing object/ and flat/")
    parser.add_argument(
        "--skip-check",
        action="store_true",
        help="Skip the JPG/HTML QA generation step after the reduction finishes.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    try:
        output_folder = process_day_folder(os.path.abspath(args.ob_folder))
        run_python_reduction(output_folder)
        if not args.skip_check:
            from check import run_check

            run_check(output_folder)
    except RuntimeError as exc:
        print(f"Error: {exc}")
        return 1
    except Exception as exc:
        print(f"Error: {exc}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
