#!/usr/bin/env python3
"""
PyEmir-based EMIR imaging reduction pipeline.

This script follows the main workflow described in the official PyEmir
imaging tutorial:

1. Prepare a tutorial-style workspace with `control.yaml` and `data/`.
2. Run `STARE_IMAGE` on each individual dithered exposure.
3. Run `FULL_DITHERED_IMAGE` for a quick first combination.
4. Optionally rerun refined combinations with improved offsets, masked sky
   estimation, detector-specific ad hoc sky correction, and the superflat
   "doughnut" fit.
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from astropy.io import fits


ORIGINAL_EMIR_UUID = "225fcaf2-7f6f-49cc-972a-70fd0aee8e96"
H2RG_EMIR_UUID = "443fc0d1-e09a-48cc-a0fd-02be6f399da2"

DEFAULT_RESULTS_SUBDIR = "reduced_pyemir"
DEFAULT_CALIBRATION_ROOT = Path("src/default_EMIR_files")

TUTORIAL_STRATEGIES = {
    "quick": ["ini", "v0"],
    "refine": ["ini", "v0", "v1", "v4"],
    "old": ["ini", "v0", "v1", "v4", "v5", "v6"],
    "h2rg": ["ini", "v0", "v1", "v4", "v5_h2rg"],
}


@dataclass
class ImagingSequence:
    sn_name: str
    band: str
    detector_kind: str
    date_obs: str
    files: list[Path]
    object_label: str


def normalize_object_name(obj: str) -> str:
    cleaned = str(obj).strip()
    match = re.search(
        r"\b(?:SN\s*)?(20\d{2})[_\s-]*([A-Za-z]+)(?=[^A-Za-z]|$)",
        cleaned,
        flags=re.IGNORECASE,
    )
    if not match:
        return cleaned
    year, suffix = match.groups()
    normalized = f"{year}{suffix.lower()}"
    return re.sub(
        r"\b(?:SN\s*)?(20\d{2})[_\s-]*([A-Za-z]+)(?=[^A-Za-z]|$)",
        normalized,
        cleaned,
        count=1,
        flags=re.IGNORECASE,
    )


def extract_sn_name(obj: str) -> str | None:
    match = re.search(r"(20\d{2}[A-Za-z]+)", normalize_object_name(obj))
    if not match:
        return None
    return match.group(1)


def normalize_band_name(band: str) -> str | None:
    cleaned = str(band).strip().replace(" ", "")
    return cleaned or None


def safe_name(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._+-]+", "_", str(value).strip())
    return cleaned.strip("_") or "UNKNOWN"


def final_object_basename(sequence: ImagingSequence) -> str:
    candidate = normalize_object_name(sequence.object_label or sequence.sn_name)
    if re.fullmatch(r"20\d{2}[A-Za-z]+", candidate, flags=re.IGNORECASE):
        candidate = f"SN{candidate}"
    return safe_name(candidate)


def final_product_filename(sequence: ImagingSequence) -> str:
    return f"{final_object_basename(sequence)}_{sequence.date_obs}_{sequence.band}.fits"


def read_single_fits_header(path: Path):
    try:
        with fits.open(path) as hdul:
            return hdul[0].header.copy()
    except Exception as exc:
        raise RuntimeError(f"Failed to read FITS header from {path}: {exc}") from exc


def classify_detector(header) -> str:
    detver = str(header.get("DETVER", "")).strip().upper()
    if "H2RG" in detver:
        return "h2rg"

    date_obs = str(header.get("DATE-OBS", "")).strip()
    if date_obs and date_obs[:10] >= "2023-07-01":
        return "h2rg"
    return "original"


def discover_imaging_sequences(ob_folder: Path) -> list[ImagingSequence]:
    object_dir = ob_folder / "object"
    if not object_dir.is_dir():
        raise RuntimeError(f"Missing object directory: {object_dir}")

    grouped: dict[tuple[str, str], dict[str, object]] = {}
    for path in sorted(object_dir.iterdir()):
        if path.suffix.lower() != ".fits":
            continue
        header = read_single_fits_header(path)
        obsmode = str(header.get("OBSMODE", "")).strip()
        if obsmode and obsmode != "STARE_IMAGE":
            continue

        band = normalize_band_name(header.get("FILTER", ""))
        if band is None:
            continue

        object_label = normalize_object_name(str(header.get("OBJECT", "")).strip())
        sn_name = extract_sn_name(object_label) or safe_name(object_label) or ob_folder.name
        key = (sn_name, band)

        info = grouped.setdefault(
            key,
            {
                "sn_name": sn_name,
                "band": band,
                "detector_kind": classify_detector(header),
                "date_obs": str(header.get("DATE-OBS", "")).split("T")[0].replace("-", ""),
                "files": [],
                "object_label": object_label or sn_name,
            },
        )
        info["files"].append(path)

    sequences = []
    for key in sorted(grouped):
        info = grouped[key]
        files = sorted(
            info["files"],
            key=lambda p: str(read_single_fits_header(p).get("DATE-OBS", p.name)),
        )
        sequences.append(
            ImagingSequence(
                sn_name=info["sn_name"],
                band=info["band"],
                detector_kind=info["detector_kind"],
                date_obs=info["date_obs"] or ob_folder.name,
                files=files,
                object_label=info["object_label"],
            )
        )
    if not sequences:
        raise RuntimeError(f"No STARE_IMAGE science frames found in {object_dir}")
    return sequences


def ensure_symlink_or_copy(src: Path, dst: Path):
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() or dst.is_symlink():
        return
    try:
        dst.symlink_to(src.resolve())
    except OSError:
        shutil.copy2(src, dst)


def find_calibration_file(calibration_root: Path, filename: str, optional: bool = False) -> Path | None:
    candidate = calibration_root / filename
    if candidate.exists():
        return candidate
    if optional:
        return None
    raise RuntimeError(f"Missing calibration file {candidate}")


def build_control_yaml(flat_original: str, flat_h2rg: str) -> str:
    rect_entries = [
        ("J", "J", "rect_wpoly_MOSlibrary_grism_J_filter_J.json", 11),
        ("H", "H", "rect_wpoly_MOSlibrary_grism_H_filter_H.json", 12),
        ("K", "Ksp", "rect_wpoly_MOSlibrary_grism_K_filter_Ksp.json", 13),
        ("LR", "YJ", "rect_wpoly_MOSlibrary_grism_LR_filter_YJ.json", 14),
        ("LR", "HK", "rect_wpoly_MOSlibrary_grism_LR_filter_HK.json", 15),
    ]

    def block_lines(config_uuid: str, flat_name: str) -> list[str]:
        lines = [f"    '{config_uuid}':"]
        lines.append("      - {id: 2, type: 'MasterBadPixelMask', tags: {}, content: 'master_bpm.fits'}")
        lines.append("      - {id: 3, type: 'MasterDark', tags: {}, content: 'master_dark_zeros.fits'}")
        lines.append(
            f"      - {{id: 4, type: 'MasterIntensityFlat', tags: {{}}, content: '{flat_name}'}}"
        )
        lines.append(
            f"      - {{id: 5, type: 'MasterSpectralFlat', tags: {{}}, content: '{flat_name}'}}"
        )
        for grism, filt, content, idx in rect_entries:
            lines.append(
                "      - "
                + "{id: %d, type: 'MasterRectWave', tags: {grism: %s, filter: %s}, content: '%s'}"
                % (idx, grism, filt, content)
            )
        return lines

    lines = [
        "version: 1",
        "products:",
        "  EMIR:",
    ]
    lines.extend(block_lines(ORIGINAL_EMIR_UUID, flat_original))
    lines.extend(block_lines(H2RG_EMIR_UUID, flat_h2rg))
    return "\n".join(lines) + "\n"


def build_stare_blocks(file_names: list[str], repeat: int, reprojection: str, enabled: bool) -> list[str]:
    blocks = []
    for start in range(0, len(file_names), repeat):
        group = file_names[start : start + repeat]
        if len(group) != repeat:
            raise RuntimeError(
                f"Found trailing group of {len(group)} file(s), but repeat={repeat}."
            )
        id_label = Path(group[0]).name[:10]
        lines = [
            f"id: _{id_label}",
            "instrument: EMIR",
            "mode: STARE_IMAGE",
            "frames:",
        ]
        lines.extend([f" - {name}" for name in group])
        lines.extend(
            [
                "requirements:",
                f"  reprojection_method: {reprojection}",
                f"enabled: {'True' if enabled else 'False'}",
            ]
        )
        blocks.append("\n".join(lines))
    return blocks


def combined_children(file_names: list[str], repeat: int) -> list[str]:
    children = []
    for start in range(0, len(file_names), repeat):
        group = file_names[start : start + repeat]
        if len(group) != repeat:
            raise RuntimeError(
                f"Found trailing group of {len(group)} file(s), but repeat={repeat}."
            )
        children.append(f"_{Path(group[0]).name[:10]}")
    return children


def stage_requirements(stage: str) -> dict[str, object]:
    if stage == "v0":
        return {"iterations": 0, "sky_images": 0, "refine_offsets": False}
    if stage == "v1":
        return {"iterations": 0, "sky_images": 0, "refine_offsets": True}
    if stage == "v2":
        return {
            "iterations": 0,
            "sky_images": 0,
            "offsets": "user_offsets.txt",
            "refine_offsets": False,
        }
    if stage == "v3":
        return {
            "iterations": 0,
            "sky_images": 0,
            "offsets": "user_offsets.txt",
            "refine_offsets": True,
        }
    if stage == "v4":
        return {"iterations": 1, "sky_images": 6, "refine_offsets": True}
    if stage == "v5":
        return {
            "iterations": 1,
            "sky_images": 6,
            "refine_offsets": True,
            "nside_adhoc_sky_correction": 10,
        }
    if stage == "v5_h2rg":
        return {
            "iterations": 1,
            "sky_images": 6,
            "refine_offsets": True,
            "adhoc_sky_correction_h2rg": True,
        }
    if stage == "v6":
        return {
            "iterations": 1,
            "sky_images": 6,
            "refine_offsets": True,
            "nside_adhoc_sky_correction": 10,
            "fit_doughnut": True,
        }
    raise ValueError(f"Unsupported stage: {stage}")


def build_combined_block(file_names: list[str], repeat: int, stage: str) -> str:
    lines = [
        f"id: _combined_{stage}",
        "instrument: EMIR",
        "mode: FULL_DITHERED_IMAGE",
        "children:",
    ]
    lines.extend([f" - {child}" for child in combined_children(file_names, repeat)])
    lines.append("requirements:")
    for key, value in stage_requirements(stage).items():
        if isinstance(value, bool):
            rendered = "True" if value else "False"
        else:
            rendered = str(value)
        lines.append(f"  {key}: {rendered}")
    lines.append("enabled: True")
    return "\n".join(lines)


def write_yaml(path: Path, blocks: list[str]):
    path.write_text("\n---\n".join(blocks) + "\n")


def determine_stage_chain(detector_kind: str, strategy: str, use_user_offsets: bool) -> list[str]:
    if strategy == "auto":
        base = list(TUTORIAL_STRATEGIES["h2rg" if detector_kind == "h2rg" else "old"])
    else:
        if strategy not in TUTORIAL_STRATEGIES:
            raise RuntimeError(f"Unknown strategy {strategy}")
        base = list(TUTORIAL_STRATEGIES[strategy])

    if use_user_offsets:
        if "v1" in base:
            base.remove("v1")
        insert_at = base.index("v4") if "v4" in base else len(base)
        base[insert_at:insert_at] = ["v2", "v3"]
    return base


def build_numina_command(python_executable: str, yaml_name: str) -> list[str]:
    entrypoint = shutil.which("numina")
    if entrypoint is not None:
        return [entrypoint, "run", yaml_name, "--link-files", "-r", "control.yaml"]
    return [
        python_executable,
        "-m",
        "numina.user.cli",
        "run",
        yaml_name,
        "--link-files",
        "-r",
        "control.yaml",
    ]


def run_command(cmd: list[str], cwd: Path):
    print("\nExecuting:")
    print("  " + " ".join(cmd))
    result = subprocess.run(cmd, cwd=str(cwd))
    if result.returncode != 0:
        raise RuntimeError(f"Command failed with exit code {result.returncode}: {' '.join(cmd)}")


def result_dir_for_stage(stage: str) -> str:
    return f"obsid_combined_{stage}_results"


def work_dir_for_stage(stage: str) -> str:
    return f"obsid_combined_{stage}_work"


def result_image_for_stage(workspace: Path, stage: str) -> Path:
    return workspace / result_dir_for_stage(stage) / "reduced_image.fits"


def copy_user_offsets(user_offsets: Path, data_dir: Path, expected_count: int):
    lines = [line.strip() for line in user_offsets.read_text().splitlines() if line.strip()]
    if len(lines) != expected_count:
        raise RuntimeError(
            f"{user_offsets} has {len(lines)} non-empty line(s), but {expected_count} "
            "input image(s) are being combined."
        )
    shutil.copy2(user_offsets, data_dir / "user_offsets.txt")


def prepare_workspace(
    workspace: Path,
    sequence: ImagingSequence,
    calibration_root: Path,
    reprojection: str,
    repeat: int,
    use_user_offsets: bool,
    user_offsets: Path | None,
):
    data_dir = workspace / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    for src in sequence.files:
        ensure_symlink_or_copy(src, data_dir / src.name)

    calib_names = [
        "master_bpm.fits",
        "master_dark_zeros.fits",
        "master_flat_ones.fits",
        "master_flat_spec.fits",
        "rect_wpoly_MOSlibrary_grism_H_filter_H.json",
        "rect_wpoly_MOSlibrary_grism_J_filter_J.json",
        "rect_wpoly_MOSlibrary_grism_K_filter_Ksp.json",
        "rect_wpoly_MOSlibrary_grism_LR_filter_HK.json",
        "rect_wpoly_MOSlibrary_grism_LR_filter_YJ.json",
    ]
    for name in calib_names:
        ensure_symlink_or_copy(find_calibration_file(calibration_root, name), data_dir / name)

    flat_h2rg_path = find_calibration_file(
        calibration_root, "master_flat_spec_H2RG.fits", optional=True
    )
    if flat_h2rg_path is not None:
        ensure_symlink_or_copy(flat_h2rg_path, data_dir / flat_h2rg_path.name)
        flat_h2rg_name = flat_h2rg_path.name
    else:
        flat_h2rg_name = "master_flat_spec.fits"

    (workspace / "control.yaml").write_text(
        build_control_yaml(
            flat_original="master_flat_spec.fits",
            flat_h2rg=flat_h2rg_name,
        )
    )

    file_names = [path.name for path in sequence.files]
    (data_dir / "list_images.txt").write_text("\n".join(file_names) + "\n")

    if use_user_offsets:
        if user_offsets is None:
            raise RuntimeError("--use-user-offsets requires --user-offsets")
        copy_user_offsets(user_offsets, data_dir, expected_count=len(file_names) // repeat)

    write_yaml(
        workspace / "dithered_ini.yaml",
        build_stare_blocks(file_names, repeat=repeat, reprojection=reprojection, enabled=True),
    )

    stare_disabled = build_stare_blocks(
        file_names, repeat=repeat, reprojection=reprojection, enabled=False
    )
    for stage in ["v0", "v1", "v2", "v3", "v4", "v5", "v5_h2rg", "v6"]:
        blocks = list(stare_disabled)
        blocks.append(build_combined_block(file_names, repeat=repeat, stage=stage))
        write_yaml(workspace / f"dithered_{stage}.yaml", blocks)

    manifest = {
        "sn_name": sequence.sn_name,
        "object_label": sequence.object_label,
        "band": sequence.band,
        "detector_kind": sequence.detector_kind,
        "date_obs": sequence.date_obs,
        "n_images": len(sequence.files),
        "repeat": repeat,
        "reprojection": reprojection,
        "calibration_root": str(calibration_root.resolve()),
        "workspace": str(workspace.resolve()),
    }
    (workspace / "pipeline_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


def run_tutorial_pipeline(
    workspace: Path,
    detector_kind: str,
    strategy: str,
    use_user_offsets: bool,
    python_executable: str,
    prepare_only: bool,
):
    stage_chain = determine_stage_chain(
        detector_kind=detector_kind,
        strategy=strategy,
        use_user_offsets=use_user_offsets,
    )
    if prepare_only:
        return stage_chain[-1]

    run_command(build_numina_command(python_executable, "dithered_ini.yaml"), cwd=workspace)
    for stage in stage_chain:
        if stage == "ini":
            continue
        run_command(build_numina_command(python_executable, f"dithered_{stage}.yaml"), cwd=workspace)
    return stage_chain[-1]


def link_or_copy_product(source: Path, destination: Path):
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() or destination.is_symlink():
        destination.unlink()
    try:
        destination.symlink_to(source.resolve())
    except OSError:
        shutil.copy2(source, destination)


def export_aligned_frames(workspace: Path, final_stage: str) -> Path:
    source_dir = workspace / work_dir_for_stage(final_stage)
    if not source_dir.exists():
        raise RuntimeError(f"Aligned-image work directory not found: {source_dir}")

    aligned_paths = sorted(source_dir.glob("*_r.fits"))
    if not aligned_paths:
        raise RuntimeError(f"No aligned '*_r.fits' files found in {source_dir}")

    output_dir = workspace.parent / "aligned_frames"
    output_dir.mkdir(parents=True, exist_ok=True)
    for path in output_dir.glob("*.fits"):
        path.unlink()

    for src in aligned_paths:
        link_or_copy_product(src, output_dir / src.name)

    return output_dir


def cleanup_intermediate_products(workspace: Path, final_stage: str, keep_aligned: bool) -> dict[str, str | int | None]:
    aligned_dir = None
    if keep_aligned:
        aligned_dir = export_aligned_frames(workspace, final_stage)

    final_result_dir = workspace / result_dir_for_stage(final_stage)
    if not final_result_dir.exists():
        raise RuntimeError(f"Final result directory not found: {final_result_dir}")

    removed_dirs = 0
    for path in sorted(workspace.iterdir()):
        if path == final_result_dir:
            continue
        if aligned_dir is not None and path == aligned_dir:
            continue
        if path.name in {
            "data",
            "control.yaml",
            "pipeline_manifest.json",
            "final_reduced_image.fits",
        }:
            continue
        if path.name.startswith("dithered_"):
            continue

        if path.is_dir() and path.name.startswith("obsid_"):
            shutil.rmtree(path)
            removed_dirs += 1

    return {
        "aligned_dir": str(aligned_dir) if aligned_dir is not None else None,
        "kept_result_dir": str(final_result_dir),
        "removed_dirs": removed_dirs,
    }


def link_final_product(workspace: Path, final_stage: str, sequence: ImagingSequence):
    source = result_image_for_stage(workspace, final_stage)
    if not source.exists():
        raise RuntimeError(f"Final reduced image not found: {source}")

    final_link = workspace / "final_reduced_image.fits"
    link_or_copy_product(source, final_link)

    named_product = workspace.parent / final_product_filename(sequence)
    link_or_copy_product(source, named_product)
    return named_product


def parse_args():
    parser = argparse.ArgumentParser(
        description="Reduce EMIR imaging data with a tutorial-style PyEmir pipeline."
    )
    parser.add_argument("ob_folder", help="OB folder containing object/ FITS images.")
    parser.add_argument(
        "--results-subdir",
        default=DEFAULT_RESULTS_SUBDIR,
        help=f"Subdirectory created inside the OB folder for the new pipeline (default: {DEFAULT_RESULTS_SUBDIR}).",
    )
    parser.add_argument(
        "--strategy",
        default="auto",
        choices=["auto", "quick", "refine", "old", "h2rg"],
        help="Tutorial stage chain to run. 'auto' chooses old/H2RG from the data.",
    )
    parser.add_argument(
        "--repeat",
        type=int,
        default=1,
        help="Number of consecutive exposures per dither position.",
    )
    parser.add_argument(
        "--reprojection",
        default="interp",
        choices=["interp", "adaptive", "exact", "none"],
        help="Reprojection method passed to the STARE_IMAGE stage.",
    )
    parser.add_argument(
        "--calibration-root",
        default=str(DEFAULT_CALIBRATION_ROOT),
        help="Folder containing the tutorial-style calibration files.",
    )
    parser.add_argument(
        "--prepare-only",
        action="store_true",
        help="Only prepare the PyEmir workspace and YAML files; do not run numina.",
    )
    parser.add_argument(
        "--python-executable",
        default=sys.executable,
        help="Python executable used when falling back to 'python -m numina.user.cli'.",
    )
    parser.add_argument(
        "--use-user-offsets",
        action="store_true",
        help="Insert the tutorial v2/v3 user-offset stages before the masked sky stages.",
    )
    parser.add_argument(
        "--user-offsets",
        help="ASCII file with one x y offset pair per dither position, as in the tutorial.",
    )
    parser.add_argument(
        "--cleanup-intermediate",
        choices=["none", "keep-aligned"],
        default="none",
        help="After a successful run, remove heavy intermediate obsid_* directories. "
        "'keep-aligned' preserves the final aligned '*_r.fits' frames in an aligned_frames/ folder.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    ob_folder = Path(args.ob_folder).resolve()
    calibration_root = Path(args.calibration_root).resolve()
    user_offsets = Path(args.user_offsets).resolve() if args.user_offsets else None

    if args.repeat < 1:
        raise RuntimeError("--repeat must be >= 1")

    sequences = discover_imaging_sequences(ob_folder)
    root_output = ob_folder / args.results_subdir
    root_output.mkdir(parents=True, exist_ok=True)

    summary = []
    for sequence in sequences:
        workspace = root_output / sequence.sn_name / sequence.band / "pyemir_workspace"
        print("\n" + "=" * 60)
        print(f"Preparing {sequence.sn_name} {sequence.band}")
        print("=" * 60)
        print(f"  detector : {sequence.detector_kind}")
        print(f"  n_images : {len(sequence.files)}")
        print(f"  workspace: {workspace}")

        prepare_workspace(
            workspace=workspace,
            sequence=sequence,
            calibration_root=calibration_root,
            reprojection=args.reprojection,
            repeat=args.repeat,
            use_user_offsets=args.use_user_offsets,
            user_offsets=user_offsets,
        )
        final_stage = run_tutorial_pipeline(
            workspace=workspace,
            detector_kind=sequence.detector_kind,
            strategy=args.strategy,
            use_user_offsets=args.use_user_offsets,
            python_executable=args.python_executable,
            prepare_only=args.prepare_only,
        )

        named_product = None
        cleanup_info = None
        if not args.prepare_only:
            named_product = link_final_product(workspace, final_stage, sequence)
            if args.cleanup_intermediate == "keep-aligned":
                cleanup_info = cleanup_intermediate_products(
                    workspace=workspace,
                    final_stage=final_stage,
                    keep_aligned=True,
                )

        summary.append(
            {
                "sn_name": sequence.sn_name,
                "band": sequence.band,
                "detector_kind": sequence.detector_kind,
                "workspace": str(workspace),
                "final_stage": final_stage,
                "named_product": str(named_product) if named_product is not None else None,
                "cleanup": cleanup_info,
            }
        )

    summary_path = root_output / "summary_reduce_imaging_2.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print("\nSummary written to", summary_path)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RuntimeError as exc:
        print(f"Error: {exc}")
        raise SystemExit(1)
