#!/usr/bin/env python3
"""Build the VinaMetal release assets from a compiled build directory.

Creates, inside the output directory:

  VinaMetal-<version>-<platform>.tar.gz
  VinaMetal-<version>-<platform>.zip
  VinaMetal-<version>-Examples.zip

The platform package contains the executables, examples, parameter data,
documentation and license files. Source code is distributed through the
repository; binaries are never committed to it.

Usage:
  scripts/package_release.py --build-dir build --platform Linux-x64
  scripts/package_release.py --build-dir build --platform Windows-x64 --version 1.0

Exit status is non-zero if an executable is missing.
"""

import argparse
import hashlib
import os
import shutil
import sys
import tarfile
import tempfile
import zipfile

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

ROOT_FILES = [
    "README.md",
    "LICENSE",
    "CITATION.cff",
    "CHANGELOG.md",
    "CONTRIBUTING.md",
]

ROOT_DIRS = [
    "docs",
    "examples",
    "data",
    "scripts",
]

IGNORE_DIRS = {"__pycache__", ".git", ".github"}
IGNORE_SUFFIXES = (".pyc", ".pyo")


def iter_files(root):
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = sorted(d for d in dirnames if d not in IGNORE_DIRS)
        for name in sorted(filenames):
            if name.endswith(IGNORE_SUFFIXES):
                continue
            full = os.path.join(dirpath, name)
            yield full, os.path.relpath(full, os.path.dirname(root))


def copy_file(src, dst, executable=False):
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copyfile(src, dst)
    if executable:
        mode = os.stat(dst).st_mode
        os.chmod(dst, mode | 0o755)


def stage_platform(stage_root, build_dir, version, windows):
    binaries = ["vina.exe", "vina_split.exe"] if windows else ["vina", "vina_split"]
    for name in binaries:
        src = os.path.join(build_dir, name)
        if not os.path.isfile(src):
            raise SystemExit(
                f"error: '{src}' not found. Build the project first "
                f"(cmake --build {build_dir})."
            )
        copy_file(src, os.path.join(stage_root, name), executable=True)

    for name in ROOT_FILES:
        src = os.path.join(REPO_ROOT, name)
        if not os.path.isfile(src):
            raise SystemExit(f"error: required file '{src}' is missing")
        copy_file(src, os.path.join(stage_root, name))

    for name in ROOT_DIRS:
        src_dir = os.path.join(REPO_ROOT, name)
        for full, rel in iter_files(src_dir):
            copy_file(full, os.path.join(stage_root, rel))

    readme = os.path.join(stage_root, "examples", "README.md")
    if not os.path.isfile(readme):
        with open(readme, "w", encoding="utf-8") as fh:
            fh.write(
                "Each example directory contains receptor.pdbqt, ligand.pdbqt "
                "and config.txt. Run from inside the example directory:\n\n"
                "    ../../vina --config config.txt --out out.pdbqt\n"
            )
    return binaries


def make_tar_gz(stage_parent, top_name, output_path):
    with tarfile.open(output_path, "w:gz") as tar:
        tar.add(os.path.join(stage_parent, top_name), arcname=top_name)


def make_zip(stage_parent, top_name, output_path):
    with zipfile.ZipFile(output_path, "w", zipfile.ZIP_DEFLATED) as zf:
        base = os.path.join(stage_parent, top_name)
        for dirpath, dirnames, filenames in os.walk(base):
            dirnames.sort()
            for name in sorted(filenames):
                full = os.path.join(dirpath, name)
                arcname = os.path.relpath(full, stage_parent)
                zf.write(full, arcname)


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--build-dir", required=True,
                        help="directory containing the built vina/vina_split")
    parser.add_argument("--platform", required=True,
                        help="asset label, e.g. Linux-x64 or Windows-x64")
    parser.add_argument("--version", default="1.0",
                        help="release version used in asset names (default: 1.0)")
    parser.add_argument("--output-dir", default=os.path.join(REPO_ROOT, "dist"),
                        help="directory for the generated assets (default: dist)")
    args = parser.parse_args()

    build_dir = os.path.abspath(args.build_dir)
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    top_name = f"VinaMetal-v{args.version}-{args.platform}"
    windows = "windows" in args.platform.lower()

    with tempfile.TemporaryDirectory(prefix="vinametal-pkg-") as tmp:
        stage_root = os.path.join(tmp, top_name)
        os.makedirs(stage_root)
        stage_platform(stage_root, build_dir, args.version, windows)

        tar_path = os.path.join(output_dir, top_name + ".tar.gz")
        zip_path = os.path.join(output_dir, top_name + ".zip")
        make_tar_gz(tmp, top_name, tar_path)
        make_zip(tmp, top_name, zip_path)

    examples_name = f"VinaMetal-v{args.version}-Examples"
    with tempfile.TemporaryDirectory(prefix="vinametal-ex-") as tmp:
        stage = os.path.join(tmp, examples_name)
        shutil.copytree(os.path.join(REPO_ROOT, "examples"), os.path.join(stage, "examples"))
        make_zip(tmp, examples_name, os.path.join(output_dir, examples_name + ".zip"))

    for name in (tar_path, zip_path, os.path.join(output_dir, examples_name + ".zip")):
        print(f"{sha256(name)}  {os.path.relpath(name, output_dir)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
