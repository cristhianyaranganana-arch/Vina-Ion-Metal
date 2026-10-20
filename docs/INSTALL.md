# Installation

VinaMetal is a C++11 project. The only external dependency is Boost (≥ 1.65)
with the `system`, `thread`, `serialization`, `filesystem` and
`program_options` libraries.

## Linux

Tested on Ubuntu 22.04 with GCC 11.4 and Boost 1.74.

    sudo apt-get install build-essential cmake libboost-all-dev
    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build -j
    ctest --test-dir build --output-on-failure

The executables are `build/vina` and `build/vina_split`. The binary reports its
version as `AutoDock Vina v1.0.0`.

If Boost is installed in a non-standard location, pass
`-DBOOST_ROOT=/path/to/boost` to CMake.

## Windows

A local build on Windows is possible with Visual Studio, CMake and a full Boost
installation, but it is **not officially supported in this release** and no
prebuilt Windows binary is provided. The CI workflow currently builds and
releases **Linux x86-64 assets only**.

If you want to experiment on Windows, install the Boost libraries your compiler
can find and run:

    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build --config Release -j
    ctest --test-dir build -C Release --output-on-failure

## Tests

    ctest --test-dir build --output-on-failure

- `metal_gradient`: finite-difference check of the analytic metal gradients.
- `cn_detection`: runs `--score_only` on the two bundled examples and checks
  the automatically adjusted CN target.

The byte-identical no-metal regression needs an unmodified AutoDock Vina 1.2.7
binary:

    tests/test_regression_no_metal.sh build/vina /path/to/vina_original_1.2.7

## Runtime data

All metal parameters are compiled into the binary. The files under
`data/ad4_reference/` are shipped for traceability and manual inspection; they
are not read at runtime. Overrides can be supplied with `--metal_params`; see
`docs/USAGE.md`.

## Building release assets

`scripts/package_release.py` assembles the Linux platform archive and the
examples archive from a build directory:

    python3 scripts/package_release.py \
        --build-dir build --platform Linux-x64 --version 1.0

It writes `dist/VinaMetal-v1.0-Linux-x64.tar.gz`,
`dist/VinaMetal-v1.0-Linux-x64.zip` and `dist/VinaMetal-v1.0-Examples.zip`, and
prints their SHA-256 hashes. The archives contain the executables, examples,
documentation, parameter data and license files; source code is distributed
through the repository.
