#!/usr/bin/env python3
"""
Simplified WABBIT unit testing framework.

This script runs all WABBIT tests defined in the `tests` list below. It supports
two types of tests with different execution models:

Test types:
- "wabbit-internal": Internal WABBIT tests that run directly in TESTING/{root_folder}/
  using a dedicated log file per test: {name}_{dim}D_{wavelet}.log
  These tests verify wavelet operations like refine/coarsen and decomposition.
  
- "simulation": Full simulation tests that run in TESTING/{root_folder}/{name}/tmp/
  with all output captured in log.txt. Input files are copied from the test directory
  to tmp/, commands are executed there, and HDF5 output files are compared against
  reference files in the test directory. The tmp/ directory is cleaned up after each test.

Directory structure:
  TESTING/
    wavelets/          - Wavelet-specific tests
    conv/             - Convection tests
    acm/              - ACM (Adaptive Cartesian Mesh) tests
    insects/          - Insect flight simulation tests

The test definitions use Python list comprehensions for wabbit-internal and adaptive
tests to avoid repetition, with explicit entries for tests requiring unique configurations.
"""

import os
import sys
import argparse
import subprocess
import shutil
import glob
import time
import contextlib
import fnmatch
import itertools

try:
    import wabbit_tools
except ImportError:
    print("ERROR: wabbit_tools module not found")
    sys.exit(881)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TEST DEFINITIONS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Each test is a dictionary with the following fields:
# - name (str): Test name, used for logging and directory structure
# - type (str): Either "wabbit-internal" or "simulation"
# - root_folder (str): Subdirectory under TESTING/ where test files are located
#
# For wabbit-internal tests:
# - wavelet (str): Wavelet type (e.g., "CDF20", "CDF40")
# - dim (int): Dimensionality (2 or 3)
# - commands (list): List of command templates to execute
#
# For simulation tests:
# - wavelet (str, optional): Wavelet type for adaptive tests
# - input_files (list): Files to copy from test directory to tmp/ before execution
# - commands (list): List of command templates to execute in tmp/
#
# Command templates support {placeholders} that are replaced at runtime:
# - {mpi_command}: MPI execution command (e.g., "nice mpirun -n 4")
# - {run_dir}: Path to WABBIT directory
# - {memory}: Memory allocation flag
# - {wavelet}: Wavelet type from test definition
# - {name}: Test name from test definition
# - {dim}: Dimensionality from test definition
#
# String separators ("---") are used to group tests in output.


tests = [
    # ===== WABBIT-INTERNAL TESTS =====
    # These tests verify internal WABBIT wavelet operations.
    # They run quickly and don't require reference HDF5 files.
    "\n--- WABBIT internal test (wavelet transform)---",

    # Equidistant refine-coarsening tests. After refining and subsequent coarsening, the original data should be recovered. 
    # If it is not, the wavelet coefficients have a serious problem.
    # 2D
    *[{"name": "equi_refineCoarsen_FWT_IWT", "type": "wabbit-internal", "root_folder": "wavelets", "wavelet": wavelet, "dim": 2,
       "groups": ["wabbit_internal", "refine_coarsen"],
       "commands": [ "{mpi_command} {run_dir}/wabbit-post --refine-coarsen-test --wavelet={wavelet} --memory={memory} --dim={dim}",
                     "{mpi_command} {run_dir}/wabbit-post --wavelet-decomposition-unit-test --wavelet={wavelet} --memory={memory} --dim={dim}" ]}
      for wavelet in ["CDF20", "CDF22", "CDF24", "CDF26", "CDF28", "CDF40", "CDF42", "CDF44", "CDF46", "CDF60", "CDF62", "CDF64", "CDF66", "CDF80", "CDF82"]],

    # 3D
    *[{"name": "equi_refineCoarsen_FWT_IWT", "type": "wabbit-internal", "root_folder": "wavelets", "wavelet": wavelet, "dim": 3,
       "groups": ["wabbit_internal", "refine_coarsen"] + (["short"] if wavelet in ["CDF20", "CDF22", "CDF40", "CDF42", "CDF60", "CDF62"] else []),
       "commands": [ "{mpi_command} {run_dir}/wabbit-post --refine-coarsen-test --wavelet={wavelet} --memory={memory} --dim={dim}",
                     "{mpi_command} {run_dir}/wabbit-post --wavelet-decomposition-unit-test --wavelet={wavelet} --memory={memory} --dim={dim}" ]}
      for wavelet in ["CDF20", "CDF22", "CDF40", "CDF42", "CDF44", "CDF60", "CDF62"]],

    # ghost_nodes tests
    # 2D
    *[{"name": "ghost_nodes", "type": "wabbit-internal", "root_folder": "wavelets", "wavelet": wavelet, "dim": 2,
       "groups": ["wabbit_internal", "ghost_nodes"],
       "commands": ["{mpi_command} {run_dir}/wabbit-post --ghost-nodes-test --wavelet={wavelet} --memory={memory} --dim={dim} --Jmax=5"]}
      for wavelet in ["CDF20", "CDF40", "CDF60", "CDF80"]],
     # 3D
    *[{"name": "ghost_nodes", "type": "wabbit-internal", "root_folder": "wavelets", "wavelet": wavelet, "dim": 3,
       "groups": ["wabbit_internal", "ghost_nodes", "short"],
       "commands": ["{mpi_command} {run_dir}/wabbit-post --ghost-nodes-test --wavelet={wavelet} --memory={memory} --dim={dim}"]}
      for wavelet in ["CDF20", "CDF40", "CDF60"]],

    # invertibility tests (2D only)
    *[{"name": "invertibility", "type": "wabbit-internal", "root_folder": "wavelets", "wavelet": wavelet, "dim": 2,
       "groups": ["wabbit_internal", "invertibility"] + (["short"] if wavelet in ["CDF20", "CDF22", "CDF40", "CDF42", "CDF60", "CDF62"] else []),
       "commands": ["{mpi_command} {run_dir}/wabbit-post --wavelet-decomposition-invertibility-test --wavelet={wavelet} --memory={memory} --dim={dim} --Jmax=7"]}
      for wavelet in ["CDF20", "CDF22", "CDF24", "CDF26", "CDF28", "CDF40", "CDF42", "CDF44", "CDF46", "CDF60", "CDF62", "CDF64", "CDF66", "CDF80", "CDF82"]],

    # ===== SIMULATION TESTS =====
    # Full simulation tests that produce HDF5 output files.
    # These tests copy input files to a tmp/ directory, run the simulation,
    # then compare the HDF5 output against reference files.
    #"---simulation---",

    # Adaptive mesh refinement tests
    # These test the adaptive refine/coarsen functionality by starting with
    # a coarse mesh (vor_000020000000.h5), refining it, then coarsening it back.
    # Wavelet is parameterized via {wavelet} in the commands.
    # Same idea as equidistant refine-coarsening test above, but on a non-equidistant grid.
    *[{"name": f"adaptive_{wavelet}", "type": "simulation", "root_folder": "wavelets", "wavelet": wavelet, "dim": 2,
       "groups": ["adaptive", "refine_coarsen", "short"],
       "input_files": ["../vor_000020000000.h5"],
       "commands": [ "{mpi_command} {run_dir}/wabbit-post --refine-everywhere vor_000020000000.h5 vor_00100.h5 --wavelet={wavelet} --time=1.0",
                     "{mpi_command} {run_dir}/wabbit-post --coarsen-everywhere vor_00100.h5 vor_00200.h5 --wavelet={wavelet} --time=2.0",
                     # remove input files after completion to not check these
                     # for one wavelet (unlifted ones), we also test the refined solution (should be the same for all wavelets with the same interpolation order X of CDFXY)
                     f"rm vor_000020000000.h5{' vor_00100.h5' if wavelet not in ['CDF20', 'CDF40', 'CDF60'] else ''}" ]}
      for wavelet in ["CDF20", "CDF22", "CDF40", "CDF42", "CDF44", "CDF60", "CDF62"]],

    "\n---Blob convection tests (convection module)---",
    # Blob convection tests
    # Test equispaced and adaptive mesh simulations with blob convection.
    # These verify basic simulation setup and execution.
    {"name": "blob_equi_2D_CDF40", "type": "simulation", "root_folder": "conv", "wavelet": 40, "dim": 2,
        "groups": ["blob", "equi"], "input_files": ["*.ini"], "commands": ["{mpi_command} {run_dir}/wabbit blob-convection.ini --memory={memory}"]},
    {"name": "blob_equi_3D_CDF20", "type": "simulation", "root_folder": "conv", "wavelet": 20, "dim": 3,
        "groups": ["blob", "equi", "short"], "input_files": ["*.ini"], "commands": ["{mpi_command} {run_dir}/wabbit blob-convection.ini --memory={memory}"]},
    {"name": "blob_equi_3D_CDF40", "type": "simulation", "root_folder": "conv", "wavelet": 40, "dim": 3,
        "groups": ["blob", "equi", "short"], "input_files": ["*.ini"], "commands": ["{mpi_command} {run_dir}/wabbit blob-convection.ini --memory={memory}"]},
    {"name": "blob_equi_avg_2D_CDF40", "type": "simulation", "root_folder": "conv", "wavelet": 40, "dim": 2,
        "groups": ["blob", "equi"], "input_files": ["*.ini"], "commands": ["{mpi_command} {run_dir}/wabbit blob-convection.ini --memory={memory}"]},

    {"name": "blob_adaptive_2D_CDF20", "type": "simulation", "root_folder": "conv", "wavelet": 20, "dim": 2,
        "groups": ["blob", "adaptive"], "input_files": ["*.ini"], "commands": ["{mpi_command} {run_dir}/wabbit blob-convection.ini --memory={memory}"]},
    {"name": "blob_adaptive_2D_CDF22", "type": "simulation", "root_folder": "conv", "wavelet": 22, "dim": 2,
        "groups": ["blob", "adaptive"], "input_files": ["*.ini"], "commands": ["{mpi_command} {run_dir}/wabbit blob-convection.ini --memory={memory}"]},
    {"name": "blob_adaptive_2D_CDF40", "type": "simulation", "root_folder": "conv", "wavelet": 40, "dim": 2,
        "groups": ["blob", "adaptive", "short"], "input_files": ["*.ini"], "commands": ["{mpi_command} {run_dir}/wabbit blob-convection.ini --memory={memory}"]},
    {"name": "blob_adaptive_2D_CDF42", "type": "simulation", "root_folder": "conv", "wavelet": 42, "dim": 2,
        "groups": ["blob", "adaptive"], "input_files": ["*.ini"], "commands": ["{mpi_command} {run_dir}/wabbit blob-convection.ini --memory={memory}"]},

    {"name": "blob_adaptive_3D_CDF40", "type": "simulation", "root_folder": "conv", "wavelet": 40, "dim": 3,
            "groups": ["blob", "adaptive", "short"], "input_files": ["*.ini"], "commands": ["{mpi_command} {run_dir}/wabbit blob-convection.ini --memory={memory}"]},
    {"name": "blob_adaptive_3D_CDF44", "type": "simulation", "root_folder": "conv", "wavelet": 44, "dim": 3,
        "groups": ["blob", "adaptive", "short"], "input_files": ["*.ini"], "commands": ["{mpi_command} {run_dir}/wabbit blob-convection.ini --memory={memory}"]},

    "\n---ACM tests, with and without penalization---",
    # ACM (Adaptive Cartesian Mesh) cylinder tests
    # Test cylindrical geometry simulations with different wavelet types.
    # acm tests
    {"name": "flowPastCylinder_FD4_CDF40", "type": "simulation", "root_folder": "acm/flowPastCylinder", "wavelet": 40, "dim": 2,
        "groups": ["acm", "adaptive", "cylinder", "short"], "input_files": ["PARAMS_flowPastCylinder.ini"],
        "commands": ["{mpi_command} {run_dir}/wabbit PARAMS_flowPastCylinder.ini --memory={memory}"]},
    {"name": "flowPastCylinder_FD4_CDF44", "type": "simulation", "root_folder": "acm/flowPastCylinder", "wavelet": 44, "dim": 2,
        "groups": ["acm", "adaptive", "cylinder", "short"], "input_files": ["PARAMS_flowPastCylinder.ini"],
        "commands": ["{mpi_command} {run_dir}/wabbit PARAMS_flowPastCylinder.ini --memory={memory}"]},
    {"name": "flowPastCylinder_threshold_norm", "type": "simulation", "root_folder": "acm/flowPastCylinder", "wavelet": 44, "dim": 2,
        "groups": ["acm", "adaptive", "cylinder"], "input_files": ["PARAMS_flowPastCylinder.ini"],
        "commands": ["{mpi_command} {run_dir}/wabbit PARAMS_flowPastCylinder.ini --memory={memory}"]},
    {"name": "flowPastCylinder_significant_refinement", "type": "simulation", "root_folder": "acm/flowPastCylinder", "wavelet": 44, "dim": 2,
        "groups": ["acm", "adaptive", "cylinder"], "input_files": ["PARAMS_flowPastCylinder.ini"],
        "commands": ["{mpi_command} {run_dir}/wabbit PARAMS_flowPastCylinder.ini --memory={memory}"]},
    {"name": "flowPastCylinder_nonquadratic_domain", "type": "simulation", "root_folder": "acm/flowPastCylinder", "wavelet": 40, "dim": 2,
        "groups": ["acm", "adaptive", "cylinder"], "input_files": ["PARAMS_flowPastCylinder.ini"],
        "commands": ["{mpi_command} {run_dir}/wabbit PARAMS_flowPastCylinder.ini --memory={memory}"]},

    # Three vortices tests
    # Test vortex interaction simulations with different finite difference orders (FD2, FD4, FD6)
    # and both equispaced and adaptive meshes.
    # 3vortices tests
    {"name": "3vorticesEquiFD2_CDF20", "type": "simulation", "root_folder": "acm/3vortices", "wavelet": 20, "dim": 2,
        "groups": ["acm", "3vortices", "equi"], "input_files": ["*_000010000000.h5", "*.ini"],
        "commands": ["{mpi_command} {run_dir}/wabbit PARAMS_3vortices.ini --memory={memory}"]},
    {"name": "3vorticesEquiFD4_CDF40", "type": "simulation", "root_folder": "acm/3vortices", "wavelet": 40, "dim": 2,
        "groups": ["acm", "3vortices", "equi"], "input_files": ["*_000010000000.h5", "*.ini"],
        "commands": ["{mpi_command} {run_dir}/wabbit PARAMS_3vortices.ini --memory={memory}"]},
    {"name": "3vorticesEquiFD6_CDF60", "type": "simulation", "root_folder": "acm/3vortices", "wavelet": 60, "dim": 2,
        "groups": ["acm", "3vortices", "equi"], "input_files": ["*_000010000000.h5", "*.ini"],
        "commands": ["{mpi_command} {run_dir}/wabbit PARAMS_3vortices.ini --memory={memory}"]},

    {"name": "3vorticesAdaptFD2_CDF20", "type": "simulation", "root_folder": "acm/3vortices", "wavelet": 20, "dim": 2,
        "groups": ["acm", "3vortices", "adaptive"], "input_files": ["*_000010000000.h5", "*.ini"],
        "commands": ["{mpi_command} {run_dir}/wabbit PARAMS_3vortices.ini --memory={memory}"]},
    {"name": "3vorticesAdaptFD2_CDF22", "type": "simulation", "root_folder": "acm/3vortices", "wavelet": 22, "dim": 2,
        "groups": ["acm", "3vortices", "adaptive"], "input_files": ["*_000010000000.h5", "*.ini"],
        "commands": ["{mpi_command} {run_dir}/wabbit PARAMS_3vortices.ini --memory={memory}"]},
    {"name": "3vorticesAdaptFD4_CDF40", "type": "simulation", "root_folder": "acm/3vortices", "wavelet": 40, "dim": 2,
        "groups": ["acm", "3vortices", "adaptive"], "input_files": ["*_000010000000.h5", "*.ini"],
        "commands": ["{mpi_command} {run_dir}/wabbit PARAMS_3vortices.ini --memory={memory}"]},
    {"name": "3vorticesAdaptFD4_CDF42", "type": "simulation", "root_folder": "acm/3vortices", "wavelet": 42, "dim": 2,
        "groups": ["acm", "3vortices", "adaptive"], "input_files": ["*_000010000000.h5", "*.ini"],
        "commands": ["{mpi_command} {run_dir}/wabbit PARAMS_3vortices.ini --memory={memory}"]},
    {"name": "3vorticesAdaptFD6_CDF60", "type": "simulation", "root_folder": "acm/3vortices", "wavelet": 60, "dim": 2,
        "groups": ["acm", "3vortices", "adaptive"], "input_files": ["*_000010000000.h5", "*.ini"],
        "commands": ["{mpi_command} {run_dir}/wabbit PARAMS_3vortices.ini --memory={memory}"]},
    {"name": "3vorticesAdaptFD6_CDF62", "type": "simulation", "root_folder": "acm/3vortices", "wavelet": 62, "dim": 2,
        "groups": ["acm", "3vortices", "adaptive"], "input_files": ["*_000010000000.h5", "*.ini"],
        "commands": ["{mpi_command} {run_dir}/wabbit PARAMS_3vortices.ini --memory={memory}"]},

    # Taylor-Green vortex tests
    # Classic fluid dynamics benchmark: decaying vortex flow.
    # Tests different finite difference orders with equispaced meshes.
    # taylorGreen tests
    {"name": "taylorGreenEqui_FD2_CDF20", "type": "simulation", "root_folder": "acm/taylorGreen", "wavelet": 20, "dim": 3,
        "groups": ["acm", "taylorGreen", "equi"], "input_files": ["*.ini"],
        "commands": ["{mpi_command} {run_dir}/wabbit PARAMS_taylor_green.ini --memory={memory}"]},
    {"name": "taylorGreenEqui_FD4_CDF40", "type": "simulation", "root_folder": "acm/taylorGreen", "wavelet": 40, "dim": 3,
        "groups": ["acm", "taylorGreen", "equi"], "input_files": ["*.ini"],
        "commands": ["{mpi_command} {run_dir}/wabbit PARAMS_taylor_green.ini --memory={memory}"]},
    {"name": "taylorGreenEqui_FD6_CDF60", "type": "simulation", "root_folder": "acm/taylorGreen", "wavelet": 60, "dim": 3,
        "groups": ["acm", "taylorGreen", "equi"], "input_files": ["*.ini"],
        "commands": ["{mpi_command} {run_dir}/wabbit PARAMS_taylor_green.ini --memory={memory}"]},

    # Bumblebee flow test
    # Insect-scale fluid flow simulation with kinematics.
    # bumblebeeFlowEquiFD4
    {"name": "bumblebeeFlowEquiFD4_CDF40", "type": "simulation", "root_folder": "acm", "wavelet": 40, "dim": 3,
        "groups": ["acm", "insects"], "input_files": ["*.ini"],
        "commands": ["{mpi_command} {run_dir}/wabbit PARAMS.ini --memory={memory}"]},

    "\n---Mask generation tests (dry runs)---",
    # Dry run tests (insect flight simulations)
    # These are quick validation tests that don't produce full simulation output.
    # They verify the simulation setup without running the full computation.
    # Used for testing various insect geometries and configurations.
    {"name": "dry_fractal_tree_CDF22", "type": "simulation", "root_folder": "insects", "wavelet": 22, "dim": 3,
        "groups": ["mask"], "input_files": ["*.ini"],
        "commands": ["{mpi_command} {run_dir}/wabbit-post --dry-run PARAMS_dry_run.ini --memory={memory} --pruned"]},
    {"name": "dry_bumblebee_CDF22", "type": "simulation", "root_folder": "insects", "wavelet": 22, "dim": 3,
        "groups": ["insects", "mask", "short"], "input_files": ["*.ini"],
        "commands": ["{mpi_command} {run_dir}/wabbit-post --dry-run PARAMS_dry_run.ini --memory={memory} --pruned"]},
    {"name": "dry_emundus_4wings_CDF22", "type": "simulation", "root_folder": "insects", "wavelet": 22, "dim": 3,
        "groups": ["insects", "mask"], "input_files": ["*.ini"],
        "commands": ["{mpi_command} {run_dir}/wabbit-post --dry-run PARAMS_dry_run.ini --memory={memory} --pruned --save-us"]},
    {"name": "dry_muscaComplete_CDF22", "type": "simulation", "root_folder": "insects", "wavelet": 22, "dim": 3,
        "groups": ["insects", "mask"], "input_files": ["*.ini"],
        "commands": ["{mpi_command} {run_dir}/wabbit-post --dry-run PARAMS_dry_run.ini --memory={memory} --pruned --save-us"]},
    {"name": "dry_dipteraFourier_CDF22", "type": "simulation", "root_folder": "insects", "wavelet": 22, "dim": 3,
        "groups": ["insects", "mask"], "input_files": ["*.ini"],
        "commands": ["{mpi_command} {run_dir}/wabbit-post --dry-run PARAMS_dry_run.ini --memory={memory} --pruned --save-us"]},
    {"name": "dry_dipteraHermite_CDF22", "type": "simulation", "root_folder": "insects", "wavelet": 22, "dim": 3,
        "groups": ["insects", "mask"], "input_files": ["*.ini"],
        "commands": ["{mpi_command} {run_dir}/wabbit-post --dry-run PARAMS_dry_run.ini --memory={memory} --pruned --save-us"]},
    {"name": "dry_dipteraBodyRotation_CDF22", "type": "simulation", "root_folder": "insects", "wavelet": 22, "dim": 3,
        "groups": ["insects", "mask"], "input_files": ["*.ini"],
        "commands": ["{mpi_command} {run_dir}/wabbit-post --dry-run PARAMS_dry_run.ini --memory={memory} --pruned --save-us"]},
    {"name": "dry_paratuposaComplete_CDF22", "type": "simulation", "root_folder": "insects", "wavelet": 22, "dim": 3,
        "groups": ["insects", "mask"], "input_files": ["*.ini"],
        "commands": ["{mpi_command} {run_dir}/wabbit-post --dry-run PARAMS_dry_run.ini --memory={memory} --pruned --save-us"]},
    {"name": "dry_butterflyKineloaderV2_CDF22", "type": "simulation", "root_folder": "insects", "wavelet": 22, "dim": 3,
        "groups": ["insects", "mask"], "input_files": ["*.ini", "*.kineloader", "*.superstl"],
        "commands": ["{mpi_command} {run_dir}/wabbit-post --dry-run PARAMS_dry_run.ini --memory={memory} --pruned --save-us"]},
    {"name": "dry_3Dbristles_CDF22", "type": "simulation", "root_folder": "insects", "wavelet": 22, "dim": 3,
        "groups": ["insects", "mask"], "input_files": ["*.ini"],
        "commands": ["{mpi_command} {run_dir}/wabbit-post --dry-run PARAMS_dry_run.ini --memory={memory} --pruned"]},
    # wing geometry models
    {"name": "dry_Insects-Wing-Fourier", "type": "simulation", "root_folder": "insects", "wavelet": 22, "dim": 3,
        "groups": ["insects", "mask"], "input_files": ["*.ini"],
        "commands": ["{mpi_command} {run_dir}/wabbit-post --dry-run PARAMS.ini --memory={memory} --pruned"]},
    {"name": "dry_Insects-Wing-Linear", "type": "simulation", "root_folder": "insects", "wavelet": 22, "dim": 3,
        "groups": ["insects", "mask"], "input_files": ["*.ini"],
        "commands": ["{mpi_command} {run_dir}/wabbit-post --dry-run PARAMS.ini --memory={memory} --pruned"]},
    {"name": "dry_Insects-Wing-Polygon", "type": "simulation", "root_folder": "insects", "wavelet": 22, "dim": 3,
        "groups": ["insects", "mask"], "input_files": ["*.ini"],
        "commands": ["{mpi_command} {run_dir}/wabbit-post --dry-run PARAMS.ini --memory={memory} --pruned"]},
    {"name": "dry_Insects-Wing-BumblebeeHardcoded", "type": "simulation", "root_folder": "insects", "wavelet": 22, "dim": 3,
        "groups": ["insects", "mask"], "input_files": ["*.ini"],
        "commands": ["{mpi_command} {run_dir}/wabbit-post --dry-run PARAMS.ini --memory={memory} --pruned"]},
    {"name": "dry_Insects-Wing-FourierCorrugated", "type": "simulation", "root_folder": "insects", "wavelet": 22, "dim": 3,
        "groups": ["insects", "mask"], "input_files": ["*.ini"],
        "commands": ["{mpi_command} {run_dir}/wabbit-post --dry-run PARAMS.ini --memory={memory} --pruned"]},
    {"name": "dry_Insects-Wing-FourierDamaged", "type": "simulation", "root_folder": "insects", "wavelet": 22, "dim": 3,
        "groups": ["insects", "mask"], "input_files": ["*.ini"],
        "commands": ["{mpi_command} {run_dir}/wabbit-post --dry-run PARAMS.ini --memory={memory} --pruned"]},
    # kinematics
    {"name": "dry_Insects-Kinematics-Fourier", "type": "simulation", "root_folder": "insects", "wavelet": 22, "dim": 3,
        "groups": ["insects", "mask"], "input_files": ["*.ini"],
        "commands": ["{mpi_command} {run_dir}/wabbit-post --dry-run PARAMS.ini --memory={memory} --pruned --save-us"]},
    {"name": "dry_Insects-Kinematics-Hermite", "type": "simulation", "root_folder": "insects", "wavelet": 22, "dim": 3,
        "groups": ["insects", "mask"], "input_files": ["*.ini"],
        "commands": ["{mpi_command} {run_dir}/wabbit-post --dry-run PARAMS.ini --memory={memory} --pruned --save-us"]},
    {"name": "dry_Insects-Kinematics-SuzukiHardcoded", "type": "simulation", "root_folder": "insects", "wavelet": 22, "dim": 3,
        "groups": ["insects", "mask"], "input_files": ["*.ini"],
        "commands": ["{mpi_command} {run_dir}/wabbit-post --dry-run PARAMS.ini --memory={memory} --pruned --save-us"]},

    {"name": "dry_Insects-CompleteModel-Dragonfly", "type": "simulation", "root_folder": "insects", "wavelet": 22, "dim": 3,
        "groups": ["insects", "mask"], "input_files": ["*.ini","*.sstl"],
        "commands": ["{mpi_command} {run_dir}/wabbit-post --dry-run PARAMS.ini --memory={memory} --pruned --save-us"]},
    {"name": "dry_snowman_2D_CDF22", "type": "simulation", "root_folder": "insects", "wavelet": 22, "dim": 2,
        "groups": ["mask"], "input_files": ["*.ini"],
        "commands": ["{mpi_command} {run_dir}/wabbit-post --dry-run PARAMS_dry_run.ini --memory={memory} --pruned"]},
    {"name": "dry_snowman_3D_CDF22", "type": "simulation", "root_folder": "insects", "wavelet": 22, "dim": 3,
        "groups": ["mask"], "input_files": ["*.ini"],
        "commands": ["{mpi_command} {run_dir}/wabbit-post --dry-run PARAMS_dry_run.ini --memory={memory} --pruned"]},
]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EXECUTION
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

FAIL_COLOR = '\033[31;1m'
PASS_COLOR = '\033[92;1m'
END_COLOR = '\033[0m'


def format_string(string, **kwargs):
    """
    Safely format a string with the given keyword arguments.
    
    This function attempts to format the string using the provided keyword arguments.
    If a KeyError occurs (i.e., a placeholder in the string doesn't have a matching
    keyword argument), it returns the original string unchanged instead of raising
    an exception. This allows command templates to work even when some placeholders
    are not applicable to all test types.
    
    Args:
        string: The string to format, potentially containing {placeholders}
        **kwargs: Keyword arguments to substitute into the string
    
    Returns:
        The formatted string, or the original string if formatting fails
    """
    try:
        return string.format(**kwargs)
    except KeyError:
        return string


def run_wabbit_internal(test, mpi_command, memory, run_dir, verbose=False):
    """
    Run a wabbit-internal test.
    
    Wabbit-internal tests verify internal WABBIT wavelet operations without
    running full simulations. They execute directly in the TESTING/{root_folder}/
    directory and write output to a dedicated log file per test:
    {name}_{dim}D_{wavelet}.log
    
    These tests are fast and don't require reference HDF5 files for comparison.
    They typically test wavelet decomposition, reconstruction, and other core
    wavelet operations.
    
    Args:
        test: Test definition dictionary with name, type, root_folder, wavelet, dim, commands
        mpi_command: MPI execution command (e.g., "nice mpirun -n 4")
        memory: Memory allocation flag (e.g., "8.0GB")
        run_dir: Path to the WABBIT directory containing the wabbit executable
        verbose: If True, print command output to screen in addition to log file
    
    Returns:
        subprocess.CompletedProcess: Result of the last command executed, or None if test dir not found
    """
    test_dir = os.path.join(run_dir, "TESTING", test["root_folder"])
    log_file = os.path.join(test_dir, f'{test["name"]}_{test["dim"]}D_{test["wavelet"]}.log')

    if not os.path.isdir(test_dir):
        print(f"    Skip: {test_dir} not found")
        return None

    with open(log_file, 'w') as log_f:
        for command in test.get("commands", [test.get("command")]):
            formatted_cmd = format_string(
                command,
                mpi_command=mpi_command,
                memory=memory,
                run_dir=run_dir,
                wavelet=test.get("wavelet", ""),
                name=test.get("name", ""),
                dim=test.get("dim", "")
            )
            if verbose:
                print(f"    {formatted_cmd}")
            result = subprocess.run(formatted_cmd, shell=True, capture_output=True, text=True, cwd=test_dir)
            if result.stdout:
                log_f.write(result.stdout)
            if result.stderr:
                log_f.write(result.stderr)
            if verbose:
                print(result.stdout + result.stderr, end='')
            if result.returncode != 0:
                return result
    return result


def run_simulation(test, mpi_command, memory, run_dir, verbose=False, write_diff=False, keep_tmp=False):
    """
    Run a simulation test.
    
    Simulation tests run full WABBIT simulations that produce HDF5 output files.
    The workflow is:
    1. Remove existing tmp/ directory to discard old data
    2. Create fresh tmp/ directory under TESTING/{root_folder}/{name}/
    3. Copy input files from test directory to tmp/
    4. Change to tmp/ directory
    5. Log tmp/ directory contents
    6. Execute all commands in tmp/
    7. Compare HDF5 output files against reference files in test directory
    8. Clean up tmp/ directory (unless keep_tmp=True or test failed)
    
    All command output and HDF5 comparison output is captured in log.txt.
    
    Args:
        test: Test definition dictionary with name, type, root_folder, input_files, commands
        mpi_command: MPI execution command
        memory: Memory allocation flag
        run_dir: Path to WABBIT directory
        verbose: If True, print command output to screen
        write_diff: If True, write difference files for failing HDF5 comparisons
        keep_tmp: If True, keep tmp directory after test (also kept on failure)
    
    Returns:
        subprocess.CompletedProcess: Result of the last command, or None if test dir not found
    """
    # Build paths: test directory, tmp directory, and log file
    test_dir = os.path.join(run_dir, "TESTING", test["root_folder"], test["name"])
    tmp_dir = os.path.join(test_dir, "tmp")
    log_file = os.path.join(test_dir, "log.txt")

    # Skip if test directory doesn't exist
    if not os.path.isdir(test_dir):
        print(f"    Skip: {test_dir} not found")
        return None

    # Remove existing tmp directory to start fresh (discard old data)
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir, ignore_errors=True)
    
    # Create tmp directory for this test run
    os.makedirs(tmp_dir, exist_ok=True)
    original_dir = os.getcwd()

    try:
        # Copy all input files from test directory to tmp/
        # Supports glob patterns like "*.kineloader" in input_files
        for pattern in test.get("input_files", []):
            source_pattern = os.path.join(test_dir, pattern)
            for source_file in glob.glob(source_pattern):
                if os.path.isfile(source_file):
                    destination_file = os.path.join(tmp_dir, os.path.basename(source_file))
                    shutil.copy2(source_file, destination_file)

        # Run commands in tmp/ directory
        os.chdir(tmp_dir)

        # Open log file and write tmp/ contents for debugging
        with open(log_file, 'w') as log_f:
            log_f.write("tmp/ folder contents:\n")
            for item in sorted(os.listdir(tmp_dir)):
                log_f.write(f"  {item}\n")
            log_f.write("\n")

            for command in test.get("commands", [test.get("command")]):
                formatted_cmd = format_string(
                    command,
                    mpi_command=mpi_command,
                    memory=memory,
                    run_dir=run_dir,
                    wavelet=test.get("wavelet", ""),
                    name=test.get("name", ""),
                    dim=test.get("dim", "")
                )
                if verbose:
                    print(f"    {formatted_cmd}")
                result = subprocess.run(formatted_cmd, shell=True, capture_output=True, text=True)
                if result.stdout:
                    log_f.write(result.stdout)
                if result.stderr:
                    log_f.write(result.stderr)
                if verbose:
                    print(result.stdout + result.stderr, end='')
                if result.returncode != 0:
                    break

            # Compare HDF5 files if commands succeeded
            if result and result.returncode == 0:
                try:
                    log_f.write("\nHDF5 file comparison:\n")
                    # Always use verbose=True for comparison so output goes to log
                    if not compare_hdf5_files(test_dir, tmp_dir, True, write_diff, log_f):
                        result.returncode = 1
                except Exception as error:
                    log_f.write(f"    Compare error: {error}\n")
                    if verbose:
                        print(f"    Compare error: {error}")
                    result.returncode = 1

        return result
    finally:
        # Always restore original directory
        os.chdir(original_dir)
        # Keep tmp dir if keep_tmp flag is set or if test failed
        if not (keep_tmp or (result and result.returncode != 0)):
            if os.path.exists(tmp_dir):
                shutil.rmtree(tmp_dir, ignore_errors=True)


def compare_hdf5_files(test_dir, tmp_dir, verbose=False, write_diff=False, log_f=None):
    """
    Compare all HDF5 output files against reference files.
    
    This function compares each HDF5 file produced in tmp_dir with its corresponding
    reference file in test_dir. It uses the wabbit_tools.WabbitHDF5file class
    to read both files and check if they are "close" (within numerical tolerance).
    
    If write_diff is True and files differ, a difference file (diff-{filename})
    is written to the test directory for debugging.
    
    All comparison output (including any print statements from wabbit_tools)
    is redirected to log_f if provided, otherwise printed to stdout/stderr.
    
    Args:
        test_dir: Path to the test directory containing reference HDF5 files
        tmp_dir: Path to tmp/ directory containing newly generated HDF5 files
        verbose: If True, wabbit_tools functions produce verbose output
        write_diff: If True, write difference files for mismatched comparisons
        log_f: File handle to write comparison output to (optional)
    
    Returns:
        bool: True if all files match, False if any comparison failed
    """
    all_ok = True
    # Get all HDF5 files produced by the simulation
    tmp_files = sorted(glob.glob(os.path.join(tmp_dir, "*.h5")))

    # check if there are any HDF5 files to compare, otherwise something is wrong
    if not tmp_files:
        msg = "    No HDF5 output files found in tmp/ directory\n"
        if log_f:
            log_f.write(msg)
        else:
            print(msg, end='')
        return False

    for tmp_file in tmp_files:
        filename = os.path.basename(tmp_file)
        # Skip diff and new files (these are outputs from comparison, not simulation)
        if filename.startswith("diff-") or filename.startswith("new-"):
            continue
        # Reference file should be in the test directory with same name
        ref_file = os.path.join(test_dir, filename)

        # Check if reference file exists
        if not os.path.isfile(ref_file):
            msg = f"    No reference file for {filename}\n"
            if log_f:
                log_f.write(msg)
            else:
                print(msg, end='')
            all_ok = False
            continue

        try:
            # Create WabbitHDF5file objects for reference and new output
            state_ref = wabbit_tools.WabbitHDF5file()
            state_new = wabbit_tools.WabbitHDF5file()
            
            # Read and compare files, redirecting all output to log if provided
            if log_f:
                with contextlib.redirect_stdout(log_f), contextlib.redirect_stderr(log_f):
                    # Read both files - these may produce verbose output
                    state_ref.read(ref_file, verbose=verbose)
                    state_new.read(tmp_file, verbose=verbose)
                    # Check if states are close (within numerical tolerance)
                    if not state_ref.isClose(state_new, verbose=verbose):
                        all_ok = False
                        # Write difference file if requested
                        if write_diff:
                            state_diff = state_ref - state_new
                            state_diff.write(os.path.join(test_dir, f"diff-{filename}"))
            else:
                state_ref.read(ref_file, verbose=verbose)
                state_new.read(tmp_file, verbose=verbose)
                if not state_ref.isClose(state_new, verbose=verbose):
                    all_ok = False
                    if write_diff:
                        state_diff = state_ref - state_new
                        state_diff.write(os.path.join(test_dir, f"diff-{filename}"))
        except Exception as error:
            msg = f"    Error comparing {filename}: {error}\n"
            if log_f:
                log_f.write(msg)
            else:
                print(msg, end='')
            all_ok = False

    # check if any hdf5 file from the reference is missing
    ref_files = sorted(glob.glob(os.path.join(test_dir, "*.h5")))
    for ref_file in ref_files:
        # Skip diff and new files (these are outputs from comparison, not simulation)
        if os.path.basename(ref_file).startswith(("diff-", "new-")):
            continue

        # Check if the corresponding output file exists in tmp
        tmp_file = os.path.join(tmp_dir, os.path.basename(ref_file))
        if not os.path.isfile(tmp_file):
            msg = f"    Missing output file for reference {os.path.basename(ref_file)}\n"
            if log_f:
                log_f.write(msg)
            else:
                print(msg, end='')
            all_ok = False

    return all_ok


def iter_tests_with_headers(tests):
    """Yield (section_header, test_dict) pairs, skipping the header strings themselves."""
    header = None
    for item in tests:
        if isinstance(item, str):
            header = item.strip()
            continue
        yield header, item


def test_matches(test, name_patterns, group_patterns):
    """Check whether a test matches the requested -t/-g filters (wildcards via fnmatch)."""
    if name_patterns and not any(fnmatch.fnmatch(test["name"], pat) for pat in name_patterns):
        return False
    if group_patterns:
        test_groups = test.get("groups", [])
        if not any(fnmatch.fnmatch(group, pat) for group in test_groups for pat in group_patterns):
            return False
    return True


def print_test_catalog(tests):
    """Print all available tests grouped by section header, and all known groups."""
    printed_header = None
    all_groups = set()

    grouped_by_name = itertools.groupby(iter_tests_with_headers(tests), key=lambda h_t: (h_t[0], h_t[1]["name"]))
    for (header, name), entries in grouped_by_name:
        entries = list(entries)
        test_groups = entries[0][1].get("groups", [])
        all_groups.update(test_groups)

        if header != printed_header:
            print(f"\n{header}")
            printed_header = header

        variant_info = f"  ({len(entries)} variants)" if len(entries) > 1 else ""
        groups_info = f"  [{', '.join(test_groups)}]" if test_groups else ""
        print(f"  {name}{variant_info}{groups_info}")

    print("\nAvailable groups:")
    for group in sorted(all_groups):
        print(f"  {group}")


def main():
    """
    Main entry point for running WABBIT tests.
    
    This function:
    1. Parses command-line arguments for MPI configuration, memory, etc.
    2. Validates that the wabbit executable exists
    3. Iterates through all test definitions and runs them
    4. Prints colored PASS/FAIL/SKIP status for each test
    5. Prints summary statistics at the end
    
    Test execution order follows the order in the `tests` list.
    String separators (starting with "---") are printed but not executed.
    
    Color coding:
    - Green: PASS (test completed successfully)
    - Red: FAIL (test failed or comparison mismatch)
    - Red: SKIP (test directory not found)
    
    For failing tests, the path to the log file is printed to help with debugging.
    """
    parser = argparse.ArgumentParser(description='Run WABBIT unit tests.')
    parser.add_argument('-v', '--verbose', action='store_true', help='Print test output to screen')
    mpi_group = parser.add_mutually_exclusive_group()
    mpi_group.add_argument('--nprocs', type=int, default=4, help='Number of processors (default: 4)')
    mpi_group.add_argument('--mpi_command', type=str, default=None, help='MPI command (default: "nice mpirun -n {nprocs}")')
    parser.add_argument('--memory', type=str, default='8.0GB', help='Memory flag (default: "8.0GB")')
    parser.add_argument('--write-diff', action='store_true', help='Write difference files for failing tests')
    parser.add_argument('--wabbit-dir', type=str, default=None, help='Directory containing wabbit executable')
    parser.add_argument('--keep-tmp', action='store_true', help='Keep tmp directories for all simulation tests')
    parser.add_argument('--update-failed-tests', action='store_true', help='Update reference data for failed simulation tests')
    parser.add_argument('--list', action='store_true', help='List all available tests and groups, then exit')
    parser.add_argument('-t', '--tests', nargs='+', default=None, metavar='NAME',
                         help='Run only tests matching these names (wildcards allowed, e.g. "blob_*")')
    parser.add_argument('-g', '--groups', nargs='+', default=None, metavar='GROUP',
                         help='Run only tests belonging to these groups (see --list, wildcards allowed)')
    args = parser.parse_args()

    if args.list:
        print_test_catalog(tests)
        return

    mpi_command = args.mpi_command or f"nice mpirun -n {args.nprocs}"
    memory = args.memory
    run_dir = os.path.abspath(args.wabbit_dir) if args.wabbit_dir else os.getcwd()

    # Validate wabbit executable exists
    if not os.path.isfile(os.path.join(run_dir, "wabbit")):
        print(f"ERROR: wabbit executable not found in {run_dir}")
        sys.exit(1)

    print(f"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print(f"         WABBIT unit testing ")
    print(f"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

    # Print configuration for verification
    print(f"Configuration:")
    print(f"  MPI command: {mpi_command}")
    print(f"  Memory: {memory}")
    print(f"  WABBIT directory: {run_dir}")
    print()

    # Select tests according to -t/--tests and -g/--groups (both default to "run everything")
    selected_tests = [(header, test_item) for header, test_item in iter_tests_with_headers(tests)
                       if test_matches(test_item, args.tests, args.groups)]

    if not selected_tests:
        print("No tests matched the given --tests/--groups filters.")
        return

    happy_count = 0  # Number of passed tests
    sad_count = 0     # Number of failed or skipped tests
    start_time = time.time()
    printed_header = None

    # Iterate through the selected test definitions
    for header, test_item in selected_tests:
        # Print the section header once, the first time a selected test appears under it
        if header != printed_header:
            print(header)
            printed_header = header

        # Extract test metadata for display
        test_name = test_item.get("name", "?")
        test_type = test_item.get("type", "?")
        test_wavelet = test_item.get("wavelet", "-----")
        test_dim = test_item.get("dim", "-")
        
        # Print test info: name (of specific length) wavelet=... dim=...D
        print(f"  {test_name:<40} wavelet={test_wavelet} dim={test_dim}D ", end="")

        test_start = time.time()
        
        # Dispatch to appropriate test runner based on type
        if test_type == "wabbit-internal":
            result = run_wabbit_internal(test_item, mpi_command, memory, run_dir, verbose=args.verbose)
        elif test_type == "simulation":
            result = run_simulation(test_item, mpi_command, memory, run_dir, verbose=args.verbose, write_diff=args.write_diff, keep_tmp=args.keep_tmp)
        else:
            result = None  # Unknown test type - will be counted as SKIP

        elapsed = time.time() - test_start

        # Classify and count test results
        if result is None:
            # Test directory not found - test was skipped
            print(f"{FAIL_COLOR}SKIP{END_COLOR} ({elapsed:.1f}s)")
            sad_count += 1
        elif result.returncode == 0:
            # All commands succeeded and HDF5 comparison passed (if applicable)
            print(f"{PASS_COLOR}PASS{END_COLOR} ({elapsed:.1f}s)")
            happy_count += 1
        else:
            # Command failed or HDF5 comparison failed
            print(f"{FAIL_COLOR}FAIL{END_COLOR} ({elapsed:.1f}s)")
            sad_count += 1
            # Print log file path for debugging
            if test_type == "simulation":
                print(f"    Log: {os.path.join(run_dir, 'TESTING', test_item['root_folder'], test_item['name'], 'log.txt')}")
            else:
                log_name = f"{test_item['name']}_{test_item['dim']}D_{test_item['wavelet']}.log"
                print(f"    Log: {os.path.join(run_dir, 'TESTING', test_item['root_folder'], log_name)}")

    # Print summary
    total_time = time.time() - start_time
    print()
    print(f"Results: {PASS_COLOR}{happy_count} passed{END_COLOR}, {FAIL_COLOR}{sad_count} failed{END_COLOR} in {total_time:.1f}s")

    # Handle --update-failed-tests: update reference data for failed simulation tests
    if args.update_failed_tests:
        print("\n" + "=" * 70)
        print("UPDATE FAILED TESTS MODE")
        print("=" * 70)
        print("The following simulation tests failed and have tmp/ directories remaining.")
        print("You can update their reference data with the output from this run.")
        print()
        
        # Collect failed simulation tests with tmp dirs
        failed_sim_tests = []
        for header, test_item in selected_tests:
            if isinstance(test_item, str):
                continue
            if test_item.get("type") == "simulation":
                test_dir = os.path.join(run_dir, "TESTING", test_item["root_folder"], test_item["name"])
                tmp_dir = os.path.join(test_dir, "tmp")
                if os.path.exists(tmp_dir):
                    failed_sim_tests.append(test_item)
        
        if not failed_sim_tests:
            print("No failed simulation tests with tmp/ directories found.")
        else:
            print(f"Found {len(failed_sim_tests)} failed simulation test(s) with tmp/ directories:")
            print()
            
            for test_item in failed_sim_tests:
                test_name = test_item.get("name", "?")
                test_dir = os.path.join(run_dir, "TESTING", test_item["root_folder"], test_item["name"])
                tmp_dir = os.path.join(test_dir, "tmp")
                
                print(f"  Test: {test_name}")
                print(f"    tmp dir: {tmp_dir}")
                print(f"    Reference dir: {test_dir}")
                
                # Only simulation-output HDF5 files are reference data.  Do
                # not update or delete diagnostic diff-/new- files.
                h5_files = [
                    path for path in sorted(glob.glob(os.path.join(tmp_dir, "*.h5")))
                    if not os.path.basename(path).startswith(("diff-", "new-"))
                ]
                reference_h5_files = [
                    path for path in sorted(glob.glob(os.path.join(test_dir, "*.h5")))
                    if not os.path.basename(path).startswith(("diff-", "new-"))
                ]
                output_names = {os.path.basename(path) for path in h5_files}
                stale_reference_files = [
                    path for path in reference_h5_files
                    if os.path.basename(path) not in output_names
                ]
                if h5_files:
                    print(f"    HDF5 files to update: {len(h5_files)}")
                    for f in h5_files:
                        print(f"      - {os.path.basename(f)}")
                else:
                    print(f"    No HDF5 files found in tmp/")

                if stale_reference_files:
                    print(f"    Stale reference HDF5 files to delete: {len(stale_reference_files)}")
                    for f in stale_reference_files:
                        print(f"      - {os.path.basename(f)}")
                
                # Ask for confirmation before overwriting or deleting reference data.
                response = input("    Replace the listed references and delete the listed stale references? [yes/no]: ").strip().lower()
                if response in ["yes", "y"]:
                    updated_count = 0
                    for h5_file in h5_files:
                        dest = os.path.join(test_dir, os.path.basename(h5_file))
                        shutil.copy2(h5_file, dest)
                        updated_count += 1
                    deleted_count = 0
                    for reference_file in stale_reference_files:
                        os.remove(reference_file)
                        deleted_count += 1
                    print(f"    Updated {updated_count} reference file(s); deleted {deleted_count} stale reference file(s)")
                    # Clean up tmp dir after updating
                    if os.path.exists(tmp_dir):
                        shutil.rmtree(tmp_dir, ignore_errors=True)
                else:
                    print(f"    Skipped updating {test_name}")
                print()


if __name__ == "__main__":
    main()
