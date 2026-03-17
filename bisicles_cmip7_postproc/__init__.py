"""
bisicles_cmip7_postproc
====================

Post-processing package for BISICLES ice sheet model output from UKESM
simulations.  Converts BISICLES plot HDF5 files into CF/CMIP7-compliant
NetCDF files by wrapping the compiled BISICLES file tools.

Two main workflows are provided:

**Scalar diagnostics timeseries**
    Use the ``diagnostics`` file tool to compute integrated quantities
    (ice volume, mass, area, mass balance fluxes) and write a CF-compliant
    scalar timeseries NetCDF covering a simulation period.

    .. code-block:: python

        from bisicles_cmip7_postproc.diagnostics import process_directory

        process_directory(
            "/run/output/",
            "GrIS_diagnostics.nc",
            ice_sheet="GrIS",
            experiment="historical",
        )

**2D spatial fields**
    Use the ``flatten`` file tool to project the multi-level AMR grid onto
    a single uniform grid and write a CMIP7/CF-compliant 2D NetCDF.
    Variables are renamed to CMIP7 standard names and unit-converted.
    Derived fields (sftgrf, sftflf) are computed from the BISICLES mask.

    .. code-block:: python

        from bisicles_cmip7_postproc.flatten import process_plotfile, process_directory

        # Single file
        process_plotfile(
            "plot.000100.2d.hdf5",
            "plot.000100_cmip7.nc",
            level=2,
            epsg_code=3413,
            ice_sheet="GrIS",
        )

        # Full run directory
        process_directory(
            "/run/output/",
            output_dir="/run/postproc/",
            level=2,
            epsg_code=3413,
            ice_sheet="GrIS",
        )

Command-line tools
------------------
After installation (``pip install -e .``):

    bisicles-diagnostics --input /run/output/ --output timeseries.nc --ice-sheet GrIS
    bisicles-flatten     --input /run/output/ --output-dir /run/postproc/ --epsg 3413

Modules
-------
cmip7_vars
    CMIP7/CF variable metadata tables and BISICLES->CMIP7 field mapping.
cf_utils
    Utilities for CF metadata (global attributes, time encoding, CRS variables).
diagnostics
    Wrapper for the BISICLES diagnostics executable; scalar timeseries writer.
flatten
    Wrapper for the BISICLES flatten executable; 2D CMIP7 NetCDF writer.
cli
    Argparse-based command-line entry points.
"""

from .diagnostics import (
    run_diagnostics,
    parse_diagnostics_csv,
    process_single_file as process_diagnostics_single,
    process_directory as process_diagnostics_directory,
)
from .flatten import (
    run_flatten,
    process_plotfile,
    process_directory as process_flatten_directory,
)
from .filename_parser import parse_bisicles_filename, BISICLESFileInfo
from .config import load_config, run_from_config

__all__ = [
    "run_diagnostics",
    "parse_diagnostics_csv",
    "process_diagnostics_single",
    "process_diagnostics_directory",
    "run_flatten",
    "process_plotfile",
    "process_flatten_directory",
    "parse_bisicles_filename",
    "BISICLESFileInfo",
    "load_config",
    "run_from_config",
]
