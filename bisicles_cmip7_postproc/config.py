"""
Config-file driver for bisicles_cmip7_postproc.

Allows all options that can be passed on the command line to be specified in a
YAML or JSON config file instead.  YAML is the preferred format; JSON is a
built-in fallback that requires no extra dependencies.

Config file structure
---------------------
The top-level key ``tool`` selects the workflow (``"flatten"`` or
``"diagnostics"``).  All other keys correspond directly to the Python keyword
arguments of the underlying processing functions, using **underscore-separated**
names (e.g. ``ice_density`` not ``--ice-density``).

Example flatten config (YAML)::

    tool: flatten

    # Input/output
    input: /scratch/cx209c/run/output/
    output_dir: /scratch/cx209c/postproc/

    # Grid
    epsg: 3413
    level: 2

    # Time
    reference_year: 1850
    calendar: 360_day   # UKESM coupled runs use a 360-day calendar

    # Physical constants
    ice_density: 918.0
    water_density: 1028.0
    h_min: 1.0
    cmip7_only: true

    # CF/CMIP7 metadata
    institution: University of Bristol
    source: BISICLES-UKESM1.0
    experiment: historical
    variant_label: r1i1p1f3
    ice_sheet: GrIS

Example diagnostics config (YAML)::

    tool: diagnostics

    input: /scratch/dx030c/run/output/
    output_dir: /scratch/dx030c/postproc/

    reference_year: 1850
    calendar: 360_day

    ice_density: 918.0
    water_density: 1028.0
    gravity: 9.81
    h_min: 1.0

    institution: University of Bristol
    experiment: historical
    variant_label: r1i1p1f3
    ice_sheet: GrIS

The same options can also be expressed in JSON::

    {
        "tool": "flatten",
        "input": "/scratch/cx209c/run/output/",
        "output_dir": "/scratch/cx209c/postproc/",
        "epsg": 3413,
        "calendar": "360_day",
        "ice_sheet": "GrIS"
    }

Programmatic use
----------------
::

    from bisicles_cmip7_postproc.config import run_from_config

    # Run directly from a config file
    run_from_config("my_run.yaml")

    # Override individual settings at call time
    run_from_config("my_run.yaml", overrides={"verbose": False})
"""

import json
from pathlib import Path


# ---------------------------------------------------------------------------
# Config loading
# ---------------------------------------------------------------------------

def load_config(path):
    """
    Load a YAML or JSON config file and return a dict.

    YAML is tried first (requires ``pyyaml``).  If pyyaml is not installed the
    file is parsed as JSON.  Files with a ``.json`` extension are always parsed
    as JSON.

    Parameters
    ----------
    path : str or Path
        Path to the config file.

    Returns
    -------
    dict
        Configuration key-value pairs.

    Raises
    ------
    FileNotFoundError
        If the config file does not exist.
    ValueError
        If the file cannot be parsed as YAML or JSON.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")

    text = path.read_text(encoding="utf-8")

    if path.suffix.lower() == ".json":
        return json.loads(text)

    # Try YAML first, fall back to JSON
    try:
        import yaml  # pyyaml
        return yaml.safe_load(text)
    except ImportError:
        pass

    # pyyaml not available – try JSON
    try:
        return json.loads(text)
    except json.JSONDecodeError as e:
        raise ValueError(
            f"Could not parse '{path}' as JSON (pyyaml is not installed so "
            f"YAML files are not supported). Install pyyaml with:\n"
            f"    pip install pyyaml\n"
            f"JSON parse error: {e}"
        ) from e


# ---------------------------------------------------------------------------
# Flatten runner
# ---------------------------------------------------------------------------

# Keys that are valid for process_plotfile / process_directory
_FLATTEN_KEYS = {
    "input", "output_dir", "plot_pattern",
    "exe_path",                          # canonical executable path key
    "level", "epsg_code", "epsg",        # accept both forms
    "x0", "y0",
    "reference_year", "calendar",
    "cmip7_only", "keep_intermediate",
    "ice_density", "water_density", "h_min",
    "x_min", "x_max", "nx", "y_min", "y_max", "ny", "grid_spec",
    "verbose",
    # nc_kwargs (passed through to write_cmip7_per_variable_netcdfs)
    "institution", "source", "experiment", "variant_label", "ice_sheet",
    "extra_attrs", "frequency",
    "ismip7_mode", "model_id", "member_id",
    "source_id", "ism_id", "ism_member_id",
    "esm_id", "forcing_member_id", "set_counter",
    "group", "contact_name", "contact_email",
}

_DIAG_KEYS = {
    "input", "output_dir",
    "plot_pattern", "exe_path",          # canonical executable path key
    "ice_density", "water_density", "gravity", "h_min",
    "mask_file", "mask_no_start", "mask_no_end",
    "reference_year", "calendar", "verbose",
    # nc_kwargs
    "institution", "source", "experiment", "variant_label", "ice_sheet",
    "extra_attrs", "frequency",
    "ismip7_mode", "model_id", "member_id",
    "source_id", "ism_id", "ism_member_id",
    "esm_id", "forcing_member_id", "set_counter",
    "group", "contact_name", "contact_email",
}

_NC_KWARGS_KEYS = {"institution", "source", "experiment", "variant_label",
                   "ice_sheet", "extra_attrs", "frequency",
                   "ismip7_mode", "model_id", "member_id",
                   # ISMIP7 DRS filename components
                   "source_id", "ism_id", "ism_member_id",
                   "esm_id", "forcing_member_id", "set_counter",
                   # ISMIP7 mandatory global attributes
                   "group", "contact_name", "contact_email"}


def _normalise_flatten_cfg(cfg):
    """Normalise config dict keys to the canonical names expected by flatten."""
    # Allow 'epsg' as shorthand for 'epsg_code'
    if "epsg" in cfg and "epsg_code" not in cfg:
        cfg["epsg_code"] = cfg.pop("epsg")
    else:
        cfg.pop("epsg", None)

    # Ensure boolean flags
    for key in ("cmip7_only", "keep_intermediate", "verbose", "ismip7_mode"):
        if key in cfg:
            cfg[key] = bool(cfg[key])

    # Assemble grid_spec tuple from individual x_min/x_max/nx/y_min/y_max/ny keys
    # (these come from YAML config; CLI assembles grid_spec before calling process_*)
    _gs_keys = ("x_min", "x_max", "nx", "y_min", "y_max", "ny")
    _gs_parts = tuple(cfg.pop(k, None) for k in _gs_keys)
    if any(v is not None for v in _gs_parts):
        if not all(v is not None for v in _gs_parts):
            raise ValueError(
                "Must specify all of x_min, x_max, nx, y_min, y_max, ny, or none of them"
            )
        cfg["grid_spec"] = _gs_parts

    return cfg


def _normalise_diag_cfg(cfg):
    """Normalise config dict keys to the canonical names expected by diagnostics."""
    for key in ("verbose", "ismip7_mode"):
        if key in cfg:
            cfg[key] = bool(cfg[key])

    return cfg


def run_flatten_from_config(cfg, overrides=None):
    """
    Run the flatten workflow from a config dict.

    Parameters
    ----------
    cfg : dict
        Config as returned by :func:`load_config` (with ``tool`` key already
        stripped, or still present — it is ignored).
    overrides : dict, optional
        Key-value pairs that override values from *cfg*.
    """
    from .flatten import process_plotfile, process_directory

    cfg = dict(cfg)
    cfg.pop("tool", None)
    if overrides:
        cfg.update(overrides)

    _normalise_flatten_cfg(cfg)

    input_path = Path(cfg.pop("input"))
    output_dir = cfg.pop("output_dir", None)
    plot_pattern = cfg.pop("plot_pattern", "plot.*.2d.hdf5")
    verbose = cfg.pop("verbose", True)

    # Default output directory to the location of the input file(s)
    if output_dir is None:
        output_dir = str(input_path) if input_path.is_dir() else str(input_path.parent)

    # Split nc_kwargs from processing kwargs
    nc_kwargs = {k: cfg.pop(k) for k in list(cfg) if k in _NC_KWARGS_KEYS}

    if input_path.is_dir():
        cfg.pop("keep_intermediate", None)  # only valid for single-file mode
        process_directory(
            directory=input_path,
            output_dir=output_dir,
            plot_pattern=plot_pattern,
            verbose=verbose,
            **cfg,
            **nc_kwargs,
        )
    elif input_path.is_file():
        process_plotfile(
            plot_file=input_path,
            output_dir=output_dir,
            verbose=verbose,
            **cfg,
            **nc_kwargs,
        )
    else:
        raise FileNotFoundError(f"Input path does not exist: {input_path}")


def run_diagnostics_from_config(cfg, overrides=None):
    """
    Run the diagnostics workflow from a config dict.

    Parameters
    ----------
    cfg : dict
        Config as returned by :func:`load_config`.
    overrides : dict, optional
        Key-value pairs that override values from *cfg*.
    """
    from .diagnostics import process_single_file, process_directory

    cfg = dict(cfg)
    cfg.pop("tool", None)
    if overrides:
        cfg.update(overrides)

    _normalise_diag_cfg(cfg)

    input_path = Path(cfg.pop("input"))
    output_dir = cfg.pop("output_dir", None)
    plot_pattern = cfg.pop("plot_pattern", "plot.*.2d.hdf5")
    verbose = cfg.pop("verbose", True)

    # Default output directory to the location of the input file(s)
    if output_dir is None:
        output_dir = str(input_path) if input_path.is_dir() else str(input_path.parent)

    nc_kwargs = {k: cfg.pop(k) for k in list(cfg) if k in _NC_KWARGS_KEYS}

    if input_path.is_dir():
        process_directory(
            directory=input_path,
            output_dir=output_dir,
            plot_pattern=plot_pattern,
            verbose=verbose,
            **cfg,
            **nc_kwargs,
        )
    elif input_path.is_file():
        process_single_file(
            plot_file=input_path,
            output_dir=output_dir,
            **cfg,
            **nc_kwargs,
        )
    else:
        raise FileNotFoundError(f"Input path does not exist: {input_path}")


# ---------------------------------------------------------------------------
# Top-level entry point
# ---------------------------------------------------------------------------

def run_from_config(config_file, overrides=None):
    """
    Load a config file and run the appropriate BISICLES post-processing workflow.

    Parameters
    ----------
    config_file : str or Path
        Path to a YAML or JSON config file.  Must contain a ``tool`` key set to
        either ``"flatten"`` or ``"diagnostics"``.
    overrides : dict, optional
        Key-value pairs that override values from the config file.  Useful for
        making small adjustments without editing the file (e.g. changing the
        output path in a batch loop).

    Raises
    ------
    ValueError
        If the ``tool`` key is missing or unrecognised.

    Examples
    --------
    ::

        from bisicles_cmip7_postproc.config import run_from_config

        run_from_config("GrIS_flatten.yaml")
        run_from_config("GrIS_diag.yaml", overrides={"calendar": "360_day"})
    """
    cfg = load_config(config_file)
    tool = cfg.get("tool", "").lower()

    if tool == "flatten":
        run_flatten_from_config(cfg, overrides=overrides)
    elif tool == "diagnostics":
        run_diagnostics_from_config(cfg, overrides=overrides)
    else:
        raise ValueError(
            f"Config file must contain 'tool: flatten' or 'tool: diagnostics'. "
            f"Got: {cfg.get('tool')!r}"
        )
