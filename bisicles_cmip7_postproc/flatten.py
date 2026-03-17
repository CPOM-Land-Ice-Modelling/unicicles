"""
Wrapper for the BISICLES flatten file tool and post-processing utilities for
writing CMIP7/CF-compliant 2D spatial NetCDF files from BISICLES plot HDF5 files.

Workflow
--------
1. Call the flatten executable to collapse the multi-level AMR hierarchy onto a
   single uniform grid and write a preliminary NetCDF file.
2. Open that NetCDF file in Python, apply unit conversions, rename variables to
   CMIP7 names, compute derived fields (sftgrf, sftflf from mask + iceFrac),
   and write a new, fully CF-compliant NetCDF file with all required metadata.

The flatten tool usage::

    flatten2d.*.ex <input.2d.hdf5> <output.nc> <level> [x0 [y0]]

where ``level`` is the refinement level to flatten onto (0 = coarsest grid,
1 = first refined level, etc.).

Grounded/floating mask fallback
---------------------------------
The derived fields ``sftgrf`` and ``sftflf`` require knowledge of which cells
contain grounded versus floating ice.  BISICLES plot files may or may not
include an explicit ``mask`` variable or ``iceFrac`` field.  When these are
absent, the module falls back to computing them from ice geometry:

* The **flotation criterion** determines whether ice at a given cell is
  grounded or floating based solely on thickness (*h*) and bed topography
  (*b*)::

      flotation_thickness = max(0, -b * rho_w / rho_i)
      floating  if  h > h_min  AND  b < 0  AND  h <= flotation_thickness
      grounded  if  h > h_min  AND  NOT floating

* When ``iceFrac`` is absent a binary approximation is used:
  ``iceFrac = 1`` where ``h > h_min``, else ``0``.  This overestimates ice
  extent at sub-grid margins; the ``bisicles_name`` attribute on the written
  variable records how the field was derived.

The function :func:`compute_flotation_mask` is also available as a public API
for use in other contexts.
"""

import os
import subprocess
import tempfile
from pathlib import Path

import numpy as np

from .cmip7_vars import (
    FIELD_MAPPING,
    CF_FIELD_MAPPING,
    DERIVED_FIELDS,
    GROUNDED_MASK_VAL,
    FLOATING_MASK_VAL,
    OPEN_SEA_MASK_VAL,
    OPEN_LAND_MASK_VAL,
)
from .cf_utils import (
    FILL_VALUE,
    UKESM_GRID_ORIGINS,
    get_global_attributes,
    add_time_variable,
    add_time_bounds,
    add_xy_variables,
    add_crs_variable,
)
from .filename_parser import parse_bisicles_filename


# ---------------------------------------------------------------------------
# Locating the flatten executable
# ---------------------------------------------------------------------------

_DEFAULT_EXE_GLOB = "flatten2d*.ex"

_FILETOOLS_REL = Path(__file__).parent.parent.parent / "bisicles-uob" / "code" / "filetools"


def find_flatten_exe(exe_path=None):
    """
    Locate the BISICLES flatten executable.

    Parameters
    ----------
    exe_path : str or Path, optional
        Explicit path to the executable.

    Returns
    -------
    str
        Absolute path to the executable.

    Raises
    ------
    FileNotFoundError
    """
    if exe_path is not None:
        p = Path(exe_path)
        if p.is_dir():
            raise ValueError(
                f"exe_path points to a directory, not an executable: {exe_path}\n"
                "Provide the full path to the executable file, e.g.:\n"
                "  exe_path: /path/to/filetools/flatten2d.Linux.64.g++.gfortran.DEBUG.OPT.ex\n"
                "or leave exe_path unset for auto-detection."
            )
        if p.is_file():
            return str(p)
        raise FileNotFoundError(f"Flatten executable not found at {exe_path}")

    # Use os.listdir (requires only read permission on the directory).
    try:
        entries = os.listdir(_FILETOOLS_REL)
    except (PermissionError, FileNotFoundError, OSError):
        entries = []

    import fnmatch
    for name in sorted(entries):
        if fnmatch.fnmatch(name, _DEFAULT_EXE_GLOB):
            return str(_FILETOOLS_REL / name)

    raise FileNotFoundError(
        "Could not locate the BISICLES flatten executable. "
        f"Searched for '{_DEFAULT_EXE_GLOB}' in:\n"
        f"  {_FILETOOLS_REL}\n"
        "Pass exe_path explicitly, e.g.:\n"
        "  process_plotfile(..., exe_path='/path/to/flatten2d.*.ex')\n"
        "or set 'exe_path' in your config file."
    )


# ---------------------------------------------------------------------------
# Running the flatten tool
# ---------------------------------------------------------------------------

def run_flatten(
    input_hdf5,
    output_nc,
    level=0,
    exe_path=None,
    x0=None,
    y0=None,
):
    """
    Run the BISICLES flatten tool to project an AMR hierarchy onto a uniform grid.

    Parameters
    ----------
    input_hdf5 : str or Path
        Path to the BISICLES plot HDF5 file.
    output_nc : str or Path
        Path for the output NetCDF file (must end in .nc).
    level : int
        AMR level to flatten onto.  0 = coarsest grid; higher values give finer
        resolution.
    exe_path : str or Path, optional
        Path to the flatten executable (auto-detected if not given).
    x0 : float, optional
        X-coordinate origin override (metres).
    y0 : float, optional
        Y-coordinate origin override (metres).

    Raises
    ------
    RuntimeError
        If the flatten tool exits with a non-zero return code.
    """
    if level < 0:
        raise ValueError(
            f"level must be >= 0 (0 = coarsest grid). Got: {level}"
        )
    exe = find_flatten_exe(exe_path)
    args = [exe, str(input_hdf5), str(output_nc), str(level)]
    if x0 is not None and y0 is not None:
        args += [str(x0), str(y0)]

    result = subprocess.run(args, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"flatten tool failed (return code {result.returncode}):\n"
            f"stdout: {result.stdout}\n"
            f"stderr: {result.stderr}"
        )


# ---------------------------------------------------------------------------
# Flotation-based mask computation
# ---------------------------------------------------------------------------

def compute_flotation_mask(
    thickness,
    topography,
    ice_density=918.0,
    water_density=1028.0,
    h_min=1.0,
):
    """
    Compute an integer ice-type mask from ice geometry using the flotation criterion.

    A cell is classified as floating when the ice is thin enough that the
    ocean can support it hydrostatically:

    .. code-block:: text

        flotation_thickness = max(0, -topography * rho_w / rho_i)
        floating  if  thickness > h_min  AND  topography < 0
                                         AND  thickness <= flotation_thickness
        grounded  if  thickness > h_min  AND  NOT floating

    This replicates the logic in BISICLES ``SigmaCSF.ChF`` / ``LevelSigmaCS``.

    Parameters
    ----------
    thickness : ndarray
        Ice thickness in metres.
    topography : ndarray
        Bed topography (positive up) in metres.  Negative values are below
        present-day sea level.
    ice_density : float
        Ice density in kg m-3 (default 918.0).
    water_density : float
        Ocean water density in kg m-3 (default 1028.0).
    h_min : float
        Minimum ice thickness threshold in metres below which a cell is
        treated as ice-free (default 1.0 m).

    Returns
    -------
    mask : ndarray of int, same shape as *thickness*
        Integer mask with values:

        * ``GROUNDED_MASK_VAL`` (1) – grounded ice
        * ``FLOATING_MASK_VAL`` (2) – floating ice
        * ``OPEN_SEA_MASK_VAL`` (4) – ice-free ocean (bed < 0)
        * ``OPEN_LAND_MASK_VAL`` (8) – ice-free land (bed >= 0)
    """
    # Thickness of ice needed just to ground: if bed is above sea level,
    # ice is always grounded regardless of thickness.
    flotation_thickness = np.maximum(
        0.0, -topography * (water_density / ice_density)
    )

    has_ice = thickness > h_min
    below_sea = topography < 0.0

    is_floating = has_ice & below_sea & (thickness <= flotation_thickness)
    is_grounded = has_ice & ~is_floating
    is_open_sea = ~has_ice & below_sea

    # Start with everything as open land, then overwrite
    mask = np.full(thickness.shape, OPEN_LAND_MASK_VAL, dtype=np.int32)
    mask = np.where(is_open_sea, OPEN_SEA_MASK_VAL, mask)
    mask = np.where(is_floating, FLOATING_MASK_VAL, mask)
    mask = np.where(is_grounded, GROUNDED_MASK_VAL, mask)
    return mask


# ---------------------------------------------------------------------------
# Reading the flatten tool output and applying CMIP7 conversions
# ---------------------------------------------------------------------------

def _read_flatten_nc(nc_path):
    """
    Read a NetCDF file produced by the flatten tool and return its contents.

    Returns
    -------
    dict with keys:
      'x', 'y'       : 1-D coordinate arrays (m)
      'time'         : scalar simulation time (years)
      'epsg'         : EPSG code (int) or None
      'variables'    : dict of {bisicles_name: ndarray (ny, nx)}
      'raw_attrs'    : dict of global attributes from the flatten NC
    """
    from netCDF4 import Dataset

    result = {"variables": {}, "x": None, "y": None, "time": None, "epsg": None}

    with Dataset(str(nc_path), "r") as ds:
        result["raw_attrs"] = {k: ds.getncattr(k) for k in ds.ncattrs()}

        # Coordinate variables
        if "x" in ds.variables:
            result["x"] = ds.variables["x"][:].data
        if "y" in ds.variables:
            result["y"] = ds.variables["y"][:].data
        if "time" in ds.variables:
            tv = ds.variables["time"]
            t_data = tv[:]
            # flatten tool may write a 1-D array; take the first (or only) value
            result["time"] = float(np.asarray(t_data).flat[0])
            # Try to extract simulation time in years from units attribute
            if hasattr(tv, "units"):
                result["time_units"] = tv.units
        if "crs" in ds.variables:
            crs_v = ds.variables["crs"]
            if hasattr(crs_v, "epsg_code"):
                try:
                    result["epsg"] = int(str(crs_v.epsg_code).replace("EPSG:", ""))
                except ValueError:
                    pass

        # Data variables (skip coordinate/metadata variables)
        skip = {"x", "y", "time", "crs", "level"}
        for name in ds.variables:
            if name in skip:
                continue
            v = ds.variables[name]
            if v.ndim < 2:
                continue
            arr = v[:]
            if hasattr(arr, "data"):
                arr = arr.data.astype(float)
                if hasattr(v, "_FillValue"):
                    arr[arr == float(v._FillValue)] = np.nan
            result["variables"][name] = arr

    return result


def _compute_derived_fields(
    variables,
    ice_density=918.0,
    water_density=1028.0,
    h_min=1.0,
):
    """
    Compute derived CMIP7 area-fraction fields (sftgrf, sftflf) in-place.

    Prioritises model-output fields where available, falling back to
    geometry-based computation when they are absent:

    * **mask** – used directly if present; otherwise computed from
      ``thickness`` and ``bedTopography`` via the flotation criterion
      (:func:`compute_flotation_mask`).
    * **iceFrac** – used directly if present; otherwise approximated as a
      binary field (1 where ``thickness > h_min``, else 0).

    The fallback paths are recorded on the returned ``_derived_sources`` dict
    which ``write_cmip7_netcdf`` uses to annotate the output variables.

    Parameters
    ----------
    variables : dict
        Mutable dict of {bisicles_name: ndarray}.  Modified in-place.
    ice_density : float
        Ice density in kg m-3, used only when computing the flotation mask.
    water_density : float
        Ocean water density in kg m-3, used only when computing the flotation mask.
    h_min : float
        Ice thickness threshold in metres, used when deriving iceFrac and
        when computing the flotation mask.

    Returns
    -------
    sources : dict
        Maps derived variable name -> provenance string describing how the
        field was computed.  Empty if no derived fields could be produced.
    """
    sources = {}

    # ------------------------------------------------------------------
    # Resolve ice fraction
    # ------------------------------------------------------------------
    if "iceFrac" in variables:
        ice_frac = variables["iceFrac"]
        icefrac_source = "iceFrac from plot file"
    elif "thickness" in variables:
        ice_frac = np.where(variables["thickness"] > h_min, 1.0, 0.0)
        icefrac_source = (
            f"binary approximation (1 where thickness > {h_min} m, else 0); "
            "sub-grid margin fractions are not represented"
        )
    else:
        return sources  # Cannot proceed without at least thickness

    # ------------------------------------------------------------------
    # Resolve grounded/floating mask
    # ------------------------------------------------------------------
    if "mask" in variables:
        mask = variables["mask"]
        mask_source = "mask from plot file"
    elif "thickness" in variables and "bedTopography" in variables:
        mask = compute_flotation_mask(
            variables["thickness"],
            variables["bedTopography"],
            ice_density=ice_density,
            water_density=water_density,
            h_min=h_min,
        )
        mask_source = (
            f"computed from flotation criterion using thickness and bedTopography "
            f"(rho_ice={ice_density} kg/m3, rho_water={water_density} kg/m3, "
            f"h_min={h_min} m)"
        )
    else:
        return sources  # Cannot determine grounded/floating without geometry

    # ------------------------------------------------------------------
    # Compute sftgrf and sftflf
    # ------------------------------------------------------------------
    ice_frac_clean = np.where(np.isnan(ice_frac), 0.0, ice_frac)
    mask_int = np.round(mask).astype(np.int32)

    variables["sftgrf"] = np.where(mask_int == GROUNDED_MASK_VAL, ice_frac_clean, 0.0)
    variables["sftflf"] = np.where(mask_int == FLOATING_MASK_VAL, ice_frac_clean, 0.0)

    provenance = f"ice_frac: {icefrac_source}; mask: {mask_source}"
    sources["sftgrf"] = provenance
    sources["sftflf"] = provenance
    return sources


def write_cmip7_netcdf(
    flatten_data,
    output_nc,
    epsg_code=None,
    x0=None,
    y0=None,
    reference_year=1850,
    calendar="gregorian",
    cmip7_only=False,
    ice_density=918.0,
    water_density=1028.0,
    h_min=1.0,
    institution="",
    source="BISICLES adaptive mesh refinement ice sheet model",
    experiment="",
    variant_label="",
    ice_sheet="",
    extra_attrs=None,
    file_info=None,
):
    """
    Write a CMIP7/CF-compliant 2D spatial NetCDF from flatten tool output.

    Parameters
    ----------
    flatten_data : dict
        Dictionary as returned by :func:`_read_flatten_nc`.
    output_nc : str or Path
        Path for the output NetCDF file.
    epsg_code : int, optional
        EPSG code for the projection.  Overrides the value read from the flatten
        NetCDF if supplied.  If None, the value from flatten_data is used.
    reference_year : int
        Reference year for the time axis (default 1850).
    calendar : str
        CF calendar name for the time axis: ``"gregorian"`` (default, 365.25
        days/year as used by BISICLES internally) or ``"360_day"`` (360
        days/year as used by UKESM).
    cmip7_only : bool
        If True, only write variables present in :data:`FIELD_MAPPING` and the
        derived fields.  If False (default), also write any unmapped BISICLES
        variables with their original names so no data is lost.
    ice_density : float
        Ice density in kg m-3 used when deriving the grounded/floating mask
        from the flotation criterion (default 918.0).  Ignored if the plot
        file contains an explicit ``mask`` variable.
    water_density : float
        Ocean water density in kg m-3, used with *ice_density* for the
        flotation criterion (default 1028.0).
    h_min : float
        Minimum ice thickness in metres below which a cell is treated as
        ice-free when computing the fallback mask or iceFrac (default 1.0 m).
    institution : str
        Written to the global 'institution' attribute.
    source : str
        Written to the global 'source' attribute.
    experiment : str
        Experiment identifier.
    variant_label : str
        CMIP variant label.
    ice_sheet : str
        Ice sheet identifier (e.g. 'GrIS', 'AIS').
    extra_attrs : dict, optional
        Additional global attributes.
    file_info : BISICLESFileInfo, optional
        Parsed filename metadata from :func:`~.filename_parser.parse_bisicles_filename`.
        When supplied the simulation time, time bounds, and ``cell_methods``
        time component are derived from the filename rather than the HDF5
        internal time (which is always 0 in UKESM-coupled runs).
    """
    from netCDF4 import Dataset

    x = flatten_data.get("x")
    y = flatten_data.get("y")
    variables = dict(flatten_data.get("variables", {}))
    epsg = epsg_code if epsg_code is not None else flatten_data.get("epsg")

    # Determine time from filename metadata (preferred) or HDF5 internal value
    if file_info is not None:
        time_years = file_info.time_years
        is_time_mean = file_info.is_time_mean
        time_start_years = file_info.start_time_years
        time_end_years = file_info.end_time_years
        time_cell_method = file_info.cell_methods_time
    else:
        time_years = flatten_data.get("time")
        is_time_mean = False
        time_start_years = None
        time_end_years = None
        time_cell_method = "time: point"

    if x is None or y is None:
        raise ValueError("Flatten data is missing x or y coordinate arrays.")

    # Compute derived fields before we start writing; capture provenance
    derived_sources = _compute_derived_fields(
        variables,
        ice_density=ice_density,
        water_density=water_density,
        h_min=h_min,
    )

    global_attrs = get_global_attributes(
        institution=institution,
        source=source,
        experiment=experiment,
        variant_label=variant_label,
        ice_sheet=ice_sheet,
    )
    if extra_attrs:
        global_attrs.update(extra_attrs)

    ny, nx = len(y), len(x)
    crs_name = "crs"

    with Dataset(str(output_nc), "w", format="NETCDF4") as ds:
        ds.setncatts(global_attrs)

        # Dimensions
        ds.createDimension("time", 1)
        ds.createDimension("y", ny)
        ds.createDimension("x", nx)

        # Coordinate variables
        add_xy_variables(ds, x, y, epsg_code=epsg)
        if time_years is not None:
            add_time_variable(ds, [time_years], reference_year=reference_year, calendar=calendar)
            if is_time_mean and time_start_years is not None:
                add_time_bounds(
                    ds,
                    [time_start_years],
                    [time_end_years],
                    reference_year=reference_year,
                    calendar=calendar,
                )

        # CRS (grid_mapping) variable
        if epsg is not None:
            add_crs_variable(ds, epsg, x0=x0, y0=y0)
        grid_mapping = crs_name if epsg is not None else None

        # ------------------------------------------------------------------
        # Write data variables
        # Track which input variable names have been written.
        written_bisicles_names = set()   # BISICLES internal names (FIELD_MAPPING)
        written_cf_names = set()         # CF/CMIP7 output names (CF_FIELD_MAPPING)

        def _write_var(ds, out_name, arr_raw, mapping, time_cell_method,
                       grid_mapping, bisicles_name=None):
            """Write one 2-D data variable with full CF/CMIP7 metadata."""
            arr = arr_raw * mapping["conversion_factor"]
            arr = np.where(np.isnan(arr), FILL_VALUE, arr)
            arr = arr[np.newaxis, :, :]  # add time dimension
            # cell_methods is stored in CMIP7 format with "time: mean" embedded.
            # For snapshot files substitute "time: point" in place of "time: mean".
            cell_methods = mapping["cell_methods"]
            if time_cell_method == "time: point":
                cell_methods = cell_methods.replace("time: mean", "time: point")
            var = ds.createVariable(
                out_name, "f4", ("time", "y", "x"), fill_value=FILL_VALUE
            )
            var[:] = arr
            var.standard_name = mapping["standard_name"]
            var.long_name = mapping["long_name"]
            var.units = mapping["cmip7_units"]
            var.cell_methods = cell_methods
            var.modeling_realm = mapping["modeling_realm"]
            var.coordinates = "time y x"
            if grid_mapping:
                var.grid_mapping = grid_mapping
            if "comment" in mapping:
                var.comment = mapping["comment"]
            if bisicles_name is not None:
                var.bisicles_name = bisicles_name
            var.bisicles_units = mapping["bisicles_units"]

        # 1. BISICLES internal name -> CMIP7 (standard plot files)
        for bisicles_name, mapping in FIELD_MAPPING.items():
            if bisicles_name not in variables:
                continue
            _write_var(ds, mapping["cmip7_name"], variables[bisicles_name],
                       mapping, time_cell_method, grid_mapping,
                       bisicles_name=bisicles_name)
            written_bisicles_names.add(bisicles_name)

        # 1b. BISICLES CF output name -> CMIP7 (CF-plot files: plot.CF-*.hdf5)
        # Variable names in the file are already CMIP7-compatible but need
        # unit conversion (per-year -> per-second) and full CF metadata.
        for cf_name, mapping in CF_FIELD_MAPPING.items():
            if cf_name not in variables:
                continue
            if cf_name in written_cf_names:
                continue  # already written (shouldn't happen, but be safe)
            _write_var(ds, mapping["cmip7_name"], variables[cf_name],
                       mapping, time_cell_method, grid_mapping)
            written_cf_names.add(cf_name)

        # 2. Derived fields (sftgrf, sftflf computed from iceFrac + mask).
        # Skip any that were already written directly from the CF-plot file.
        for derived_name, dmeta in DERIVED_FIELDS.items():
            if derived_name in written_cf_names:
                continue  # CF file provided these directly; don't overwrite
            if derived_name not in variables:
                continue
            arr = variables[derived_name] * dmeta["conversion_factor"]
            arr = np.where(np.isnan(arr), FILL_VALUE, arr)
            arr = arr[np.newaxis, :, :]
            cell_methods = dmeta["cell_methods"]
            if time_cell_method == "time: point":
                cell_methods = cell_methods.replace("time: mean", "time: point")
            var = ds.createVariable(
                derived_name, "f4", ("time", "y", "x"), fill_value=FILL_VALUE
            )
            var[:] = arr
            var.standard_name = dmeta["standard_name"]
            var.long_name = dmeta["long_name"]
            var.units = dmeta["cmip7_units"]
            var.cell_methods = cell_methods
            var.modeling_realm = dmeta["modeling_realm"]
            var.coordinates = "time y x"
            if grid_mapping:
                var.grid_mapping = grid_mapping
            if "comment" in dmeta:
                var.comment = dmeta["comment"]
            if derived_name in derived_sources:
                var.derived_from = derived_sources[derived_name]

        # 3. Truly unmapped variables (preserved unless cmip7_only=True)
        if not cmip7_only:
            for bname, arr in variables.items():
                if bname in written_bisicles_names:
                    continue
                if bname in written_cf_names:
                    continue
                if bname in DERIVED_FIELDS:
                    continue  # derived, already handled above
                if arr.ndim < 2:
                    continue
                safe_name = bname.replace("/", "_")
                arr_out = np.where(np.isnan(arr), FILL_VALUE, arr)
                arr_out = arr_out[np.newaxis, :, :]
                var = ds.createVariable(
                    safe_name, "f4", ("time", "y", "x"), fill_value=FILL_VALUE
                )
                var[:] = arr_out
                var.long_name = bname
                var.coordinates = "time y x"
                var.note = "Unmapped BISICLES field; units as in original plot file."
                if grid_mapping:
                    var.grid_mapping = grid_mapping


# ---------------------------------------------------------------------------
# High-level convenience functions
# ---------------------------------------------------------------------------

def process_plotfile(
    plot_file,
    output_nc,
    exe_path=None,
    level=0,
    epsg_code=None,
    x0=None,
    y0=None,
    reference_year=1850,
    calendar="gregorian",
    cmip7_only=False,
    ice_density=918.0,
    water_density=1028.0,
    h_min=1.0,
    keep_intermediate=False,
    verbose=True,
    **nc_kwargs,
):
    """
    Flatten a BISICLES plot HDF5 file and write a CMIP7/CF-compliant NetCDF.

    This is the main entry point for processing a single plot file.

    Parameters
    ----------
    plot_file : str or Path
        Input BISICLES plot HDF5 file.
    output_nc : str or Path
        Output CF/CMIP7 NetCDF file.
    exe_path : str or Path, optional
        Path to the flatten executable (auto-detected if not given).
    level : int
        AMR level to flatten onto (0 = coarsest grid).
    epsg_code : int, optional
        EPSG code for the projection.  Overrides value from HDF5 metadata.
    x0 : float, optional
        X-coordinate of the lower-left corner of the domain (metres).
        Defaults to the UKESM standard value for the given ``epsg_code``
        (see :data:`cf_utils.UKESM_GRID_ORIGINS`).
    y0 : float, optional
        Y-coordinate of the lower-left corner of the domain (metres).
        Defaults to the UKESM standard value for the given ``epsg_code``.
    reference_year : int
        Reference year for the CF time axis.
    calendar : str
        CF calendar name: ``"gregorian"`` (default) or ``"360_day"`` for
        UKESM-coupled runs.
    cmip7_only : bool
        If True, only write CMIP7-named variables.
    ice_density : float
        Ice density in kg m-3 for the flotation-based mask fallback (default 918.0).
    water_density : float
        Ocean water density in kg m-3 for the flotation-based mask fallback
        (default 1028.0).
    h_min : float
        Minimum ice thickness threshold in metres for the mask/iceFrac fallback
        (default 1.0 m).
    keep_intermediate : bool
        If True, keep the intermediate flatten NetCDF file (saved next to
        output_nc with suffix ``_flatten_raw.nc``).
    verbose : bool
        Print progress messages.
    **nc_kwargs
        Keyword arguments for :func:`write_cmip7_netcdf`
        (institution, source, experiment, variant_label, ice_sheet, extra_attrs).
    """
    plot_file = Path(plot_file)
    output_nc = Path(output_nc)

    # Intermediate file location
    if keep_intermediate:
        intermediate_nc = output_nc.with_name(
            output_nc.stem + "_flatten_raw" + output_nc.suffix
        )
    else:
        fd, _tmp = tempfile.mkstemp(suffix=".nc", prefix="bisicles_flatten_")
        os.close(fd)
        intermediate_nc = Path(_tmp)

    # Parse simulation time from the filename (UKESM BISICLES internal
    # time is always 0, so the filename is the authoritative source).
    file_info = parse_bisicles_filename(plot_file)
    if file_info is not None and verbose:
        kind = "time-mean" if file_info.is_time_mean else "snapshot"
        print(
            f"  Parsed filename: ice_sheet={file_info.ice_sheet}, "
            f"time={file_info.time_years:.2f} yr ({kind})"
        )

    # Override ice_sheet from filename when not set explicitly
    if file_info is not None and not nc_kwargs.get("ice_sheet"):
        nc_kwargs["ice_sheet"] = file_info.ice_sheet

    # Fall back to UKESM standard grid origins when x0/y0 not supplied
    if (x0 is None or y0 is None) and epsg_code in UKESM_GRID_ORIGINS:
        origin = UKESM_GRID_ORIGINS[epsg_code]
        if x0 is None:
            x0 = origin["x0"]
        if y0 is None:
            y0 = origin["y0"]
        if verbose:
            print(f"  Using UKESM default grid origin: x0={x0}, y0={y0}")

    try:
        if verbose:
            print(f"Flattening {plot_file.name} onto level {level}...")
        run_flatten(
            plot_file,
            intermediate_nc,
            level=level,
            exe_path=exe_path,
            x0=x0,
            y0=y0,
        )

        if verbose:
            print(f"Reading flattened data from {intermediate_nc.name}...")
        flatten_data = _read_flatten_nc(intermediate_nc)

        if verbose:
            bisicles_vars = list(flatten_data["variables"].keys())
            print(f"  Variables found: {bisicles_vars}")

        if verbose:
            print(f"Writing CMIP7 NetCDF to {output_nc}...")
        write_cmip7_netcdf(
            flatten_data,
            output_nc,
            epsg_code=epsg_code,
            x0=x0,
            y0=y0,
            reference_year=reference_year,
            calendar=calendar,
            cmip7_only=cmip7_only,
            ice_density=ice_density,
            water_density=water_density,
            h_min=h_min,
            file_info=file_info,
            **nc_kwargs,
        )
    finally:
        if not keep_intermediate and intermediate_nc.exists():
            os.unlink(intermediate_nc)

    if verbose:
        print("Done.")


def process_directory(
    directory,
    output_dir=None,
    plot_pattern="plot.*.2d.hdf5",
    exe_path=None,
    level=0,
    epsg_code=None,
    x0=None,
    y0=None,
    reference_year=1850,
    calendar="gregorian",
    cmip7_only=False,
    ice_density=918.0,
    water_density=1028.0,
    h_min=1.0,
    verbose=True,
    **nc_kwargs,
):
    """
    Flatten all BISICLES plot files in a directory, writing one CMIP7 NetCDF
    per plot file.

    Parameters
    ----------
    directory : str or Path
        Directory containing BISICLES plot HDF5 files.
    output_dir : str or Path, optional
        Directory to write output NetCDF files.  Defaults to the same directory
        as the input files.
    plot_pattern : str
        Glob pattern for plot files.
    exe_path : str or Path, optional
        Path to the flatten executable.
    level : int
        AMR level to flatten onto.
    epsg_code : int, optional
        EPSG code for the projection.
    x0 : float, optional
        X-coordinate of the lower-left corner of the domain (metres).
        Defaults to the UKESM standard value for the given ``epsg_code``.
    y0 : float, optional
        Y-coordinate of the lower-left corner of the domain (metres).
        Defaults to the UKESM standard value for the given ``epsg_code``.
    reference_year : int
        Reference year for the CF time axis.
    calendar : str
        CF calendar name: ``"gregorian"`` (default) or ``"360_day"`` for
        UKESM-coupled runs.
    cmip7_only : bool
        If True, only write CMIP7-named variables.
    ice_density : float
        Ice density in kg m-3 for the flotation-based mask fallback (default 918.0).
    water_density : float
        Ocean water density in kg m-3 for the flotation-based mask fallback
        (default 1028.0).
    h_min : float
        Minimum ice thickness threshold in metres for the mask/iceFrac fallback
        (default 1.0 m).
    verbose : bool
        Print progress messages.
    **nc_kwargs
        Passed to :func:`write_cmip7_netcdf`.

    Returns
    -------
    list of Path
        Paths of output NetCDF files written.
    """
    from .diagnostics import find_plot_files

    directory = Path(directory)
    if output_dir is None:
        output_dir = directory
    else:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

    plot_files = find_plot_files(directory, pattern=plot_pattern)
    if verbose:
        print(f"Found {len(plot_files)} plot files in {directory}")

    output_files = []
    for i, pf in enumerate(plot_files):
        out_nc = output_dir / (pf.stem.replace(".2d", "") + "_cmip7.nc")
        if verbose:
            print(f"  [{i+1}/{len(plot_files)}] {pf.name} -> {out_nc.name}")
        process_plotfile(
            pf,
            out_nc,
            exe_path=exe_path,
            level=level,
            epsg_code=epsg_code,
            x0=x0,
            y0=y0,
            reference_year=reference_year,
            calendar=calendar,
            cmip7_only=cmip7_only,
            ice_density=ice_density,
            water_density=water_density,
            h_min=h_min,
            verbose=False,
            **nc_kwargs,
        )
        output_files.append(out_nc)

    if verbose:
        print(f"Done. Wrote {len(output_files)} files.")
    return output_files
