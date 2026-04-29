"""
Utilities for creating CF-compliant NetCDF metadata.

Provides helpers for:
  - Global attribute dictionaries (Conventions, history, institution, etc.)
  - CF time coordinate encoding (BISICLES years -> CF days since reference)
  - Grid mapping (CRS) variable attributes for common BISICLES projections
  - Adding CF coordinate attributes to existing variables
  - Computing 2-D lat/lon auxiliary coordinates from projected grids (requires pyproj)
"""

import datetime
import numpy as np

CF_CONVENTIONS = "CF-1.12 CMIP-7.0"
FILL_VALUE = 1.0e20

# UKESM-BISICLES domain origins (lower-left corner of the model grid, metres).
# These are the x0/y0 values passed to the flatten file tool for each ice sheet.
UKESM_GRID_ORIGINS = {
    3413: {"x0": -654650.0,  "y0": -3385950.0},  # GrIS  (EPSG:3413)
    3031: {"x0": -3072000.0, "y0": -3072000.0},  # AIS   (EPSG:3031)
}


def get_global_attributes(
    institution="",
    source="BISICLES adaptive mesh refinement ice sheet model", # should this be UKESM?
    experiment="",
    variant_label="",
    ice_sheet="",
    references="",
    extra_history="",
    grid_label="gn",
    grid="",
    nominal_resolution="",
    source_files=None,
    **kwargs,
):
    """
    Return a dict of CF/CMIP7 global attributes.

    Parameters
    ----------
    institution : str
        Name of the institution that produced the data.
    source : str
        Model description string.
    experiment : str
        Experiment identifier (e.g. 'historical', 'ssp585').
    variant_label : str
        Variant label (e.g. 'r1i1p1f3').
    ice_sheet : str
        Ice sheet identifier, e.g. 'GrIS' or 'AIS'.
    references : str
        Relevant references or DOIs.
    extra_history : str
        Additional history text prepended before the auto-generated entry.
    grid_label : str
        CMIP grid label.  ``"gn"`` (native grid, default) means the data are
        on the model's native projected grid rather than regridded.
    grid : str
        Human-readable description of the grid, e.g.
        ``"Antarctic Polar Stereographic (EPSG:3031)"``.
    nominal_resolution : str
        Approximate grid spacing as a CMIP string, e.g. ``"5 km"``.
    **kwargs
        Any additional key-value pairs to include as global attributes.
    """
    now = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
    history = f"{now}: Created by bisicles_cmip7_postproc"
    if extra_history:
        history = f"{extra_history}; {history}"

    attrs = {
        "Conventions": CF_CONVENTIONS,
        "source": source,
        "history": history,
        "institution": institution,
        "experiment": experiment,
        "variant_label": variant_label,
        "ice_sheet": ice_sheet,
        "references": references,
        "grid_label": grid_label,
        "grid": grid,
        "nominal_resolution": nominal_resolution,
    }
    if source_files:
        attrs["source_file"] = " ".join(source_files) if isinstance(source_files, list) else source_files
    # Add any extra kwargs
    attrs.update(kwargs)
    # Remove empty strings to keep the output clean
    return {k: v for k, v in attrs.items() if v != ""}


_DAYS_PER_YEAR = {
    "gregorian": 365.25,  # tropical year used by BISICLES internally
    "360_day":   360.0,   # 360-day calendar used by UKESM
}

def years_to_days(time_years, reference_year=1850, calendar="gregorian"):
    """
    Convert BISICLES simulation time (in years) to days since a reference year.

    BISICLES uses a simple year counter.  Two calendars are supported:

    ``"gregorian"``
        365.25 days/year — matches BISICLES internal time (tropical
        year approximation).
    ``"360_day"``
        360 days/year — matches the UKESM calendar. Use this when 
        processing output from UKESM-coupled BISICLES runs.

    Parameters
    ----------
    time_years : float or array-like
        Time in simulation years.
    reference_year : int
        The year of the reference epoch for the time axis (default 1850).
    calendar : str
        CF calendar name: ``"gregorian"`` (default) or ``"360_day"``.

    Returns
    -------
    days : ndarray
        Time in days since reference_year-01-01.
    units : str
        CF-compliant time units string.
    calendar : str
        Calendar name (echoed back for use as a NetCDF attribute).
    """
    if calendar not in _DAYS_PER_YEAR:
        raise ValueError(
            f"Unsupported calendar '{calendar}'. "
            f"Choose from: {list(_DAYS_PER_YEAR)}"
        )
    days_per_year = _DAYS_PER_YEAR[calendar]
    days = (np.asarray(time_years, dtype=float) - reference_year) * days_per_year
    units = f"days since {reference_year:04d}-01-01 00:00:00"
    return days, units, calendar


def get_crs_variable_attrs(epsg_code, x0=None, y0=None):
    """
    Return CF grid_mapping variable attributes for common BISICLES projections.

    Parameters
    ----------
    epsg_code : int
        EPSG code of the coordinate reference system used in the BISICLES run.
    x0 : float, optional
        X coordinate of the lower-left corner of the BISICLES domain (metres).
        Stored as a ``bisicles_domain_x0`` attribute for reference — this is
        the domain origin, *not* the projection false easting (which is always
        0 for EPSG:3031 and EPSG:3413 and must not be changed).
    y0 : float, optional
        Y coordinate of the lower-left corner of the BISICLES domain (metres).
        Stored as a ``bisicles_domain_y0`` attribute for reference.

    Returns
    -------
    dict
        Attributes to set on the grid_mapping variable in the NetCDF file.
        Returns an empty dict if the EPSG code is not recognised.
    """
    _crs_table = {
        # Antarctic Polar Stereographic (standard BISICLES AIS grid)
        # UKESM domain origin: x0=-3072000, y0=-3072000
        3031: {
            "grid_mapping_name": "polar_stereographic",
            "straight_vertical_longitude_from_pole": 0.0,
            "latitude_of_projection_origin": -90.0,
            "standard_parallel": -71.0,
            "false_easting": 0.0,
            "false_northing": 0.0,
            "semi_major_axis": 6378137.0,
            "inverse_flattening": 298.257223563,
            "reference_ellipsoid_name": "WGS84",
            "horizontal_datum_name": "World Geodetic System 1984",
            "prime_meridian_name": "Greenwich",
            "crs_wkt": (
                'PROJCS["WGS 84 / Antarctic Polar Stereographic",'
                'GEOGCS["WGS 84",DATUM["WGS_1984",'
                'SPHEROID["WGS 84",6378137,298.257223563]],'
                'PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],'
                'PROJECTION["Polar_Stereographic"],'
                'PARAMETER["latitude_of_origin",-71],'
                'PARAMETER["central_meridian",0],'
                'PARAMETER["false_easting",0],'
                'PARAMETER["false_northing",0],'
                'UNIT["metre",1]]'
            ),
            "proj4_params": (
                "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 "
                "+k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
            ),
            "epsg_code": "EPSG:3031",
        },
        # NSIDC Sea Ice Polar Stereographic North (standard BISICLES GrIS grid)
        # UKESM domain origin: x0=-654650, y0=-3385950
        3413: {
            "grid_mapping_name": "polar_stereographic",
            "straight_vertical_longitude_from_pole": -45.0,
            "latitude_of_projection_origin": 90.0,
            "standard_parallel": 70.0,
            "false_easting": 0.0,
            "false_northing": 0.0,
            "semi_major_axis": 6378137.0,
            "inverse_flattening": 298.257223563,
            "reference_ellipsoid_name": "WGS84",
            "horizontal_datum_name": "World Geodetic System 1984",
            "prime_meridian_name": "Greenwich",
            "crs_wkt": (
                'PROJCS["WGS 84 / NSIDC Sea Ice Polar Stereographic North",'
                'GEOGCS["WGS 84",DATUM["WGS_1984",'
                'SPHEROID["WGS 84",6378137,298.257223563]],'
                'PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],'
                'PROJECTION["Polar_Stereographic"],'
                'PARAMETER["latitude_of_origin",70],'
                'PARAMETER["central_meridian",-45],'
                'PARAMETER["false_easting",0],'
                'PARAMETER["false_northing",0],'
                'UNIT["metre",1]]'
            ),
            "proj4_params": (
                "+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 "
                "+k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
            ),
            "epsg_code": "EPSG:3413",
        },
        # WGS84 geographic (used for some regional BISICLES setups)
        4326: {
            "grid_mapping_name": "latitude_longitude",
            "longitude_of_prime_meridian": 0.0,
            "semi_major_axis": 6378137.0,
            "inverse_flattening": 298.257223563,
            "epsg_code": "EPSG:4326",
        },
    }
    attrs = dict(_crs_table.get(epsg_code, {"epsg_code": f"EPSG:{epsg_code}"}))
    # Record the BISICLES domain origin as informational attributes.
    # Note: these are NOT the projection false_easting/false_northing — the
    # proj4/WKT parameters correctly remain 0 for EPSG:3031 and EPSG:3413.
    # The domain origin is already encoded in the x/y coordinate arrays.
    if x0 is not None:
        attrs["bisicles_domain_x0"] = float(x0)
    if y0 is not None:
        attrs["bisicles_domain_y0"] = float(y0)
    return attrs


def add_time_variable(ds, time_years, reference_year=1850, calendar="gregorian"):
    """
    Add a CF-compliant time variable to an open NetCDF4 Dataset.

    Parameters
    ----------
    ds : netCDF4.Dataset
        Open dataset to write into (must already have a 'time' dimension).
    time_years : array-like
        Time values in simulation years.
    reference_year : int
        Reference epoch year.
    calendar : str
        CF calendar name: ``"gregorian"`` (default) or ``"360_day"``.

    Returns
    -------
    time_var : netCDF4.Variable
        The created time variable.
    """
    days, units, calendar = years_to_days(time_years, reference_year, calendar)
    time_var = ds.createVariable("time", "f8", ("time",))
    time_var[:] = days
    time_var.standard_name = "time"
    time_var.long_name = "time"
    time_var.units = units
    time_var.calendar = calendar
    time_var.axis = "T"
    return time_var


def add_time_bounds(ds, start_years, end_years, reference_year=1850, calendar="gregorian"):
    """
    Add a ``time_bnds`` variable to an open NetCDF4 Dataset.

    This is required by CF-1.12 for time-averaged data.  The ``time``
    variable's ``bounds`` attribute is also set to ``'time_bnds'``.

    Parameters
    ----------
    ds : netCDF4.Dataset
        Open dataset that already has a ``time`` dimension and variable.
    start_years : array-like
        Start of each averaging period in simulation years.
    end_years : array-like
        End of each averaging period in simulation years.
    reference_year : int
        Reference epoch year (must match the ``time`` variable).
    calendar : str
        CF calendar name: ``"gregorian"`` (default) or ``"360_day"``.
        Must match the calendar used when the ``time`` variable was created.
    """
    import numpy as np

    start_days, _, _ = years_to_days(np.asarray(start_years, dtype=float), reference_year, calendar)
    end_days, _, _ = years_to_days(np.asarray(end_years, dtype=float), reference_year, calendar)

    if "bnds" not in ds.dimensions:
        ds.createDimension("bnds", 2)

    bnds_var = ds.createVariable("time_bnds", "f8", ("time", "bnds"))
    bnds_var[:, 0] = start_days
    bnds_var[:, 1] = end_days

    # Link the time variable to its bounds
    ds.variables["time"].bounds = "time_bnds"
    return bnds_var


def add_xy_variables(ds, x_values, y_values, epsg_code=None):
    """
    Add CF-compliant x/y coordinate variables to an open NetCDF4 Dataset.

    Parameters
    ----------
    ds : netCDF4.Dataset
        Open dataset (must have 'x' and 'y' dimensions).
    x_values : array-like
        X coordinate values in metres.
    y_values : array-like
        Y coordinate values in metres.
    epsg_code : int or None
        EPSG code used to determine axis names. For geographic coordinates
        (EPSG:4326) the standard names are longitude/latitude; for projected
        coordinates they are projection_x_coordinate / projection_y_coordinate.
    """
    if epsg_code == 4326:
        x_std = "longitude"
        y_std = "latitude"
        x_long = "Longitude"
        y_long = "Latitude"
        x_units = "degrees_east"
        y_units = "degrees_north"
    else:
        x_std = "projection_x_coordinate"
        y_std = "projection_y_coordinate"
        x_long = "X Coordinate of Projection"
        y_long = "Y Coordinate of Projection"
        x_units = "m"
        y_units = "m"

    x_var = ds.createVariable("x", "f8", ("x",))
    x_var[:] = x_values
    x_var.standard_name = x_std
    x_var.long_name = x_long
    x_var.units = x_units
    x_var.axis = "X"

    y_var = ds.createVariable("y", "f8", ("y",))
    y_var[:] = y_values
    y_var.standard_name = y_std
    y_var.long_name = y_long
    y_var.units = y_units
    y_var.axis = "Y"

    return x_var, y_var


def add_crs_variable(ds, epsg_code, x0=None, y0=None):
    """
    Add a scalar grid_mapping variable for the coordinate reference system.

    Parameters
    ----------
    ds : netCDF4.Dataset
        Open dataset.
    epsg_code : int
        EPSG code of the projection.
    x0 : float, optional
        X coordinate of the lower-left corner of the BISICLES domain (metres).
        Recorded as a ``bisicles_domain_x0`` attribute for reference.
    y0 : float, optional
        Y coordinate of the lower-left corner of the BISICLES domain (metres).
        Recorded as a ``bisicles_domain_y0`` attribute for reference.

    Returns
    -------
    crs_var : netCDF4.Variable or None
        The created CRS variable, or None if no attributes are available.
    """
    attrs = get_crs_variable_attrs(epsg_code, x0=x0, y0=y0)
    if not attrs:
        return None
    crs_var = ds.createVariable("crs", "i4")
    crs_var[:] = 1
    crs_var.setncatts(attrs)
    return crs_var


def compute_latlon_arrays(x, y, epsg_code):
    """
    Compute 2-D latitude and longitude arrays from projected x/y coordinates.

    Requires the ``pyproj`` package (``pip install pyproj`` or
    ``pip install bisicles-cmip7-postproc[geo]``).

    Parameters
    ----------
    x : array-like, shape (nx,)
        Projected x coordinates in metres.
    y : array-like, shape (ny,)
        Projected y coordinates in metres.
    epsg_code : int
        EPSG code of the source projection (e.g. 3031 or 3413).

    Returns
    -------
    lat_2d : ndarray, shape (ny, nx)
        Latitude values in degrees north.
    lon_2d : ndarray, shape (ny, nx)
        Longitude values in degrees east.

    Raises
    ------
    ImportError
        If ``pyproj`` is not installed.
    """
    try:
        from pyproj import Transformer
    except ImportError:
        raise ImportError(
            "pyproj is required to compute lat/lon coordinate arrays. "
            "Install it with:\n"
            "    pip install pyproj\n"
            "or:\n"
            "    pip install 'bisicles-cmip7-postproc[geo]'"
        )
    transformer = Transformer.from_crs(
        f"EPSG:{epsg_code}", "EPSG:4326", always_xy=True
    )
    xx, yy = np.meshgrid(x, y)
    lon_2d, lat_2d = transformer.transform(xx, yy)
    return lat_2d, lon_2d


def add_latlon_variables(ds, lat_2d, lon_2d):
    """
    Add CF-compliant 2-D ``lat(y, x)`` and ``lon(y, x)`` auxiliary coordinate
    variables to an open NetCDF4 Dataset.

    These are *auxiliary* coordinate variables (not dimension coordinates) and
    must be listed in the ``coordinates`` attribute of data variables, as
    required by CF-1.12 for projected grids.

    Parameters
    ----------
    ds : netCDF4.Dataset
        Open dataset (must already have ``y`` and ``x`` dimensions).
    lat_2d : ndarray, shape (ny, nx)
        Latitude values in degrees north.
    lon_2d : ndarray, shape (ny, nx)
        Longitude values in degrees east.

    Returns
    -------
    lat_var, lon_var : netCDF4.Variable
        The created latitude and longitude variables.
    """
    lat_var = ds.createVariable("lat", "f8", ("y", "x"))
    lat_var[:] = lat_2d
    lat_var.standard_name = "latitude"
    lat_var.long_name = "Latitude"
    lat_var.units = "degrees_north"

    lon_var = ds.createVariable("lon", "f8", ("y", "x"))
    lon_var[:] = lon_2d
    lon_var.standard_name = "longitude"
    lon_var.long_name = "Longitude"
    lon_var.units = "degrees_east"

    return lat_var, lon_var


