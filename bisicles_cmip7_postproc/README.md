# bisicles_cmip7_postproc

Post-processing of BISICLES ice sheet model output into CF-1.12 /
CMIP7-compliant NetCDF files.

Two independent workflows are provided:

  flatten      -- 2D spatial fields on a uniform grid, one file per variable
  diagnostics  -- Scalar timeseries (volume, area, fluxes, ...), one per variable

Both workflows accept a single BISICLES HDF5 plot file or a directory of plot
files. In directory mode all timesteps are combined onto a single time axis.

Commands installed after setup:

  bike-cmip7-postproc-flatten
  bike-cmip7-postproc-diagnostics
  bike-cmip7-postproc-run

## Installation

Requires Python >= 3.8, numpy, and netCDF4.

```bash
# From the repo root — editable mode means source changes take effect immediately
pip install -e bisicles_cmip7_postproc/

# With optional extras (recommended)
pip install -e "bisicles_cmip7_postproc/[yaml,geo]"
```

Optional extras:

  yaml  -- adds pyyaml; needed for YAML config files (JSON works without it)
  geo   -- adds pyproj; needed for 2-D lat/lon coordinates (required for full
           CF/CMOR compliance)

## Prerequisites

Both workflows call external BISICLES executables that must be compiled
separately from the bisicles-uob source tree:

```text
flatten      -- flatten2d.PLATFORM.ex
diagnostics  -- diagnostics2d.PLATFORM.ex
```

The path to the relevant executable must be supplied via `--exe-path` on the
command line or `exe_path` in a config file.

## Quick start

### Flatten workflow — 2D spatial fields

Single file:

```bash
bike-cmip7-postproc-flatten \
    --input  /data/run/bisicles_dx030c_1y_18880101-18890101_plot-AIS.hdf5 \
    --output-dir /data/postproc/ \
    --exe-path /opt/bisicles/flatten2d.Linux.64.g++.gfortran.OPT.ex \
    --epsg 3031 \
    --ice-sheet AIS \
    --calendar 360_day \
    --experiment historical \
    --variant-label r1i1p1f3
```

Directory of plot files (all timesteps combined into one output per variable):

```bash
bike-cmip7-postproc-flatten \
    --input  /data/run/ \
    --output-dir /data/postproc/ \
    --exe-path /opt/bisicles/flatten2d.Linux.64.g++.gfortran.OPT.ex \
    --epsg 3031 \
    --level 2 \
    --ice-sheet AIS \
    --calendar 360_day \
    --experiment historical \
    --variant-label r1i1p1f3 \
    --institution "University of Bristol"
```

Output files are named `{cmip7_name}.nc`, one per variable.

### Diagnostics workflow — scalar timeseries

```bash
bike-cmip7-postproc-diagnostics \
    --input  /data/run/ \
    --output-dir /data/postproc/ \
    --exe-path /opt/bisicles/diagnostics2d.Linux.64.g++.gfortran.OPT.ex \
    --ice-sheet AIS \
    --calendar 360_day \
    --experiment historical \
    --variant-label r1i1p1f3
```

Output files are named `{cmip7_name}.nc`, or `{cmip7_name}_mask{N}.nc` for
regional mask variants.

## Config file usage

All options can be supplied via a YAML (or JSON) config file instead of on the
command line. This is the recommended approach for batch processing.

```bash
# Run from a config file
bike-cmip7-postproc-run my_config.yaml

# Config file with command-line overrides
bike-cmip7-postproc-run my_config.yaml --set calendar 360_day --set level 2

# Pass a config file to the workflow-specific commands
bike-cmip7-postproc-flatten     --config my_flatten_config.yaml
bike-cmip7-postproc-diagnostics --config my_diag_config.yaml
```

Template config files with all options documented are in the `config/` directory:

  config/flatten_defaults.yaml      -- flatten workflow
  config/diagnostics_defaults.yaml  -- diagnostics workflow

Copy one of these, fill in the required fields, and run.

## Configuration reference

Options marked [required] have no default and must always be supplied.
Options marked [UKESM] are most relevant for UKESM-coupled runs.

All CLI flags have a matching config-file key: replace hyphens with underscores
and drop the leading `--` (e.g. `--output-dir` becomes `output_dir`).

### Shared options (both workflows)

**--input** / `input`  [required]
    Path to a BISICLES HDF5 plot file, or a directory of plot files.

**--output-dir** / `output_dir`  [default: same directory as input]
    Directory where per-variable NetCDF output files are written.

**--plot-pattern** / `plot_pattern`  [default: "plot.*.2d.hdf5"]
    Glob pattern used to find plot files when --input is a directory.

**--exe-path** / `exe_path`  [required]
    Full path to the BISICLES flatten or diagnostics executable.

**--reference-year** / `reference_year`  [default: 1850]
    Reference year for the CF time axis (days since YYYY-01-01).

**--calendar** / `calendar`  [default: gregorian]  [UKESM]
    CF calendar for converting simulation years to days.
      gregorian  -- 365.25 days/year (standard BISICLES)
      360_day    -- 360 days/year (UKESM-coupled runs)

**--ice-sheet** / `ice_sheet`  [default: auto-detected from filename]
    Ice sheet identifier written to NetCDF global attributes.
      GrIS  -- Greenland Ice Sheet
      AIS   -- Antarctic Ice Sheet
    Parsed automatically from UKESM-style filenames when not set.

**--institution** / `institution`  [default: ""]
    Institution name written to NetCDF global attributes.

**--source** / `source`  [default: "BISICLES adaptive mesh refinement..."]
    Model description string written to NetCDF global attributes.

**--experiment** / `experiment`  [default: ""]
    Experiment identifier, e.g. "historical" or "ssp585".

**--variant-label** / `variant_label`  [default: ""]
    CMIP variant label, e.g. "r1i1p1f3".

`extra_attrs`  (config file only)  [default: null]
    Dict of additional key-value pairs written as NetCDF global attributes.

```yaml
extra_attrs:
  contact: your.name@institution.ac.uk
  references: doi:10.5194/gmd-xx-xxxx-xxxx
```

**--quiet** / `verbose: false`
    Suppress progress messages to stdout.

### Flatten-only options

**--epsg** / `epsg`  [default: read from HDF5 file]
    EPSG code for the output coordinate reference system.
      3413  -- NSIDC Polar Stereographic North (GrIS)
      3031  -- Antarctic Polar Stereographic (AIS)

**--level** / `level`  [default: 0]
    AMR refinement level to flatten onto. 0 = coarsest grid; higher values
    give finer resolution. Must be >= 0.

**--x0** / `x0`  [default: UKESM standard for the given EPSG]
**--y0** / `y0`  [default: UKESM standard for the given EPSG]
    X/Y coordinates of the lower-left corner of the domain in metres.
    When not set and epsg is a known UKESM grid, the standard origin is used:

```text
Ice sheet  EPSG   x0            y0
---------  ----   ----------    -----------
GrIS       3413   -654650.0     -3385950.0
AIS        3031   -3072000.0    -3072000.0
```

**--cmip7-only** / `cmip7_only`  [default: false]
    Only write CMIP7-standard variables. When false, unmapped BISICLES
    fields are also written with their original names.

**--keep-intermediate** / `keep_intermediate`  [default: false]
    Keep the raw flatten NetCDF written by the executable (single-file mode
    only; saved alongside the output as `*_flatten_raw.nc`).

**--ice-density** / `ice_density`  [default: 918.0 kg/m3]
**--water-density** / `water_density`  [default: 1028.0 kg/m3]
**--h-min** / `h_min`  [default: 1.0 m]
    Physical constants used for the flotation mask fallback. Applied when the
    plot file contains no explicit mask or iceFrac variable; the grounded/
    floating partition is derived from the flotation criterion using ice
    thickness and bed topography.

### Diagnostics-only options

**--ice-density** / `ice_density`  [default: 918.0 kg/m3]
**--water-density** / `water_density`  [default: 1028.0 kg/m3]
**--gravity** / `gravity`  [default: 9.81 m/s2]
**--h-min** / `h_min`  [default: 1.0 m]
    Physical constants passed to the diagnostics executable.

**--mask-file** / `mask_file`  [default: null]
    Path to an HDF5 mask file for regional diagnostics (e.g. drainage basins).
    When not set only whole-ice-sheet diagnostics are produced.

**--mask-no-start** / `mask_no_start`  [default: 0]
**--mask-no-end** / `mask_no_end`  [default: 0]
    First and last mask region indices to process (inclusive). Only used when
    mask_file is set.

## Example config file (flatten)

```yaml
tool: flatten

input: /data/run/
output_dir: /data/postproc/

exe_path: /opt/bisicles/flatten2d.Linux.64.g++.gfortran.OPT.ex

epsg: 3031
level: 2
calendar: 360_day

institution: "University of Bristol"
experiment: historical
variant_label: r1i1p1f3
ice_sheet: AIS

extra_attrs:
  contact: your.name@institution.ac.uk
  references: doi:10.5194/gmd-xx-xxxx-xxxx
```

## Output file structure

Each output NetCDF file contains:

- A single data variable with a `time` dimension (and `y`, `x` for 2D fields).
- A `time` coordinate in CF `days since YYYY-01-01` units, with `time_bnds`
  added for time-mean files.
- `x` and `y` projected coordinates, plus 2-D `lat`/`lon` auxiliary coordinates
  (flatten workflow only; requires pyproj).
- A `crs` grid-mapping variable with full projection metadata (flatten only).
- Global attributes: `Conventions`, `history`, `institution`, `experiment`,
  `variant_label`, `ice_sheet`, and `source_file` (the basename(s) of the
  input HDF5 file(s) that produced this output).

## Standalone BISICLES / ISMIP7

Three additional commands are provided for postprocessing **standalone** BISICLES
simulations (not coupled to UKESM) and producing ISMIP7-compliant output:

  bike-ismip7-postproc-flatten
  bike-ismip7-postproc-diagnostics
  bike-ismip7-postproc-run

The key differences from the CMIP7-coupled workflow are:

- **Output filenames** follow the ISMIP7 DRS convention:
  `{varname}_{icesheet}_{exp}_{model}_{member}_{freq}_{startyr}-{endyr}.nc`
  e.g. `lithk_AIS_ctrl_proj_BISICLES_r1_yr_2015-2100.nc`
- **`Conventions`** global attribute is set to `"CF-1.12"` (not `"CF-1.12 CMIP-7.0"`).
- **`model_id`** and **`member_id`** are written as global attributes.
- **Simulation time** is read directly from the BISICLES HDF5 internal time
  counter (not from the filename), which is the authoritative source for
  standalone runs.
- **Calendar** defaults to `gregorian` (365.25 days/year).
- **Variable names and metadata** are identical to the CMIP7 workflow.

### ISMIP7 flatten example

```bash
bike-ismip7-postproc-flatten \
    --input  /data/standalone/output/ \
    --output-dir /data/postproc/ \
    --exe-path /opt/bisicles/flatten2d.Linux.64.g++.gfortran.OPT.ex \
    --epsg 3031 \
    --level 0 \
    --ice-sheet AIS \
    --experiment ctrl_proj \
    --model-id BISICLES \
    --member-id r1 \
    --frequency yr \
    --reference-year 1950
```

### ISMIP7 diagnostics example

```bash
bike-ismip7-postproc-diagnostics \
    --input  /data/standalone/output/ \
    --output-dir /data/postproc/ \
    --exe-path /opt/bisicles/diagnostics2d.Linux.64.g++.gfortran.OPT.ex \
    --ice-sheet AIS \
    --experiment ctrl_proj \
    --model-id BISICLES \
    --member-id r1 \
    --frequency yr \
    --reference-year 1950
```

### ISMIP7 config file

Config template files are in the `config/` directory:

  config/ismip7_flatten_defaults.yaml      -- flatten workflow
  config/ismip7_diagnostics_defaults.yaml  -- diagnostics workflow

```bash
bike-ismip7-postproc-run my_ismip7_flatten.yaml
```

The config file uses `tool: flatten` or `tool: diagnostics` as usual; the
`bike-ismip7-postproc-run` command automatically sets `ismip7_mode: true`.
The ISMIP7 mode can also be enabled explicitly in any config file by adding
`ismip7_mode: true`, which then works with `bike-cmip7-postproc-run` too.

### New CLI options (ISMIP7 commands only)

**--model-id** / `model_id`  [default: "BISICLES"]
    ISMIP7 model identifier written in the DRS filename and as a global attribute.

**--member-id** / `member_id`  [default: "r1"]
    ISMIP7 member identifier written in the DRS filename and as a global attribute.

### ISMIP7 DRS filename components

| Component  | Source                        | Example        |
|------------|-------------------------------|----------------|
| varname    | CMIP7/ISMIP7 variable name    | `lithk`        |
| icesheet   | `--ice-sheet`                 | `AIS`          |
| exp        | `--experiment`                | `ctrl_proj`    |
| model      | `--model-id`                  | `BISICLES`     |
| member     | `--member-id`                 | `r1`           |
| freq       | `--frequency` (default `yr`)  | `yr`           |
| startyr    | Minimum year in time axis     | `2015`         |
| endyr      | Maximum year in time axis     | `2100`         |

## UKESM-specific notes

- Always use `--calendar 360_day` for UKESM-coupled runs. BISICLES internal
  time is always 0 in UKESM coupling; simulation time is read from the
  filename instead.

- UKESM BISICLES filenames follow the pattern
  `bisicles_SUITE_FREQ_START-END_plot-ICESHEET.hdf5`,
  e.g. `bisicles_dx030c_1y_18880101-18890101_plot-AIS.hdf5`.
  The ice sheet, simulation time, and time-mean status are all parsed
  automatically from this pattern.

- The `ice_sheet` global attribute is set automatically from the filename
  when `--ice-sheet` is not provided.
