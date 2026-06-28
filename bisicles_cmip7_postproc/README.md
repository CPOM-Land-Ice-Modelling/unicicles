# bisicles_cmip7_postproc

Post-processing of BISICLES ice sheet model output into CF-1.12-compliant
NetCDF files.

This tool supports two target workflows:

- **UKESM-coupled runs → CMIP7 submission** using the `bike-cmip7-postproc-*`
  commands.
- **Standalone BISICLES runs → ISMIP7 submission** using the
  `bike-ismip7-postproc-*` commands.

Two independent processing workflows are provided for both cases:

| Workflow      | Input file type                            | Output                                      |
|---------------|--------------------------------------------|---------------------------------------------|
| `flatten`     | Any BISICLES HDF5 plot file                | 2D spatial fields on a uniform grid         |
| `diagnostics` | Instantaneous snapshot plot files only     | Scalar timeseries (volume, area, fluxes, …) |

---

## Contents

1. [Installation](#installation)
2. [Prerequisites](#prerequisites)
3. [Quick start — config file](#quick-start--config-file)
4. [Quick start — command line](#quick-start--command-line)
5. [Configuration reference](#configuration-reference)
6. [Output file structure](#output-file-structure)
7. [Standalone BISICLES / ISMIP7](#standalone-bisicles--ismip7)
8. [UKESM-specific notes](#ukesm-specific-notes)

---

## Installation

Requires Python >= 3.8, numpy, and netCDF4.

```bash
# From the repo root — editable mode means source changes take effect immediately
pip install -e bisicles_cmip7_postproc/

# With optional extras (recommended)
pip install -e "bisicles_cmip7_postproc/[yaml,geo]"

# With regridding support (needed for --x-min/--x-max/etc.)
pip install -e "bisicles_cmip7_postproc/[yaml,geo,regrid]"
```

Optional extras:

| Extra    | Adds    | Needed for                                                       |
|----------|---------|------------------------------------------------------------------|
| `yaml`   | pyyaml  | YAML config files (JSON config works without it)                |
| `geo`    | pyproj  | 2-D lat/lon coordinates (required for full CF/CMOR compliance)   |
| `regrid` | scipy   | Regridding onto a target grid (`x_min`/`x_max`/`nx`/… options)  |

Commands installed after setup:

```text
# UKESM-coupled / CMIP7
bike-cmip7-postproc-flatten
bike-cmip7-postproc-diagnostics
bike-cmip7-postproc-run

# Standalone / ISMIP7
bike-ismip7-postproc-flatten
bike-ismip7-postproc-diagnostics
bike-ismip7-postproc-run
```

---

## Prerequisites

Both workflows call external BISICLES executables that must be compiled
separately from the bisicles-uob source tree:

```text
flatten      -- flatten2d.PLATFORM.ex
diagnostics  -- diagnostics2d.PLATFORM.ex
```

The path to the relevant executable must be supplied via `exe_path` in a
config file, or `--exe-path` on the command line.

---

## Quick start — config file

Config files are the recommended way to run these tools.  They keep all
settings in one place, are easy to re-run and share, and avoid long
command-line strings.

### Flatten workflow (2D spatial fields)

Copy the template and fill in the required fields:

```bash
cp config/flatten_defaults.yaml my_flatten.yaml
# Edit my_flatten.yaml, then:
bike-cmip7-postproc-run my_flatten.yaml
```

Minimal working example — `my_flatten.yaml`:

```yaml
tool: flatten

input:       /data/run/
output_dir:  /data/postproc/
exe_path:    /opt/bisicles/flatten2d.Linux.64.g++.gfortran.OPT.ex

epsg:        3031
level:       2
calendar:    360_day

ice_sheet:   AIS
institution: "University of Bristol"
experiment:  historical
variant_label: r1i1p1f3
```

### Diagnostics workflow (scalar timeseries)

> **Note:** time-mean CF output files (`plot.CF-*.hdf5`) are not supported
> by the diagnostics workflow — use `flatten` instead.  See
> [UKESM-specific notes](#ukesm-specific-notes) for details.

```bash
cp config/diagnostics_defaults.yaml my_diagnostics.yaml
# Edit my_diagnostics.yaml, then:
bike-cmip7-postproc-run my_diagnostics.yaml
```

Minimal working example — `my_diagnostics.yaml`:

```yaml
tool: diagnostics

input:       /data/run/
output_dir:  /data/postproc/
exe_path:    /opt/bisicles/diagnostics2d.Linux.64.g++.gfortran.OPT.ex

calendar:    360_day
ice_sheet:   AIS
experiment:  historical
variant_label: r1i1p1f3
```

### Overriding individual values from the command line

Any config-file value can be overridden on the command line without editing
the file:

```bash
bike-cmip7-postproc-run my_flatten.yaml --set calendar 360_day --set level 2
```

### Passing a config file to a workflow-specific command

The workflow-specific commands also accept a `--config` flag if you prefer
not to use `bike-cmip7-postproc-run`:

```bash
bike-cmip7-postproc-flatten     --config my_flatten.yaml
bike-cmip7-postproc-diagnostics --config my_diagnostics.yaml
```

Template config files with all options documented are in the `config/`
directory:

```text
config/flatten_defaults.yaml      -- flatten workflow
config/diagnostics_defaults.yaml  -- diagnostics workflow
```

---

## Quick start — command line

For one-off runs or scripting, all options can be passed directly on the
command line.  Config files are generally easier for repeated use.

### Flatten — single file

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

### Flatten — directory of plot files

All timesteps are combined into one output file per variable:

```bash
bike-cmip7-postproc-flatten \
    --input  /data/run/ \
    --output-dir /data/postproc/ \
    --exe-path /opt/bisicles/flatten2d.Linux.64.g++.gfortran.OPT.ex \
    --epsg 3031 --level 2 \
    --ice-sheet AIS \
    --calendar 360_day \
    --experiment historical \
    --variant-label r1i1p1f3 \
    --institution "University of Bristol"
```

### Diagnostics — directory of snapshot files

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

---

## Configuration reference

Options marked **[required]** have no default and must always be supplied.
Options marked **[UKESM]** are most relevant for UKESM-coupled runs.

All CLI flags have a matching config-file key: replace hyphens with
underscores and drop the leading `--`
(e.g. `--output-dir` becomes `output_dir`).

### Shared options (both workflows)

**--input** / `input`  [required]
  Path to a BISICLES HDF5 plot file, or a directory of plot files.

**--output-dir** / `output_dir`  [default: same directory as input]
  Directory where per-variable NetCDF output files are written.

**--plot-pattern** / `plot_pattern`  [default: `"plot.*.2d.hdf5"`]
  Glob pattern used to find plot files when `--input` is a directory.

**--exe-path** / `exe_path`  [required]
  Full path to the BISICLES flatten or diagnostics executable.

**--reference-year** / `reference_year`  [default: 1850]
  Reference year for the CF time axis (`days since YYYY-01-01`).

**--calendar** / `calendar`  [default: gregorian]  [UKESM]
  CF calendar for converting simulation years to days.
  - `gregorian` / `standard` — 365.25 days/year (standalone BISICLES / CMIP7)
  - `360_day` — 360 days/year (UKESM-coupled runs)

  `standard` is accepted as a synonym for `gregorian`.  The ISMIP7 commands
  (`bike-ismip7-*`) default to `standard`; both produce identical data but
  ISMIP7 requires the attribute to be written as `"standard"`.

**--ice-sheet** / `ice_sheet`  [default: auto-detected from filename]
  Ice sheet identifier written to NetCDF global attributes.
  - `GrIS` — Greenland Ice Sheet
  - `AIS`  — Antarctic Ice Sheet

  Parsed automatically from UKESM-style filenames when not set.

**--institution** / `institution`  [default: `""`]
  Institution name written to NetCDF global attributes.

**--source** / `source`  [default: `"BISICLES adaptive mesh refinement…"`]
  Model description string written to NetCDF global attributes.

**--experiment** / `experiment`  [default: `""`]
  Experiment identifier, e.g. `"historical"` or `"ssp585"`.

**--variant-label** / `variant_label`  [default: `""`]
  CMIP variant label, e.g. `"r1i1p1f3"`.

`extra_attrs`  (config file only)  [default: null]
  Dict of additional key-value pairs written as NetCDF global attributes:

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
  - `3413` — NSIDC Polar Stereographic North (GrIS)
  - `3031` — Antarctic Polar Stereographic (AIS)

**--level** / `level`  [default: 0]
  AMR refinement level to flatten onto. 0 = coarsest grid; higher values
  give finer resolution.

**--x0** / `x0`  and  **--y0** / `y0`  [default: see below]
  X/Y coordinates of the lower-left corner of the domain in metres.
  These must match the grid origin the BISICLES simulation was run with.

  Typical BISICLES runs (UKESM-coupled and standalone) use the standard
  BISICLES grid origins, which are the defaults when `x0`/`y0` are not set:

  ```text
  Ice sheet  EPSG   x0 (standard BISICLES)  y0 (standard BISICLES)
  ---------  ----   ----------------------  ----------------------
  GrIS       3413   -654650.0               -3385950.0
  AIS        3031   -3072000.0              -3072000.0
  ```

  Note: the ISMIP7 standard grid uses different origins.  Do **not** pass
  ISMIP7 origins here unless your BISICLES run was configured to use them —
  doing so would shift all data to the wrong geographic location.

**--x-min** / `x_min`,  **--x-max** / `x_max`,  **--nx** / `nx`  [default: null]
**--y-min** / `y_min`,  **--y-max** / `y_max`,  **--ny** / `ny`  [default: null]
  Regrid output onto a regular nx × ny grid spanning x_min..x_max and
  y_min..y_max (in the projection given by `--epsg`).  All six must be
  supplied together, or none at all.  Requires `scipy`
  (`pip install "bisicles_cmip7_postproc[regrid]"`).

  When not set, output is written on the standard BISICLES grid (no
  regridding).  To produce output on the ISMIP7 standard grid, use these
  parameters rather than running a separate regridding step:

  ```text
  Domain  EPSG   x_min         x_max        nx    y_min         y_max        ny
  ------  ----   -----------   ----------   ---   -----------   ----------   ---
  AIS     3031   -3040000.0    3040000.0    761   -3040000.0    3040000.0    761
  GrIS    3413   -720000.0     960000.0     1681  -3450000.0    -570000.0    2881
  ```

  Note: GrIS parameters should be verified against ISMIP7 documentation
  before use.

**--cmip7-only** / `cmip7_only`  [default: false]
  Only write CMIP7-standard variables.  When false, unmapped BISICLES
  fields are also written with their original names.

**--keep-intermediate** / `keep_intermediate`  [default: false]
  Keep the raw flatten NetCDF written by the executable (single-file mode
  only; saved alongside the output as `*_flatten_raw.nc`).

**--ice-density** / `ice_density`  [default: 918.0 kg m⁻³]
**--water-density** / `water_density`  [default: 1028.0 kg m⁻³]
**--h-min** / `h_min`  [default: 1.0 m]
  Physical constants for the flotation mask fallback.  Applied when the
  plot file contains no explicit mask or `iceFrac` variable; the grounded/
  floating partition is derived from the flotation criterion using ice
  thickness and bed topography.

### Diagnostics-only options

**--ice-density** / `ice_density`  [default: 918.0 kg m⁻³]
**--water-density** / `water_density`  [default: 1028.0 kg m⁻³]
**--gravity** / `gravity`  [default: 9.81 m s⁻²]
**--h-min** / `h_min`  [default: 1.0 m]
  Physical constants passed to the diagnostics executable.

**--mask-file** / `mask_file`  [default: null]
  Path to an HDF5 mask file for regional diagnostics (e.g. drainage
  basins).  When not set only whole-ice-sheet diagnostics are produced.

**--mask-no-start** / `mask_no_start`  [default: 0]
**--mask-no-end** / `mask_no_end`  [default: 0]
  First and last mask region indices to process (inclusive).  Only used
  when `mask_file` is set.

---

## Output file structure

Each output NetCDF file contains:

- A single data variable (2D spatial fields written as **f4**; scalar
  timeseries written as **f8** for CMIP7/UKESM, **f4** for ISMIP7).
- A `time` dimension (unlimited record dimension) with a CF time coordinate
  (`days since YYYY-01-01`).  Time coordinate written as **f8** for
  CMIP7/UKESM, **f4** for ISMIP7.  Time-mean files include a `time_bnds`
  variable.
- `x` and `y` projected coordinates (1-D), plus 2-D `lat`/`lon` auxiliary
  coordinates (flatten workflow only; requires pyproj).
- A `crs` grid-mapping variable with full projection metadata (flatten only).
- A `_FillValue` of `1.0e20` for CMIP7/UKESM output, or
  `9.969209968386869e+36` (netCDF4 default for f4) for ISMIP7 output.
- Global attributes: `Conventions`, `history`, `institution`, `experiment`,
  `variant_label`, `ice_sheet`, and `source_file` (the basename(s) of the
  input HDF5 file(s) that produced this output).

---

## Standalone BISICLES / ISMIP7

Three additional commands are provided for postprocessing **standalone**
BISICLES simulations (not coupled to UKESM) and producing ISMIP7-compliant
output:

```text
bike-ismip7-postproc-flatten
bike-ismip7-postproc-diagnostics
bike-ismip7-postproc-run
```

The key differences from the CMIP7-coupled workflow are:

- **Output filenames** follow the full ISMIP7 DRS convention (see below).
- **`Conventions`** global attribute is set to `"CF-1.12"`.
- **ISMIP7 mandatory global attributes** (`group`, `model`, `contact_name`,
  `contact_email`, `crs`) are written automatically.
- **Calendar** defaults to `standard` (365.25 days/year).
- **Reference year** must be 1850 (`time:units = "days since 1850-01-01"`).

The same restriction on input file types applies: the diagnostics workflow
requires instantaneous snapshot files; CF time-mean output files must be
processed with the flatten workflow.

### ISMIP7 DRS filename format

```text
{variable}_{domain_id}_{source_id}_{ism_id}_{ism_member_id}_{esm_id}_{forcing_member_id}_{experiment}_{set_counter}_{startyr}-{endyr}.nc
```

Example:

```text
lithk_AIS_BristolGlaciology_BISICLES_m001_standalone_f001_ctrl_proj_C001_2015-2100.nc
```

| Component           | Config key          | Example              |
|---------------------|---------------------|----------------------|
| `variable`          | (from variable table) | `lithk`            |
| `domain_id`         | `ice_sheet`         | `AIS`                |
| `source_id`         | `source_id`         | `BristolGlaciology`  |
| `ism_id`            | `ism_id`            | `BISICLES`           |
| `ism_member_id`     | `ism_member_id`     | `m001`               |
| `esm_id`            | `esm_id`            | `standalone`         |
| `forcing_member_id` | `forcing_member_id` | `f001`               |
| `experiment`        | `experiment`        | `ctrl_proj`          |
| `set_counter`       | `set_counter`       | `C001`               |
| `startyr-endyr`     | (from time axis)    | `2015-2100`          |

Note: frequency is **not** part of the ISMIP7 ISM submission filename.

### Quick start — ISMIP7 config file

Config templates are in the `config/` directory:

```text
config/ismip7_flatten_defaults.yaml      -- flatten workflow
config/ismip7_diagnostics_defaults.yaml  -- diagnostics workflow
```

```bash
bike-ismip7-postproc-run my_ismip7_flatten.yaml
```

The config file uses `tool: flatten` or `tool: diagnostics`; the
`bike-ismip7-postproc-run` command sets `ismip7_mode: true` automatically.
`ismip7_mode: true` can also be added to any config file and used with
`bike-cmip7-postproc-run`.

Minimal ISMIP7 config example:

```yaml
tool: flatten

input:       /data/standalone/output/
output_dir:  /data/ismip7/
exe_path:    /opt/bisicles/flatten2d.Linux.64.g++.gfortran.OPT.ex

epsg:        3031
level:       0

ice_sheet:   AIS
experiment:  ctrl_proj

source_id:   BristolGlaciology
ism_id:      BISICLES
ism_member_id: m001
esm_id:      standalone
forcing_member_id: f001
set_counter: C001

contact_name:  "Your Name"
contact_email: your.name@institution.ac.uk
reference_year: 1850
```

### ISMIP7 CLI examples

Flatten:

```bash
bike-ismip7-postproc-flatten \
    --input  /data/standalone/output/ \
    --output-dir /data/ismip7/ \
    --exe-path /opt/bisicles/flatten2d.Linux.64.g++.gfortran.OPT.ex \
    --epsg 3031 \
    --x0 -3072000.0 --y0 -3072000.0 \
    --level 0 \
    --ice-sheet AIS \
    --experiment ctrl_proj \
    --source-id BristolGlaciology \
    --ism-id BISICLES \
    --ism-member-id m001 \
    --esm-id standalone \
    --forcing-member-id f001 \
    --set-counter C001 \
    --contact-name "Your Name" \
    --contact-email your.name@institution.ac.uk \
    --reference-year 1850
```

Diagnostics:

```bash
bike-ismip7-postproc-diagnostics \
    --input  /data/standalone/output/ \
    --output-dir /data/ismip7/ \
    --exe-path /opt/bisicles/diagnostics2d.Linux.64.g++.gfortran.OPT.ex \
    --ice-sheet AIS \
    --experiment ctrl_proj \
    --source-id BristolGlaciology \
    --ism-id BISICLES \
    --ism-member-id m001 \
    --esm-id standalone \
    --forcing-member-id f001 \
    --set-counter C001 \
    --contact-name "Your Name" \
    --contact-email your.name@institution.ac.uk \
    --reference-year 1850
```

### ISMIP7 configuration options

**`source_id`**  [required]
  Modelling group name.  No underscores, dots or special characters.
  E.g. `"BristolGlaciology"`.

**`ism_id`**  [default: `"BISICLES"`]
  ISM name and version.  No underscores, dots or special characters.
  E.g. `"BISICLES"` or `"BISICLESv3-2"`.

**`ism_member_id`**  [default: `"m001"`]
  ISM choice variant (mNNN format).  Increment for each new initial state
  or ice sheet model parameter choice.

**`esm_id`**  [default: `"standalone"`]
  CMIP ESM used to produce the forcing.  Use `"standalone"` for idealised
  or observed forcing.  E.g. `"CESM2-WACCM"`, `"MRI-ESM2-0"`, `"ERA5"`.

**`forcing_member_id`**  [default: `"f001"`]
  Forcing choice variant (fNNN format).  Increment for each new downscaling
  method or melt parameterisation choice.

**`set_counter`**  [default: `"C001"`]
  Set counter linking all files for a single run in the submission
  spreadsheet.  Use `Cnnn` for CORE, `Ennn` for ESM, `Pnnn` for PPE.

**`group`**  [default: same as `source_id`]
  Modelling group name written as the ISMIP7 mandatory global attribute
  `group`.

**`contact_name`**  [required for submission]
  Contact person name(s), written as global attribute `contact_name`.

**`contact_email`**  [required for submission]
  Contact person email(s), written as global attribute `contact_email`.

### ISMIP7 grid origins and regridding

Output coordinates must match the grid origin used when running BISICLES.
Typical standalone BISICLES simulations use the **standard BISICLES grid**,
so those origins must be used:

```text
Ice sheet  EPSG   x0 (standard BISICLES)  y0 (standard BISICLES)
---------  ----   ----------------------  ----------------------
GrIS       3413   -654650.0               -3385950.0
AIS        3031   -3072000.0              -3072000.0
```

These are the defaults when `--x0`/`--y0` are not supplied, so in most
cases you do not need to pass them explicitly.

**To submit to ISMIP7 the output must be on the ISMIP7 standard grid**,
which has different origins and resolution.  Regrid directly within this
tool using the `x_min`/`x_max`/`nx`/`y_min`/`y_max`/`ny` options (requires
`scipy`):

```text
Domain  EPSG   x_min         x_max        nx    y_min         y_max        ny
------  ----   -----------   ----------   ---   -----------   ----------   ---
AIS     3031   -3040000.0    3040000.0    761   -3040000.0    3040000.0    761
GrIS    3413   -720000.0     960000.0     1681  -3450000.0    -570000.0    2881
```

Do **not** pass ISMIP7 origins to `x0`/`y0` — those control the input grid
origin (which must match your BISICLES run), not the output grid.

---

## UKESM-specific notes

- Always use `calendar: 360_day` for UKESM-coupled runs.  BISICLES internal
  time is always 0 in UKESM coupling; simulation time is read from the
  filename instead.

- UKESM BISICLES filenames follow the pattern
  `bisicles_SUITE_FREQ_START-END_plot-ICESHEET.hdf5`,
  e.g. `bisicles_dx030c_1y_18880101-18890101_plot-AIS.hdf5`.
  The ice sheet, simulation time, and time-mean status are all parsed
  automatically from this pattern.

- Time-mean CF output files follow the pattern
  `bisicles_SUITE_FREQ_START-END_plot.CF-ICESHEET.hdf5`.
  These must be processed with the `flatten` workflow; the `diagnostics`
  workflow will raise an error if given one of these files.
