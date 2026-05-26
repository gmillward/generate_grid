# GIP Apex Coordinate Grid Generation

## Overview

This repository generates the GIP static coordinate file
(`GIP_apex_coords_etc.YYYY.0.format`) used by the GT-GIP coupled
ionosphere-thermosphere model. It defines every flux tube in the GIP
computational grid — their geographic positions, magnetic field vectors,
and coordinate transforms between geographic and apex-magnetic frames.

### What "apex coordinates" means

A simple dipole magnetic field has the analytic field-line equation
`r = L · R_earth · cos²(λ)`, where λ is magnetic latitude and L is the
field-line label. Richmond's apex coordinate system moves beyond this to
use the **real IGRF magnetic field**: instead of an analytic formula, each
field line is traced numerically to its highest point — the *apex* — and
that apex location becomes the coordinate reference. This is the key
contribution of the Richmond (1995) library (`apxntrpb4lf.f`, `apex.f`,
`ggrid.f`, `magfld.f`): it precomputes the IGRF-based apex transforms on a
global geographic grid (`Apex_grid_data`) so that subsequent coordinate
conversions can be done by fast interpolation rather than tracing field
lines from scratch each time.

**Implication for the visualisation scripts:** `plot_tubes_mp1.py` uses the
dipole formula to draw the arch shape of each tube (latitude vs height).
The L-values themselves come from the real IGRF-computed
`tiegcm_defined_apex_heights`, so the tube *boundaries* are accurate — but
the *shape of each arch as drawn* is a dipole approximation. A fully
accurate geographic-latitude plot would require calling APXQ2G for every
grid point (i.e. reading the actual lat/lon arrays from the GIP output
file).

## Two-program pipeline

Controlled by `runscript.sh`:

```
apex2000_prog  <  input_date  >  outfile        # Stage 1
apex_prog      <  input_date  >> outfile        # Stage 2
```

### Stage 1 — `apex2000_prog` (`apex2000.f`)

Evaluates the IGRF magnetic field on a dense global geographic grid and
writes the intermediate file `Apex_grid_data` (~31 MB). This is a
pre-computed lookup table of apex coordinate quantities at every
geographic location — the Richmond apex-coordinate library
(`apxntrpb4lf.f`, `ggrid.f`, `magfld.f`) is used internally.

**Known limitation:** the IGRF coefficients in `magfld.f` only extend to
epoch ~2000 (NGRF=8 epochs, 1965–2005 per the comment in
`generate_apex_coordinates.f`). Generating grids for post-2005 dates
requires updating `magfld.f` with IGRF-13 coefficients (which cover up
to 2025).

### Stage 2 — `apex_prog` (`generate_apex_coordinates.f` + supporting files)

Reads `Apex_grid_data` and `tiegcm_defined_apex_heights`, then builds
the full 3D GIP grid tube by tube. Calls `calc_apex_params_2d_2.f90`
to compute derived field quantities and write the final output file.

Supporting routines: `apex.f`, `apxntrpb4lf.f`, `divve.f`, `ggrid.f`,
`magfld.f`.

## Grid structure

### Latitudinal shells (`nlp = 67`) — **fixed by electrodynamics coupling**

Defined by `tiegcm_defined_apex_heights` — a 97-row table originating
from NCAR's TIEGCM model. Each row gives a magnetic colatitude and the
corresponding apex height (km) and apex radius (km) for a field line
rooted at 90 km altitude at that latitude.

**These magnetic latitudes are not freely adjustable.** The
electrodynamics solver works on the 2D (mp, lp) magnetic
longitude/latitude grid: it integrates physical parameters (conductivity,
currents etc.) along each field line to produce 2D field-line-integrated
quantities, then solves the electrodynamics on that 2D grid. It is the
2D (mp, lp) grid that couples to TIEGCM — not the individual along-tube
points. Changing the latitudinal shells would break that coupling.
`tiegcm_defined_apex_heights` should be treated as a fixed external
constraint, not a redesign target.

**The along-tube points are freely redesignable.** They are used only
by GIP's plasma diffusion solver (the tridiagonal O+/H+ solvers) and the
field-line integration that feeds the 2D electrodynamics grid. More
points, or differently spaced points, along a tube do not affect the 2D
grid structure — only the accuracy and stability of the along-tube
computation.

The code reads the first 48 rows (southern hemisphere) and selects all
tubes with L-value between 1 and 4 (L = apex_radius / R_earth ≈
1 + apex_height_km / 6371.2). It then doubles the count by inserting
interpolated L-values between each pair, giving `nlp=67` latitudinal
shells:

- `lp=1` — outermost tube, L ≈ 3.5, apex at ~16,000 km
- `lp=34` — mid-latitude, apex at a few hundred km
- `lp=67` — innermost tube, L ≈ 1 (apex just above 90 km, equatorial
  E-region)

### Magnetic longitude (`nmp = 80`)

80 evenly-spaced magnetic meridians, one every 4.5°.

### Along-tube points (`npts = 583` intermediate; `npts2 = 13813` in GIP)

**`npts = 583` is not the number of points along a single tube.** It is
the size of the intermediate 3D generation array `(npts, nmp, nlp)` used
during grid construction, with the apex of every tube centred at
`n_mid_point = (npts+1)/2 = 292`. The actual number of points per tube
varies by `lp` — from many hundreds for the outermost large tubes down
to just 1–2 for the innermost near-equatorial tubes.

**GIP uses a packed 2D representation `(npts2, nmp)`.** For each
magnetic longitude `mp`, all 67 tubes are laid out end-to-end in a
single 1D array of total length `npts2 = 13813`:

```
mp column:  [ tube lp=1 (outermost) | tube lp=2 | ... | tube lp=67 (innermost) ]
              IN(mp,1)   IS(mp,1)                        IS(mp,67) = 13813
```

`IN(mp,lp)` and `IS(mp,lp)` are the start and end indices of tube `lp`
within this packed column. Since per-tube point counts depend only on
apex height (i.e. `lp`), not on longitude, `IN` and `IS` are the same
for all `mp`.

This avoids the waste of a full 3D array `param(npts_max, nmp, nlp)`
where most entries would be empty — large outer tubes have many more
points than small inner ones, so a uniform per-tube dimension would waste
enormous memory for the innermost tubes.

## Along-tube height spacing formula (current)

Points are placed by iterating from the footpoint (90 km) toward the
apex:

```fortran
height = height + (iht-1) * sqrt(HA_metres - height) * factor
```

where `factor = 0.2` for `lp < 11` (high-latitude tubes) and `0.4`
elsewhere. The iteration stops when `height > HA` (the apex altitude).

**Key property:** step size starts at zero (iht=1) and grows with each
step, so the grid is *coarsest near the footpoints and densest near the
apex*. This is the inverse of what is physically needed — the sharpest
O+ density gradients occur in the E-F transition region (90–200 km),
close to the footpoints.

### Actual point counts (confirmed by `analyse_grid.py`)

Point counts per hemisphere and average step sizes for selected shells:

| lp  | apex (km) | pts/hemi | avg dh (km) | note                         |
|-----|-----------|----------|-------------|------------------------------|
|   1 |   18784   |   292    |   64.2      | outermost; factor=0.2        |
|  10 |    7304   |   230    |   31.5      | last lp with factor=0.2      |
|  11 |    6541   |   158    |   41.1      | factor switches to 0.4 here  |
|  20 |    2504   |   123    |   19.8      |                              |
|  40 |     271   |    64    |    2.9      |                              |
|  47 |     158   |    50    |    1.4      |                              |
|  49 |     143   |    47    |    1.2      | GT-GIP problem zone          |
|  55 |     117   |    40    |    0.7      |                              |
|  67 |      94   |    24    |    0.15     | innermost; 3.5 km range      |

**Important correction**: The sqrt formula does NOT produce "very few
points" for the low-apex tubes. lp=49 has 47 pts/hemi with ~1 km average
spacing. The issue is the **spacing distribution**, not the count.

**The formula concentrates points at both the footpoint and the apex** —
step size starts at 0 (iht=1 adds zero), grows to a peak in the middle of
the tube, then shrinks back to 0 near the apex (since `sqrt(HA-height)→0`).

For high-apex tubes (lp=1, apex 18784 km), the step size in the 90–200 km
E-F transition grows from 0 to ~1–2 km within a few steps and then
continues growing — reaching ~60–80 km by 1000 km altitude. This gives
coarse resolution throughout the plasmasphere.

The GT-GIP O+ solver instabilities in the lp≈47–50 band at the
dawn/dusk terminator (see `../gt-gip/CLAUDE.md`) are caused by large
Peclet numbers and steep O+ gradients at the E-F boundary — not by
insufficient grid points. The Scharfetter-Gummel implementation in
GT-GIP addresses this, with 11 residual failures in the hardest cases.

## Key files

| File | Role |
|------|------|
| `npts.h` | Grid dimensions: `npts=583`, `nmp=80`, `nlp=67` |
| `tiegcm_defined_apex_heights` | 97-row table defining latitudinal shells by apex height/L-value |
| `input_date` | Epoch for grid generation (currently `2000.0`) |
| `generate_apex_coordinates.f` | Main Stage 2 program; contains height-spacing formula |
| `calc_apex_params_2d_2.f90` | Computes derived field quantities; packs to `npts2=13813` |
| `apex2000.f` | Stage 1; contains IGRF coefficients (currently limited to ~2005) |
| `apex.f`, `apxntrpb4lf.f`, `divve.f`, `ggrid.f`, `magfld.f` | Richmond apex-coordinate library |

## Compiler

Both programs now build with `gfortran -std=legacy -ffixed-line-length-132 -w -O2`
(converted from ifort on the `development` branch). Three Fortran 77
string-literal continuation incompatibilities were fixed with labeled FORMAT
statements (two in `apex2000.f`, one in `apxntrpb4lf.f`). Both stages run
successfully and produce the correct output (npts2=13813).

**Note on Apex_grid_data**: APXWRA opens the file with `STATUS='unknown'`
(does not truncate). If an old multi-epoch file exists, re-running stage 1
only overwrites the first epoch's records, leaving stale data that confuses
APXRDA. Always delete `Apex_grid_data` before regenerating.

## Planned 2026 redesign

The original grid was designed ~30 years ago under tight CPU and memory
constraints. With 2026 hardware the grid can be substantially refined.
Priority changes:

### 1. Replace the along-tube height formula

The sqrt formula must be replaced with one that concentrates points
where the physics demands it — the E-F transition region (90–200 km).
Candidates:

- **Geometric progression from footpoint:** `dh(i) = dh_0 * r^i` with
  `dh_0` set to ~1–2 km at 90 km, growing toward the apex. This gives
  fine E-region resolution while remaining manageable near the apex.
- **Piecewise:** fixed fine spacing (e.g. 2 km) up to ~300 km, then
  coarser above. Simple and predictable.
- **Density-weighted:** target a fixed number of grid points per
  scale-height of O+. Requires an a-priori density profile but gives
  physically optimal resolution.

The formula is in `generate_apex_coordinates.f` at the `do iht = 1,1000`
loop (~line 398), replicated for the northern and southern hemisphere
traversals.

### 2. Increase `npts2` (the packed array size)

Finer along-tube spacing means more total points across all tubes for
each longitude column. The packed array size `npts2 = 13813` in
`calc_apex_params_2d_2.f90` must be increased to accommodate this, as
must the matching `NPTS` parameter in GT-GIP's
`GIP_ionosphere_plasmasphere.f90`. The intermediate generation array
size `npts = 583` in `npts.h` may also need to increase if the largest
tubes (lp=1, apex ~16,000 km) get more points.

The new `npts2` can be estimated by summing the new per-tube point
counts across all 67 lp shells once the new height formula is defined.

### 3. Update IGRF to IGRF-13

`magfld.f` contains IGRF coefficients through epoch ~2005. IGRF-13
(released 2019) extends coverage to 2025 with 5-year model updates.
Required to generate grids for present-day simulations.

### 4. Latitudinal range — treat with caution

The current outer boundary is L=4 (~16,000 km apex). The magnetic
latitudes in `tiegcm_defined_apex_heights` are coupled to GIP's
electrodynamics solver (inherited from NCAR's TIEGCM) and must not be
changed without also updating that solver. Any extension of the
latitudinal range (e.g. to L=8 for better plasmasphere coverage) would
require coordinated changes across both the grid generation code and the
GIP electrodynamics module.

### 5. gfortran compatibility

Both programs need to be converted from ifort to gfortran (updating
`Makefile`, `runscript.sh`, and any ifort-specific syntax in the source).
