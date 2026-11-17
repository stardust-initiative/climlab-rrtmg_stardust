# Offline direct/diffuse scattering example

This small NumPy example composes layered scattering responses and checks them
against an independently assembled system of interface flux equations. It
requires no atmospheric files, optical tables, compiled backend, or network access.

From this directory:

```bash
python -m pip install numpy
python offline_scattering.py
python -m unittest -v test_offline_scattering.py
```

The JSON output reports numerical agreement and a controlled demonstration of
when the separate direct-beam transmission affects the reflected flux.

## Coefficients and boundary conditions

Layers are ordered from top to bottom. Each `Slab` has seven dimensionless
coefficients mapping incident downward direct beam `B`, downward diffuse flux
`D`, and upward diffuse flux `U` to its outgoing fluxes:

```text
U_top    = reflect_beam*B + reflect_down*D + transmit_up*U
D_bottom = source_beam*B + transmit_down*D + reflect_up*U
B_bottom = transmit_beam*B
```

`Slab.symmetric_cell(r_direct, r_diffuse, t_direct, t_diffuse, beam)` constructs
a vertically symmetric cell. `t_direct` is total downward transmission under
direct illumination, including both the beam and forward diffuse flux; its
forward diffuse source is therefore `t_direct - beam`. An aggregated slab can
have unequal upward and downward diffuse responses, so it retains all seven
coefficients.

`join` eliminates the internal diffuse interface fluxes when stacking two slabs.
`stack` repeats this operation. `solve_interfaces` independently solves all
interface equations with prescribed illumination at the top and bottom; it does
not call `join` or `stack`. All inputs in the example are explicit synthetic
coefficients, chosen to obey energy bounds. Coefficients are not inferred from
an atmospheric state. The solver rejects nonfinite coefficients and undefined
reflection feedback; it does not clip, regularize, or prove arbitrary supplied
coefficients physically admissible.

`reflected` couples an aggregate slab to a lower region characterized by separate
direct and diffuse reflectances. Feedback between the slab and that lower region
is retained. A Lambertian lower boundary is the special case where the two
reflectances are equal.

## What these checks establish

For one symmetric cell above a Lambertian boundary of reflectance `A`, the
reflected fraction is

```text
R = r_direct + t_diffuse * t_direct * A / (1 - r_diffuse * A).
```

The separate beam fraction cancels from this expression. The example changes
that fraction while holding total direct transmission fixed: the Lambertian
answer stays unchanged, but a lower region with different direct and diffuse
responses gives a different answer. Thus a single Lambertian-layer comparison
cannot test a general atmosphere's direct/diffuse partition. Its mathematical
limit remains useful; applying it as a complete atmospheric validation would
omit vertical structure, directional lower responses and surrounding transfer.

Other tests cover transparent and absorbing limits, energy closure in a
conservative synthetic stack, composition associativity, and agreement with
the independent dense solve under illumination from both boundaries.

These are **algebraic checks for supplied coefficients**, not independent
physical validation of optical coefficients, atmospheric radiative transfer,
or an aerosol parameterization. No optical-depth continuation of clipped
coefficients is used. Native and offline agreement with matched coefficients
would likewise validate the transfer implementation rather than the physical
accuracy of the coefficients.

## Optional comparison with native adding

With `gfortran` on PATH, run from this directory:

```bash
python compare_native.py
```

This compiles the repository's existing `parkind.f90` and
`rrtmg_sw_vrtqdr.f90` sources together with a small test-only C bridge, inside
a temporary directory. It does not rebuild or replace the installed package.
No upstream Fortran implementation is copied into the example. Existing source
notices remain intact. The script reports source hashes and compares every
interface's upward and total downward flux with the independent dense solve,
and top reflection with the generalized slab composition. Five tiny synthetic
cases include a two-layer scattering column with distinct direct/diffuse lower
reflectances, absorption, transparency, a Lambertian boundary, and a black
surface. Maximum expected errors are at floating-point roundoff.

This bridges the offline example to the actual native adding routine, but
does not test how `spcvmc` constructs its input coefficients or demonstrate the
beam patch itself. The full native regression below covers that separate step.
The optional build currently targets a Unix-like system with GNU Fortran; it
adds no runtime dependency or native API to the package.

## Native implementation regression

The repository's separate native regression exercises the compiled modified
RRTMG pathway. After installing the package and its test dependencies, run from
outside the source checkout:

```bash
python -m pytest --pyargs climlab_rrtmg.tests.test_modified_direct_beam
```

That test complements this standalone example. Running only this NumPy script
does not test the installed backend or compare a real atmospheric column.
The installed native API remains unchanged; the opt-in source-copy exporter
below is a separate diagnostic workflow.

The code uses the repository's MIT license; existing upstream notices remain
in their source files.

## Optional atmospheric coefficient export

For a matched-coefficient comparison with a finite-pressure atmosphere, the
opt-in diagnostic below constructs a six-layer **synthetic clear-sky** column
with water vapour, ozone, CO2, CH4, N2O, O2, native Rayleigh scattering, a surface,
and two prescribed effective aerosol layers. It exports native per-g-point
operators and derives the lower atmospheric direct and diffuse responses
`A_b` and `A_d`. These are atmospheric subcolumn responses at the aerosol-region
boundary, **not physical surface albedos**. This is reusable reader machinery,
not a reproduction of a particular paper figure or its atmospheric inputs.

Use an isolated development environment with GNU Fortran, a C compiler, Python,
NumPy, Meson (>=0.63.3), meson-python and Ninja. The example currently supports
little-endian Unix-like systems with GNU Fortran. From the repository root:

```bash
# Requires an activated, isolated development environment; not a production env.
python -m pip install numpy 'meson>=0.63.3' meson-python ninja
DEMO_WORK=$(mktemp -d)
python examples/offline_scattering/prepare_diagnostic.py "$DEMO_WORK/source"
python -m pip install --no-deps --no-build-isolation \
  --config-settings=compile-args=-j2 --target "$DEMO_WORK/install" "$DEMO_WORK/source"
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONPATH="$DEMO_WORK/install" \
  python examples/offline_scattering/atmospheric_diagnostic.py \
  --build-manifest "$DEMO_WORK/source/diagnostic_build.json" \
  --output "$DEMO_WORK/results"
```

Source preparation requires a new destination outside the source repository. It
copies the source and adds diagnostic stream writes only to that copy, retaining
existing notices. The separately installed diagnostic package is selected only
for this subprocess via `PYTHONPATH`; the normal installed backend and public
API stay unchanged. `OMP_NUM_THREADS` must equal 1: concurrent writes to this
simple diagnostic stream are unsupported. The fixture rejects cloudy inputs;
although native clear-sky coefficients exist inside cloudy calculations, such
calls are outside this workflow's supported scope.

The script checks the input/exported effective-layer indexing explicitly,
derives the occupied aerosol-region boundaries, and uses the exported corrected
beam directly. It does not repair or multiply the exported beam a second time.
It checks that this beam equals clean atmospheric beam times the effective
transmission. The lower region is composed independently; its `A_b`, `A_d` are
checked against native lower-subcolumn responses. Full up/down interface fluxes
are compared to the independent dense solve; spectral weights restore W/m² and
are checked against the public wrapper's broadband clear-sky outputs. A paired
identity-aerosol calculation checks zero forcing, and clean/aerosol lower-region
responses must match because the aerosol lies above that boundary.

Outputs include:

- `summary.json`: source-preparation manifest, backend/fixture hashes, pressure
  interfaces, aerosol-region indices, errors and paired aerosol forcing.
- `clean.npz`, `aerosol.npz`, `identity.npz`: lower direct/diffuse responses,
  spectral bands and weights, native layer operators and aerosol coefficients,
  offline/native interface fluxes and broadband fluxes.
- Corresponding `.bin` streams: versioned raw diagnostic exports, with hashes
  recorded in the summary. The reader rejects unknown versions, malformed
  dimensions, invalid spectral ordering, nonfinite values and truncated files.

The NPZ interface order is top-to-bottom. Coefficients and normalized interface
fluxes are dimensionless; spectral weights and weighted broadband fluxes use
W/m². The native atmospheric input grid is bottom-to-top while the effective
layer coefficients are supplied top-to-bottom. Pressure boundaries are saved
from the fixture's actual grid. Tests for the parser require no compiled backend:

```bash
cd examples/offline_scattering
python -m unittest -v test_offline_scattering.py test_atmospheric_diagnostic.py
```

The exporter is deliberately separate from the production wrapper. Agreement
checks its matched transfer algebra; it does not independently establish the
physical accuracy of gas or aerosol coefficients, clouds, or optical closures.
