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
No new native diagnostic API or coefficient-export interface is introduced.

The code uses the repository's MIT license; existing upstream notices remain
in their source files.
