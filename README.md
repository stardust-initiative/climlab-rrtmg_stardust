# climlab-rrtmg (Stardust fork)

Brian Rose, University at Albany
Stardust modifications by [stardust-initiative](https://github.com/stardust-initiative)

## About

This is a fork of [climlab-rrtmg](https://github.com/climlab/climlab-rrtmg), a stand-alone Python wrapper for the [RRTMG](http://rtweb.aer.com/rrtm_frame.html) radiation modules.

The primary use-case is to drive the RRTMG radiation processes in [climlab](https://climlab.readthedocs.io/),
but it can be used as a stand-alone radiation model if you are familiar with the
RRTMG Fortran interface. This is a lightweight wrapper that emulates the Fortran
interface as closely as possible.

Currently we are wrapping RRTMG_LW v4.85 and RRTMG_SW v4.0. The Fortran source code
is included in this repository. The latest versions 5.0 of the RRTMG source code
are available on GitHub [here](https://github.com/AER-RC)

This wrapper includes the following modifications:

- RRTMG_LW is modified to allow reporting of OLR components in spectral bands,
 as illustrated using climlab
[here](https://climlab.readthedocs.io/en/latest/courseware/Spectral_OLR_with_RRTMG.html).
This modification is strictly diagnostic does not change any other behavior of RRTMG_LW.
- RRTMG_SW is modified to allow gridpoint-specific values of the input parameter `adjes` 
which is the adjustment to total solar irradiance. 

The parameter `adjes` is used within RRTMG_SW to account for time-of-year adjustments 
to Earth-Sun distance. Latest versions of climlab compute non-uniform values of 
`irradiance_factor` to compensate for [different diurnal averages of the solar zenith angle](https://climlab.readthedocs.io/en/latest/api/climlab.radiation.insolation.html). 
climlab's `irradiance_factor` is mapped to `adjes` when this RRTMG driver is called from climlab.


## Stardust modifications

This fork extends the upstream [climlab-rrtmg](https://github.com/climlab/climlab-rrtmg)
with the following features for stratospheric aerosol climate modeling:

### OpenMP parallelization of column radiation calls

The LW and SW Fortran drivers (`Driver.f90`) are modified to parallelize
the radiation calculation over independent atmospheric columns using
OpenMP `parallel do` directives. This provides significant speedup on
multi-core machines when the model has many latitude points, since each
column's radiation calculation is independent. Set `OMP_NUM_THREADS`
to control the number of threads.

### Ensemble cloud sampling

Support for stochastic cloud overlap sampling with configurable
`n_rrtmg_repeat` samples per column per radiation timestep. At each
timestep, multiple independent cloud realizations are drawn and the
radiation fluxes are averaged over the ensemble. This reduces noise
from the McICA (Monte Carlo Independent Column Approximation) cloud
overlap scheme used in RRTMG.

The `do_col_by_col` flag enables column-by-column processing where
each column can have its own independent cloud random seed, and
`do_seed_permutation` controls whether the seed changes between
timesteps.

### Spectral flux exposure and thin-layer aerosol injection (SW)

The SW driver exposes spectral (per-band) flux outputs and supports
injection of aerosol optical properties into individual model layers.
This enables computation of spectrally-resolved shortwave radiative
forcing from stratospheric aerosol layers.

## Installation

Pre-built binaries for many platforms are available from [conda-forge](https://conda-forge.org).

To install in the current environment:
```
conda install climlab-rrtmg --channel conda-forge
```
or create a self-contained environment:
```
conda create --name my_env python=3.13 climlab-rrtmg --channel conda-forge
conda activate my_env
```

See below for instructions on how to build from source.

## Example usage

You can import the modules into a Python session with
```
from climlab_rrtmg import rrtmg_lw, rrtmg_sw
```

The main RRTMG drivers are exposed through
```
rrtmg_lw.climlab_rrtmg_lw()
```
and
```
rrtmg_sw.climlab_rrtmg_sw()
```

Please see the directory `climlab_rrtmg/tests/` directory in this repository
for working examples that set up all the necessary input arrays and call the drivers.

## Building from source

Here are instructions to create a build environment (including Fortran compiler)
with conda and build using f2py.

To build:
```
conda env create --file ./ci/requirements.yml
conda activate rrtmg_build_env
python -m pip install . --no-deps -vv
```

To run tests, do this from any directory other than the climlab-rrtmg repo:
```
pytest -v --pyargs climlab_rrtmg
```

## Version history

- Version 0.5.0 (May 2026) is the first Stardust-modifications release. Adds OpenMP column-radiation parallelization (LW and SW), ensemble cloud sampling for McICA noise reduction, and spectral-band SW flux exposure with thin-layer aerosol injection. See the *Stardust modifications* section above. Incorporates upstream changes through climlab-rrtmg v0.4.2.
- Version 0.4.1 (late Feburary 2025) is a bug fix, now properly exposing some constants defined in the Fortran code to the Python wrapper.
- Version 0.4 (released February 2025) is a major refactor of the build system so this package can run on Python 3.12 and above. The build now uses [meson](https://mesonbuild.com/).
- Version 0.3 (released January 2025) includes the modification for grid-variable `adjes` described above. It is designed to work with climlab v0.9 and higher.
- Version 0.2 is the first public release (April 2022).
The Python wrapper code has been extracted from
[climlab v0.7.13](https://github.com/brian-rose/climlab/releases/tag/v0.7.13).