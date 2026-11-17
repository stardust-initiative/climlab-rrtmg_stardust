"""Native regressions for effective-layer beam propagation.

These tests call the compiled public wrapper, not a Python imitation of spcvmc.
The nearly empty atmosphere isolates prescribed layer coefficients. A dense
boundary-value solve supplies independent interface-flux expectations.
"""
import numpy as np
import pytest
from climlab_rrtmg import rrtmg_sw


def _column(layers, *, enabled=True, cloud=0.0, pressure=1.e-6, mu=0.5, cloud_tau=0., cpdair=1004.):
    n = len(layers)
    shape = (1, n)
    zero = np.zeros(shape)
    bands = np.zeros((1, n, 14))
    cloudbands = np.zeros((14, 1, n))
    bounds = np.linspace(pressure, pressure / 100., n + 1)[None, :]
    coefficients = [np.broadcast_to(np.array([x[i] for x in layers])[None, :, None],
                                   (1, n, 14)).copy() for i in range(4)]
    args = dict(ncol=1, nlay=n, permuteseed=150, irng=1, n_ens=1,
                col_by_col=1, do_seed_permutation=0, icld=1 if np.any(cloud) else 0,
                ispec=1, iaer=0, play=(bounds[:, :-1]+bounds[:, 1:])/2,
                plev=bounds, tlay=np.full(shape, 250.),
                tlev=np.full((1, n+1), 250.), tsfc=np.array([250.]),
                cldfrac=np.full(shape, cloud), ciwp=zero, clwp=zero,
                reic=np.full(shape, 30.), relq=np.full(shape, 10.),
                tauc=cloudbands+cloud_tau, ssac=cloudbands+.999, asmc=cloudbands+.6,
                fsfc=cloudbands+.36, h2ovmr=zero, o3vmr=zero, co2vmr=zero,
                ch4vmr=zero, n2ovmr=zero, o2vmr=zero,
                asdir=np.array([.12]), aldir=np.array([.12]),
                asdif=np.array([.35]), aldif=np.array([.35]),
                kmodts=2, coszen=np.array([mu]), adjes=np.ones(1),
                dyofyr=0, scon=1361., isolvar=-1,
                inflgsw=0, iceflgsw=1, liqflgsw=1,
                tauaer=bands, ssaaer=bands, asmaer=bands,
                ecaer=np.zeros((1,n,6)), bndsolvar=np.ones(14),
                indsolvar=np.ones(2), solcycfrac=1.,
                add_aero_layer=int(enabled),
                r_mu=coefficients[0], t_mu=coefficients[1],
                r_bar=coefficients[2], t_bar=coefficients[3])
    rrtmg_sw.climlab_rrtmg_sw_ini(cpdair)
    return rrtmg_sw.climlab_rrtmg_sw_ensemble(**args)


def _boundary_solution(layers):
    # Effective t_mu propagates the direct component in this model. All
    # coefficients here are top-down; native atmospheric/output grids reverse.
    n = len(layers)
    beam = np.r_[1., np.cumprod([x[1] for x in layers])]
    matrix = np.zeros((2*(n+1), 2*(n+1)))
    rhs = np.zeros(2*(n+1))
    matrix[0, 0] = 1.  # no diffuse light incident from space
    matrix[1, -1] = 1.
    matrix[1, n] = -.35
    rhs[1] = .12 * beam[-1]
    row = 2
    for i, (r, t, rd, td) in enumerate(layers):
        matrix[row, i+1] = 1.
        matrix[row, i] = -td
        matrix[row, n+2+i] = -rd
        row += 1
        matrix[row, n+1+i] = 1.
        matrix[row, i] = -rd
        matrix[row, n+2+i] = -td
        rhs[row] = r*beam[i]
        row += 1
    solution = np.linalg.solve(matrix, rhs)
    return solution[n+1:], solution[:n+1]+beam


@pytest.mark.parametrize('layers', [
    [(0.2,0.8,0.3,0.7), (0.1,0.9,0.4,0.6)],
    [(0.05,0.55,0.1,0.5), (0.1,0.65,0.2,0.55)],
])
@pytest.mark.parametrize('cloud', [0., .5, 1., [0.,1.], [1.,0.]])
def test_effective_beam_matches_boundary_solution(layers, cloud):
    result = _column(layers, cloud=cloud)
    up, down = _boundary_solution(layers)
    incident = result[1][0,-1]
    np.testing.assert_allclose(result[0][0,::-1]/incident, up, atol=2e-7, rtol=2e-7)
    np.testing.assert_allclose(result[1][0,::-1]/incident, down, atol=2e-7, rtol=2e-7)
    absorption = np.diff((result[1]-result[0])[0])
    assert absorption.min() >= -1e-6
    if all(abs(r+t-1.)<1e-12 and abs(rd+td-1.)<1e-12 for r,t,rd,td in layers):
        assert abs(absorption).max() < 1e-5
    else:
        assert absorption.min() > 0.
    for total, spectral in [(0,6),(1,7),(3,8),(4,9)]:
        assert np.isfinite(result[spectral]).all()
        assert result[spectral].min() >= -1e-8
        np.testing.assert_allclose(result[spectral].sum(axis=-1), result[total], atol=1e-9)


@pytest.mark.parametrize('cloud', [0., .5, 1.])
def test_transparent_effective_layers_are_identity(cloud):
    layers=[(0.,1.,0.,1.)]*3
    baseline=_column(layers,enabled=False,cloud=cloud,pressure=1000.)
    modified=_column(layers,cloud=cloud,pressure=1000.)
    for old,new in zip(baseline,modified):
        np.testing.assert_allclose(new,old,rtol=0.,atol=5e-10)


def test_disabled_coefficients_cannot_change_native_fluxes():
    original=_column([(0.,1.,0.,1.)]*2,enabled=False,pressure=1000.)
    unused=_column([(.2,.8,.3,.7),(.1,.9,.4,.6)],enabled=False,pressure=1000.)
    for a,b in zip(original,unused):
        np.testing.assert_array_equal(a,b)


@pytest.mark.parametrize('cloud', [[0.,1.], [.5,.5], [1.,0.]])
def test_cloud_overlap_preserves_clear_flux_and_spectral_consistency(cloud):
    layers=[(.2,.8,.3,.7),(.05,.55,.1,.5)]
    clear=_column(layers,pressure=1000.)
    cloudy=_column(layers,cloud=cloud,pressure=1000.,cloud_tau=2.)
    # Clear-sky diagnostics must not inherit cloudy masks or cloud optical depth.
    for field in [3,4,5,8,9]:
        np.testing.assert_array_equal(clear[field],cloudy[field])
    assert not np.allclose(clear[0],cloudy[0])
    for total,spectral in [(0,6),(1,7),(3,8),(4,9)]:
        assert np.isfinite(cloudy[spectral]).all()
        assert cloudy[spectral].min() >= -1e-8
        np.testing.assert_allclose(cloudy[spectral].sum(axis=-1),cloudy[total],atol=1e-9)
    assert np.diff((cloudy[1]-cloudy[0])[0]).min() >= -1e-6


def test_specific_heat_initializer_has_finite_heating_and_correct_scaling():
    layers = [(.05, .85, .1, .8), (.04, .9, .08, .85)]
    try:
        normal = _column(layers, pressure=1000., cpdair=1004.)
        doubled = _column(layers, pressure=1000., cpdair=2008.)
    finally:
        rrtmg_sw.climlab_rrtmg_sw_ini(1004.)
    for value in normal + doubled:
        assert np.isfinite(value).all()
    for index in (0, 1, 3, 4, 6, 7, 8, 9):
        np.testing.assert_array_equal(normal[index], doubled[index])
    for index in (2, 5):
        assert np.max(np.abs(normal[index])) > 0.
        np.testing.assert_allclose(doubled[index], normal[index]/2., rtol=2e-14, atol=1e-14)
