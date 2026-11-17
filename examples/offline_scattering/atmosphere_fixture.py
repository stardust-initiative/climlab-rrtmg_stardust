"""Synthetic finite-pressure atmosphere, not a measured or manuscript profile.
Pressure is hPa; temperatures K; gas fields are volume mixing ratios.
Atmospheric arrays are bottom-to-top; effective aerosol arrays are top-to-bottom.
No clouds are supported by this diagnostic fixture.
"""
import numpy as np
from climlab_rrtmg import rrtmg_sw

def pressure_bounds(n, pressure=1000.):
    return np.geomspace(pressure, .1, n + 1)


def column(layers, *, enabled=True, cloud=0.0, pressure=1000., mu=0.5, cloud_tau=0.):
    if np.any(cloud) or cloud_tau != 0:
        raise ValueError("This diagnostic fixture supports clear sky only")
    if pressure <= .1:
        raise ValueError("Surface pressure must exceed .1 hPa")
    n = len(layers)
    shape = (1, n)
    zero = np.zeros(shape)
    bands = np.zeros((1, n, 14))
    cloudbands = np.zeros((14, 1, n))
    bounds = pressure_bounds(n, pressure)[None, :]
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
                fsfc=cloudbands+.36, h2ovmr=np.geomspace(.01, 1e-6, n)[None, :], o3vmr=np.full(shape, 5e-6), co2vmr=np.full(shape, 420e-6),
                ch4vmr=np.full(shape, 1.8e-6), n2ovmr=np.full(shape, .32e-6), o2vmr=np.full(shape, .2095),
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
    rrtmg_sw.climlab_rrtmg_sw_ini(1004.)
    return rrtmg_sw.climlab_rrtmg_sw_ensemble(**args)
