"""Compare offline slab algebra with the repository's actual Fortran adding.

Copyright (c) 2026 Stardust Labs. MIT license; see repository LICENSE.
Requires gfortran and NumPy. Builds only a temporary test library, not RRTMG.
"""
import ctypes
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

import numpy as np

from offline_scattering import Slab, reflected, solve_interfaces, stack


# (r_direct, r_diffuse, t_direct, t_diffuse, beam), top to bottom.
CASES = {
    "lambertian_cell": ([(.2, .3, .8, .7, .8)], (.25, .25)),
    "two_scattering_layers": ([(.2, .3, .8, .7, .8), (.1, .4, .9, .6, .65)], (.2, .35)),
    "absorbing_layers": ([(0., 0., .6, .5, .6), (.1, .2, .65, .55, .4)], (.2, .3)),
    "transparent_layer": ([(0., 0., 1., 1., 1.), (.1, .4, .9, .6, .65)], (.25, .25)),
    "black_surface": ([(.2, .3, .8, .7, .8), (.1, .4, .9, .6, .65)], (0., 0.)),
}


def main():
    compiler = shutil.which("gfortran")
    if compiler is None:
        raise SystemExit("Optional native comparison requires gfortran on PATH")
    here = Path(__file__).resolve().parent
    source = here.parents[1] / "climlab_rrtmg/rrtmg_sw/rrtmg_sw_v4.0/gcm_model"
    sources = [source / "modules/parkind.f90", source / "src/rrtmg_sw_vrtqdr.f90", here / "native_shim.f90"]
    rows = []
    with tempfile.TemporaryDirectory(prefix="offline-scattering-") as build:
        library = Path(build) / "adding.so"
        subprocess.run([compiler, "-O0", "-fPIC", "-shared", "-fcheck=all",
                        *(str(p) for p in sources), "-o", str(library)], cwd=build, check=True)
        dll = ctypes.CDLL(str(library))
        evaluate = dll.evaluate_adding
        array = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS")
        evaluate.argtypes = [ctypes.c_int] + [array]*7
        evaluate.restype = None
        for name, (coefficients, lower) in CASES.items():
            n = len(coefficients)
            inputs = [np.array([row[i] for row in coefficients] + [0.], dtype=np.float64) for i in range(5)]
            inputs[0][-1], inputs[1][-1] = lower
            down, up = np.zeros(n+1), np.zeros(n+1)
            evaluate(n, *inputs, down, up)
            layers = [Slab.symmetric_cell(*row) for row in coefficients]
            # An opaque final slab implements the lower direct/diffuse boundary.
            boundary = Slab(lower[0], lower[1], 0., 0., 0., 0., 0.)
            beam, diffuse, expected_up = solve_interfaces(layers + [boundary])
            expected_down = (beam + diffuse)[:-1]
            expected_up = expected_up[:-1]
            error = float(max(np.max(abs(down-expected_down)), np.max(abs(up-expected_up))))
            top_error = abs(reflected(stack(layers), *lower)-up[0])
            if error > 2e-13 or top_error > 2e-13:
                raise AssertionError(f"{name}: interface error={error}, composed error={top_error}")
            if min(down.min(), up.min()) < -2e-13:
                raise AssertionError(f"{name}: negative flux for admissible synthetic inputs")
            rows.append(dict(case=name, interface_max_error=error, composed_top_error=float(top_error),
                             native_down=down.tolist(), native_up=up.tolist()))
    print(json.dumps({"scope": "matched-coefficient adding algebra, not full atmospheric optics",
                      "source_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
                      "cases": rows}, indent=2))


if __name__ == "__main__":
    main()
