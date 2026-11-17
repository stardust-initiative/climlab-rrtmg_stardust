"""Compose direct/diffuse slab responses independently of native RRTMG.

Copyright (c) 2026 Stardust Labs. Distributed under the repository MIT license.
All coefficients and fluxes in this example are dimensionless.
"""

from dataclasses import dataclass, fields
import json
import math
from typing import Sequence

import numpy as np


@dataclass(frozen=True)
class Slab:
    """Response to incident downward beam B, downward diffuse D, upward U.

    U_top = reflect_beam*B + reflect_down*D + transmit_up*U
    D_bottom = source_beam*B + transmit_down*D + reflect_up*U
    B_bottom = transmit_beam*B

    Aggregating layers generally makes upward and downward responses unequal.
    This class represents coefficients; it does not calculate physical optics.
    """

    reflect_beam: float
    reflect_down: float
    transmit_up: float
    source_beam: float
    transmit_down: float
    reflect_up: float
    transmit_beam: float

    def __post_init__(self):
        if not all(math.isfinite(getattr(self, f.name)) for f in fields(self)):
            raise ValueError("Slab coefficients must be finite")

    @classmethod
    def identity(cls):
        """Transparent layer, with no scattering or absorption."""
        return cls(0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0)

    @classmethod
    def symmetric_cell(cls, r_direct, r_diffuse, t_direct, t_diffuse, beam):
        """Build a cell; t_direct includes both beam and diffuse transmission."""
        return cls(r_direct, r_diffuse, t_diffuse, t_direct - beam,
                   t_diffuse, r_diffuse, beam)


def denominator(value):
    """Require finite reflection feedback; never clip or regularize it."""
    if not math.isfinite(value) or value <= 0:
        raise ValueError("Reflection feedback requires a positive denominator")
    return value


def join(upper: Slab, lower: Slab) -> Slab:
    """Eliminate the two diffuse fluxes at a shared internal interface."""
    a, b = upper, lower
    den = denominator(1.0 - a.reflect_up * b.reflect_down)
    xb = (a.source_beam + a.reflect_up * b.reflect_beam * a.transmit_beam) / den
    xd = a.transmit_down / den
    xu = a.reflect_up * b.transmit_up / den
    return Slab(
        a.reflect_beam + a.transmit_up * (b.reflect_beam * a.transmit_beam + b.reflect_down * xb),
        a.reflect_down + a.transmit_up * b.reflect_down * xd,
        a.transmit_up * (b.transmit_up + b.reflect_down * xu),
        b.source_beam * a.transmit_beam + b.transmit_down * xb,
        b.transmit_down * xd,
        b.reflect_up + b.transmit_down * xu,
        a.transmit_beam * b.transmit_beam,
    )


def stack(layers: Sequence[Slab]) -> Slab:
    """Aggregate layers ordered from top to bottom."""
    result = Slab.identity()
    for layer in layers:
        result = join(result, layer)
    return result


def reflected(slab: Slab, lower_direct: float, lower_diffuse: float) -> float:
    """Upward top flux for unit incident beam and a two-response lower region.

    The lower region returns lower_direct*B + lower_diffuse*D. A Lambertian
    boundary is the special case lower_direct == lower_diffuse.
    """
    if not all(math.isfinite(x) and 0 <= x <= 1 for x in (lower_direct, lower_diffuse)):
        raise ValueError("Lower reflectances must be finite and within [0, 1]")
    den = denominator(1.0 - slab.reflect_up * lower_diffuse)
    down = (slab.source_beam + slab.reflect_up * lower_direct * slab.transmit_beam) / den
    return slab.reflect_beam + slab.transmit_up * (
        lower_direct * slab.transmit_beam + lower_diffuse * down
    )


def solve_interfaces(layers: Sequence[Slab], beam_top=1.0, diffuse_top=0.0, diffuse_bottom=0.0):
    """Independent dense solve for every diffuse interface flux.

    Boundary illumination is prescribed at both ends. This does not call join
    or stack. Unknowns are D[1:n+1], U[0:n], with D[0] and U[n] prescribed.
    Returned B, D, U each have n+1 entries ordered top to bottom.
    """
    if not all(math.isfinite(x) for x in (beam_top, diffuse_top, diffuse_bottom)):
        raise ValueError("Boundary fluxes must be finite")
    n = len(layers)
    beam = np.empty(n + 1)
    beam[0] = beam_top
    for i, layer in enumerate(layers):
        beam[i + 1] = beam[i] * layer.transmit_beam
    if n == 0:
        return beam, np.array([diffuse_top]), np.array([diffuse_bottom])
    matrix, rhs = np.zeros((2*n, 2*n)), np.zeros(2*n)
    for i, layer in enumerate(layers):
        matrix[2*i, n+i] = 1.0
        matrix[2*i+1, i] = 1.0
        rhs[2*i] = layer.reflect_beam * beam[i]
        rhs[2*i+1] = layer.source_beam * beam[i]
        if i == 0:
            rhs[2*i] += layer.reflect_down * diffuse_top
            rhs[2*i+1] += layer.transmit_down * diffuse_top
        else:
            matrix[2*i, i-1] = -layer.reflect_down
            matrix[2*i+1, i-1] = -layer.transmit_down
        if i == n-1:
            rhs[2*i] += layer.transmit_up * diffuse_bottom
            rhs[2*i+1] += layer.reflect_up * diffuse_bottom
        else:
            matrix[2*i, n+i+1] = -layer.transmit_up
            matrix[2*i+1, n+i+1] = -layer.reflect_up
    solution = np.linalg.solve(matrix, rhs)
    return beam, np.r_[diffuse_top, solution[:n]], np.r_[solution[n:], diffuse_bottom]


def demonstration():
    """Tiny deterministic examples; no external data, native code or downloads."""
    layers = [Slab(.03, .12, .7, .1, .65, .08, .5),
              Slab(.04, .09, .75, .12, .7, .06, .6)]
    total = stack(layers)
    _, down, up = solve_interfaces(layers, 1.0, .2, .15)
    composed_up = total.reflect_beam + total.reflect_down*.2 + total.transmit_up*.15
    composed_down = total.source_beam + total.transmit_down*.2 + total.reflect_up*.15
    # Hold total direct transmission fixed while changing its beam/diffuse split.
    first = Slab.symmetric_cell(.08, .04, .9, .94, .8)
    second = Slab.symmetric_cell(.08, .04, .9, .94, .5)
    expected = .08 + .94*.9*.3/(1-.04*.3)
    return {
        "interface_solve_max_error": max(abs(up[0]-composed_up), abs(down[-1]-composed_down)),
        "lambertian_formula_error": abs(reflected(first, .3, .3)-expected),
        "lambertian_beam_split_difference": reflected(second, .3, .3)-reflected(first, .3, .3),
        "two_response_beam_split_difference": reflected(second, .2, .4)-reflected(first, .2, .4),
    }


if __name__ == "__main__":
    print(json.dumps(demonstration(), indent=2))
